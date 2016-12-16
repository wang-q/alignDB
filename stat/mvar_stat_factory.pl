#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::ToXLSX;

use List::Util;
use List::MoreUtils::PP;

use lib "$FindBin::Bin/../lib";
use AlignDB::Position;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

mvar_stat_factory.pl - Variable lists for alignDB

=head1 SYNOPSIS

    perl mvar_stat_factory.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL IP/Domain
        --port          INT     MySQL port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --outfile   -o  STR     outfile filename
        --run       -r  STR     run special analysis
        --index                 add an index sheet

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server   = $conf->{database}{server} ),
    'port=i'       => \( my $port     = $conf->{database}{port} ),
    'db|d=s'       => \( my $db       = $conf->{database}{db} ),
    'username|u=s' => \( my $username = $conf->{database}{username} ),
    'password|p=s' => \( my $password = $conf->{database}{password} ),
    'output|o=s'   => \( my $outfile ),
    'run|r=s'      => \( my $run      = $conf->{stat}{run} ),
    'index'        => \( my $add_index_sheet, ),
) or Getopt::Long::HelpMessage(1);

my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $db, $server, $port;

$outfile = "$db.mvar.xlsx" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
}
else {
    $run =~ s/\"\'//s;
    my $set = AlignDB::IntSpan->new;
    if ( AlignDB::IntSpan->valid($run) ) {
        $set   = $set->add($run);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $run;
        $set->add(@tasks);
    }

    my $runlist = $set->runlist;
    $outfile =~ s/(\.xlsx)$/.$runlist$1/;
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Do stat for $db...");

my DBI $dbh = DBI->connect( $dsn, $username, $password )
    or die $DBI::errstr;
my $toxlsx = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
);
my $pos_obj = AlignDB::Position->new( dbh => $dbh );

#----------------------------#
# count freq
#----------------------------#
my $all_freq;
{
    my DBI $sth = $dbh->prepare(
        q{
        SELECT DISTINCT COUNT(s.seq_id) + 1
        FROM  sequence s
        WHERE 1 = 1
        AND s.seq_role = "Q"
        GROUP BY s.align_id
        }
    );

    my @counts;
    $sth->execute;
    while ( my ($count) = $sth->fetchrow_array ) {
        push @counts, $count;
    }
    if ( scalar @counts > 1 ) {
        die "Database corrupts, freqs are not consistent\n";
    }

    $all_freq = $counts[0];
}

#----------------------------------------------------------#
# worksheet -- indel_basic
#----------------------------------------------------------#
my $indel_basic = sub {
    my $sheet_name = 'indel_basic';
    my $sheet;
    $toxlsx->row(0);
    $toxlsx->column(1);

    my @names = qw{Type_name AVG_length STD_length COUNT};
    {    # header
        $sheet = $toxlsx->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    if ( $toxlsx->check_column( 'indel', 'indel_coding' ) ) {
        my $query_name = 'Coding_type';
        my $sql_query  = q{
            SELECT CASE i.indel_coding
                     WHEN 0 THEN 'non_coding'
                     WHEN 1 THEN 'coding'
                     WHEN NULL THEN 'NULL'
                     ELSE 'semi_coding'
                   END `coding`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY coding
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    if ( $toxlsx->check_column( 'indel', 'indel_repeats' ) ) {
        my $query_name = 'Repeat_type';
        my $sql_query  = q{
            SELECT CASE i.indel_repeats
                     WHEN 0 THEN 'non_repeat'
                     WHEN 1 THEN 'repeat' 
                     WHEN NULL THEN 'NULL'
                     ELSE 'semi_repeat'
                   END `repeat`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY `repeat`
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    if ( $toxlsx->check_column( 'indel', 'indel_slippage' ) ) {
        my $query_name = 'Slippage_type';
        my $sql_query  = q{
            SELECT CASE i.indel_slippage
                     WHEN 0 THEN 'non_slippage'
                     WHEN 1 THEN 'slippage' 
                     WHEN NULL THEN 'NULL'
                     ELSE 'semi_slippage'
                   END `slippage`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY slippage
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    if ( $toxlsx->check_column( 'indel', 'indel_type' ) ) {
        my $query_name = 'Indel_type';
        my $sql_query  = q{
            SELECT CASE i.indel_type
                     WHEN 'C' THEN 'Complex'
                     WHEN 'D' THEN 'Deletion' 
                     WHEN 'I' THEN 'Insertion'
                     ELSE 'WRONG'
                   END `type`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY type
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    if ( $toxlsx->check_column( 'indel', 'indel_freq' ) ) {
        my $query_name = 'Indel_freq';
        my $sql_query  = q{
            SELECT indel_freq,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY indel_freq
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    {
        my $query_name = 'Length_group';
        my $sql_query  = q{
            SELECT CASE
                     WHEN indel_length BETWEEN 1 AND 5 THEN '1[1-5]'
                     WHEN indel_length BETWEEN 6 AND 10 THEN '2[6-10]'
                     WHEN indel_length BETWEEN 11 AND 50 THEN '3[11-50]'
                     WHEN indel_length BETWEEN 51 AND 300 THEN '4[51-300]'
                     ELSE 'WRONG'
                   END `length`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY length
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_basic
#----------------------------------------------------------#
my $snp_basic = sub {
    my $sheet_name = 'snp_basic';
    my $sheet;
    $toxlsx->row(0);
    $toxlsx->column(1);

    my @names = qw{Type_name COUNT};
    {    # header
        $sheet = $toxlsx->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    if ( $toxlsx->check_column( 'snp', 'snp_coding' ) ) {
        my $query_name = 'Coding_type';
        my $sql_query  = q{
            SELECT CASE s.snp_coding
                     WHEN 0 THEN 'non_coding'
                     WHEN 1 THEN 'coding'
                     WHEN NULL THEN 'NULL'
                     ELSE 'WRONG'
                   END `coding`,
                   COUNT(*) `COUNT`
            FROM snp s
            GROUP BY coding
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    if ( $toxlsx->check_column( 'snp', 'snp_repeats' ) ) {
        my $query_name = 'Repeat_type';
        my $sql_query  = q{
            SELECT CASE s.snp_repeats
                     WHEN 0 THEN 'non_repeat'
                     WHEN 1 THEN 'repeat' 
                     WHEN NULL THEN 'NULL'
                     ELSE 'WRONG'
                   END `repeat`,
                   COUNT(*) `COUNT`
            FROM snp s
            GROUP BY `repeat`
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    if ( $toxlsx->check_column( 'snp', 'snp_cpg' ) ) {
        my $query_name = 'CpG_type';
        my $sql_query  = q{
            SELECT CASE s.snp_cpg
                     WHEN 0 THEN 'non_CpG'
                     WHEN 1 THEN 'CpG' 
                     WHEN NULL THEN 'NULL'
                     ELSE 'WRONG'
                   END `cpg`,
                   COUNT(*) `COUNT`
            FROM snp s
            GROUP BY cpg
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    # contents
    if ( $toxlsx->check_column( 'snp', 'snp_freq' ) ) {
        my $query_name = 'SNP_freq';
        my $sql_query  = q{
            SELECT snp_freq,
                   COUNT(*) `COUNT`
            FROM snp s
            GROUP BY snp_freq
        };
        $toxlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_list
#----------------------------------------------------------#
my $indel_list = sub {
    my $sheet_name = 'indel_list';
    my $sheet;
    $toxlsx->row(0);
    $toxlsx->column(0);

    my @names = qw{indel_id common_name chr_name indel_start indel_end
        indel_length indel_seq indel_gc indel_freq indel_occured indel_type
        indel_slippage indel_coding indel_repeats};
    {    # header
        $sheet = $toxlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # contents
        my $sql_query = q{
            SELECT 
                i.align_id, i.indel_id, s.common_name, s.chr_name, i.indel_start,
                i.indel_end, i.indel_length, i.indel_seq, i.indel_gc,
                i.indel_freq, i.indel_occured, i.indel_type, i.indel_slippage,
                i.indel_coding, i.indel_repeats
            FROM
                indel i
                    INNER JOIN
                sequence s ON i.align_id = s.align_id
            ORDER BY s.chr_name, s.chr_start
        };
        my DBI $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i ( 3, 4 ) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            $toxlsx->write_row( $sheet, { row => \@row, } );
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_list
#----------------------------------------------------------#
my $snp_list = sub {
    my $sheet_name = 'snp_list';
    my $sheet;
    $toxlsx->row(0);
    $toxlsx->column(0);

    my @names = qw{snp_id common_name chr_name snp_pos isw_distance
        snp_mutant_to snp_freq snp_occured
        snp_coding snp_repeats snp_cpg};
    {    # header
        $sheet = $toxlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # contents
        my $sql_query = q{
            SELECT 
                s.align_id, s.snp_id, se.common_name, se.chr_name, s.snp_pos,
                i.isw_distance, s.snp_mutant_to, s.snp_freq, s.snp_occured,
                s.snp_coding, s.snp_repeats, s.snp_cpg
            FROM
                snp s
                    INNER JOIN
                sequence se ON s.align_id = se.align_id
                    LEFT JOIN
                isw i ON s.isw_id = i.isw_id
            ORDER BY se.chr_name, se.chr_start
        };
        my DBI $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i (3) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            $toxlsx->write_row( $sheet, { row => \@row, } );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_codon_list
#----------------------------------------------------------#
my $snp_codon_list = sub {

    unless ( $toxlsx->check_column( 'snp', 'snp_codon_pos' ) ) {
        return;
    }

    my $sheet_name = 'snp_codon_list';
    my $sheet;
    $toxlsx->row(0);
    $toxlsx->column(0);

    my @names = qw{snp_id name mutant_to freq occured target codon_pos syn nsy };
    {    # header
        $sheet = $toxlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # contents
        my $sql_query = q{
            SELECT 
                s.align_id, s.snp_id, se.chr_name, s.snp_pos, s.snp_mutant_to,
                s.snp_freq, s.snp_occured, s.snp_codon_pos, s.snp_syn, s.snp_nsy
            FROM
                snp s
                    INNER JOIN
                sequence se ON s.align_id = se.align_id
            ORDER BY se.chr_name, se.chr_start
        };
        my DBI $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i (2) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            my $name = "$row[1]:$row[2]";
            my $target = substr $row[5], 0, 1;

            $toxlsx->write_row(
                $sheet,
                {   row => [
                        $row[0], $name,   $row[3], $row[4], $row[5],
                        $target, $row[6], $row[7], $row[8],
                    ],
                }
            );
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gene_list
#----------------------------------------------------------#
my $gene_list = sub {

    unless ( $toxlsx->check_column( 'gene', 'gene_syn' ) ) {
        return;
    }

    my $sheet_name = 'gene_list';
    my $sheet;
    $toxlsx->row(0);
    $toxlsx->column(0);

    my @names = qw{ gene_id chr_name gene_start gene_end
        gene_stable_id gene_external_name gene_biotype gene_strand
        gene_is_full gene_multitrans gene_multiexons gene_subs gene_indel
        gene_pi gene_syn gene_nsy };
    {    # header
        $sheet = $toxlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # contents
        my $sql_query = q{
            SELECT 
                w.align_id, g.gene_id, s.chr_name, w.window_start, w.window_end,
                g.gene_stable_id, g.gene_external_name, g.gene_biotype,
                g.gene_strand, g.gene_is_full, g.gene_multitrans,
                g.gene_multiexons, w.window_differences, w.window_indel,
                w.window_pi, g.gene_syn, g.gene_nsy
            FROM
                gene g
                    INNER JOIN
                window w ON g.window_id = w.window_id
                    INNER JOIN
                sequence s ON w.align_id = s.align_id
            ORDER BY s.chr_name, s.chr_start
        };
        my DBI $sth = $dbh->prepare($sql_query);
        $sth->execute();

        # write full genes
        my @data;
        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i ( 2, 3 ) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            if ( $row[8] == 0 ) {
                push @data, \@row;
            }
            else {
                $toxlsx->write_row( $sheet, { row => \@row, } );
            }
        }

        # combine partial genes
        my @ids = map { $_->[4] } @data;
        @ids = List::MoreUtils::PP::uniq(@ids);
        printf " " x 4 . "%d partial gene records\n", scalar @data;
        printf " " x 4 . "%d partial genes\n",        scalar @ids;

        for my $id (@ids) {
            my @records = grep { $_->[4] eq $id } @data;

            my $gene_id    = join ",", ( map { $_->[0] } @records );
            my $chr_name   = $records[0]->[1];
            my $gene_start = List::Util::min( map { $_->[2] } @records );
            my $gene_end   = List::Util::max( map { $_->[3] } @records );
            my $gene_subs  = List::Util::sum( map { $_->[11] } @records );
            my $gene_indel = List::Util::sum( map { $_->[12] } @records );

            my ($gene_pi);
            my @records_pi = grep { defined $_->[13] } @records;
            my $total_length = List::Util::sum( map { $_->[3] - $_->[2] + 1 } @records_pi );
            for my $record (@records_pi) {
                my $partial_length = $record->[3] - $record->[2] + 1;
                $gene_pi += $record->[13] * $partial_length / $total_length;
            }

            my ( $gene_syn, $gene_nsy );
            my @records_syn_nsy = grep { defined $_->[14] } @records;
            my $effective_length
                = List::Util::sum( map { $_->[3] - $_->[2] + 1 } @records_syn_nsy );
            if ($effective_length) {
                for my $record (@records_syn_nsy) {
                    my $partial_length = $record->[3] - $record->[2] + 1;
                    $gene_syn += $record->[14] * $partial_length / $total_length;
                    $gene_nsy += $record->[15] * $partial_length / $total_length;
                }
            }

            $toxlsx->write_row(
                $sheet,
                {   row => [
                        $gene_id,         $chr_name,        $gene_start,       $gene_end,
                        $records[0]->[4], $records[0]->[5], $records[0]->[6],  $records[0]->[7],
                        $records[0]->[8], $records[0]->[9], $records[0]->[10], $gene_subs,
                        $gene_indel,      $gene_pi,         $gene_syn,         $gene_nsy,
                    ],
                }
            );
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- strain_list
#----------------------------------------------------------#
my $strain_list = sub {
    my $sheet_name = 'strain_list';
    my $sheet;
    $toxlsx->row(0);
    $toxlsx->column(0);

    my @names = qw{strain ins del snp };
    {    # header
        $sheet = $toxlsx->write_header( $sheet_name, { header => \@names } );
    }

    my $count_of = {
        I => [],
        D => [],
        S => [],
    };
    {    # count indels
        my $sql_query = q{
            SELECT  i.indel_id, i.indel_occured
            FROM    indel i
            WHERE   i.indel_occured <> 'unknown'
            AND     i.indel_type = ?
        };

        for my $type ( 'I', 'D' ) {
            my DBI $sth = $dbh->prepare($sql_query);
            $sth->execute($type);

            while ( my ( $id, $string ) = $sth->fetchrow_array ) {
                my @chars = split //, $string;
                if ( $all_freq != scalar @chars ) {
                    warn "indel_id [$id] occured string errors\n";
                }

                for my $i ( 0 .. $#chars ) {
                    if ( $chars[$i] eq '1' ) {
                        $count_of->{$type}[$i]++;
                    }
                    elsif ( $chars[$i] eq '0' ) {    # OK, do nothing
                    }
                    else {
                        warn "indel_id [$id] occured string errors\n";
                    }
                }
            }
        }
    }

    {                                                # count snps
        my $sql_query = q{
            SELECT  s.snp_id, s.snp_occured
            FROM    snp s
            WHERE   s.snp_occured <> 'unknown'
        };
        my DBI $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my ( $id, $string ) = $sth->fetchrow_array ) {
            my @chars = split //, $string;
            if ( $all_freq != scalar @chars ) {
                warn "snp_id [$id] occured string errors\n";
            }

            for my $i ( 0 .. $#chars ) {
                if ( $chars[$i] eq '1' ) {
                    $count_of->{S}[$i]++;
                }
                elsif ( $chars[$i] eq '0' ) {    # OK, do nothing
                }
                else {
                    warn "indel_id [$id] occured string errors\n";
                }
            }
        }
    }

    {                                            # contents
        for my $i ( 1 .. $all_freq ) {
            $toxlsx->write_row(
                $sheet,
                {   row => [
                        $i,
                        $count_of->{I}[ $i - 1 ],
                        $count_of->{D}[ $i - 1 ],
                        $count_of->{S}[ $i - 1 ],
                    ],
                }
            );
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

for my $n (@tasks) {
    if ( $n == 1 )  { &$indel_basic;    next; }
    if ( $n == 2 )  { &$snp_basic;      next; }
    if ( $n == 8 )  { &$gene_list;      next; }
    if ( $n == 10 ) { &$strain_list;    next; }
    if ( $n == 51 ) { &$indel_list;     next; }
    if ( $n == 52 ) { &$snp_list;       next; }
    if ( $n == 53 ) { &$snp_codon_list; next; }
}

if ($add_index_sheet) {
    $toxlsx->add_index_sheet;
    print "Sheet [INDEX] has been generated.\n";
}

$stopwatch->end_message;
exit;

sub mean {
    @_ = grep { defined $_ } @_;
    return unless @_;
    return $_[0] unless @_ > 1;
    return List::Util::sum(@_) / scalar(@_);
}

__END__
