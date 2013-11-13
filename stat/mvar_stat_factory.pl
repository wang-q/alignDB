#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::WriteExcel;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);
use AlignDB::Position;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# stat parameter
my $run = $Config->{stat}{run};

my $outfile;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=s'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'output=s'   => \$outfile,
    'run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

my $dbh = $write_obj->dbh;
my $pos_obj = AlignDB::Position->new( dbh => $dbh );

#----------------------------#
# count freq
#----------------------------#
my $all_freq;
{

    my $sql_query = q{
            SELECT DISTINCT COUNT(q.query_id) + 1
            FROM  query q, sequence s
            WHERE q.seq_id = s.seq_id
            GROUP BY s.align_id
        };
    my $sth = $dbh->prepare($sql_query);

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
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{Type_name AVG_length STD_length COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            query_name => $query_name,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'indel', 'indel_coding' ) ) {
        my $query_name = 'Coding_type';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'indel', 'indel_repeats' ) ) {
        my $query_name = 'Repeat_type';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'indel', 'indel_slippage' ) ) {
        my $query_name = 'Slippage_type';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'indel', 'indel_type' ) ) {
        my $query_name = 'Indel_type';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'indel', 'indel_freq' ) ) {
        my $query_name = 'Indel_freq';
        $sheet_row++;
        my $sql_query = q{
            SELECT indel_freq,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY indel_freq
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    {
        my $query_name = 'Length_group';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_basic
#----------------------------------------------------------#
my $snp_basic = sub {
    my $sheet_name = 'snp_basic';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{Type_name COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            query_name => $query_name,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'snp', 'snp_coding' ) ) {
        my $query_name = 'Coding_type';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'snp', 'snp_repeats' ) ) {
        my $query_name = 'Repeat_type';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'snp', 'snp_cpg' ) ) {
        my $query_name = 'CpG_type';
        $sheet_row++;
        my $sql_query = q{
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
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    if ( $write_obj->check_column( 'snp', 'snp_freq' ) ) {
        my $query_name = 'SNP_freq';
        $sheet_row++;
        my $sql_query = q{
            SELECT snp_freq,
                   COUNT(*) `COUNT`
            FROM snp s
            GROUP BY snp_freq
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_list
#----------------------------------------------------------#
my $indel_list = sub {
    my $sheet_name = 'indel_list';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{indel_id taxon_id chr_name indel_start indel_end
            indel_length indel_seq indel_gc indel_freq indel_occured indel_type
            indel_slippage indel_coding indel_repeats};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $sql_query = q{
            SELECT 
                i.align_id, i.indel_id, a.taxon_id, a.chr_name, i.indel_start,
                i.indel_end, i.indel_length, i.indel_seq, i.indel_gc,
                i.indel_freq, i.indel_occured, i.indel_type, i.indel_slippage,
                i.indel_coding, i.indel_repeats
            FROM
                indel i
                    INNER JOIN
                (SELECT 
                    se.align_id, ta.taxon_id, c.chr_name, se.chr_start
                FROM
                    sequence se, target t, chromosome c, taxon ta
                WHERE
                    se.seq_id = t.seq_id
                        AND se.chr_id = c.chr_id
                        AND c.taxon_id = ta.taxon_id) a ON i.align_id = a.align_id
            ORDER BY a.chr_name, a.chr_start
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i ( 3, 4 ) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            ($sheet_row) = $write_obj->write_row_direct(
                $sheet,
                {   row       => \@row,
                    sheet_row => $sheet_row,
                    sheet_col => $sheet_col,
                }
            );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_list
#----------------------------------------------------------#
my $snp_list = sub {
    my $sheet_name = 'snp_list';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{snp_id taxon_id chr_name snp_pos isw_distance
            mutant_to snp_freq snp_occured
            snp_coding snp_repeats snp_cpg};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $sql_query = q{
            SELECT 
                s.align_id, s.snp_id, a.taxon_id, a.chr_name, s.snp_pos,
                i.isw_distance, s.mutant_to, s.snp_freq, s.snp_occured,
                s.snp_coding, s.snp_repeats, s.snp_cpg
            FROM
                snp s
                    INNER JOIN
                (SELECT 
                    se.align_id, ta.taxon_id, c.chr_name, se.chr_start
                FROM
                    sequence se, target t, chromosome c, taxon ta
                WHERE
                    se.seq_id = t.seq_id
                        AND se.chr_id = c.chr_id
                        AND c.taxon_id = ta.taxon_id) a ON s.align_id = a.align_id
                    LEFT JOIN
                isw i ON s.isw_id = i.isw_id
            ORDER BY a.chr_name, a.chr_start
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i (3) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            ($sheet_row) = $write_obj->write_row_direct(
                $sheet,
                {   row       => \@row,
                    sheet_row => $sheet_row,
                    sheet_col => $sheet_col,
                }
            );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_codon_list
#----------------------------------------------------------#
my $snp_codon_list = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'snp', 'snp_codon_pos' ) ) {
        return;
    }

    my $sheet_name = 'snp_codon_list';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{snp_id name mutant_to freq occured target codon_pos syn nsy };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $sql_query = q{
            SELECT 
                s.align_id, s.snp_id, a.chr_name, s.snp_pos, s.mutant_to,
                s.snp_freq, s.snp_occured, s.snp_codon_pos, s.snp_syn, s.snp_nsy
            FROM
                snp s
                    INNER JOIN
                (SELECT 
                    se.align_id, ta.taxon_id, c.chr_name, se.chr_start
                FROM
                    sequence se, target t, chromosome c, taxon ta
                WHERE
                    se.seq_id = t.seq_id
                        AND se.chr_id = c.chr_id
                        AND c.taxon_id = ta.taxon_id) a ON s.align_id = a.align_id
            WHERE s.exon_id IS NOT NULL
            ORDER BY a.chr_name, a.chr_start
        };
        my $sth = $dbh->prepare($sql_query);
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

            ($sheet_row) = $write_obj->write_row_direct(
                $sheet,
                {   row => [
                        $row[0], $name,   $row[3], $row[4], $row[5],
                        $target, $row[6], $row[7], $row[8],
                    ],
                    sheet_row => $sheet_row,
                    sheet_col => $sheet_col,
                }
            );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gene_list
#----------------------------------------------------------#
my $gene_list = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'gene', 'gene_syn' ) ) {
        return;
    }

    my $sheet_name = 'gene_list';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{ gene_id chr_name gene_start gene_end
            gene_stable_id gene_external_name gene_biotype gene_strand
            gene_is_full gene_multitrans gene_multiexons gene_syn gene_nsy };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $sql_query = q{
            SELECT 
                w.align_id, g.gene_id, a.chr_name, w.window_start,
                w.window_end, g.gene_stable_id, g.gene_external_name,
                g.gene_biotype, g.gene_strand, g.gene_is_full,
                g.gene_multitrans, g.gene_multiexons, g.gene_syn, g.gene_nsy
            FROM
                gene g
                    INNER JOIN
                window w ON g.window_id = w.window_id
                    INNER JOIN
                (SELECT 
                    se.align_id, ta.taxon_id, c.chr_name, se.chr_start
                FROM
                    sequence se, target t, chromosome c, taxon ta
                WHERE
                    se.seq_id = t.seq_id
                        AND se.chr_id = c.chr_id
                        AND c.taxon_id = ta.taxon_id) a ON w.align_id = a.align_id
            ORDER BY a.chr_name , a.chr_start
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i ( 2, 3 ) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            ($sheet_row) = $write_obj->write_row_direct(
                $sheet,
                {   row       => \@row,
                    sheet_row => $sheet_row,
                    sheet_col => $sheet_col,
                }
            );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";

};

#----------------------------------------------------------#
# worksheet -- strain_list
#----------------------------------------------------------#
my $strain_list = sub {
    my $sheet_name = 'strain_list';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{strain ins del snp };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
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
            my $sth = $dbh->prepare($sql_query);
            $sth->execute($type);

            while ( my ( $id, $string ) = $sth->fetchrow_array ) {
                my @chars = split //, $string;
                if ( $all_freq != scalar @chars ) {
                    warn "indel_id [$id] occured string errors\n";
                }

                for my $i ( 0 .. $#chars ) {
                    if ( $chars[$i] eq 'o' ) {
                        $count_of->{$type}[$i]++;
                    }
                    elsif ( $chars[$i] eq 'x' ) {

                        # OK, do nothing
                    }
                    else {
                        warn "indel_id [$id] occured string errors\n";
                    }
                }
            }
        }
    }

    {    # count snps
        my $sql_query = q{
            SELECT  s.snp_id, s.snp_occured
            FROM    snp s
            WHERE   s.snp_occured <> 'unknown'
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my ( $id, $string ) = $sth->fetchrow_array ) {
            my @chars = split //, $string;
            if ( $all_freq != scalar @chars ) {
                warn "snp_id [$id] occured string errors\n";
            }

            for my $i ( 0 .. $#chars ) {
                if ( $chars[$i] eq 'o' ) {
                    $count_of->{S}[$i]++;
                }
                elsif ( $chars[$i] eq 'x' ) {

                    # OK, do nothing
                }
                else {
                    warn "indel_id [$id] occured string errors\n";
                }
            }
        }
    }

    {    # write contents
        for my $i ( 1 .. $all_freq ) {
            ($sheet_row) = $write_obj->write_row_direct(
                $sheet,
                {   row => [
                        $i,
                        $count_of->{I}[ $i - 1 ],
                        $count_of->{D}[ $i - 1 ],
                        $count_of->{S}[ $i - 1 ],
                    ],
                    sheet_row => $sheet_row,
                    sheet_col => $sheet_col,
                }
            );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

foreach my $n (@tasks) {
    if ( $n == 1 )  { &$indel_basic;    next; }
    if ( $n == 2 )  { &$snp_basic;      next; }
    if ( $n == 3 )  { &$indel_list;     next; }
    if ( $n == 4 )  { &$snp_list;       next; }
    if ( $n == 5 )  { &$snp_codon_list; next; }
    if ( $n == 8 )  { &$gene_list;      next; }
    if ( $n == 10 ) { &$strain_list;    next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    mvar_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    mvar_stat_factory.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --output          output filename
       --run             run special analysis

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
