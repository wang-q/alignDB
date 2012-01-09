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

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# stat parameter
my $run = $Config->{stat}->{run};

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
            SELECT  i.indel_id, ta.taxon_id, c.chr_name,
                    i.indel_start, i.indel_end, i.indel_length, i.indel_seq,
                    i.indel_gc, i.indel_freq, i.indel_occured, i.indel_type,
                    i.indel_slippage, i.indel_coding, i.indel_repeats
            FROM    indel i, align a, sequence s,
                    target t, chromosome c, taxon ta
            WHERE   a.align_id = i.align_id AND
                    a.align_id = s.align_id AND
                    s.seq_id = t.seq_id AND
                    s.chr_id = c.chr_id AND
                    c.taxon_id = ta.taxon_id
            ORDER BY c.chr_name, i.indel_start
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
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
            SELECT  s.snp_id, ta.taxon_id, c.chr_name,
                    s.snp_pos, i.isw_distance,
                    s.mutant_to, s.snp_freq, s.snp_occured,
                    s.snp_coding, s.snp_repeats, s.snp_cpg
            FROM    snp s, align a, isw i, sequence se,
                    target t, chromosome c, taxon ta
            WHERE   a.align_id = s.align_id AND
                    s.isw_id = i.isw_id AND
                    a.align_id = se.align_id AND
                    se.seq_id = t.seq_id AND
                    se.chr_id = c.chr_id AND
                    c.taxon_id = ta.taxon_id
            ORDER BY c.chr_name, s.snp_pos
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$indel_basic; next; }
    if ( $n == 2 ) { &$snp_basic;   next; }
    if ( $n == 3 ) { &$indel_list;  next; }
    if ( $n == 4 ) { &$snp_list;    next; }
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
