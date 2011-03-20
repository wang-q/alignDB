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
my $run           = $Config->{stat}->{run};
my $sum_threshold = $Config->{stat}->{sum_threshold};
my $outfile       = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'server=s'    => \$server,
    'port=s'      => \$port,
    'db=s'        => \$db,
    'username=s'  => \$username,
    'password=s'  => \$password,
    'output=s'    => \$outfile,
    'run=s'       => \$run,
    'threshold=i' => \$sum_threshold,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.var.xlsx" unless $outfile;

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

    # write contents -- coding
    if ( $write_obj->check_column( 'indel_extra', 'indel_feature1' ) ) {
        my $query_name = 'Coding_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT CASE e.indel_feature1
                     WHEN 0 THEN 'non_coding'
                     WHEN 1 THEN 'coding'
                     WHEN NULL THEN 'NULL'
                     ELSE 'semi_coding'
                   END `coding`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i,
                 indel_extra e
            WHERE i.indel_id = e.indel_id
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

    # write contents -- repeats
    if ( $write_obj->check_column( 'indel_extra', 'indel_feature2' ) ) {
        my $query_name = 'Repeat_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT CASE e.indel_feature2
                     WHEN 0 THEN 'non_repeat'
                     WHEN 1 THEN 'repeat' 
                     WHEN NULL THEN 'NULL'
                     ELSE 'semi_repeat'
                   END `repeat`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i,
                 indel_extra e
            WHERE i.indel_id = e.indel_id
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

    # write contents -- slippage
    if ( $write_obj->check_column( 'indel_extra', 'indel_feature3' ) ) {
        my $query_name = 'Slippage_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT CASE e.indel_feature3
                     WHEN 0 THEN 'non_slippage'
                     WHEN 1 THEN 'slippage' 
                     WHEN NULL THEN 'NULL'
                     ELSE 'semi_slippage'
                   END `slippage`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i,
                 indel_extra e
            WHERE i.indel_id = e.indel_id
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

    # write contents -- occured
    if ( $write_obj->check_column( 'indel', 'indel_occured' ) ) {
        my $query_name = 'Occured_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT CASE i.indel_occured
                     WHEN 'N' THEN 'Not determined'
                     WHEN 'T' THEN 'Target' 
                     WHEN 'Q' THEN 'Query'
                     ELSE 'WRONG'
                   END `occured`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY occured
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents -- insertion or deletion
    if ( $write_obj->check_column( 'indel', 'indel_type' ) ) {
        my $query_name = 'Indel_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT CASE i.indel_type
                     WHEN 'C' THEN 'Complex'
                     WHEN 'I' THEN 'Insertion' 
                     WHEN 'D' THEN 'Deletion'
                     ELSE 'WRONG'
                   END `occured`,
                   AVG (i.indel_length) `AVG_length`,
                   STD(i.indel_length) `STD_length`,
                   COUNT(*) `COUNT`
            FROM indel i
            GROUP BY occured
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents -- length
    {
        my $query_name = 'Length_type';
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

    # write contents -- length - non-slippage
    if ( $write_obj->check_column( 'indel_extra', 'indel_feature3' ) ) {
        $sheet_row++;
        my $query_name = 'Length_type';
        for ( 1 .. 20 ) {
            my $sql_query = qq{
                SELECT '1-$_',
                       AVG (i.indel_length) `AVG_length`,
                       STD(i.indel_length) `STD_length`,
                       COUNT(*) `COUNT`
                FROM indel i, indel_extra e
                WHERE i.indel_id = e.indel_id
                AND e.indel_feature3 = 0
                AND indel_length BETWEEN 1 AND $_
            };
            my %option = (
                query_name => $query_name,
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }
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
        my @headers    = qw{Type_name  COUNT};
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

    # write contents -- coding
    if ( $write_obj->check_column( 'snp', 'snp_coding' ) ) {
        my $query_name = 'Coding_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT  CASE s.snp_coding
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

    # write contents -- repeats
    if ( $write_obj->check_column( 'snp', 'snp_repeats' ) ) {
        my $query_name = 'Repeat_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT  CASE s.snp_repeats
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

    # write contents -- CpG
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

    # write contents -- occured
    if ( $write_obj->check_column( 'snp', 'snp_occured' ) ) {
        my $query_name = 'Occured_type';
        $sheet_row++;
        my $sql_query = q{
            SELECT CASE s.snp_occured
                     WHEN 'N' THEN 'Not determined'
                     WHEN 'T' THEN 'Target' 
                     WHEN 'Q' THEN 'Query'
                     ELSE 'WRONG'
                   END `occured`,
                   COUNT(*) `COUNT`
            FROM snp s
            GROUP BY occured
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
# worksheet -- occured_type
#----------------------------------------------------------#
my $occured_type = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my $sheet_name = 'occured_type';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'Item';
        my $sql_query  = q{
            # header of Table summray
            SELECT  'Indel_occured', 'slippage-like',
                    'AVG_lenght', 'STD_length', 'COUNT'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'Occured_not_change';
        my $sql_query  = q{
            SELECT  CONCAT(i.indel_occured, "->", e.indel_feature4) `change`,
                    e.indel_feature3 `slippage-like`,
                    AVG(i.indel_length) AVG_length,
                    STD(i.indel_length) STD_length,
                    COUNT(*) COUNT
            FROM indel i, indel_extra e
            WHERE i.indel_id = e.indel_id
            AND e.indel_feature4 IS NOT NULL
            AND i.indel_occured = e.indel_feature4
            AND i.indel_occured IN ('T', 'Q')
            GROUP BY CONCAT(i.indel_occured, "->", e.indel_feature4),
                    e.indel_feature3
            ORDER BY CONCAT(i.indel_occured, "->", e.indel_feature4) DESC
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
        $sheet_row++;
        my $query_name = 'Occured_change';
        my $sql_query  = q{
            SELECT  CONCAT(i.indel_occured, "->", e.indel_feature4) `change`, 
                    e.indel_feature3 `slippage-like`,
                    AVG(i.indel_length) AVG_length,
                    STD(i.indel_length) STD_length,
                    COUNT(*) COUNT
            FROM indel i, indel_extra e
            WHERE i.indel_id = e.indel_id
            AND e.indel_feature4 IS NOT NULL
            AND i.indel_occured != e.indel_feature4
            AND i.indel_occured IN ('T', 'Q')
            AND e.indel_feature4 IN ('T', 'Q')
            GROUP BY CONCAT(i.indel_occured, "->", e.indel_feature4),
                    e.indel_feature3
            ORDER BY CONCAT(i.indel_occured, "->", e.indel_feature4) DESC
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

foreach my $n (@tasks) {
    if ( $n == 6 ) { &$indel_basic; next; }
    if ( $n == 8 ) { &$snp_basic;   next; }

    # useless on this time when AlignDB::Multi exists
    #if ( $n == 18 ) { &$occured_type; next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    var_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    three_stat_factory.pl [options]
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
