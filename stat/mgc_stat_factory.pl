#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::WriteExcel;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# stat parameter
my $run           = $Config->{stat}->{run};
my $sum_threshold = 100_000;
my $outfile;

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

$outfile = "$db.mgc.xlsx" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 20 );
}
else {
    $run =~ s/\"\'//s;
    my $set = AlignDB::IntSpan->new();
    if ( AlignDB::IntSpan->valid($run) ) {
        $set   = $set->add($run);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $run;
        $set->add(@tasks);
    }

    my $runlist = $set->runlist();
    $outfile =~ s/(\.xlsx)$/.$runlist$1/;
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Do stat for $db...");

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

# auto detect threshold
if ( $sum_threshold == 0 ) {
    my $dbh = $write_obj->dbh;

    my $sql_query = q{
        select sum(align_length)
        from align
    };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute;
    my ($total_length) = $sth->fetchrow_array;

    if ( $total_length <= 1_000_000 ) {
        $sum_threshold = int( $total_length / 10 );
    }
    elsif ( $total_length <= 10_000_000 ) {
        $sum_threshold = int( $total_length / 10 );
    }
    elsif ( $total_length <= 100_000_000 ) {
        $sum_threshold = int( $total_length / 20 );
    }
    elsif ( $total_length <= 1_000_000_000 ) {
        $sum_threshold = int( $total_length / 50 );
    }
    else {
        $sum_threshold = int( $total_length / 100 );
    }
}

#----------------------------------------------------------#
# worksheet -- summary
#----------------------------------------------------------#
my $summary = sub {
    my $sheet_name = 'summary';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{
            TYPE COUNT AVG_length SUM_length
            indel INDEL/100bp SNP SNP/100bp
        };
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
    {
        my $query_name = 'align';
        my $sql_query  = q{
            SELECT 'All' TYPE, 
                   COUNT(DISTINCT a.align_id) COUNT,
                   AVG(a.align_comparables) AVG_length, 
                   SUM(a.align_comparables) SUM_length,
                   SUM(a.align_indels) indel,
                   SUM(a.align_indels) / SUM(a.align_comparables) * 100 `INDEL/100bp`,
                   SUM(a.align_differences) SNP,
                   SUM(a.align_differences) / SUM(a.align_comparables) * 100 `SNP/100bp`
            FROM   align a
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
# worksheet -- segment_gc_indel
#----------------------------------------------------------#
my $segment_gc_indel = sub {

    my @segment_levels = ( 'A', 0 .. 3 );

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_gc_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # create temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_average_gc `gc`,
                           s.segment_gc_CV `cv`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY gc DESC
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],

            );
            $write_obj->excute_sql( \%option );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                # segment_sum
                SELECT t_id, length
                FROM tmp_group
            };
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
        }

        # write header
        {
            my $sql_query = q{
                # header of Table group_density
                SELECT 'AVG_gc', 'AVG_pi', 'AVG_Indel/100bp', 'AVG_CV',
                       'AVG_length', 'COUNT', 'SUM_length'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # query
        {
            my $sql_query = q{
                SELECT AVG(t.gc) `AVG_gc`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.cv) `AVG_CV`,
                       AVG(t.length) `AVG_length`,
                       COUNT(*) COUNT,
                       SUM(t.length) `SUM_length`
                FROM tmp_group t
                WHERE t_id IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@combined_segment,
            );

            ($sheet_row)
                = $write_obj->write_content_group( $sheet, \%option );
        }

        # drop temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@segment_levels) {
        &$write_sheet($_);
    }

};

#----------------------------------------------------------#
# worksheet -- segment_std_indel
#----------------------------------------------------------#
my $segment_std_indel = sub {

    my @segment_levels = ( 'A', 0 .. 3 );

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_std_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # create temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_average_gc `gc`,
                           s.segment_gc_std `std`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY std DESC
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],

            );
            $write_obj->excute_sql( \%option );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                # segment_sum
                SELECT t_id, length
                FROM tmp_group
            };
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
        }

        # write header
        {
            my $sql_query = q{
                # header of Table group_density
                SELECT 'AVG_std', 'AVG_pi', 'AVG_Indel/100bp', 'AVG_gc',
                       'AVG_length', 'COUNT', 'SUM_length'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # query
        {
            my $sql_query = q{
                SELECT AVG(t.std) `AVG_std`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.length) `AVG_length`,
                       COUNT(*) COUNT,
                       SUM(t.length) `SUM_length`
                FROM tmp_group t
                WHERE t_id IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@combined_segment,
            );

            ($sheet_row)
                = $write_obj->write_content_group( $sheet, \%option );
        }

        # drop temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@segment_levels) {
        &$write_sheet($_);
    }

};

#----------------------------------------------------------#
# worksheet -- segment_cv_indel
#----------------------------------------------------------#
my $segment_cv_indel = sub {

    my @segment_levels = ( 'A', 0 .. 3 );

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_cv_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # create temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_average_gc `gc`,
                           s.segment_gc_CV `cv`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY cv DESC
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],

            );
            $write_obj->excute_sql( \%option );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                # segment_sum
                SELECT t_id, length
                FROM tmp_group
            };
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
        }

        # write header
        {
            my $sql_query = q{
                # header of Table group_density
                SELECT 'AVG_CV', 'AVG_pi', 'AVG_Indel/100bp', 'AVG_gc',
                       'AVG_length', 'COUNT', 'SUM_length'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # query
        {
            my $sql_query = q{
                SELECT AVG(t.CV) `AVG_CV`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.length) `AVG_length`,
                       COUNT(*) COUNT,
                       SUM(t.length) `SUM_length`
                FROM tmp_group t
                WHERE t_id IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@combined_segment,
            );

            ($sheet_row)
                = $write_obj->write_content_group( $sheet, \%option );
        }

        # drop temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@segment_levels) {
        &$write_sheet($_);
    }

};

#----------------------------------------------------------#
# worksheet -- segment_mdcw_indel
#----------------------------------------------------------#
my $segment_mdcw_indel = sub {

    my @segment_levels = ( 'A', 0 .. 3 );

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_mdcw_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # create temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_average_gc `gc`,
                           s.segment_gc_mdcw `mdcw`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY mdcw DESC
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],

            );
            $write_obj->excute_sql( \%option );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                SELECT t_id, length
                FROM tmp_group
            };
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
        }

        # write header
        {
            my $sql_query = q{
                # header of Table group_density
                SELECT 'AVG_mdcw', 'AVG_pi', 'AVG_Indel/100bp', 'AVG_gc',
                       'AVG_length', 'COUNT', 'SUM_length'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # query
        {
            my $sql_query = q{
                SELECT AVG(t.mdcw) `AVG_mdcw`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.length) `AVG_length`,
                       COUNT(*) COUNT,
                       SUM(t.length) `SUM_length`
                FROM tmp_group t
                WHERE t_id IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@combined_segment,
            );

            ($sheet_row)
                = $write_obj->write_content_group( $sheet, \%option );
        }

        # drop temporary table
        {
            my $sql_query = q{
                DROP TABLE IF EXISTS tmp_group
            };
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@segment_levels) {
        &$write_sheet($_);
    }

};

foreach my $n (@tasks) {
    if ( $n == 1 )  { &$summary;              next; }
    if ( $n == 6 )  { &$segment_gc_indel;     next; }
    if ( $n == 7 )  { &$segment_std_indel;    next; }
    if ( $n == 8 )  { &$segment_cv_indel;     next; }
    if ( $n == 9 )  { &$segment_mdcw_indel;   next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    gc_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    gc_stat_factory.pl [options]
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
