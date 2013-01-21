#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::IntSpan;
use AlignDB::WriteExcel;
use AlignDB::Stopwatch;

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
my $run     = $Config->{stat}{run};
my $outfile = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=s'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'o|output=s'   => \$outfile,
    'r|run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.ofg.xlsx" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 20 );
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
# worksheet -- summary
#----------------------------------------------------------#
#
my $summary_gene = sub {
    my $sheet_name = 'summary';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my $sql_query  = q{
            SELECT 'Type', 'COUNT', 'AVG_length', 'SUM_length'
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

    {    # write contents
        my $query_name = 'ofg count';
        my $sql_query  = q{
            SELECT  CONCAT(o.ofg_tag, "_", o.ofg_type) Type,
                    COUNT(*) COUNT
              FROM ofg o
            GROUP BY o.ofg_tag, o.ofg_type
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    {    # write contents
        my $query_name = 'ofg';
        my $sql_query  = q{
            SELECT CONCAT(o.ofg_tag, "_", o.ofg_type) Type,
                   COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length
              FROM ofg o, window w
             WHERE w.window_id = o.window_id
            GROUP BY o.ofg_tag, o.ofg_type
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    {    # write contents
        my $query_name = 'ofgsw_outside';
        my $sql_query  = q{
            SELECT CONCAT(o.ofg_tag, "_", o.ofg_type) Type,
                   COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length
              FROM ofg o, window w, ofgsw s
             WHERE w.window_id = s.window_id AND
                   o.ofg_id = s.ofg_id
            GROUP BY o.ofg_tag, o.ofg_type
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- ofg_all
#----------------------------------------------------------#
my $ofg_all = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my $sheet_name = "ofg_all";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_pi STD_pi AVG_indel STD_indel AVG_gc
            STD_gc AVG_cv STD_cv AVG_repeats STD_repeats COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # query
        my $sql_query = q{
            SELECT s.ofgsw_distance `distance`,
                   AVG(w.window_pi) `avg_pi`,
                   STD(w.window_pi) `avg_pi`,
                   AVG(w.window_indel / w.window_length * 100) `avg_indel`,
                   STD(w.window_indel / w.window_length * 100) `avg_indel`,
                   AVG(w.window_target_gc) `avg_gc`,
                   STD(w.window_target_gc) `avg_gc`,
                   AVG(s.ofgsw_cv) `avg_cv`,
                   STD(s.ofgsw_cv) `avg_cv`,
                   AVG(w.window_repeats) `avg_repeats`,
                   STD(w.window_repeats) `avg_repeats`,
                   COUNT(*) count
              FROM ofg o,
                   ofgsw s,
                   window w
             WHERE o.ofg_id = s.ofg_id AND
                   s.window_id = w.window_id
            GROUP BY s.ofgsw_distance
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
# worksheet -- ofg_coding
#----------------------------------------------------------#
my $ofg_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my @levels = ( [ "coding", 1 ], [ "noncoding", 0 ], );

    my $write_sheet = sub {
        my ( $order, $value ) = @_;
        my $sheet_name = "ofg_$order";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_pi STD_pi AVG_indel STD_indel AVG_gc
                STD_gc AVG_cv STD_cv AVG_repeats STD_repeats COUNT};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT s.ofgsw_distance `distance`,
                   AVG(w.window_pi) `avg_pi`,
                   STD(w.window_pi) `avg_pi`,
                   AVG(w.window_indel / w.window_length * 100) `avg_indel`,
                   STD(w.window_indel / w.window_length * 100) `avg_indel`,
                   AVG(w.window_target_gc) `avg_gc`,
                   STD(w.window_target_gc) `avg_gc`,
                   AVG(s.ofgsw_cv) `avg_cv`,
                   STD(s.ofgsw_cv) `avg_cv`,
                   AVG(w.window_repeats) `avg_repeats`,
                   STD(w.window_repeats) `avg_repeats`,
                   COUNT(*) count
                FROM ofgsw s,
                     window w,
                     (
                      SELECT s.ofgsw_id
                      FROM ofgsw s,
                            ofg o,
                           window w
                      WHERE o.window_id = w.window_id AND
                            w.window_coding = ? AND
                            s.ofg_id = o.ofg_id
                     ) sw
                WHERE s.window_id = w.window_id AND
                      s.ofgsw_id = sw.ofgsw_id
                GROUP BY s.ofgsw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $value, ],
            );

            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for (@levels) {
        $write_sheet->(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- ofg_coding_pure
#----------------------------------------------------------#
my $ofg_coding_pure = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my @levels = ( [ "coding", 1 ], [ "noncoding", 0 ], );

    my $write_sheet = sub {
        my ( $order, $value ) = @_;
        my $sheet_name = "ofg_$order" . "_pure";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_pi STD_pi AVG_indel STD_indel AVG_gc
                STD_gc AVG_cv STD_cv AVG_repeats STD_repeats COUNT};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT s.ofgsw_distance `distance`,
                   AVG(w.window_pi) `avg_pi`,
                   STD(w.window_pi) `avg_pi`,
                   AVG(w.window_indel / w.window_length * 100) `avg_indel`,
                   STD(w.window_indel / w.window_length * 100) `avg_indel`,
                   AVG(w.window_target_gc) `avg_gc`,
                   STD(w.window_target_gc) `avg_gc`,
                   AVG(s.ofgsw_cv) `avg_cv`,
                   STD(s.ofgsw_cv) `avg_cv`,
                   AVG(w.window_repeats) `avg_repeats`,
                   STD(w.window_repeats) `avg_repeats`,
                   COUNT(*) count
                FROM ofgsw s,
                     window w,
                     (
                      SELECT s.ofgsw_id
                      FROM ofgsw s,
                            ofg o,
                           window w
                      WHERE o.window_id = w.window_id AND
                            w.window_coding = ? AND
                            s.ofg_id = o.ofg_id
                     ) sw
                WHERE s.window_id = w.window_id AND
                      s.ofgsw_id = sw.ofgsw_id AND
                      w.window_coding = ?
                GROUP BY s.ofgsw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $value, $value, ],
            );

            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for (@levels) {
        $write_sheet->(@$_);
    }
};


#----------------------------------------------------------#
# worksheet -- ofg_dG
#----------------------------------------------------------#
my $ofg_dG = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_dG' ) ) {
        return;
    }

    my $sheet_name = "ofg_dG";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_pi STD_pi AVG_indel STD_indel AVG_gc
            STD_gc AVG_cv STD_cv AVG_dG STD_dG COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # query
        my $sql_query = q{
            SELECT s.ofgsw_distance `distance`,
                   AVG(w.window_pi) `avg_pi`,
                   STD(w.window_pi) `std_pi`,
                   AVG(w.window_indel / w.window_length * 100) `avg_indel`,
                   STD(w.window_indel / w.window_length * 100) `std_indel`,
                   AVG(w.window_target_gc) `avg_gc`,
                   STD(w.window_target_gc) `std_gc`,
                   AVG(s.ofgsw_cv) `avg_cv`,
                   STD(s.ofgsw_cv) `std_cv`,
                   AVG(s.ofgsw_dG) `avg_dG`,
                   STD(s.ofgsw_dG) `std_dG`,
                   COUNT(*) count
              FROM ofg o,
                   ofgsw s,
                   window w
            WHERE     1 = 1
            AND o.ofg_id = s.ofg_id
            AND s.window_id = w.window_id
            AND s.ofgsw_distance != 0
            GROUP BY s.ofgsw_distance
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
    if ( $n == 1 ) { &$summary_gene;    next; }
    if ( $n == 2 ) { &$ofg_all;         next; }
    if ( $n == 3 ) { &$ofg_coding;      next; }
    if ( $n == 4 ) { &$ofg_coding_pure; next; }
    if ( $n == 5 ) { &$ofg_dG; next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    ofg_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    ofg_stat_factory.pl [options]
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
