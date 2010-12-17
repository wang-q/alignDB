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
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# stat parameter
my $run     = $Config->{stat}->{run};
my $outfile = "";

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

$outfile = "$db.gce.xls" unless $outfile;

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
    $outfile =~ s/(\.xls)$/.$runlist$1/;
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Do stat for $db...");

my $write_excel_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

#----------------------------------------------------------#
# worksheet -- summary_gene
#----------------------------------------------------------#
#
my $summary_gene = sub {
    my $sheet_name = 'summary';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my $sql_query  = q{
            SELECT 'Type', 'COUNT', 'AVG_length', 'SUM_length', 'AVG_pi',
                   'indel', 'INDEL/100bp', 'ns_indel', 'ns_INDEL/100bp'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_excel_obj->write_header_sql( $sheet_name, \%option );
    }

    {    # write contents
        my $query_name = 'gce count';
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
        ($sheet_row)
            = $write_excel_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    {    # write contents
        my $query_name = 'gce';
        my $sql_query  = q{
            SELECT CONCAT(o.ofg_tag, "_", o.ofg_type) Type,
                   COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_feature3) `ns_indel`,
                   SUM(w.window_feature3) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
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
        ($sheet_row)
            = $write_excel_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    {    # write contents
        my $query_name = 'ofgsw_outside';
        my $sql_query  = q{
            SELECT CONCAT(o.ofg_tag, "_", o.ofg_type) Type,
                   COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_feature3) `ns_indel`,
                   SUM(w.window_feature3) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
              FROM ofg o, window w, ofgsw s
             WHERE w.window_id = o.window_id AND
                   o.ofg_id = s.ofg_id AND
                   ASCII (s.ofgsw_type) IN (ASCII ('L'), ASCII ('R'))
            GROUP BY o.ofg_tag, o.ofg_type
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row)
            = $write_excel_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    {    # write contents
        my $query_name = 'ofgsw_inside';
        my $sql_query  = q{
            SELECT CONCAT(o.ofg_tag, "_", o.ofg_type) Type,
                   COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_feature3) `ns_indel`,
                   SUM(w.window_feature3) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
              FROM ofg o, window w, ofgsw s
             WHERE w.window_id = o.window_id AND
                   o.ofg_id = s.ofg_id AND
                   ASCII (s.ofgsw_type) IN (ASCII ('l'), ASCII ('r'))
            GROUP BY o.ofg_tag, o.ofg_type
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row)
            = $write_excel_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gce_all
#----------------------------------------------------------#
my $gce_all = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_excel_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my $sheet_name = "gce_all";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $sql_query = q{
            SELECT  'distance',
                    'AVG_pi',
                    'AVG_Indel/100bp',
                    'AVG_gc',
                    'COUNT'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_excel_obj->write_header_sql( $sheet_name, \%option );
    }

    {    # query
        my $sql_query = q{
            SELECT s.ofgsw_distance `distance`,
                   AVG (w.window_pi) `avg_pi`,
                   AVG (w.window_indel / w.window_length * 100) `avg_indel/100bp`,
                   AVG (w.window_target_gc) `avg_gc`,
                   COUNT(*) count
              FROM ofg o,
                   ofgsw s,
                   window w
             WHERE o.ofg_id = s.ofg_id AND
                   s.window_id = w.window_id AND
                   o.ofg_tag = 'all'
            GROUP BY s.ofgsw_distance
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );

        ($sheet_row)
            = $write_excel_obj->write_content_direct( $sheet, \%option );
    }

    $sheet->set_zoom(75);

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gce_coldspot
#----------------------------------------------------------#
my $gce_type = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_excel_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my @gce_types = qw{C CC M X total CO NCO};

    my $write_sheet = sub {
        my ($gce_tag) = @_;
        my $sheet_name = "gce" . "_" . $gce_tag;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my $sql_query = q{
                SELECT  'distance',
                        'AVG_pi',
                        'AVG_Indel/100bp',
                        'AVG_gc',
                        'COUNT'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_excel_obj->write_header_sql( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT s.ofgsw_distance `distance`,
                       AVG (w.window_pi) `avg_pi`,
                       AVG (w.window_indel / w.window_length * 100) `avg_indel/100bp`,
                       AVG (w.window_target_gc) `avg_gc`,
                       COUNT(*) count
                FROM ofg g,
                     ofgsw s,
                     window w
                WHERE g.ofg_id = s.ofg_id AND
                      s.window_id = w.window_id AND
                      g.ofg_type = ?
                GROUP BY s.ofgsw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [$gce_tag],
            );

            ($sheet_row)
                = $write_excel_obj->write_content_direct( $sheet, \%option );
        }

        $sheet->set_zoom(75);

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@gce_types) {
        &$write_sheet($_);
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$summary_gene; next; }
    if ( $n == 2 ) { &$gce_all;      next; }
    if ( $n == 3 ) { &$gce_type;     next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    gce_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    gce_stat_factory.pl [options]
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
