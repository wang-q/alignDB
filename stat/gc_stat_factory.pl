#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Statistics::R;
use Tie::IxHash;
use Number::Format qw(round);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::ToXLSX;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

gc_stat_factory.pl - GC stats for alignDB

=head1 SYNOPSIS

    perl gc_stat_factory.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --outfile   -o  STR     outfile filename
        --run       -r  STR     run special analysis
        --combine       INT     
        --piece         INT
        --alt_level             use 200 .. 900, 1k, 2k, 3k, 4k, 5k segment levels
        --replace       STR=STR replace strings in axis names
        --index                 add an index sheet
        --chart                 add charts

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'output|o=s'   => \( my $outfile ),
    'run|r=s'      => \( my $run      = $Config->{stat}{run} ),
    'combine=i'    => \( my $combine  = 0 ),
    'piece=i'      => \( my $piece    = 0 ),
    'alt_level'    => \my $alt_level,
    'replace=s'    => \my %replace,
    'index'        => \( my $add_index_sheet, ),
    'chart'        => \( my $add_chart, ),
) or HelpMessage(1);

$outfile = "$db.gc.xlsx" unless $outfile;

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

my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password )
    or die "Cannot connect to MySQL database at $db:$server";
my $write_obj = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
    replace => \%replace,
);

my $sql_file = AlignDB::SQL::Library->new( lib => "$FindBin::Bin/sql.lib" );

# auto detect combine threshold
if ( $combine == 0 ) {
    ($combine) = $write_obj->calc_threshold;
}

# auto detect combine threshold
if ( $piece == 0 ) {
    ( undef, $piece ) = $write_obj->calc_threshold;
}

#----------------------------------------------------------#
# chart -- wave
#----------------------------------------------------------#
my $chart_wave_distance = sub {
    my $sheet   = shift;
    my $data    = shift;
    my $x_title = shift;

    my %opt = (
        x_column    => 0,
        y_column    => 1,
        first_row   => 1,
        last_row    => 16,
        x_max_scale => 15,
        y_data      => $data->[1],
        x_title     => $x_title,
        y_title     => "Nucleotide diversity",
        top         => 1,
        left        => 10,
    );
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 3;
    $opt{y_data}   = $data->[3];
    $opt{y_title}  = "Indel per 100 bp";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 5;
    $opt{y_data}   = $data->[5];
    $opt{y_title}  = "Window CV";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );
};

my $chart_series = sub {
    my $sheet   = shift;
    my $data_of = shift;

    # write charting data
    my @keys = keys %{$data_of};
    $write_obj->row(2);
    $write_obj->column(7);

    $write_obj->write_column( $sheet, { column => $data_of->{ $keys[0] }[0], } );
    for my $key (@keys) {
        $write_obj->write_column(
            $sheet,
            {   query_name => $key,
                column     => $data_of->{$key}[1],
            }
        );
    }

    my %opt = (
        x_column      => 7,
        y_column      => 8,
        y_last_column => 8 + @keys - 1,
        first_row     => 2,
        last_row      => 12,
        x_min_scale   => 0,
        x_max_scale   => 10,
        y_data        => [ map { $data_of->{$_}[1] } @keys ],
        x_title       => "Distance to GC trough",
        y_title       => "Indel per 100 bp",
        top           => 14,
        left          => 7,
        height        => 480,
        width         => 480,
    );
    $write_obj->draw_dd( $sheet, \%opt );
};

#----------------------------------------------------------#
# worksheet -- summary
#----------------------------------------------------------#
my $summary = sub {
    my $sheet_name = 'summary';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{TYPE COUNT AVG_length SUM_length indel INDEL/100bp };
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # contents
        my $query_name = 'crest_trough';
        my $sql_query  = q{
            SELECT e.extreme_type TYPE,
                   COUNT(e.window_id) COUNT, 
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`
            FROM extreme e, window w, (SELECT SUM(align_comparables) sum_length
                                       FROM align) a,
                                      (SELECT COUNT(indel_id) sum_indel
                                       FROM indel) i 
            WHERE w.window_id = e.window_id
            GROUP BY e.extreme_type
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'gsw';
        my $sql_query  = q{
            SELECT g.gsw_type TYPE, 
                   COUNT(g.window_id) COUNT, 
                   AVG(w.window_length) AVG_length, 
                   SUM(w.window_length) SUM_length, 
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`
            FROM gsw g, window w, (SELECT SUM(align_comparables) sum_length
                                   FROM align) a,
                                  (SELECT COUNT(indel_id) sum_indel
                                   FROM indel) i 
            WHERE w.window_id = g.window_id
            GROUP BY g.gsw_type
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'align';
        my $sql_query  = q{
            SELECT 'All' TYPE, 
                   COUNT(DISTINCT a.align_id) COUNT,
                   AVG(a.align_comparables) AVG_length, 
                   SUM(a.align_comparables) SUM_length,
                   SUM(i.indel) indel,
                   SUM(i.indel) / SUM(a.align_comparables) * 100 `INDEL/100bp`
            FROM    align a,
                    (SELECT a.align_id,
                            COUNT(i.indel_id) indel
                    FROM align a, indel i
                    WHERE a.align_id = i.align_id
                    GROUP BY a.align_id) i
            WHERE a.align_id = i.align_id
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

my $segment_summary = sub {
    my $sheet_name = 'segment_summary';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{TYPE AVG MIN MAX STD COUNT SUM};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $column_stat = sub {
        my $query_name = shift;
        my $column     = shift;

        my $sql_query = q{
            # summray stat of _COLUMN_
            SELECT  segment_type TYPE,
                    AVG(_COLUMN_) AVG,
                    MIN(_COLUMN_) MIN,
                    MAX(_COLUMN_) MAX,
                    STD(_COLUMN_) STD,
                    COUNT(_COLUMN_) COUNT,
                    SUM(_COLUMN_) SUM
            FROM    segment
            GROUP BY segment_type
        };

        $sql_query =~ s/_COLUMN_/$column/g;

        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    };

    {    # contents
        $column_stat->( 'mean', 'segment_gc_mean' );
        $column_stat->( 'std',  'segment_gc_std' );
        $column_stat->( 'cv',   'segment_gc_cv' );
        $column_stat->( 'mdcw', 'segment_gc_mdcw' );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_to_trough
#----------------------------------------------------------#
my $distance_to_trough = sub {
    my $sheet_name = 'distance_to_trough';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combined = $write_obj->make_combine(
        {   sql_query  => $sql_file->retrieve('gc-wave_combine-0')->as_sql,
            threshold  => $combine,
            standalone => [0],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');

    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    if ($add_chart) {    # chart
        $chart_wave_distance->( $sheet, $data, "Distance to GC trough" );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_to_crest
#----------------------------------------------------------#
my $distance_to_crest = sub {
    my $sheet_name = 'distance_to_crest';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combine_sql = $sql_file->retrieve('gc-wave_combine-0');
    $combine_sql->replace( { gsw_distance => 'gsw_distance_crest' } );
    my $combined = $write_obj->make_combine(
        {   sql_query  => $combine_sql->as_sql,
            threshold  => $combine,
            standalone => [0],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
    $thaw_sql->replace( { gsw_distance => 'gsw_distance_crest' } );
    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );
        $thaw_sql->replace( { gsw_distance => 'gsw_distance_crest' } );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    if ($add_chart) {    # chart
        $chart_wave_distance->( $sheet, $data, "Distance to GC crest" );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- wave_length
#----------------------------------------------------------#
my $wave_length = sub {
    my $sheet_name = 'wave_length';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combine_sql = $sql_file->retrieve('gc-wave_combine-0');
    $combine_sql->replace( { gsw_distance => 'FLOOR(gsw_wave_length / 100)' } );
    my $combined = $write_obj->make_combine(
        {   sql_query  => $combine_sql->as_sql,
            threshold  => $combine,
            standalone => [],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
    $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_wave_length / 100)' } );
    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );
        $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_wave_length / 100)' } );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    if ($add_chart) {    # chart
        $chart_wave_distance->( $sheet, $data, "Distance to GC crest" );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- amplitude
#----------------------------------------------------------#
my $amplitude = sub {
    my $sheet_name = 'amplitude';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combine_sql = $sql_file->retrieve('gc-wave_combine-0');
    $combine_sql->add_where( 'gsw_distance' => \'>= 10' );
    $combine_sql->replace( { gsw_distance => 'FLOOR(gsw_amplitude / 0.01)' } );
    my $combined = $write_obj->make_combine(
        {   sql_query  => $combine_sql->as_sql,
            threshold  => $combine,
            standalone => [],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
    $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_wave_length / 100)' } );
    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );
        $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_wave_length / 100)' } );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gradient
#----------------------------------------------------------#
my $gradient = sub {
    my $sheet_name = 'gradient';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combine_sql = $sql_file->retrieve('gc-wave_combine-0');
    $combine_sql->replace( { gsw_distance => 'gsw_gradient' } );
    my $combined = $write_obj->make_combine(
        {   sql_query  => $combine_sql->as_sql,
            threshold  => $combine,
            standalone => [],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
    $thaw_sql->replace( { gsw_distance => 'gsw_gradient' } );
    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );
        $thaw_sql->replace( { gsw_distance => 'gsw_gradient' } );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- trough_gc
#----------------------------------------------------------#
my $trough_gc = sub {
    my $sheet_name = 'trough_gc';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combine_sql = $sql_file->retrieve('gc-wave_combine-0');
    $combine_sql->replace( { gsw_distance => 'FLOOR(gsw_trough_gc / 0.01)' } );
    my $combined = $write_obj->make_combine(
        {   sql_query  => $combine_sql->as_sql,
            threshold  => $combine,
            standalone => [],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
    $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_trough_gc / 0.01)' } );
    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );
        $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_trough_gc / 0.01)' } );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";

};

#----------------------------------------------------------#
# worksheet -- crest_gc
#----------------------------------------------------------#
my $crest_gc = sub {
    my $sheet_name = 'crest_gc';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combine_sql = $sql_file->retrieve('gc-wave_combine-0');
    $combine_sql->replace( { gsw_distance => 'FLOOR(gsw_crest_gc / 0.01)' } );
    my $combined = $write_obj->make_combine(
        {   sql_query  => $combine_sql->as_sql,
            threshold  => $combine,
            standalone => [],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
    $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_crest_gc / 0.01)' } );
    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );
        $thaw_sql->replace( { gsw_distance => 'FLOOR(gsw_crest_gc / 0.01)' } );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- window_gc
#----------------------------------------------------------#
my $window_gc = sub {
    my $sheet_name = 'window_gc';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    # make combine
    my $combine_sql = q{
        SELECT
          FLOOR(window_average_gc / 0.01),
          COUNT(*)
        FROM gsw
          INNER JOIN window ON
            gsw.window_id = window.window_id
        WHERE (window_average_gc IS NOT NULL)
        GROUP BY
          FLOOR(window_average_gc / 0.01)
    };
    my $combined = $write_obj->make_combine(
        {   sql_query  => $combine_sql,
            threshold  => $combine,
            standalone => [],
            merge_last => 1,
        }
    );

    my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
    $thaw_sql->replace( { gsw_distance => 'FLOOR(window_average_gc / 0.01)' } );
    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('gc-wave_comb_pi_indel_cv-0');
        $thaw_sql->add_where( 'gsw_distance' => $comb );
        $thaw_sql->replace( { gsw_distance => 'FLOOR(window_average_gc / 0.01)' } );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_wave_length_series
#----------------------------------------------------------#
my $d_wave_length_series = sub {
    my $sheet_name = 'd_wave_length_series';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @levels = ( [ 5, 10 ], [ 11, 20 ], [ 21, 30 ], [ 31, 999 ], );

    my $sql_query = q{
        SELECT gsw_distance distance,
               AVG(w.window_indel / w.window_length * 100) AVG_indel,
               COUNT(w.window_indel) COUNT,
               STD(w.window_indel / w.window_length * 100) STD_indel
        FROM gsw g, window w
        WHERE g.window_id = w.window_id
        AND g.gsw_distance <= 10
        AND FLOOR(gsw_wave_length / 100) BETWEEN ? AND ?
        GROUP BY g.gsw_distance
        ORDER BY g.gsw_distance ASC
    };

    my @names = $write_obj->sql2names( $sql_query, { bind_value => $levels[0] } );
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    for (@levels) {
        my $group_name = $_->[0] . "--" . $_->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $sql_query,
                query_name => $group_name,
                bind_value => $_,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_amplitude_series
#----------------------------------------------------------#
my $d_amplitude_series = sub {
    my $sheet_name = 'd_amplitude_series';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @levels = ( [ 10, 20 ], [ 21, 30 ], [ 31, 100 ] );

    my $sql_query = q{
        SELECT  g.gsw_distance gsw_distance,
                AVG(w.window_indel / w.window_length * 100) AVG_indel,
                COUNT(*) COUNT,
                STD(w.window_indel / w.window_length * 100) STD_indel
        FROM    gsw g,
                window w
        WHERE g.window_id = w.window_id
        AND g.gsw_distance <= 10
        AND FLOOR(gsw_amplitude / 0.01) BETWEEN ? AND ?
        GROUP BY g.gsw_distance
        ORDER BY g.gsw_distance ASC
    };

    my @names = $write_obj->sql2names( $sql_query, { bind_value => $levels[0] } );
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    for (@levels) {
        my $group_name = $_->[0] . "--" . $_->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $sql_query,
                query_name => $group_name,
                bind_value => $_,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_gradient_series
#----------------------------------------------------------#
#my $d_gradient_series = sub {
#
#    # find quartiles
#    my @levels;
#    {
#        my $sql_query = q{
#            SELECT g.gsw_gradient
#            FROM gsw g
#            WHERE 1 = 1
#        };
#        my %option = ( sql_query => $sql_query, );
#        my $quartiles = $write_obj->quantile_sql( \%option, 4 );
#        $_ = round( $_, 5 ) for @{$quartiles};
#        @levels = (
#            [ $quartiles->[0], $quartiles->[1] ],    # 1/4
#            [ $quartiles->[1], $quartiles->[2] ],    # 2/4
#            [ $quartiles->[2], $quartiles->[3] ],    # 3/4
#            [ $quartiles->[3], $quartiles->[4] ],    # 4/4
#        );
#    }
#
#    my $sheet_name = 'd_gradient_series';
#    my $sheet;
#    my ( $sheet_row, $sheet_col );
#
#    {                                                # header
#        my $query_name = 'd_gradient_series';
#        my @headers    = qw{gsw_distance AVG_indel COUNT STD_indel};
#        ( $sheet_row, $sheet_col ) = ( 0, 1 );
#        my %option = (
#            sheet_row  => $sheet_row,
#            sheet_col  => $sheet_col,
#            header     => \@headers,
#            query_name => $query_name,
#        );
#        ( $sheet, $sheet_row )
#            = $write_obj->write_header_direct( $sheet_name, \%option );
#    }
#
#    {    # contents
#        my $sql_query = q{
#            SELECT  g.gsw_distance gsw_distance,
#                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
#                    COUNT(*) COUNT,
#                    STD(w.window_indel / w.window_length * 100) STD_indel
#            FROM    gsw g,
#                    window w
#            WHERE g.window_id = w.window_id
#            AND g.gsw_distance <= 10
#            AND g.gsw_gradient BETWEEN ? AND ?
#            GROUP BY g.gsw_distance
#            ORDER BY g.gsw_distance ASC
#        };
#        my %option = (
#            sql_query => $sql_query,
#            sheet_row => $sheet_row,
#            sheet_col => $sheet_col,
#            group     => \@levels,
#        );
#        ($sheet_row) = $write_obj->write_content_series( $sheet, \%option );
#    }
#
#    print "Sheet \"$sheet_name\" has been generated.\n";
#};

#----------------------------------------------------------#
# worksheet -- d_gc_series
#----------------------------------------------------------#
my $d_gc_series = sub {
    my $sheet_name = 'd_gc_series';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # find quartiles
    my @levels;
    {
        my $sql_query = q{
            SELECT w.window_average_gc
            FROM gsw g, window w
            WHERE 1 = 1
            AND g.window_id = w.window_id
            AND g.gsw_distance > 0
            AND g.gsw_distance <= 10
            AND w.window_average_gc IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        my $quartiles = $write_obj->quantile_sql( \%option, 4 );
        $_ = round( $_, 3 ) for @{$quartiles};
        @levels = (
            [ $quartiles->[0], $quartiles->[1] ],    # 1/4
            [ $quartiles->[1], $quartiles->[2] ],    # 2/4
            [ $quartiles->[2], $quartiles->[3] ],    # 3/4
            [ $quartiles->[3], $quartiles->[4] ],    # 4/4
        );
    }

    my $sql_query = q{
        SELECT  g.gsw_distance gsw_distance,
                AVG(w.window_indel / w.window_length * 100) AVG_indel,
                COUNT(*) COUNT,
                STD(w.window_indel / w.window_length * 100) STD_indel
        FROM    gsw g,
                window w
        WHERE g.window_id = w.window_id
        AND g.gsw_distance > 0
        AND g.gsw_distance <= 10
        AND w.window_average_gc BETWEEN ? AND ?
        GROUP BY g.gsw_distance 
        ORDER BY g.gsw_distance ASC
    };

    my @names = $write_obj->sql2names( $sql_query, { bind_value => $levels[0] } );
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    for (@levels) {
        my $group_name = $_->[0] . "--" . $_->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $sql_query,
                query_name => $group_name,
                bind_value => $_,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_trough_gc_series
#----------------------------------------------------------#
my $d_trough_gc_series = sub {
    my $sheet_name = 'd_trough_gc_series';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # find quartiles
    my @levels;
    {
        my $sql_query = q{
            SELECT gsw_trough_gc
            FROM gsw
            WHERE 1 = 1
            AND gsw_trough_gc IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        my $quartiles = $write_obj->quantile_sql( \%option, 4 );
        $_ = round( $_, 4 ) for @{$quartiles};
        @levels = (
            [ $quartiles->[0], $quartiles->[1] ],    # 1/4
            [ $quartiles->[1], $quartiles->[2] ],    # 2/4
            [ $quartiles->[2], $quartiles->[3] ],    # 3/4
            [ $quartiles->[3], $quartiles->[4] ],    # 4/4
        );
    }

    my $sql_query = q{
        SELECT  g.gsw_distance gsw_distance,
                AVG(w.window_indel / w.window_length * 100) AVG_indel,
                COUNT(*) COUNT,
                STD(w.window_indel / w.window_length * 100) STD_indel
        FROM    gsw g,
                window w
        WHERE g.window_id = w.window_id
        AND g.gsw_distance <= 10
        AND g.gsw_trough_gc BETWEEN ? AND ?
        GROUP BY g.gsw_distance 
        ORDER BY g.gsw_distance ASC
    };

    my @names = $write_obj->sql2names( $sql_query, { bind_value => $levels[0] } );
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    for (@levels) {
        my $group_name = $_->[0] . "--" . $_->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $sql_query,
                query_name => $group_name,
                bind_value => $_,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_crest_gc_series
#----------------------------------------------------------#
my $d_crest_gc_series = sub {
    my $sheet_name = 'd_crest_gc_series';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # find quartiles
    my @levels;
    {
        my $sql_query = q{
            SELECT gsw_crest_gc
            FROM gsw
            WHERE 1 = 1
            AND gsw_crest_gc IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        my $quartiles = $write_obj->quantile_sql( \%option, 4 );
        $_ = round( $_, 4 ) for @{$quartiles};
        @levels = (
            [ $quartiles->[0], $quartiles->[1] ],    # 1/4
            [ $quartiles->[1], $quartiles->[2] ],    # 2/4
            [ $quartiles->[2], $quartiles->[3] ],    # 3/4
            [ $quartiles->[3], $quartiles->[4] ],    # 4/4
        );
    }

    my $sql_query = q{
        SELECT  g.gsw_distance gsw_distance,
                AVG(w.window_indel / w.window_length * 100) AVG_indel,
                COUNT(*) COUNT,
                STD(w.window_indel / w.window_length * 100) STD_indel
        FROM    gsw g,
                window w
        WHERE g.window_id = w.window_id
        AND g.gsw_distance <= 10
        AND g.gsw_crest_gc BETWEEN ? AND ?
        GROUP BY g.gsw_distance 
        ORDER BY g.gsw_distance ASC
    };

    my @names = $write_obj->sql2names( $sql_query, { bind_value => $levels[0] } );
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    for (@levels) {
        my $group_name = $_->[0] . "--" . $_->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $sql_query,
                query_name => $group_name,
                bind_value => $_,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- bed_count_trough
#----------------------------------------------------------#
my $bed_count_trough = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'gsw', 'gsw_bed_count' ) ) {
        return;
    }

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT gsw_distance gsw_distance,
                   COUNT(*) COUNT
            FROM gsw g
            WHERE 1 = 1
            GROUP BY gsw_distance
        };
        my $standalone = [];
        my %option     = (
            sql_query  => $sql_query,
            threshold  => $combine,
            standalone => $standalone,
            merge_last => 1,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    my $sheet_name = 'bed_count_trough';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # header
        my @headers
            = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv AVG_bed STD_bed COUNT };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # contents
        my $sql_query = q{
            SELECT  AVG(gsw_distance) distance_to_trough,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    AVG(g.gsw_bed_count) AVG_bed,
                    STD(g.gsw_bed_count) STD_bed,
                    COUNT(w.window_id) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND gsw_distance IN
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            combined  => \@combined,
        );
        ($sheet_row) = $write_obj->write_content_combine( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- bed_count_crest
#----------------------------------------------------------#
my $bed_count_crest = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'gsw', 'gsw_bed_count' ) ) {
        return;
    }

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT gsw_distance_crest gsw_distance_crest,
                   COUNT(*) COUNT
            FROM gsw g
            WHERE 1 = 1
            GROUP BY gsw_distance_crest
        };
        my $standalone = [];
        my %option     = (
            sql_query  => $sql_query,
            threshold  => $combine,
            standalone => $standalone,
            merge_last => 1,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    my $sheet_name = 'bed_count_crest';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # header
        my @headers
            = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv AVG_bed STD_bed COUNT };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # contents
        my $sql_query = q{
            SELECT  AVG(gsw_distance_crest) distance_to_crest,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    AVG(g.gsw_bed_count) AVG_bed,
                    STD(g.gsw_bed_count) STD_bed,
                    COUNT(w.window_id) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND gsw_distance_crest IN
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            combined  => \@combined,
        );
        ($sheet_row) = $write_obj->write_content_combine( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- segment_gc_indel
#----------------------------------------------------------#
my $segment_gc_indel = sub {

    my @segment_levels = ( 1 .. 3 );
    if ($alt_level) {
        @segment_levels = ( 2 .. 10, 20, 30, 40, 50 );
    }

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_gc_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
                           s.segment_gc_cv `cv`,
                           w.window_coding `coding`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY gc DESC, pi, indel
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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        {    # header
            my @headers
                = qw{AVG_gc AVG_pi AVG_Indel/100bp AVG_CV AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.gc) `AVG_gc`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.cv) `AVG_CV`,
                       AVG(t.coding) `AVG_coding`,
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

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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

    my @segment_levels = ( 1 .. 3 );
    if ($alt_level) {
        @segment_levels = ( 2 .. 10, 20, 30, 40, 50 );
    }

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_std_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
                           w.window_coding `coding`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY std DESC, pi, indel
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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        {    # header
            my @headers
                = qw{AVG_std AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.std) `AVG_std`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.coding) `AVG_coding`,
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

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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

    my @segment_levels = ( 1 .. 3 );
    if ($alt_level) {
        @segment_levels = ( 2 .. 10, 20, 30, 40, 50 );
    }

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_cv_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
                           s.segment_gc_cv `cv`,
                           w.window_coding `coding`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY cv DESC, pi, indel
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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        {    # header
            my @headers
                = qw{AVG_CV AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length Range_gc};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.CV) `AVG_CV`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.coding) `AVG_coding`,
                       AVG(t.length) `AVG_length`,
                       COUNT(*) COUNT,
                       SUM(t.length) `SUM_length`,
                       MAX(t.gc) - MIN(t.gc) `Range_gc`
                FROM tmp_group t
                WHERE t_id IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@combined_segment,
            );

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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

    my @segment_levels = ( 1 .. 3 );
    if ($alt_level) {
        @segment_levels = ( 2 .. 10, 20, 30, 40, 50 );
    }

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_mdcw_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
                           w.window_coding `coding`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY mdcw DESC, pi, indel
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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        {    # header
            my @headers
                = qw{AVG_mdcw AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.mdcw) `AVG_mdcw`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.coding) `AVG_coding`,
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

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
# worksheet -- segment_coding_indel
#----------------------------------------------------------#
my $segment_coding_indel = sub {
    unless ( $write_obj->check_column( 'window', 'window_coding' ) ) {
        return;
    }

    my @segment_levels = ( 1 .. 3 );
    if ($alt_level) {
        @segment_levels = ( 2 .. 10, 20, 30, 40, 50 );
    }

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_coding_indel' . "_$segment_type";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
                           w.window_coding `coding`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY coding DESC, pi, indel
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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        {    # header
            my @headers
                = qw{AVG_coding AVG_pi AVG_Indel/100bp AVG_gc AVG_CV AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.coding) `AVG_coding`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.CV) `AVG_CV`,
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

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@segment_levels) {
        &$write_sheet($_);
    }

};

# coding; repeat; non-repeat & non-coding
my $segment_gc_indel_cr = sub {

    my @segment_levels = (3);
    my @feature_levels = ( [ 1, 0 ], [ 0, 1 ], [ 0, 0 ] );

    my $write_sheet = sub {
        my ( $segment_type, $feature_types ) = @_;
        my $sheet_name
            = 'segment_gc_indel_cr'
            . "_$segment_type" . '_'
            . $feature_types->[0]
            . $feature_types->[1];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
                           w.window_coding `coding`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    AND w.window_coding = ?
                    AND w.window_repeats = ?
                    ORDER BY cv DESC, pi, indel
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [ $segment_type, $feature_types->[0], $feature_types->[1] ],

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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        {    # header
            my @headers
                = qw{AVG_gc AVG_pi AVG_Indel/100bp AVG_CV AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.gc) `AVG_gc`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.cv) `AVG_CV`,
                       AVG(t.coding) `AVG_coding`,
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

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for my $i (@segment_levels) {
        for my $j (@feature_levels) {
            &$write_sheet( $i, $j );
        }
    }

};

# coding; repeat; non-repeat & non-coding
my $segment_cv_indel_cr = sub {

    my @segment_levels = (3);
    my @feature_levels = ( [ 1, 0 ], [ 0, 1 ], [ 0, 0 ] );

    my $write_sheet = sub {
        my ( $segment_type, $feature_types ) = @_;
        my $sheet_name
            = 'segment_cv_indel_cr'
            . "_$segment_type" . '_'
            . $feature_types->[0]
            . $feature_types->[1];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
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
                           w.window_coding `coding`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    AND w.window_coding = ?
                    AND w.window_repeats = ?
                    ORDER BY cv DESC, pi, indel
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [ $segment_type, $feature_types->[0], $feature_types->[1] ],

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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        {    # header
            my @headers
                = qw{AVG_CV AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.gc) `AVG_CV`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.cv) `AVG_gc`,
                       AVG(t.coding) `AVG_coding`,
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

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for my $i (@segment_levels) {
        for my $j (@feature_levels) {
            &$write_sheet( $i, $j );
        }
    }

};

for my $n (@tasks) {
    if ( $n == 1 ) { &$summary;              &$segment_summary;    next; }
    if ( $n == 2 ) { &$distance_to_trough;   &$distance_to_crest;  next; }
    if ( $n == 3 ) { &$wave_length;          &$amplitude;          &$gradient; next; }
    if ( $n == 4 ) { &$trough_gc;            &$crest_gc;           &$window_gc; next; }
    if ( $n == 6 ) { &$d_wave_length_series; &$d_amplitude_series; next; }                #
    if ( $n == 7 ) { &$d_gc_series;        next; }
    if ( $n == 8 ) { &$d_trough_gc_series; &$d_crest_gc_series; next; }
    if ( $n == 9 ) { &$bed_count_trough;   &$bed_count_crest; next; }

    if ( $n == 10 ) { &$segment_gc_indel;     next; }
    if ( $n == 11 ) { &$segment_std_indel;    next; }
    if ( $n == 12 ) { &$segment_cv_indel;     next; }
    if ( $n == 13 ) { &$segment_mdcw_indel;   next; }
    if ( $n == 14 ) { &$segment_coding_indel; next; }

    if ( $n == 30 ) { &$segment_gc_indel_cr; next; }
    if ( $n == 31 ) { &$segment_cv_indel_cr; next; }
}

$stopwatch->end_message;
exit;

__END__
