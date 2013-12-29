#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Number::Format qw(round);

use AlignDB::IntSpan;
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::WriteExcel;

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
my $run     = 'all';
my $combine = 0;
my $piece   = 0;
my $outfile;

# use 200 .. 900, 1k, 2k, 3k, 4k, 5k segment levels
my $alt_level;

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
    'cb|combine=i' => \$combine,
    'pc|piece=i'   => \$piece,
    'alt_level'    => \$alt_level,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
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
            indel INDEL/100bp ns_indel ns_INDEL/100bp
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

    {    # write contents
        my $query_name = 'crest_trough';
        my $sql_query  = q{
            SELECT e.extreme_type TYPE,
                   COUNT(e.window_id) COUNT, 
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100 `ns_INDEL/100bp`
            FROM extreme e, window w, (SELECT SUM(align_comparables) sum_length
                                       FROM align) a,
                                      (SELECT COUNT(indel_id) sum_indel
                                       FROM indel) i 
            WHERE w.window_id = e.window_id
            GROUP BY e.extreme_type
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'gsw';
        my $sql_query  = q{
            SELECT g.gsw_type TYPE, 
                   COUNT(g.window_id) COUNT, 
                   AVG(w.window_length) AVG_length, 
                   SUM(w.window_length) SUM_length, 
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100 `ns_INDEL/100bp`
            FROM gsw g, window w, (SELECT SUM(align_comparables) sum_length
                                   FROM align) a,
                                  (SELECT COUNT(indel_id) sum_indel
                                   FROM indel) i 
            WHERE w.window_id = g.window_id
            GROUP BY g.gsw_type
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write contents
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
# worksheet -- distance_to_trough
#----------------------------------------------------------#
my $distance_to_trough = sub {

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

    my $sheet_name = 'distance_to_trough';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_distance_to_trough AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(gsw_distance) distance_to_trough,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
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
# worksheet -- distance_to_crest
#----------------------------------------------------------#
my $distance_to_crest = sub {

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

    my $sheet_name = 'distance_to_crest';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_distance_to_crest AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(gsw_distance_crest) distance_to_crest,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
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

    {    # write header
        my @headers
            = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv AVG_bed STD_bed COUNT };
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

    {    # write header
        my @headers
            = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv AVG_bed STD_bed COUNT };
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
# worksheet -- wave_length
#----------------------------------------------------------#
my $wave_length = sub {

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT FLOOR(gsw_wave_length / 100) wave_length,
                   COUNT(*) COUNT
            FROM gsw g
            WHERE 1 = 1
            GROUP BY FLOOR(gsw_wave_length / 100)
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

    my $sheet_name = 'wave_length';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_wave_length AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(FLOOR(gsw_wave_length / 100)) wave_length,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    COUNT(w.window_indel) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND FLOOR(gsw_wave_length / 100) IN
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
# worksheet -- amplitude
#----------------------------------------------------------#
my $amplitude = sub {

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT FLOOR(gsw_amplitude / 0.01) amplitude,
                   COUNT(*) COUNT
            FROM gsw g
            WHERE 1 = 1
            AND FLOOR(gsw_amplitude / 0.01) >= 10
            GROUP BY FLOOR(gsw_amplitude / 0.01)
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

    my $sheet_name = 'amplitude';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_amplitude AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(FLOOR(gsw_amplitude / 0.01)) AVG_amplitude,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    COUNT(w.window_indel) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND FLOOR(gsw_amplitude / 0.01) >= 10
            AND FLOOR(gsw_amplitude / 0.01) IN 
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
# worksheet -- trough_gc
#----------------------------------------------------------#
my $trough_gc = sub {

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT FLOOR(gsw_trough_gc / 0.01) trough_gc,
                   COUNT(*) COUNT
            FROM gsw g
            GROUP BY FLOOR(gsw_trough_gc / 0.01)
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

    my $sheet_name = 'trough_gc';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_trough_gc AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(FLOOR(gsw_trough_gc / 0.01)) AVG_trough_gc,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    COUNT(w.window_indel) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND FLOOR(gsw_trough_gc / 0.01) IN 
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
# worksheet -- crest_gc
#----------------------------------------------------------#
my $crest_gc = sub {

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT FLOOR(gsw_crest_gc / 0.01) crest_gc,
                   COUNT(*) COUNT
            FROM gsw g
            GROUP BY FLOOR(gsw_crest_gc / 0.01)
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

    my $sheet_name = 'crest_gc';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_crest_gc AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(FLOOR(gsw_crest_gc / 0.01)) AVG_crest_gc,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    COUNT(w.window_indel) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND FLOOR(gsw_crest_gc / 0.01) IN 
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
# worksheet -- window_gc
#----------------------------------------------------------#
my $window_gc = sub {

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT FLOOR(w.window_average_gc / 0.01) window_gc,
                   COUNT(*) COUNT
            FROM gsw g, window w
            WHERE 1 = 1
            AND g.window_id = w.window_id
            GROUP BY FLOOR(w.window_average_gc / 0.01)
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

    my $sheet_name = 'window_gc';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_window_gc AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(FLOOR(w.window_average_gc / 0.01)) AVG_window_gc,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    COUNT(w.window_indel) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND FLOOR(w.window_average_gc / 0.01) IN 
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
# worksheet -- gradient
#----------------------------------------------------------#
my $gradient = sub {

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT FLOOR(gsw_gradient / 0.00001) gradient,
                   COUNT(*) COUNT
            FROM gsw g
            GROUP BY FLOOR(gsw_gradient / 0.00001)
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

    my $sheet_name = 'gradient';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ AVG_gradient AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            SELECT  AVG(FLOOR(gsw_gradient / 0.00001)) AVG_gradient,
                    AVG(w.window_pi) AVG_pi,
                    STD(w.window_pi) STD_pi,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    STD(w.window_indel / w.window_length * 100) STD_indel,
                    AVG(g.gsw_cv) AVG_cv,
                    STD(g.gsw_cv) STD_cv,
                    COUNT(w.window_indel) COUNT
            FROM gsw g, window w
            WHERE g.window_id = w.window_id
            AND FLOOR(gsw_gradient / 0.00001) IN
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
# worksheet -- d_wave_length_series
#----------------------------------------------------------#
my $d_wave_length_series = sub {

    my $sheet_name = 'd_wave_length_series';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'd_wave_length_series';
        my @headers    = qw{distance AVG_indel COUNT STD_indel};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @levels = ( [ 5, 14 ], [ 15, 24 ], [ 25, 34 ], [ 35, 999 ], );

    {    # write contents
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
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group     => \@levels,
        );
        ($sheet_row) = $write_obj->write_content_series( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_amplitude_series
#----------------------------------------------------------#
my $d_amplitude_series = sub {

    my $sheet_name = 'd_amplitude_series';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'd_amplitude_series';
        my @headers    = qw{gsw_distance AVG_indel COUNT STD_indel};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @levels = ( [ 10, 19 ], [ 20, 29 ], [ 30, 100 ] );

    {    # write contents
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
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group     => \@levels,
        );
        ($sheet_row) = $write_obj->write_content_series( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
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
#    {                                                # write header
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
#    {    # write contents
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
        };
        my %option = ( sql_query => $sql_query, );
        my $quartiles = $write_obj->quantile_sql( \%option, 10 );
        $_ = round( $_, 4 ) for @{$quartiles};
        @levels = (

            [ $quartiles->[0], $quartiles->[1] ],    # 1/10
            [ $quartiles->[1], $quartiles->[2] ],    # 2/10
            [ $quartiles->[2], $quartiles->[3] ],    # 3/10
            [ $quartiles->[3], $quartiles->[4] ],    # 4/10
            [ $quartiles->[5], $quartiles->[6] ],    # 5/10
            [ $quartiles->[6], $quartiles->[7] ],    # 6/10
            [ $quartiles->[7], $quartiles->[8] ],    # 7/10
            [ $quartiles->[8], $quartiles->[9] ],    # 8/10
                 #[ $quartiles->[9],  $quartiles->[10] ],    # 9/10
                 #[ $quartiles->[10], $quartiles->[11] ],    # 10/10
        );
    }

    my $sheet_name = 'd_gc_series';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {            # write header
        my $query_name = 'd_gc_series';
        my @headers    = qw{gsw_distance AVG_indel COUNT STD_indel};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
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
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group     => \@levels,
        );
        ($sheet_row) = $write_obj->write_content_series( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_trough_gc_series
#----------------------------------------------------------#
my $d_trough_gc_series = sub {

    # find quartiles
    my @levels;
    {
        my $sql_query = q{
            SELECT g.gsw_trough_gc
            FROM gsw g
            WHERE 1 = 1
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

    my $sheet_name = 'd_trough_gc_series';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {                                                # write header
        my $query_name = 'd_trough_gc_series';
        my @headers    = qw{gsw_trough_gc AVG_indel COUNT STD_indel};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
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
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group     => \@levels,
        );
        ($sheet_row) = $write_obj->write_content_series( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_crest_gc_series
#----------------------------------------------------------#
my $d_crest_gc_series = sub {

    # find quartiles
    my @levels;
    {
        my $sql_query = q{
            SELECT g.gsw_crest_gc
            FROM gsw g
            WHERE 1 = 1
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

    my $sheet_name = 'd_crest_gc_series';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {                                                # write header
        my $query_name = 'd_crest_gc_series';
        my @headers    = qw{gsw_crest_gc AVG_indel COUNT STD_indel};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
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
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group     => \@levels,
        );
        ($sheet_row) = $write_obj->write_content_series( $sheet, \%option );
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

        {    # write header
            my @headers
                = qw{AVG_gc AVG_pi AVG_Indel/100bp AVG_CV AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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

        {    # write header
            my @headers
                = qw{AVG_std AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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

        {    # write header
            my @headers
                = qw{AVG_CV AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length Range_gc};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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

        {    # write header
            my @headers
                = qw{AVG_mdcw AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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

        {    # write header
            my @headers
                = qw{AVG_coding AVG_pi AVG_Indel/100bp AVG_gc AVG_CV AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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
                sql_query => $sql_query,
                bind_value =>
                    [ $segment_type, $feature_types->[0], $feature_types->[1] ],

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

        {    # write header
            my @headers
                = qw{AVG_gc AVG_pi AVG_Indel/100bp AVG_CV AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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
                sql_query => $sql_query,
                bind_value =>
                    [ $segment_type, $feature_types->[0], $feature_types->[1] ],

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

        {    # write header
            my @headers
                = qw{AVG_CV AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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

my $segment_summary = sub {
    my $sheet_name = 'segment_summary';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{TYPE AVG MIN MAX STD COUNT SUM};
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

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );

        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    };

    {    # write contents
        $column_stat->( 'mean', 'segment_gc_mean' );
        $column_stat->( 'std',  'segment_gc_std' );
        $column_stat->( 'cv',   'segment_gc_cv' );
        $column_stat->( 'mdcw', 'segment_gc_mdcw' );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$summary;            &$segment_summary;   next; }
    if ( $n == 2 ) { &$distance_to_trough; &$distance_to_crest; next; }
    if ( $n == 3 ) { &$wave_length; &$amplitude; &$gradient; next; }
    if ( $n == 4 ) { &$trough_gc;   &$crest_gc;  next; }
    if ( $n == 5 ) { &$window_gc;            next; }
    if ( $n == 6 ) { &$d_wave_length_series; &$d_amplitude_series; next; }
    if ( $n == 7 ) { &$d_gc_series;          next; }
    if ( $n == 8 ) { &$d_trough_gc_series;   &$d_crest_gc_series; next; }
    if ( $n == 9 ) { &$bed_count_trough;     &$bed_count_crest; next; }

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
