#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
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
my $run           = $Config->{stat}{run};
my $sum_threshold = $Config->{stat}{sum_threshold};
my $outfile;

# use 100 .. 900 segment levels
my $alt_level;

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
    'alt_level'   => \$alt_level,
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

my $lib = "$FindBin::Bin/sql.lib";
my $sql_file = AlignDB::SQL::Library->new( lib => $lib );

# auto detect threshold
if ( $sum_threshold == 0 ) {
    my $dbh = $write_obj->dbh;

    my $sql_query = q{
        SELECT SUM(align_length)
        FROM align
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
#
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

my $combined_distance = sub {

    # make combine
    my $sql_query = q{
        SELECT gsw_distance distance,
               COUNT(*) COUNT
        FROM gsw g
        GROUP BY gsw_distance
    };
    my $threshold  = 1000;
    my $standalone = [ 0, 1 ];
    my %option     = (
        sql_query  => $sql_query,
        threshold  => $threshold,
        standalone => $standalone,
    );
    my @combined = @{ $write_obj->make_combine( \%option ) };

    #----------------------------------------------------------#
    # worksheet -- combined_distance
    #----------------------------------------------------------#
    {
        my $sheet_name = 'combined_distance';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            $sql_query = q{
                SELECT  AVG(gsw_distance) AVG_distance,
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
            %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                combined  => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $combined_density = sub {

    # make combine
    my $sql_query = q{
        SELECT gsw_density density,
               COUNT(*) COUNT
        FROM gsw g
        WHERE gsw_density != 0
        GROUP BY gsw_density
    };
    my $threshold  = 1000;
    my $standalone = [ 0, 1 ];
    my %option     = (
        sql_query  => $sql_query,
        threshold  => $threshold,
        standalone => $standalone,
    );
    my @combined = @{ $write_obj->make_combine( \%option ) };

    #----------------------------------------------------------#
    # worksheet -- combined_density
    #----------------------------------------------------------#
    {
        my $sheet_name = 'combined_density';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            $sql_query = q{
                SELECT  AVG(gsw_density) AVG_density,
                        AVG(w.window_pi) AVG_pi,
                        STD(w.window_pi) STD_pi,
                        AVG(w.window_indel / w.window_length * 100) AVG_indel,
                        STD(w.window_indel / w.window_length * 100) STD_indel,
                        AVG(g.gsw_cv) AVG_cv,
                        STD(g.gsw_cv) STD_cv,
                        COUNT(w.window_indel) COUNT
                FROM gsw g, window w
                WHERE g.window_id = w.window_id
                AND gsw_density IN
            };
            %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                combined  => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $combined_amplitude = sub {

    # make combine
    my $sql_query = q{
        SELECT gsw_amplitude amplitude,
               COUNT(*) COUNT
        FROM gsw g
        WHERE gsw_amplitude >= 10
        GROUP BY gsw_amplitude
    };
    my $threshold  = 1000;
    my $standalone = [ 0 .. 9 ];
    my %option     = (
        sql_query  => $sql_query,
        threshold  => $threshold,
        standalone => $standalone,
    );
    my @combined = @{ $write_obj->make_combine( \%option ) };

    #----------------------------------------------------------#
    # worksheet -- combined_amplitude
    #----------------------------------------------------------#
    {
        my $sheet_name = 'combined_amplitude';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            $sql_query = q{
                SELECT  AVG(gsw_amplitude) AVG_amplitude,
                        AVG(w.window_pi) AVG_pi,
                        STD(w.window_pi) STD_pi,
                        AVG(w.window_indel / w.window_length * 100) AVG_indel,
                        STD(w.window_indel / w.window_length * 100) STD_indel,
                        AVG(g.gsw_cv) AVG_cv,
                        STD(g.gsw_cv) STD_cv,
                        COUNT(w.window_indel) COUNT
                FROM gsw g, window w
                WHERE g.window_id = w.window_id
                AND gsw_amplitude IN
            };
            %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                combined  => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $combined_a2d = sub {

    # make combine
    my $sql_query = q{
        SELECT FLOOR(gsw_amplitude / gsw_density) a2d,
               COUNT(*) COUNT
        FROM gsw g
        WHERE gsw_amplitude >= 10
        AND FLOOR(gsw_amplitude / gsw_density) >= 0
        GROUP BY a2d
    };
    my $threshold  = 1000;
    my $standalone = [ 0, 1 ];
    my %option     = (
        sql_query  => $sql_query,
        threshold  => $threshold,
        standalone => $standalone,
    );
    my @combined = @{ $write_obj->make_combine( \%option ) };

    #----------------------------------------------------------#
    # worksheet -- combined_a2d
    #----------------------------------------------------------#
    {
        my $sheet_name = 'combined_a2d';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ AVG_distance AVG_pi STD_pi AVG_indel STD_indel AVG_cv STD_cv COUNT };
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
            $sql_query = q{
                SELECT  AVG(FLOOR(gsw_amplitude / gsw_density)) AVG_a2d,
                        AVG(w.window_pi) AVG_pi,
                        STD(w.window_pi) STD_pi,
                        AVG(w.window_indel / w.window_length * 100) AVG_indel,
                        STD(w.window_indel / w.window_length * 100) STD_indel,
                        AVG(g.gsw_cv) AVG_cv,
                        STD(g.gsw_cv) STD_cv,
                        COUNT(w.window_indel) COUNT
                FROM gsw g, window w
                WHERE g.window_id = w.window_id
                AND FLOOR(gsw_amplitude / gsw_density) IN
            };
            %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                combined  => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $dd_group = sub {

    #----------------------------------------------------------#
    # worksheet -- dd_group
    #----------------------------------------------------------#
    {
        my $sheet_name = 'dd_group';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my $query_name = 'dd_group';
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

        my @dd_density_group
            = ( [ 1, 2 ], [ 3, 6 ], [ 7, 10 ], [ 11, 14 ], [ 15, 999 ], );

        # write contents
        {
            my $sql_query = q{
                SELECT CONCAT(gsw_type, gsw_distance) distance,
                       AVG(w.window_indel / w.window_length * 100) AVG_indel,
                       COUNT(w.window_indel) COUNT,
                       STD(w.window_indel / w.window_length * 100) STD_indel
                FROM gsw g, window w
                WHERE g.gsw_type = ?
                AND g.window_id = w.window_id
                AND g.gsw_density BETWEEN ? AND ?
                AND g.gsw_distance <= ?
                GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
                ORDER BY g.gsw_distance
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@dd_density_group,
            );
            ($sheet_row) = $write_obj->write_content_dd_gc( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- da_group
#----------------------------------------------------------#
#
my $da_group = sub {
    my $sheet_name = 'da_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'da_group';
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

    my @amplitude_group = ( [ 10, 14 ], [ 15, 19 ], [ 20, 100 ], );

    # write contents
    {
        my $sql_query_D = q{
            # amplitude_group for D windows
            SELECT  CONCAT(g.gsw_type, g.gsw_distance) gsw_distance,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    COUNT(*) COUNT,
                    STD(w.window_indel / w.window_length * 100) STD_indel
            FROM    gsw g,
                    window w
            WHERE g.window_id = w.window_id
            AND g.gsw_type = 'L'
            AND g.gsw_density > 5
            AND g.gsw_distance <= 5
            AND g.gsw_amplitude BETWEEN ? AND ?
            GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
            ORDER BY g.gsw_distance DESC
        };
        my $sql_query_A = q{
            # amplitude_group for A windows
            SELECT  CONCAT(g.gsw_type, g.gsw_distance) gsw_distance,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    COUNT(*) COUNT,
                    STD(w.window_indel / w.window_length * 100) STD_indel
            FROM    gsw g,
                    window w
            WHERE g.window_id = w.window_id
            AND g.gsw_type = 'R'
            AND g.gsw_density > 5
            AND g.gsw_distance <= 5
            AND g.gsw_amplitude BETWEEN ? AND ?
            GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
        };
        my %option = (
            sql_query_1 => $sql_query_D,
            sql_query_2 => $sql_query_A,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@amplitude_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- da2d_group
#----------------------------------------------------------#
#
my $da2d_group = sub {
    my $sheet_name = 'da2d_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'da2d_group';
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

    my @amplitude_group = ( [ 0, 0 ], [ 1, 1 ], [ 2, 2 ], [ 3, 100 ], );

    # write contents
    {
        my $sql_query_D = q{
            # amplitude_group for D windows
            SELECT  CONCAT(g.gsw_type, g.gsw_distance) gsw_distance,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    COUNT(*) COUNT,
                    STD(w.window_indel / w.window_length * 100) STD_indel
            FROM    gsw g,
                    window w
            WHERE g.window_id = w.window_id
            AND g.gsw_type = 'L'
            AND g.gsw_density > 5
            AND g.gsw_distance <= 5
            AND FLOOR(g.gsw_amplitude / g.gsw_density) BETWEEN ? AND ?
            GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
            ORDER BY g.gsw_distance DESC
        };
        my $sql_query_A = q{
            # amplitude_group for A windows
            SELECT  CONCAT(g.gsw_type, g.gsw_distance) gsw_distance,
                    AVG(w.window_indel / w.window_length * 100) AVG_indel,
                    COUNT(*) COUNT,
                    STD(w.window_indel / w.window_length * 100) STD_indel
            FROM    gsw g,
                    window w
            WHERE g.window_id = w.window_id
            AND g.gsw_type = 'R'
            AND g.gsw_density > 5
            AND g.gsw_distance <= 5
            AND FLOOR(g.gsw_amplitude / g.gsw_density) BETWEEN ? AND ?
            GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
        };
        my %option = (
            sql_query_1 => $sql_query_D,
            sql_query_2 => $sql_query_A,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@amplitude_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- extreme_amplitude_group
#----------------------------------------------------------#
#my $extreme_amplitude_group = sub {
#    my $sheet_name = 'extreme_amplitude_group';
#    my $sheet;
#    my ($sheet_row, $sheet_col);
#
#    # write header
#    {
#        my $query_name = 'extreme_amplitude_group';
#        my $sql_query  = q{
#            # header of Table group_density
#            SELECT 'gsw_distance', 'AVG_indel', 'COUNT', 'STD_indel'
#        };
#        ($sheet_row, $sheet_col) = (0, 1);
#        my %option = (
#            sql_query  => $sql_query,
#            sheet_row  => $sheet_row,
#            sheet_col  => $sheet_col,
#            query_name => $query_name,
#        );
#        ($sheet, $sheet_row) =
#          $write_obj->write_header_sql($sheet_name, \%option);
#    }
#
#    my @amplitude_group = (
#        [ 0,    0.1499, 0,    0.1499, ],
#        [ 0.15, 1,      0,    0.1499, ],
#        [ 0,    0.1499, 0.15, 1, ],
#        [ 0.15, 1,      0.15, 1, ],
#    );
#
#    # write contents
#    {
#        my $sql_query_D = q{
#            # extreme_amplitude_group for D windows
#            SELECT  CONCAT(g.gsw_type, g.gsw_distance) gsw_distance,
#                    AVG(w.window_indel) AVG_indel,
#                    COUNT(*) COUNT,
#                    STD(w.window_indel) STD_indel
#            FROM    gsw g,
#                    window w,
#                    (SELECT extreme_id
#                     FROM extreme
#                     WHERE extreme_left_amplitude BETWEEN ? AND ?
#                     AND extreme_right_amplitude BETWEEN ? AND ?
#                    ) e
#            WHERE g.window_id = w.window_id
#            AND g.gsw_type = 'L'
#            AND g.gsw_density > 5
#            AND g.gsw_distance <= 5
#            AND g.extreme_id = e.extreme_id
#            GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
#            ORDER BY g.gsw_distance DESC
#        };
#        my $sql_query_A = q{
#            # extreme_amplitude_group for A windows
#            SELECT  CONCAT(g.gsw_type, g.gsw_distance) gsw_distance,
#                    AVG(w.window_indel) AVG_indel,
#                    COUNT(*) COUNT,
#                    STD(w.window_indel) STD_indel
#            FROM    gsw g,
#                    window w,
#                    (SELECT extreme_id
#                     FROM extreme
#                     WHERE extreme_left_amplitude BETWEEN ? AND ?
#                     AND extreme_right_amplitude BETWEEN ? AND ?
#                    ) e
#            WHERE g.window_id = w.window_id
#            AND g.gsw_type = 'R'
#            AND g.gsw_density > 5
#            AND g.gsw_distance <= 5
#            AND g.prev_extreme_id = e.extreme_id
#            GROUP BY CONCAT(g.gsw_type, g.gsw_distance)
#        };
#        my %option = (
#            sql_query_1 => $sql_query_D,
#            sql_query_2 => $sql_query_A,
#            sheet_row   => $sheet_row,
#            sheet_col   => $sheet_col,
#            group       => \@amplitude_group,
#        );
#        ($sheet_row) = $write_obj->write_content_indel($sheet, \%option);
#    }
#
#    print "Sheet \"$sheet_name\" has been generated.\n";
#};

#----------------------------------------------------------#
# worksheet -- segment_gc_indel
#----------------------------------------------------------#
my $segment_gc_indel = sub {

    my @segment_levels = ( 'A', 0 .. 3 );
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
                           s.segment_gc_CV `cv`,
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
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
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

    my @segment_levels = ( 'A', 0 .. 3 );
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
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
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

    my @segment_levels = ( 'A', 0 .. 3 );
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
                           s.segment_gc_CV `cv`,
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
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
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

    my @segment_levels = ( 'A', 0 .. 3 );
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
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
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

    my @segment_levels = ( 'A', 0 .. 3 );
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
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
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

#----------------------------------------------------------#
# worksheet -- segment_extreme_indel
#----------------------------------------------------------#
#my $segment_extreme_indel = sub {
#
#    my @segment_levels =
#      ([ '', 1 ], [ 10000, 2 ], [ 5000, 3 ], [ 1000, 4 ], [ 500, 5 ],);
#
#    my $write_sheet = sub {
#        my ($segment_type) = @_;
#        my $sheet_name = 'segment_extreme_indel' . "_$segment_type";
#        my $sheet;
#        my ($sheet_row, $sheet_col);
#
#        # create temporary table
#        {
#            my $sql_query = q{
#                DROP TABLE IF EXISTS extreme_group
#            };
#            my %option = (sql_query => $sql_query,);
#            $write_obj->excute_sql(\%option);
#        }
#
#        {
#            my $sql_query = q{
#                # create temporary table
#                CREATE TABLE extreme_group (e_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (e_id))
#                    ENGINE=MyISAM
#                    SELECT w.window_pi `pi`,
#                           w.window_indel `indel`,
#                           w.window_average_gc `gc`,
#                           s.segment_feature4 `extreme`,
#                           w.window_length `length`
#                    FROM segment s, window w
#                    WHERE s.window_id = w.window_id
#                    AND s.segment_type = ?
#                    ORDER BY extreme DESC
#            };
#            my %option = (
#                sql_query => $sql_query,
#                bind_value => [ $segment_type ],
#
#            );
#            $write_obj->excute_sql(\%option);
#        }
#
#        # make group
#        my @combined_segment;
#        {
#            my $sql_query = q{
#                # segment_sum
#                SELECT e_id, length
#                FROM extreme_group
#            };
#            my $threshold = $sum_threshold;
#            my $merge_last = 1;
#            my %option     = (
#                sql_query  => $sql_query,
#                threshold  => $threshold,
#                merge_last => $merge_last,
#            );
#            @combined_segment = @{ $write_obj->make_combine(\%option) };
#        }
#
#        # write header
#        {
#            my $sql_query = q{
#                # header of Table group_density
#                SELECT 'AVG_extreme', 'AVG_pi', 'AVG_Indel/100bp', 'AVG_gc',
#                       'AVG_length', 'COUNT', 'SUM_length'
#            };
#            ($sheet_row, $sheet_col) = (0, 1);
#            my %option = (
#                sql_query => $sql_query,
#                sheet_row => $sheet_row,
#                sheet_col => $sheet_col,
#            );
#            ($sheet, $sheet_row) =
#              $write_obj->write_header_sql($sheet_name, \%option);
#        }
#
#        # query
#        {
#            my $sql_query = q{
#                SELECT AVG(e.extreme) `AVG_extreme`,
#                       AVG(e.pi) `AVG_pi`,
#                       AVG(e.indel / e.length * 100) `AVG_Indel/100bp`,
#                       AVG(e.gc) `AVG_gc`,
#                       AVG(e.length) `AVG_length`,
#                       COUNT(*) COUNT,
#                       SUM(e.length) `SUM_length`
#                FROM extreme_group e
#                WHERE e_id IN
#            };
#            my %option = (
#                sql_query => $sql_query,
#                sheet_row => $sheet_row,
#                sheet_col => $sheet_col,
#                group     => \@combined_segment,
#            );
#
#            ($sheet_row) = $write_obj->write_content_group($sheet, \%option);
#        }
#
#        # drop temporary table
#        {
#            my $sql_query = q{
#                DROP TABLE IF EXISTS extreme_group
#            };
#            my %option = (sql_query => $sql_query,);
#            $write_obj->excute_sql(\%option);
#        }
#
#        $sheet->set_zoom(75);
#
#        print "Sheet \"$sheet_name\" has been generated.\n";
#    };
#
#    foreach (@segment_levels) {
#        my $segment_type = $_->[1];
#        &$write_sheet($segment_type);
#    }
#
#};

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
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
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
            my $threshold  = $sum_threshold;
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
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
    if ( $n == 1 ) { &$summary;            next; }
    if ( $n == 2 ) { &$combined_distance;  next; }
    if ( $n == 3 ) { &$combined_density;   next; }
    if ( $n == 4 ) { &$combined_amplitude; next; }
    if ( $n == 5 ) { &$combined_a2d;       next; }
    if ( $n == 6 ) { &$dd_group;           next; }
    if ( $n == 7 ) { &$da_group;           &$da2d_group; next; }

    if ( $n == 8 )  { &$segment_gc_indel;     next; }
    if ( $n == 9 )  { &$segment_std_indel;    next; }
    if ( $n == 10 ) { &$segment_cv_indel;     next; }
    if ( $n == 11 ) { &$segment_mdcw_indel;   next; }
    if ( $n == 12 ) { &$segment_coding_indel; next; }
    if ( $n == 13 ) { &$segment_summary;      next; }

    #if ($n == 8) { &$extreme_amplitude_group; next; }
    #if ($n == 24) { &$segment_extreme_indel;  next; }

    #if ( $n == 26 ) { &$segment_gc_indel_cr; next; }
    #if ( $n == 27 ) { &$segment_cv_indel_cr; next; }
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
