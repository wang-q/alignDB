#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::WriteExcel;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
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

$outfile = "$db.three.xlsx" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
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

my $lib = "$FindBin::Bin/sql.lib";
my $sql_file = AlignDB::SQL::Library->new( lib => $lib );

#----------------------------------------------------------#
# worksheet -- summary
#----------------------------------------------------------#
my $summary = sub {
    my $sheet_name = 'summary';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{AVG MIN MAX STD COUNT SUM};
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
        my $table      = shift;
        my $column     = shift;
        my $where      = shift;

        my $sql_query = q{
            # summray stat of _COLUMN_
            SELECT AVG(_COLUMN_) AVG,
                   MIN(_COLUMN_) MIN,
                   MAX(_COLUMN_) MAX,
                   STD(_COLUMN_) STD,
                   COUNT(_COLUMN_) COUNT,
                   SUM(_COLUMN_) SUM
            FROM _TABLE_
        };

        $sql_query =~ s/_TABLE_/$table/g;
        $sql_query =~ s/_COLUMN_/$column/g;
        $sql_query .= $where if $where;

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );

        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    };

    {    # write contents
        &$column_stat( 'align_length',       'align', 'align_length' );
        &$column_stat( 'indel_length',       'indel', 'indel_length' );
        &$column_stat( 'indel_left_extand',  'indel', 'left_extand' );
        &$column_stat( 'indel_right_extand', 'indel', 'right_extand' );
        &$column_stat( 'indel_windows', 'isw', 'isw_length',
            'WHERE isw_distance <= 0' );
        &$column_stat(
            'indel_free_windows', 'isw',
            'isw_length',         'WHERE isw_distance > 0'
        );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance
#----------------------------------------------------------#
my $distance = sub {
    my $sheet_name = 'distance';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{distance AVG_pi AVG_d_indel AVG_d_noindel AVG_d_complex COUNT Di/Dn};
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
        my $thaw_sql = $sql_file->retrieve('three-distance-0');
        my %option   = (
            sql_query => $thaw_sql->as_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write footer
        $sheet_row += 2;
        my $thaw_sql = $sql_file->retrieve('three-distance_total-0');
        my %option   = (
            sql_query      => $thaw_sql->as_sql,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_window0eq0
#----------------------------------------------------------#
my $distance_window0eq0 = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw', 'isw_indel_id' ) ) {
        return;
    }

    my $sheet_name = 'distance_window0eq0';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{distance AVG_pi AVG_d_indel AVG_d_noindel AVG_d_complex COUNT Di/Dn};
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
            SELECT w.isw_distance distance,
                   AVG(w.isw_pi) AVG_pi,
                   AVG(w.isw_d_indel) AVG_d_indel,
                   AVG(w.isw_d_noindel) AVG_d_noindel,
                   AVG(w.isw_d_complex) AVG_d_complex,
                   COUNT(*) COUNT,
                   AVG(w.isw_d_indel) / AVG(w.isw_d_noindel)  `Di/Dn`
            FROM isw w,
                (SELECT i.indel_id, COUNT(*) COUNT
                FROM indel i, isw w
                WHERE i.indel_id = w.isw_indel_id
                AND w.isw_distance = 0
                AND w.isw_pi = 0
                GROUP BY i.indel_id
                HAVING count = 2) i
            WHERE w.isw_indel_id = i.indel_id
            AND w.isw_distance >= 0
            AND isw_d_indel IS NOT NULL
            GROUP BY isw_distance
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write footer
        my $sql_query = q{
            SELECT 'Total',
                   AVG(w.isw_pi) AVG_pi,
                   AVG(w.isw_d_indel) AVG_d_indel,
                   AVG(w.isw_d_noindel) AVG_d_noindel,
                   AVG(w.isw_d_complex) AVG_d_complex,
                   COUNT(*) COUNT,
                   AVG(w.isw_d_indel) / AVG(w.isw_d_noindel)  `Di/Dn`
            FROM isw w,
                (SELECT i.indel_id, COUNT(*) COUNT
                FROM indel i, isw w
                WHERE i.indel_id = w.isw_indel_id
                AND w.isw_distance = 0
                AND w.isw_pi = 0
                GROUP BY i.indel_id
                HAVING count = 2) i
            WHERE w.isw_indel_id = i.indel_id
            AND w.isw_distance >= 0
            AND isw_d_indel IS NOT NULL
        };
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_iswpi_le_2
#----------------------------------------------------------#
my $distance_iswpi_le_2 = sub {
    my $sheet_name = 'distance_iswpi_le_2';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{distance AVG_pi AVG_d_indel AVG_d_noindel AVG_d_complex COUNT Di/Dn};
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
        my $thaw_sql = $sql_file->retrieve('three-distance-0');
        $thaw_sql->add_where( 'isw_pi' => \'<= 0.02' );
        my %option = (
            sql_query => $thaw_sql->as_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write footer
        $sheet_row += 2;
        my $thaw_sql = $sql_file->retrieve('three-distance_total-0');
        my %option   = (
            sql_query      => $thaw_sql->as_sql,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_indel_bottom90
#----------------------------------------------------------#
my $distance_indel_bottom90 = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw', 'isw_indel_id' ) ) {
        return;
    }

    # create temporary table
    {
        my $sql_query = q{
            DROP TABLE IF EXISTS indel_90
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    {

        # init objects
        my $dbh = $write_obj->dbh;

        # count indel
        my $indel_sth = $dbh->prepare(
            q{
            SELECT count(distinct i.indel_id)
            FROM indel i, isw w
            WHERE i.indel_id = w.isw_indel_id
            AND w.isw_distance = 0
            }
        );
        $indel_sth->execute;
        my ($count) = $indel_sth->fetchrow_array;
        my $count_offset = int( $count * 0.9 );

        my $sql_query = qq{
            # create temporary table
            CREATE TABLE indel_90 (indel_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (indel_id))
                ENGINE=MyISAM
                SELECT i.indel_id, avg(w.isw_pi) avg_pi
                FROM indel i, isw w
                WHERE i.indel_id = w.isw_indel_id
                AND w.isw_distance = 0
                GROUP BY i.indel_id
                order by avg_pi ASC
                limit $count_offset
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    #----------------------------#
    # distance_indel_bottom90
    #----------------------------#
    {
        my $sheet_name = 'distance_indel_bottom90';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{distance AVG_pi AVG_d_indel AVG_d_noindel AVG_d_complex COUNT Di/Dn};
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
                SELECT w.isw_distance distance,
                       AVG(w.isw_pi) AVG_pi,
                       AVG(w.isw_d_indel) AVG_d_indel,
                       AVG(w.isw_d_noindel) AVG_d_noindel,
                       AVG(w.isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG(w.isw_d_indel) / AVG(w.isw_d_noindel)  `Di/Dn`
                FROM isw w, indel_90 i
                WHERE w.isw_indel_id = i.indel_id
                AND w.isw_distance >= 0
                AND isw_d_indel IS NOT NULL
                GROUP BY isw_distance
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        {    # write footer
            my $sql_query = qq{
                SELECT 'Total',
                       AVG(w.isw_pi) AVG_pi,
                       AVG(w.isw_d_indel) AVG_d_indel,
                       AVG(w.isw_d_noindel) AVG_d_noindel,
                       AVG(w.isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG(w.isw_d_indel) / AVG(w.isw_d_noindel)  `Di/Dn`
                FROM isw w, indel_90 i
                WHERE w.isw_indel_id = i.indel_id
                AND w.isw_distance >= 0
                AND isw_d_indel IS NOT NULL
            };
            my %option = (
                sql_query      => $sql_query,
                sheet_row      => $sheet_row,
                sheet_col      => $sheet_col,
                content_format => 'TOTAL',
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    # drop temporary table
    {
        my $sql_query = q{
            DROP TABLE IF EXISTS indel_90
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }
};

my $combined_distance = sub {

    # make combine
    my @combined;
    {
        my $thaw_sql = $sql_file->retrieve('common-distance_combine-0');
        $thaw_sql->add_where( 'isw_distance' => \'>= 0' );
        $thaw_sql->add_where( 'isw_d_indel'  => \'IS NOT NULL' );
        my $threshold  = 1000;
        my $standalone = [ -1, 0 ];
        my %option     = (
            sql_query  => $thaw_sql->as_sql,
            threshold  => $threshold,
            standalone => $standalone,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    #----------------------------------------------------------#
    # worksheet -- combined_distance
    #----------------------------------------------------------#
    {
        my $sheet_name = 'combined_distance';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{distance AVG_pi AVG_d_indel AVG_d_noindel AVG_d_complex COUNT Di/Dn};
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
            my $thaw_sql = $sql_file->retrieve('three-distance_avg-0');
            my %option   = (
                sql_obj     => $thaw_sql,
                sheet_row   => $sheet_row,
                sheet_col   => $sheet_col,
                combined    => \@combined,
                combine_col => 'isw_distance',
            );
            ($sheet_row)
                = $write_obj->write_content_combine_obj( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $distance_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature1' ) ) {
        return;
    }

    # make combine
    my ( $isw_f1_1, $isw_f1_2 ) = ( 1, 1 );
    my @combined;
    {
        my $thaw_sql = $sql_file->retrieve('common-distance_coding_combine-2');
        $thaw_sql->add_where( 'isw_distance' => \'>= 0' );
        $thaw_sql->add_where( 'isw_d_indel'  => \'IS NOT NULL' );
        my $threshold  = 1000;
        my $standalone = [ -1, 0 ];
        my %option     = (
            sql_query  => $thaw_sql->as_sql,
            threshold  => $threshold,
            standalone => $standalone,
            bind_value => [ $isw_f1_1, $isw_f1_2 ],
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    #----------------------------------------------------------#
    # worksheet -- distance_coding
    #----------------------------------------------------------#
    {
        my $sheet_name = 'distance_coding';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{distance AVG_pi AVG_d_indel AVG_d_noindel AVG_d_complex COUNT Di/Dn};
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
            my $thaw_sql = $sql_file->retrieve('three-distance_coding-2');
            $thaw_sql->bind( [ $isw_f1_1, $isw_f1_2 ] );
            my %option = (
                sql_obj     => $thaw_sql,
                sheet_row   => $sheet_row,
                sheet_col   => $sheet_col,
                combined    => \@combined,
                combine_col => 'isw_distance',
            );
            ($sheet_row)
                = $write_obj->write_content_combine_obj( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $distance_non_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature1' ) ) {
        return;
    }

    # make combine
    my ( $isw_f1_1, $isw_f1_2 ) = ( 0, 0 );
    my @combined;
    {
        my $thaw_sql = $sql_file->retrieve('common-distance_coding_combine-2');
        $thaw_sql->add_where( 'isw_distance' => \'>= 0' );
        $thaw_sql->add_where( 'isw_d_indel'  => \'IS NOT NULL' );
        my $threshold  = 1000;
        my $standalone = [ -1, 0 ];
        my %option     = (
            sql_query  => $thaw_sql->as_sql,
            threshold  => $threshold,
            standalone => $standalone,
            bind_value => [ $isw_f1_1, $isw_f1_2 ],
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    #----------------------------------------------------------#
    # worksheet -- distance_non_coding
    #----------------------------------------------------------#
    {
        my $sheet_name = 'distance_non_coding';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{distance AVG_pi AVG_d_indel AVG_d_noindel AVG_d_complex COUNT Di/Dn};
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
            my $thaw_sql = $sql_file->retrieve('three-distance_coding-2');
            $thaw_sql->bind( [ $isw_f1_1, $isw_f1_2 ] );
            my %option = (
                sql_obj     => $thaw_sql,
                sheet_row   => $sheet_row,
                sheet_col   => $sheet_col,
                combined    => \@combined,
                combine_col => 'isw_distance',
            );
            ($sheet_row)
                = $write_obj->write_content_combine_obj( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- distance_dir_dnr
#----------------------------------------------------------#
#
my $distance_dir_dnr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature4' ) ) {
        return;
    }

    my $sheet_name = 'distance_dir_dnr';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_pi AVG_d_ir AVG_d_nr AVG_d_total COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q~
            # distance effect
            SELECT isw_distance distance,
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature4) AVG_d_ir,
                   AVG(isw_feature5) AVG_d_nr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw i, isw_extra e
            WHERE isw_distance >= 0
            AND isw_feature4 IS NOT NULL
            and i.isw_id = e.isw_id
            GROUP BY isw_distance
        ~;
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write footer
    {
        $sheet_row += 2;
        my $sql_query = q~
            # isw average
            SELECT 'Total',
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature4) AVG_d_ir,
                   AVG(isw_feature5) AVG_d_nr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw i, isw_extra e
            WHERE isw_distance >= 0
            AND isw_feature4 IS NOT NULL
            and i.isw_id = e.isw_id
        ~;
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_dtr_dqr
#----------------------------------------------------------#
#
my $distance_dtr_dqr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature6' ) ) {
        return;
    }

    my $sheet_name = 'distance_dtr_dqr';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_pi AVG_d_tr AVG_d_qr AVG_d_total COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q~
            # distance effect
            SELECT isw_distance distance,
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature6) AVG_d_tr,
                   AVG(isw_feature7) AVG_d_qr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw i, isw_extra e
            WHERE isw_distance >= 0
            AND isw_feature6 IS NOT NULL
            and i.isw_id = e.isw_id
            GROUP BY isw_distance
        ~;
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write footer
    {
        $sheet_row += 2;
        my $sql_query = q~
            # isw average
            SELECT 'Total',
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature6) AVG_d_tr,
                   AVG(isw_feature7) AVG_d_qr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw i, isw_extra e
            WHERE isw_distance >= 0
            AND isw_feature6 IS NOT NULL
            and i.isw_id = e.isw_id
        ~;
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- target_dtr_dqr
#----------------------------------------------------------#
#
my $target_dtr_dqr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature6' ) ) {
        return;
    }

    my $sheet_name = 'target_dtr_dqr';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_pi AVG_d_tr AVG_d_qr AVG_d_total COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q~
            # distance effect
            SELECT isw_distance distance,
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature6) AVG_d_tr,
                   AVG(isw_feature7) AVG_d_qr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw_extra e,
                 (SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "T"
                 AND isw_type = 'R'
                 AND isw.indel_id = i.indel_id
                 UNION 
                 SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "T"
                 AND isw_type = 'L'
                 AND isw.prev_indel_id = i.indel_id) i
            WHERE i.isw_distance >= 0
            AND e.isw_feature6 IS NOT NULL
            AND i.isw_id = e.isw_id
            GROUP BY isw_distance
        ~;
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write footer
    {
        $sheet_row += 2;
        my $sql_query = q~
            # isw average
            SELECT 'Total',
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature6) AVG_d_tr,
                   AVG(isw_feature7) AVG_d_qr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw_extra e,
                 (SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "T"
                 AND isw_type = 'R'
                 AND isw.indel_id = i.indel_id
                 UNION 
                 SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "T"
                 AND isw_type = 'L'
                 AND isw.prev_indel_id = i.indel_id) i
            WHERE i.isw_distance >= 0
            AND e.isw_feature6 IS NOT NULL
            AND i.isw_id = e.isw_id
        ~;
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- query_dtr_dqr
#----------------------------------------------------------#
#
my $query_dtr_dqr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature6' ) ) {
        return;
    }

    my $sheet_name = 'query_dtr_dqr';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_pi AVG_d_tr AVG_d_qr AVG_d_total COUNT};
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q~
            # distance effect
            SELECT isw_distance distance,
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature6) AVG_d_tr,
                   AVG(isw_feature7) AVG_d_qr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw_extra e,
                 (SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "Q"
                 AND isw_type = 'R'
                 AND isw.indel_id = i.indel_id
                 UNION 
                 SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "Q"
                 AND isw_type = 'L'
                 AND isw.prev_indel_id = i.indel_id) i
            WHERE i.isw_distance >= 0
            AND e.isw_feature6 IS NOT NULL
            AND i.isw_id = e.isw_id
            GROUP BY isw_distance
        ~;
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write footer
    {
        $sheet_row += 2;
        my $sql_query = q~
            # isw average
            SELECT 'Total',
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_feature6) AVG_d_tr,
                   AVG(isw_feature7) AVG_d_qr,
                   AVG(isw_feature8) AVG_d_total,
                   COUNT(*) COUNT
            FROM isw_extra e,
                 (SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "Q"
                 AND isw_type = 'R'
                 AND isw.indel_id = i.indel_id
                 UNION 
                 SELECT isw_id, isw_distance, isw_pi
                 FROM   isw,
                        indel i
                 WHERE  i.indel_occured = "Q"
                 AND isw_type = 'L'
                 AND isw.prev_indel_id = i.indel_id) i
            WHERE i.isw_distance >= 0
            AND e.isw_feature6 IS NOT NULL
            AND i.isw_id = e.isw_id
        ~;
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_distance
#----------------------------------------------------------#
#
my $snp_distance = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ssw', 'ssw_id' ) ) {
        return;
    }

    my @snp_levels = (
        [ 'all',  0,  999 ],
        [ '-1',   -1, -1 ],
        [ '0-5',  0,  5 ],
        [ '6-10', 6,  10 ],

        #[ '0-0',    0,  0 ],
        #[ '1-1',    1,  1 ],
        #[ '2-2',    2,  2 ],
        #[ '3-3',    3,  3 ],
        #[ '4-4',    4,  4 ],
        #[ '5-5',    5,  5 ],
        #[ '6-6',    6,  6 ],
        #[ '7-7',    7,  7 ],
        #[ '8-8',    8,  8 ],
        #[ '9-9',    9,  9 ],
        #[ '10-10',  10, 10 ],
        [ '11-20',  11, 20 ],
        [ '21-30',  21, 30 ],
        [ '31-40',  31, 40 ],
        [ '41-50',  41, 50 ],
        [ '51-999', 51, 999 ],
        [ '21-999', 21, 999 ],
        [ '31-999', 31, 999 ],
        [ '41-999', 41, 999 ],
    );

    my $write_sheet = sub {
        my ($snp_levels) = @_;
        my $sheet_name = 'snp_distance' . "_" . $snp_levels->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q~
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_snp', 'AVG_d_nosnp', 'AVG_d_complex',
                        'COUNT', 'Ds/Dns'
            ~;
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q~
                # distance effect
                SELECT  ssw_distance `distance`,
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM ssw s, window w, snp, isw i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                AND i.isw_distance BETWEEN ? AND ?
                GROUP BY ssw_distance
            ~;
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $snp_levels->[1], $snp_levels->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        # write footer
        {
            $sheet_row += 2;
            my $sql_query = q~
                # isw average
                SELECT  'Total',
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM ssw s, window w, snp, isw i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                AND i.isw_distance BETWEEN ? AND ?
            ~;
            my %option = (
                sql_query      => $sql_query,
                sheet_row      => $sheet_row,
                sheet_col      => $sheet_col,
                bind_value     => [ $snp_levels->[1], $snp_levels->[2] ],
                content_format => 'TOTAL',
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@snp_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- snp_distance_non_slip
#----------------------------------------------------------#
#
my $snp_distance_non_slip = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ssw', 'ssw_id' ) ) {
        return;
    }

    my @snp_levels = (
        [ 'all',  0,  999 ],
        [ '-1',   -1, -1 ],
        [ '0-5',  0,  5 ],
        [ '6-10', 6,  10 ],

        #[ '0-0',    0,  0 ],
        #[ '1-1',    1,  1 ],
        #[ '2-2',    2,  2 ],
        #[ '3-3',    3,  3 ],
        #[ '4-4',    4,  4 ],
        #[ '5-5',    5,  5 ],
        #[ '6-6',    6,  6 ],
        #[ '7-7',    7,  7 ],
        #[ '8-8',    8,  8 ],
        #[ '9-9',    9,  9 ],
        #[ '10-10',  10, 10 ],
        [ '11-20',  11, 20 ],
        [ '21-30',  21, 30 ],
        [ '31-40',  31, 40 ],
        [ '41-50',  41, 50 ],
        [ '51-999', 51, 999 ],
        [ '21-999', 21, 999 ],
        [ '31-999', 31, 999 ],
        [ '41-999', 41, 999 ],
    );

    my $write_sheet = sub {
        my ($snp_levels) = @_;
        my $sheet_name = 'snp_distance' . "_" . $snp_levels->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q~
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_snp', 'AVG_d_nosnp', 'AVG_d_complex',
                        'COUNT', 'Ds/Dns'
            ~;
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q~
                # distance effect
                SELECT  ssw_distance `distance`,
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM    ssw s, window w, snp,
                       (SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND isw_type = 'R'
                        AND isw_distance BETWEEN ? AND ?
                        AND isw.indel_id = i.indel_id
                        UNION 
                        SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND isw_type = 'L'
                        AND isw_distance BETWEEN ? AND ?
                        AND isw.prev_indel_id = i.indel_id) i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                GROUP BY ssw_distance
            ~;
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [
                    $snp_levels->[1], $snp_levels->[2],
                    $snp_levels->[1], $snp_levels->[2]
                ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        # write footer
        {
            $sheet_row += 2;
            my $sql_query = q~
                # isw average
                SELECT  'Total',
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM    ssw s, window w, snp, 
                       (SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND isw_type = 'R'
                        AND isw_distance BETWEEN ? AND ?
                        AND isw.indel_id = i.indel_id
                        UNION 
                        SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND isw_type = 'L'
                        AND isw_distance BETWEEN ? AND ?
                        AND isw.prev_indel_id = i.indel_id) i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
            ~;
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [
                    $snp_levels->[1], $snp_levels->[2],
                    $snp_levels->[1], $snp_levels->[2]
                ],
                content_format => 'TOTAL',
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@snp_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- snp_LR_distance_
#----------------------------------------------------------#
#
my $snp_LR_distance = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ssw', 'ssw_id' ) ) {
        return;
    }

    my @snp_levels = (
        [ 'R', 0, 0 ],
        [ 'L', 0, 0 ],
        [ 'R', 1, 1 ],
        [ 'L', 1, 1 ],
        [ 'R', 2, 2 ],
        [ 'L', 2, 2 ],
        [ 'R', 3, 3 ],
        [ 'L', 3, 3 ],
        [ 'R', 4, 4 ],
        [ 'L', 4, 4 ],
        [ 'R', 5, 5 ],
        [ 'L', 5, 5 ],
    );

    my $write_sheet = sub {
        my ($snp_levels) = @_;
        my $sheet_name = 'snp_distance' . "_" . join "-", @$snp_levels;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q~
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_snp', 'AVG_d_nosnp', 'AVG_d_complex',
                        'COUNT', 'Ds/Dns'
            ~;
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query_R = q~
                # distance effect
                SELECT  CONCAT(s.ssw_type, s.ssw_distance) `distance`,
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM ssw s, window w, snp, isw i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                AND i.isw_type = ?
                AND i.isw_distance BETWEEN ? AND ?
                AND s.ssw_type = 'L'
                GROUP BY CONCAT(s.ssw_type, s.ssw_distance)
                ORDER BY s.ssw_distance DESC
            ~;
            my $sql_query_L = q~
                # distance effect
                SELECT  CONCAT(s.ssw_type, s.ssw_distance) `distance`,
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM ssw s, window w, snp, isw i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                AND i.isw_type = ?
                AND i.isw_distance BETWEEN ? AND ?
                AND s.ssw_type = 'R'
                GROUP BY CONCAT(s.ssw_type, s.ssw_distance)
                ORDER BY s.ssw_distance 
            ~;

            $sheet_row++;    # add a blank line
            my %option = (
                sql_query  => $sql_query_R,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [@$snp_levels],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );

            $sheet_row++;    # add a blank line
            %option = (
                sql_query  => $sql_query_L,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [@$snp_levels],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@snp_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- snp_indel
#----------------------------------------------------------#
#
my $snp_indel = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ssw', 'ssw_id' ) ) {
        return;
    }

    my $sheet_name = 'snp_indel';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q~
            # header of Table distance
            SELECT 'distance', 'AVG_indel', 'COUNT', 'STD_indel'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q~
            # distance effect
            SELECT AVG(s.ssw_distance) distance,
                   AVG(w.window_indel) AVG_indel,
                   COUNT(w.window_indel) COUNT,
                   STD(w.window_indel) STD_indel
            FROM ssw s, window w
            WHERE s.window_id = w.window_id
            GROUP BY s.ssw_distance
        ~;
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write footer
    {
        $sheet_row += 2;
        my $sql_query = q~
            # distance effect
            SELECT AVG(s.ssw_distance) distance,
                   AVG(w.window_indel) AVG_indel,
                   COUNT(w.window_indel) COUNT,
                   STD(w.window_indel) STD_indel
            FROM ssw s, window w
            WHERE s.window_id = w.window_id
        ~;
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_occured
#----------------------------------------------------------#
my $indel_occured = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my $sheet_name = 'indel_occured';
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
        my $query_name = 'Original';
        my $sql_query  = q{
            SELECT  i.indel_occured,
                    e.indel_feature3,
                    AVG(i.indel_length) AVG_length,
                    STD(i.indel_length) STD_length,
                    COUNT(*) COUNT
            FROM indel i, indel_extra e
            WHERE i.indel_id = e.indel_id
            GROUP BY i.indel_occured, e.indel_feature3
            ORDER BY i.indel_occured DESC
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
        my $query_name = 'Another_ref';
        my $sql_query  = q{
            SELECT  e.indel_feature4,
                    e.indel_feature3,
                    AVG(i.indel_length) AVG_length,
                    STD(i.indel_length) STD_length,
                    COUNT(*) COUNT
            FROM indel i, indel_extra e
            WHERE i.indel_id = e.indel_id
            GROUP BY e.indel_feature4, e.indel_feature3
            ORDER BY e.indel_feature4 DESC
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
        my $query_name = 'Occured_vanish';
        my $sql_query  = q{
            SELECT  CONCAT(i.indel_occured, "->") `vanish`, 
                    e.indel_feature3 `slippage-like`,
                    AVG(i.indel_length) AVG_length,
                    STD(i.indel_length) STD_length,
                    COUNT(*) COUNT
            FROM indel i, indel_extra e
            WHERE i.indel_id = e.indel_id
            AND e.indel_feature4 IS NULL
            GROUP BY CONCAT(i.indel_occured, "->"),
                    e.indel_feature3
            ORDER BY CONCAT(i.indel_occured, "->") DESC
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
# worksheet -- snp_occured
#----------------------------------------------------------#
my $snp_occured = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'snp_extra', 'snp_feature4' ) ) {
        return;
    }

    my $sheet_name = 'snp_occured';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'Item';
        my $sql_query  = q{
            # header of Table
            SELECT 'SNP_occured', 'COUNT'
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
        my $query_name = 'Original';
        my $sql_query  = q{
            SELECT  s.snp_occured,
                    COUNT(*) COUNT
            FROM snp s
            GROUP BY s.snp_occured
            ORDER BY s.snp_occured DESC
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
        my $query_name = 'Another_ref';
        my $sql_query  = q{
            SELECT  e.snp_feature4,
                    COUNT(*) COUNT
            FROM snp_extra e
            GROUP BY e.snp_feature4
            ORDER BY e.snp_feature4 DESC
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
        my $query_name = 'Occured_not_change';
        my $sql_query  = q{
            SELECT  CONCAT(s.snp_occured, "->", e.snp_feature4) `change`,
                    COUNT(*) COUNT
            FROM snp s, snp_extra e
            WHERE s.snp_id = e.snp_id
            AND e.snp_feature4 IS NOT NULL
            AND s.snp_occured = e.snp_feature4
            GROUP BY CONCAT(s.snp_occured, "->", e.snp_feature4)
            ORDER BY CONCAT(s.snp_occured, "->", e.snp_feature4) DESC
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
            SELECT  CONCAT(s.snp_occured, "->", e.snp_feature4) `change`,
                    COUNT(*) COUNT
            FROM snp s, snp_extra e
            WHERE s.snp_id = e.snp_id
            AND e.snp_feature4 IS NOT NULL
            AND s.snp_occured != e.snp_feature4
            GROUP BY CONCAT(s.snp_occured, "->", e.snp_feature4)
            ORDER BY CONCAT(s.snp_occured, "->", e.snp_feature4) DESC
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
        my $query_name = 'Occured_vanish';
        my $sql_query  = q{
            SELECT  CONCAT(s.snp_occured, "->") `vanish`, 
                    COUNT(*) COUNT
            FROM snp s, snp_extra e
            WHERE s.snp_id = e.snp_id
            AND e.snp_feature4 IS NULL
            GROUP BY CONCAT(s.snp_occured, "->")
            ORDER BY CONCAT(s.snp_occured, "->") DESC
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
# worksheet -- indel_distance_slip
#----------------------------------------------------------#
#
my $indel_distance_slip = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature3' ) ) {
        return;
    }

    my @slip_levels = ( [ 'non-slip', 0, 0 ], [ 'slip', 1, 1 ], );

    my $write_sheet = sub {
        my ($slip_levels) = @_;
        my $sheet_name = 'indel_distance_' . $slip_levels->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q~
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            ~;
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q~
                # distance effect
                SELECT isw_distance distance,
                       AVG(isw_pi) AVG_pi,
                       AVG(isw_d_indel) AVG_d_indel,
                       AVG(isw_d_noindel) AVG_d_noindel,
                       AVG(isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG(isw_d_indel) / AVG(isw_d_noindel)  `Di/Dn`
                FROM    isw,
                       (SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN ? AND ?
                        AND isw_type = 'R'
                        AND isw.indel_id = i.indel_id
                        UNION 
                        SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN ? AND ?
                        AND isw_type = 'L'
                        AND isw.prev_indel_id = i.indel_id) i
                WHERE  isw_distance >= 0
                AND isw_d_indel IS NOT NULL 
                AND isw.isw_id = i.isw_id
                GROUP BY isw_distance
            ~;
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [
                    $slip_levels->[1], $slip_levels->[2],
                    $slip_levels->[1], $slip_levels->[2]
                ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        # write footer
        {
            $sheet_row += 2;
            my $sql_query = q~
                # isw average
                SELECT 'Total',
                       AVG(isw_pi) AVG_pi,
                       AVG(isw_d_indel) AVG_d_indel,
                       AVG(isw_d_noindel) AVG_d_noindel,
                       AVG(isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG(isw_d_indel) / AVG(isw_d_noindel)  `Di/Dn`
                FROM isw,
                       (SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN ? AND ?
                        AND isw_type = 'R'
                        AND isw.indel_id = i.indel_id
                        UNION 
                        SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN ? AND ?
                        AND isw_type = 'L'
                        AND isw.prev_indel_id = i.indel_id) i
                WHERE  isw_distance >= 0
                AND isw_d_indel IS NOT NULL 
                AND isw.isw_id = i.isw_id
            ~;
            my %option = (
                sql_query      => $sql_query,
                sheet_row      => $sheet_row,
                sheet_col      => $sheet_col,
                content_format => 'TOTAL',
                bind_value     => [
                    $slip_levels->[1], $slip_levels->[2],
                    $slip_levels->[1], $slip_levels->[2]
                ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@slip_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_slip_di_group
#----------------------------------------------------------#
#
my $indel_slip_di_group = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature3' ) ) {
        return;
    }

    my $sheet_name = 'indel_slip_di_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'indel_slip_di_group';
        my $sql_query  = q~
            # header
            SELECT 'isw_distance', 'AVG_d_indel', 'COUNT', 'STD_d_indel'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    my @indel_group = ( [ 0, 0 ], [ 1, 1 ], );

    # write contents
    {
        my $sql_query_R = q~
            # indel repeat effect for R windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance,
                   AVG(isw_d_indel) AVG_d_indel,
                   COUNT(*) COUNT,
                   STD(isw_d_indel) STD_d_indel
            FROM indel, isw, indel_extra
            WHERE isw.isw_type = 'R'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.indel_id = indel.indel_id
            AND indel.indel_id = indel_extra.indel_id
            AND indel_extra.indel_feature3 BETWEEN ? AND ?
            AND indel.indel_seq NOT LIKE "%N%"
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
            ORDER BY CONCAT(isw.isw_type, isw.isw_distance) DESC
        ~;
        my $sql_query_L = q~
            # indel repeat effect for L windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance,
                   AVG(isw_d_indel) AVG_d_indel,
                   COUNT(*) COUNT,
                   STD(isw_d_indel) STD_d_indel
            FROM indel, isw, indel_extra
            WHERE isw.isw_type = 'L'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.prev_indel_id = indel.indel_id
            AND indel.indel_id = indel_extra.indel_id
            AND indel_extra.indel_feature3 BETWEEN ? AND ?
            AND indel.indel_seq NOT LIKE "%N%"
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
        ~;
        my %option = (
            sql_query_1 => $sql_query_R,
            sql_query_2 => $sql_query_L,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_slip_dn_group
#----------------------------------------------------------#
#
my $indel_slip_dn_group = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature3' ) ) {
        return;
    }

    my $sheet_name = 'indel_slip_dn_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'indel_slip_dn_group';
        my $sql_query  = q~
            # header
            SELECT 'isw_distance', 'AVG_d_noindel', 'COUNT', 'STD_d_noindel'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    my @indel_group = ( [ 0, 0 ], [ 1, 1 ], );

    # write contents
    {
        my $sql_query_R = q~
            # indel repeat effect for R windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance,
                   AVG(isw_d_noindel) AVG_d_noindel,
                   COUNT(*) COUNT,
                   STD(isw_d_noindel) STD_d_noindel
            FROM indel, isw, indel_extra
            WHERE isw.isw_type = 'R'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.indel_id = indel.indel_id
            AND indel.indel_id = indel_extra.indel_id
            AND indel_extra.indel_feature3 BETWEEN ? AND ?
            AND indel.indel_seq NOT LIKE "%N%"
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
            ORDER BY CONCAT(isw.isw_type, isw.isw_distance) DESC
        ~;
        my $sql_query_L = q~
            # indel repeat effect for L windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance,
                   AVG(isw_d_noindel) AVG_d_noindel,
                   COUNT(*) COUNT,
                   STD(isw_d_noindel) STD_d_noindel
            FROM indel, isw, indel_extra
            WHERE isw.isw_type = 'L'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.prev_indel_id = indel.indel_id
            AND indel.indel_id = indel_extra.indel_id
            AND indel_extra.indel_feature3 BETWEEN ? AND ?
            AND indel.indel_seq NOT LIKE "%N%"
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
        ~;
        my %option = (
            sql_query_1 => $sql_query_R,
            sql_query_2 => $sql_query_L,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_type_gc_10
#----------------------------------------------------------#
my $indel_type_gc_10 = sub {
    my $sheet_name = 'indel_type_gc_10';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # indel_type groups
    my @indel_types
        = ( [ 'Insertion', ['I'] ], [ 'Deletion', ['D'] ], [ 'Complex', ['C'] ],
        );

    # write header
    {
        my $query_name = 'indel_type_gc_10';
        my $sql_query  = q{
            # header of Table group_density
            SELECT 'indel_length', 'AVG_gc_ratio', 'COUNT', 'STD_gc_ratio'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    foreach (@indel_types) {
        $sheet_row++;
        my $sql_query = q{
            # indel_length distribution
            SELECT  indel_length,
                    AVG(indel_gc) AVG_indel_gc,
                    COUNT(*),
                    STD(indel_gc) STD_indel_gc
            FROM indel
            WHERE indel_type = ?
            AND indel_length <= 10
            GROUP BY indel_length
        };
        my %option = (
            sql_query  => $sql_query,
            query_name => $_->[0],
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            bind_value => $_->[1],
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_type_gc_100
#----------------------------------------------------------#
my $indel_type_gc_100 = sub {
    my $sheet_name = 'indel_type_gc_100';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # indel_type groups
    my @indel_types
        = ( [ 'Insertion', ['I'] ], [ 'Deletion', ['D'] ], [ 'Complex', ['C'] ],
        );

    # write header
    {
        my $query_name = 'indel_type_gc_100';
        my $sql_query  = q{
            # header of Table group_density
            SELECT 'indel_length', 'AVG_gc_ratio', 'COUNT', 'STD_gc_ratio'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    foreach (@indel_types) {
        $sheet_row++;
        my $sql_query = q{
            # indel_length distribution
            SELECT  CEIL(indel_length / 10) * 10,
                    AVG(indel_gc) AVG_indel_gc,
                    COUNT(*),
                    STD(indel_gc) STD_indel_gc
            FROM indel
            WHERE indel_type = ?
            AND indel_length <= 100
            GROUP BY CEIL(indel_length / 10)
        };
        my %option = (
            sql_query  => $sql_query,
            query_name => $_->[0],
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            bind_value => $_->[1],
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- insertion_length
#----------------------------------------------------------#
my $insertion_length = sub {
    my $sheet_name = 'insertion_length';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $sql_query = q{
        # insertion_length distribution
        SELECT  indel_length,
                COUNT(indel_length) indel_number,
                AVG(indel_gc) AVG_indel_gc,
                SUM(indel_length) indel_sum
        FROM indel
        WHERE indel_type = "I"
        GROUP BY indel_length
    };

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_highlight( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- deletion_length
#----------------------------------------------------------#
my $deletion_length = sub {
    my $sheet_name = 'deletion_length';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $sql_query = q{
        # insertion_length distribution
        SELECT  indel_length,
                COUNT(indel_length) indel_number,
                AVG(indel_gc) AVG_indel_gc,
                SUM(indel_length) indel_sum
        FROM indel
        WHERE indel_type = "D"
        GROUP BY indel_length
    };

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_highlight( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- complex_length
#----------------------------------------------------------#
my $complex_length = sub {
    my $sheet_name = 'complex_length';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $sql_query = q{
        # insertion_length distribution
        SELECT  indel_length,
                COUNT(indel_length) indel_number,
                AVG(indel_gc) AVG_indel_gc,
                SUM(indel_length) indel_sum
        FROM indel
        WHERE indel_type = "C"
        GROUP BY indel_length
    };

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_highlight( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_base_change
#----------------------------------------------------------#
#
my $snp_base_change = sub {
    my $sheet_name = 'snp_base_change';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'snp_base_change';
        my $sql_query  = q~
            # header
            SELECT 'base_change', 'snp_number'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'Target';
        my $sql_query  = q~
        SELECT  s.base_change, s.snp_number / total.total * 100
        FROM    (SELECT count(*) total
                FROM snp
                WHERE snp_occured = "T"
                ) total,
                (SELECT CONCAT(ref_base, "->", target_base) base_change,
                        COUNT(snp_id) snp_number
                FROM snp
                WHERE snp_occured = "T"
                GROUP BY CONCAT(ref_base, "->", target_base)
                ) s
        ~;
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {
        $sheet_row++;
        my $query_name = 'Query';
        my $sql_query  = q~
        SELECT  s.base_change, s.snp_number / total.total * 100
        FROM    (SELECT count(*) total
                FROM snp
                WHERE snp_occured = "Q"
                ) total,
                (SELECT CONCAT(ref_base, "->", query_base) base_change,
                        COUNT(snp_id) snp_number
                FROM snp
                WHERE snp_occured = "Q"
                GROUP BY CONCAT(ref_base, "->", query_base)
                ) s
        ~;
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- target_snp
#----------------------------------------------------------#
#
my $target_snp = sub {
    my $sheet_name = 'target_snp';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'target_snp';
        my $sql_query  = q~
            # header
            SELECT 'base_change', 'snp_number'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'T_indel';
        my $sql_query  = q~
        SELECT  s.base_change, s.snp_number / total.total * 100
        FROM    (SELECT COUNT(*) total
                FROM indel, isw, snp
                WHERE indel.indel_occured = "T"
                AND snp.snp_occured = "T"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                ) total,
                (SELECT CONCAT(ref_base, "->", target_base) base_change,
                        COUNT(*) snp_number
                FROM indel, isw, snp
                WHERE indel.indel_occured = "T"
                AND snp.snp_occured = "T"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                GROUP BY CONCAT(ref_base, "->", target_base)
                ) s
        ~;
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {
        $sheet_row++;
        my $query_name = 'T_noindel';
        my $sql_query  = q~
        SELECT  s.base_change, s.snp_number / total.total * 100
        FROM    (SELECT COUNT(*) total
                FROM indel, isw, snp
                WHERE indel.indel_occured = "T"
                AND snp.snp_occured = "Q"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                ) total,
                (SELECT CONCAT(ref_base, "->", query_base) base_change,
                        COUNT(*) snp_number
                FROM indel, isw, snp
                WHERE indel.indel_occured = "T"
                AND snp.snp_occured = "Q"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                GROUP BY CONCAT(ref_base, "->", query_base)
                ) s
        ~;
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- query_snp
#----------------------------------------------------------#
#
my $query_snp = sub {
    my $sheet_name = 'query_snp';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'query_snp';
        my $sql_query  = q~
            # header
            SELECT 'base_change', 'snp_number'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'Q_indel';
        my $sql_query  = q~
        SELECT  s.base_change, s.snp_number / total.total * 100
        FROM    (SELECT COUNT(*) total
                FROM indel, isw, snp
                WHERE indel.indel_occured = "Q"
                AND snp.snp_occured = "Q"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                ) total,
                (SELECT CONCAT(ref_base, "->", query_base) base_change,
                        COUNT(*) snp_number
                FROM indel, isw, snp
                WHERE indel.indel_occured = "Q"
                AND snp.snp_occured = "Q"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                GROUP BY CONCAT(ref_base, "->", query_base)
                ) s
        ~;
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {
        $sheet_row++;
        my $query_name = 'Q_noindel';
        my $sql_query  = q~
        SELECT  s.base_change, s.snp_number / total.total * 100
        FROM    (SELECT COUNT(*) total
                FROM indel, isw, snp
                WHERE indel.indel_occured = "Q"
                AND snp.snp_occured = "T"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                ) total,
                (SELECT CONCAT(ref_base, "->", target_base) base_change,
                        COUNT(*) snp_number
                FROM indel, isw, snp
                WHERE indel.indel_occured = "Q"
                AND snp.snp_occured = "T"
                AND indel.indel_id = isw.indel_id
                AND isw.isw_id = snp.isw_id
                GROUP BY CONCAT(ref_base, "->", target_base)
                ) s
        ~;
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $query_name,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_snp
#----------------------------------------------------------#
my $distance_snp = sub {
    my $sheet_name = 'distance_snp';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # Six base pair groups
    my @base_pair = qw/
        A->C A->G A->T
        C->A C->G C->T
        G->A G->C G->T
        T->A T->C T->G
        /;

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ '', @base_pair ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
        ~;
        my $sql_query2 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) = ?
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
        ~;
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp2( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_tri_trv
#----------------------------------------------------------#
my $distance_tri_trv = sub {
    my $sheet_name = 'distance_tri_trv';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # Six base pair groups
    my @base_pair = (
        [ 'Tri', qw/AG CT GA TC/ ],
        [ 'Trv', qw/AC AT CA CG GC GT TA TG/ ],

        #['G->A,T->C', qw/GA TC/],
        #['A->T,T->A', qw/AT TA/],
    );
    my @headers;
    foreach (@base_pair) {
        push @headers, $_->[0];
    }

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ '', @headers ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN 5 AND 10
            
        ~;
        my $sql_query2 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN 5 AND 10
        ~;
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp3( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_snp_non_cpg
#----------------------------------------------------------#
my $distance_snp_non_cpg = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'snp', 'snp_cpg' ) ) {
        return;
    }

    my $sheet_name = 'distance_snp_non_cpg';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # Six base pair groups
    my @base_pair = qw/
        A->C A->G A->T
        C->A C->G C->T
        G->A G->C G->T
        T->A T->C T->G
        /;

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ '', @base_pair ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_cpg = 0
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
        ~;
        my $sql_query2 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_cpg = 0
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) = ?
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
        ~;
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp2( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_tri_trv_non_cpg
#----------------------------------------------------------#
my $distance_tri_trv_non_cpg = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'snp', 'snp_cpg' ) ) {
        return;
    }

    my $sheet_name = 'distance_tri_trv_non_cpg';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # Six base pair groups
    my @base_pair = (
        [ 'Tri', qw/AG CT GA TC/ ],
        [ 'Trv', qw/AC AT CA CG GC GT TA TG/ ],

        #['G->A,T->C', qw/GA TC/],
        #['A->T,T->A', qw/AT TA/],
    );
    my @headers;
    foreach (@base_pair) {
        push @headers, $_->[0];
    }

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ '', @headers ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_cpg = 0
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_cpg = 0
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN 5 AND 10
        ~;
        my $sql_query2 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_cpg = 0
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN -1 AND 10
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND s.snp_cpg = 0
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN 5 AND 10
        ~;
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp3( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_tri_trv_coding
#----------------------------------------------------------#
my $distance_tri_trv_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature1' ) ) {
        return;
    }

    my $sheet_name = 'distance_tri_trv_coding';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # Six base pair groups
    my @base_pair = (
        [ 'Tri', qw/AG CT GA TC/ ],
        [ 'Trv', qw/AC AT CA CG GC GT TA TG/ ],

        #['G->A,T->C', qw/GA TC/],
        #['A->T,T->A', qw/AT TA/],
    );
    my @headers;
    foreach (@base_pair) {
        push @headers, $_->[0];
    }

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ '', @headers ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN -1 AND 10
            AND e.isw_feature1 = 1
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN 5 AND 10
            AND e.isw_feature1 = 1
            
        ~;
        my $sql_query2 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN -1 AND 10
            AND e.isw_feature1 = 1
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN 5 AND 10
            AND e.isw_feature1 = 1
        ~;
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp3( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_tri_trv_non_coding
#----------------------------------------------------------#
my $distance_tri_trv_non_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature1' ) ) {
        return;
    }

    my $sheet_name = 'distance_tri_trv_non_coding';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # Six base pair groups
    my @base_pair = (
        [ 'Tri', qw/AG CT GA TC/ ],
        [ 'Trv', qw/AC AT CA CG GC GT TA TG/ ],

        #['G->A,T->C', qw/GA TC/],
        #['A->T,T->A', qw/AT TA/],
    );
    my @headers;
    foreach (@base_pair) {
        push @headers, $_->[0];
    }

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ '', @headers ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN -1 AND 10
            AND e.isw_feature1 = 0
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND i.isw_distance BETWEEN 5 AND 10
            AND e.isw_feature1 = 0
            
        ~;
        my $sql_query2 = q~
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN -1 AND 10
            AND e.isw_feature1 = 0
            GROUP BY i.isw_distance
            UNION
            SELECT '5-10', COUNT(s.snp_id) snp_number
            FROM snp s, isw i, isw_extra e
            WHERE s.isw_id = i.isw_id
            AND i.isw_id = e.isw_id
            AND s.snp_occured IN ("T", "Q")
            AND CONCAT(ref_base, IF(s.snp_occured = "T", target_base, query_base)) IN (in_list)
            AND i.isw_distance BETWEEN 5 AND 10
            AND e.isw_feature1 = 0
        ~;
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp3( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_cpg
#----------------------------------------------------------#
#
my $distance_cpg = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw', 'isw_cpg_pi' ) ) {
        return;
    }

    my $sheet_name = 'distance_cpg';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q~
            # header of Table distance
            SELECT  'distance', 'CpG/100bp', 'COUNT', 'STD'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q{
            SELECT isw_distance distance,
                   AVG(isw_cpg_pi) AVG_cpg,
                   COUNT(isw_cpg_pi) COUNT,
                   STD(isw_cpg_pi) STD
            FROM isw i
            GROUP BY isw_distance
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write footer
    {
        $sheet_row += 2;
        my $sql_query = q{
            SELECT 'Total',
                   AVG(isw_cpg_pi) AVG_cpg,
                   COUNT(isw_cpg_pi) COUNT,
                   STD(isw_cpg_pi) STD
            FROM isw i
        };
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- d_d_ref
#----------------------------------------------------------#
my $d_d_ref = sub {
    my $sheet_name = 'd_d_ref';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [
                qw{item AVG_pi align_comparables align_identities
                    align_differences align_gaps align_ns align_error
                    }
            ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # init objects
    my $dbh  = $write_obj->dbh;
    my $fmt  = $write_obj->format;
    my @cols = @{ $write_obj->columns };

    # alignments
    my $align_query = q~
        SELECT align_id
        FROM align 
    ~;
    my $align_sth = $dbh->prepare($align_query);

    my $seq_query = q~
        SELECT t.target_seq,
               q.query_seq,
               r.ref_seq
        FROM align a, target t, query q, reference r
        WHERE a.align_id = ?
        AND a.align_id = q.align_id
        AND a.align_id = t.align_id
        AND a.align_id = r.align_id
    ~;
    my $seq_sth = $dbh->prepare($seq_query);

    my $result_matrix = [ [] ];
    $align_sth->execute();
    while ( my @row = $align_sth->fetchrow_array ) {
        my ($align_id) = @row;
        print "$align_id\r";

        $seq_sth->execute($align_id);
        my ( $target_seq, $query_seq, $ref_seq ) = $seq_sth->fetchrow_array;

        my @result;
        $result[0] = &pair_seq_stat( $target_seq, $query_seq );
        $result[1] = &pair_seq_stat( $target_seq, $ref_seq );
        $result[2] = &pair_seq_stat( $query_seq,  $ref_seq );

        foreach my $i ( 0 .. 2 ) {
            foreach my $j ( 1 .. 6 ) {
                $result_matrix->[$i][$j] += $result[$i]->[$j];
            }
        }
    }

    # calc total pi
    my @names = qw{TQ TR QR};
    foreach my $i ( 0 .. 2 ) {
        $result_matrix->[$i][0]
            = $result_matrix->[$i][3] / $result_matrix->[$i][1];
        unshift @{ $result_matrix->[$i] }, $names[$i];
        for ( my $j = 0; $j < scalar @{ $result_matrix->[$i] }; $j++ ) {
            $sheet->write(
                $sheet_row,
                $j + $sheet_col,
                $result_matrix->[$i][$j],
                $fmt->{NORMAL}
            );
        }
        $sheet_row++;
    }

    print "Sheet \"$sheet_name\" has been generated.\n";

};

#----------------------------------------------------------#
# worksheet -- ds_dns
#----------------------------------------------------------#
my $ds_dns = sub {
    my $sheet_name = 'ds_dns';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [qw{distance AVG_d_snp AVG_d_nosnp COUNT Ds/Dns}],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # init objects
    my $dbh  = $write_obj->dbh;
    my $fmt  = $write_obj->format;
    my @cols = @{ $write_obj->columns };

    my $sql_query1 = qq~
        SELECT  i.isw_id, i.indel_id,
                i.isw_distance, i.isw_density,
                (floor((i.isw_density + 1) / 2) - i.isw_distance) reverse_distance,
                s.snp_occured, COUNT(s.snp_id) * 1.0000 / i.isw_length d
        FROM isw i
        LEFT JOIN snp s ON s.isw_id = i.isw_id
        GROUP BY i.isw_id, s.snp_occured
        HAVING ((s.snp_occured != 'N'
        OR s.snp_occured IS NULL)
        AND i.isw_distance > 30
        AND i.isw_density > 70)
    ~;

    my $sth = $dbh->prepare($sql_query1);
    $sth->execute();

    my $each_interval = {};
    while ( my @row = $sth->fetchrow_array ) {
        unless ( $row[5] ) {
            next;
        }
        $each_interval->{ $row[1] }{ $row[4] }{ $row[5] } += $row[-1];
    }

    my $sql_query2 = qq~
        SELECT  distinct i.isw_id, i.indel_id,
                i.isw_distance, i.isw_density,
                (floor((i.isw_density + 1) / 2) - i.isw_distance) reverse_distance
        FROM isw i
        LEFT JOIN snp s ON s.isw_id = i.isw_id
        GROUP BY i.isw_id, s.snp_occured
        HAVING i.isw_distance > 30
        AND i.isw_density > 70
    ~;

    $sth = $dbh->prepare($sql_query2);
    $sth->execute();

    while ( my @row = $sth->fetchrow_array ) {
        $each_interval->{ $row[1] }{ $row[4] }{count}++;
    }

    foreach my $interval ( keys %$each_interval ) {
        foreach my $rd ( keys %{ $each_interval->{$interval} } ) {
            if ( exists $each_interval->{$interval}{$rd}{T} ) {
                $each_interval->{$interval}{$rd}{T}
                    = $each_interval->{$interval}{$rd}{T}
                    / $each_interval->{$interval}{$rd}{count};
            }
            else {
                $each_interval->{$interval}{$rd}{T} = 0;
            }
            if ( exists $each_interval->{$interval}{$rd}{Q} ) {
                $each_interval->{$interval}{$rd}{Q}
                    = $each_interval->{$interval}{$rd}{Q}
                    / $each_interval->{$interval}{$rd}{count};
            }
            else {
                $each_interval->{$interval}{$rd}{Q} = 0;
            }
            delete $each_interval->{$interval}{$rd}{count};
        }
    }

    # Target SNP and Query SNP probability
    my $interval_tq = {};
    foreach my $interval ( keys %$each_interval ) {
        foreach my $rd ( keys %{ $each_interval->{$interval} } ) {
            $interval_tq->{$interval}{T} += $each_interval->{$interval}{$rd}{T};
            $interval_tq->{$interval}{Q} += $each_interval->{$interval}{$rd}{Q};
        }
    }

    my $rd_count = {};
    foreach ( 1 .. 1000 ) {
        srand( time ^ $$ );
        foreach my $interval ( keys %$each_interval ) {
            foreach my $rd ( keys %{ $each_interval->{$interval} } ) {
                my $t_frq     = $interval_tq->{$interval}{T};
                my $q_frq     = $interval_tq->{$interval}{Q};
                my $total_frq = $t_frq + $q_frq;
                unless ($total_frq) {
                    next;
                }
                my ( $Ds, $Dns );
                if ( $t_frq == 0 ) {
                    ( $Ds, $Dns ) = qw/Q T/;
                }
                elsif ( $q_frq == 0 ) {
                    ( $Ds, $Dns ) = qw/T Q/;
                }
                else {
                    my $random = rand($total_frq);
                    $Ds  = $random < $t_frq ? 'T' : 'Q';
                    $Dns = $random < $t_frq ? 'Q' : 'T';
                }
                $rd_count->{$rd}{Ds}  += $each_interval->{$interval}{$rd}{$Ds};
                $rd_count->{$rd}{Dns} += $each_interval->{$interval}{$rd}{$Dns};
                $rd_count->{$rd}{count}++;
            }
        }
    }

    foreach my $rd ( sort { $a <=> $b } keys %{$rd_count} ) {
        $rd_count->{$rd}{Ds}  = $rd_count->{$rd}{Ds} / $rd_count->{$rd}{count};
        $rd_count->{$rd}{Dns} = $rd_count->{$rd}{Dns} / $rd_count->{$rd}{count};
        $sheet->write( $sheet_row, 0 + $sheet_col, $rd, $fmt->{NORMAL} );
        $sheet->write(
            $sheet_row,
            0 + 1 + $sheet_col,
            $rd_count->{$rd}{Ds},
            $fmt->{NORMAL}
        );
        $sheet->write(
            $sheet_row,
            0 + 2 + $sheet_col,
            $rd_count->{$rd}{Dns},
            $fmt->{NORMAL}
        );
        $sheet->write(
            $sheet_row,
            0 + 3 + $sheet_col,
            $rd_count->{$rd}{count},
            $fmt->{NORMAL}
        );
        $sheet->write(
            $sheet_row,
            0 + 4 + $sheet_col,
            $rd_count->{$rd}{Dns}
            ? $rd_count->{$rd}{Ds} / $rd_count->{$rd}{Dns}
            : '#N/A',
            $fmt->{NORMAL}
        );
        $sheet_row++;
    }

    print "Sheet \"$sheet_name\" has been generated.\n";

    #DumpFile("interval.yaml", $each_interval);
    #DumpFile("interval_tq.yaml", $interval_tq);
    #DumpFile("rd.yaml", $rd_count);

};

#----------------------------------------------------------#
# worksheet -- di_dn_ttest
#----------------------------------------------------------#
#
my $di_dn_ttest = sub {
    my $sheet_name = 'di_dn_ttest';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q~
            SELECT  'AVG_distance', 'AVG_pi', 'AVG_d_indel', 'AVG_d_noindel',
                    'Di-Dn', 'Di/Dn', 'COUNT', 'P_value'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    my @group_distance = ( [0], [1], [2], [3], [4], [5], [ 5 .. 10 ], );

    my @p_values;

    # ttest
    {
        my $sql_query = qq~
            SELECT isw_d_indel,
                   isw_d_noindel
            FROM    isw,
                   (SELECT isw_id
                    FROM   isw,
                           indel i,
                           indel_extra e
                    WHERE  i.indel_id = e.indel_id
                    AND e.indel_feature3 BETWEEN 0 AND 0
                    AND isw_type = 'R'
                    AND isw.indel_id = i.indel_id
                    UNION 
                    SELECT isw_id
                    FROM   isw,
                           indel i,
                           indel_extra e
                    WHERE  i.indel_id = e.indel_id
                    AND e.indel_feature3 BETWEEN 0 AND 0
                    AND isw_type = 'L'
                    AND isw.prev_indel_id = i.indel_id) i
            WHERE  isw_distance >= 0
            AND isw_d_indel IS NOT NULL 
            AND isw.isw_id = i.isw_id
            AND isw_distance IN 
        ~;

        foreach (@group_distance) {
            my @range      = @$_;
            my $in_list    = '(' . join( ',', @range ) . ')';
            my $sql_query2 = $sql_query . $in_list;
            $sql_query2 .= " \nLIMIT 10000";    # prevent to exceed excel limit

            my %option  = ( sql_query => $sql_query2, );
            my $ttest   = $write_obj->column_ttest( \%option );
            my $p_value = $ttest->{t_prob};

            $p_value = "NA" if !defined $p_value;
            push @p_values, [$p_value];

            {
                my $ttest_sheet_name = "ttest_" . join( '-', @range[ 0, -1 ] );
                my $stat_sheet;

                my ( $stat_sheet_row, $stat_sheet_col ) = ( 0, 0 );
                my %option = (
                    sheet_row => $stat_sheet_row,
                    sheet_col => $stat_sheet_col,
                    header    => [qw{d_indel d_noindel}],
                );
                ( $stat_sheet, $stat_sheet_row )
                    = $write_obj->write_header_direct( $ttest_sheet_name,
                    \%option );

                %option = (
                    sql_query => $sql_query2,
                    sheet_row => $stat_sheet_row,
                    sheet_col => $stat_sheet_col,
                );
                $write_obj->write_content_direct( $stat_sheet, \%option );
            }
        }
    }

    {
        my $sql_query = qq~
            SELECT AVG(isw_distance) AVG_distance,
                   AVG(isw_pi) AVG_pi,
                   AVG(isw_d_indel) AVG_d_indel,
                   AVG(isw_d_noindel) AVG_d_noindel,
                   AVG(isw_d_indel) - AVG(isw_d_noindel)  `Di-Dn`,
                   AVG(isw_d_indel) / AVG(isw_d_noindel)  `Di/Dn`,
                   COUNT(*) COUNT
            FROM    isw,
                   (SELECT isw_id
                    FROM   isw,
                           indel i,
                           indel_extra e
                    WHERE  i.indel_id = e.indel_id
                    AND e.indel_feature3 BETWEEN 0 AND 0
                    AND isw_type = 'R'
                    AND isw.indel_id = i.indel_id
                    UNION 
                    SELECT isw_id
                    FROM   isw,
                           indel i,
                           indel_extra e
                    WHERE  i.indel_id = e.indel_id
                    AND e.indel_feature3 BETWEEN 0 AND 0
                    AND isw_type = 'L'
                    AND isw.prev_indel_id = i.indel_id) i
            WHERE  isw_distance >= 0
            AND isw_d_indel IS NOT NULL 
            AND isw.isw_id = i.isw_id
            AND isw_distance IN 
        ~;
        my %option = (
            sql_query     => $sql_query,
            sheet_row     => $sheet_row,
            sheet_col     => $sheet_col,
            group         => \@group_distance,
            append_column => \@p_values,
        );
        ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_distance_type
#----------------------------------------------------------#
#
my $indel_distance_type = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel', 'indel_type' ) ) {
        return;
    }

    # indel_type groups
    my @indel_types = ( [ 'insertion', 'I' ], [ 'deletion', 'D' ], );

    my $write_sheet = sub {
        my $level      = shift;
        my $sheet_name = 'indel_distance_' . $level->[0];
        my $indel_type = $level->[1];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                # distance effect
                SELECT isw_distance distance,
                       AVG(isw_pi) AVG_pi,
                       AVG(isw_d_indel) AVG_d_indel,
                       AVG(isw_d_noindel) AVG_d_noindel,
                       AVG(isw_d_complex) AVG_d_complex,
                       COUNT(* ) COUNT,
                       AVG(isw_d_indel) / AVG(isw_d_noindel) `Di/Dn`
                FROM   isw,
                       (SELECT isw_id
                        FROM   isw,
                               indel i
                        WHERE i.indel_occured IS NOT NULL
                        AND i.indel_type = ?
                        AND isw_type = 'R'
                        AND isw.indel_id = i.indel_id
                        UNION 
                        SELECT isw_id
                        FROM   isw,
                               indel i
                        WHERE i.indel_occured IS NOT NULL
                        AND i.indel_type = ?
                        AND isw_type = 'L'
                        AND isw.prev_indel_id = i.indel_id) i
                WHERE  isw_distance >= 0
                AND isw_d_indel IS NOT NULL 
                AND isw.isw_id = i.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $indel_type, $indel_type ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@indel_types) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_distance_tq
#----------------------------------------------------------#
#
my $indel_distance_tq = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel', 'indel_occured' ) ) {
        return;
    }

    my @levels = ( [ 'target', 'T', 'T' ], [ 'query', 'Q', 'Q' ], );

    my $write_sheet = sub {
        my ($levels) = @_;
        my $sheet_name = 'indel_distance_' . $levels->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q~
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            ~;
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q~
                # distance effect
                SELECT isw_distance distance,
                       AVG(isw_pi) AVG_pi,
                       AVG(isw_d_indel) AVG_d_indel,
                       AVG(isw_d_noindel) AVG_d_noindel,
                       AVG(isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG(isw_d_indel) / AVG(isw_d_noindel)  `Di/Dn`
                FROM    isw,
                        (SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND i.indel_occured = ?
                        AND isw_type = 'R'
                        AND isw.indel_id = i.indel_id
                        UNION 
                        SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND i.indel_occured = ?
                        AND isw_type = 'L'
                        AND isw.prev_indel_id = i.indel_id) i                        
                WHERE  isw_distance >= 0
                AND isw_d_indel IS NOT NULL 
                AND isw.isw_id = i.isw_id
                GROUP BY isw_distance
            ~;
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $levels->[1], $levels->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        # write footer
        {
            $sheet_row += 2;
            my $sql_query = q~
                # isw average
                SELECT 'Total',
                       AVG(isw_pi) AVG_pi,
                       AVG(isw_d_indel) AVG_d_indel,
                       AVG(isw_d_noindel) AVG_d_noindel,
                       AVG(isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG(isw_d_indel) / AVG(isw_d_noindel)  `Di/Dn`
                FROM isw,
                        (SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND i.indel_occured = ?
                        AND isw_type = 'R'
                        AND isw.indel_id = i.indel_id
                        UNION 
                        SELECT isw_id
                        FROM   isw,
                               indel i,
                               indel_extra e
                        WHERE  i.indel_id = e.indel_id
                        AND e.indel_feature3 BETWEEN 0 AND 0
                        AND i.indel_occured = ?
                        AND isw_type = 'L'
                        AND isw.prev_indel_id = i.indel_id) i     
                WHERE  isw_distance >= 0
                AND isw_d_indel IS NOT NULL 
                AND isw.isw_id = i.isw_id
            ~;
            my %option = (
                sql_query      => $sql_query,
                sheet_row      => $sheet_row,
                sheet_col      => $sheet_col,
                content_format => 'TOTAL',
                bind_value     => [ $levels->[1], $levels->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@levels) {
        &$write_sheet($_);
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 )  { &$summary;             next; }
    if ( $n == 2 )  { &$distance;            &$combined_distance; next; }
    if ( $n == 3 )  { &$distance_coding;     &$distance_non_coding; next; }
    if ( $n == 4 )  { &$snp_distance;        &$snp_indel; next; }
    if ( $n == 8 )  { &$indel_occured;       &$snp_occured; next; }
    if ( $n == 10 ) { &$indel_distance_slip; next; }
    if ( $n == 11 ) { &$indel_slip_di_group; &$indel_slip_dn_group; next; }
    if ( $n == 12 ) { &$indel_type_gc_10;    &$indel_type_gc_100; next; }
    if ( $n == 15 ) { &$insertion_length;    &$deletion_length; next; }
    if ( $n == 16 ) { &$complex_length;      next; }
    if ( $n == 17 ) { &$snp_base_change;     next; }
    if ( $n == 18 ) { &$target_snp;          &$query_snp; next; }
    if ( $n == 19 ) { &$distance_snp;        &$distance_tri_trv; next; }
    if ( $n == 20 ) {
        &$distance_snp_non_cpg;
        &$distance_tri_trv_non_cpg;
        next;
    }
    if ( $n == 21 ) {
        &$distance_tri_trv_coding;
        &$distance_tri_trv_non_coding;
        next;
    }
    if ( $n == 22 ) { &$distance_cpg;     next; }
    if ( $n == 23 ) { &$distance_dir_dnr; &$distance_dtr_dqr; next; }
    if ( $n == 24 ) { &$target_dtr_dqr;   &$query_dtr_dqr; next; }

    if ( $n == 52 ) { &$ds_dns;                  next; }
    if ( $n == 54 ) { &$snp_distance_non_slip;   next; }
    if ( $n == 56 ) { &$snp_LR_distance;         next; }
    if ( $n == 58 ) { &$di_dn_ttest;             next; }
    if ( $n == 60 ) { &$d_d_ref;                 next; }
    if ( $n == 66 ) { &$indel_distance_type;     next; }
    if ( $n == 71 ) { &$indel_distance_tq;       next; }
    if ( $n == 72 ) { &$distance_window0eq0;     next; }
    if ( $n == 73 ) { &$distance_iswpi_le_2;     next; }
    if ( $n == 74 ) { &$distance_indel_bottom90; next; }
}

$stopwatch->end_message();
exit;

__END__

=head1 NAME

    three_stat_factory.pl - Generate statistical Excel files from alignDB

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
