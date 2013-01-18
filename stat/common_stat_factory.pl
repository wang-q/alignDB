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
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;

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
my $outfile = "";

my $help = 0;
my $man  = 0;

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
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
    $outfile = "$db.common.xlsx" unless $outfile;
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

    unless ($outfile) {
        my $runlist = $set->runlist;
        $outfile = "$db.common.$runlist.xlsx";
    }
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

#----------------------------#
# count freq
#----------------------------#
my $all_freq = $write_obj->get_freq;

if ( $all_freq < 2 ) {
    die "all_freq is $all_freq, are you sure this is a AlignDB database?\n";
}

#----------------------------------------------------------#
# worksheet -- basic
#----------------------------------------------------------#
my $basic = sub {
    my $sheet_name = 'basic';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{VALUE};
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
        my $query_name = 'No. of strains';
        my %option     = (
            query_name => $query_name,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            row        => [$all_freq],
        );
        ($sheet_row) = $write_obj->write_row_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'Target length (Mb)';
        my $sql_query  = q{
            SELECT  SUM(s.seq_length) / 1000000.00
            FROM    sequence s, target t
            WHERE   s.seq_id = t.seq_id
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
        my $query_name = 'Aligned length (Mb)';
        my $sql_query  = q{
            SELECT  SUM(a.align_length) / 1000000.00
            FROM    align a
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
        my $query_name = 'Indels per 100 bp';
        my $sql_query  = q{
            SELECT  SUM(a.align_indels) / SUM(a.align_comparables) * 100.0
            FROM    align a
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
        my $query_name = 'SNVs per 100 bp';
        my $sql_query  = q{
            SELECT  SUM(a.align_differences) * 1.0 / SUM(a.align_comparables) * 100.0
            FROM    align a
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
        my $query_name = 'D on average';
        my $sql_query  = q{
            SELECT -0.75 * log2( 1 - ( 4.0 / 3.0 ) * original.Pi )
            FROM (
                SELECT AVG(a.align_pi) Pi
                FROM align a
                ) original
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
        my $query_name = 'GC-content';
        my $sql_query  = q{
            SELECT sum(a.align_length * a.align_average_gc )
                   / sum(a.align_length)
            FROM align a
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
# worksheet -- process
#----------------------------------------------------------#
my $process = sub {
    my $sheet_name = 'process';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{Order Operation Duration Cmd_line};
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
            SELECT  meta_value
            FROM    meta m
            WHERE   meta_key IN ("a_operation","d_duration", "e_cmd_line")
        };
        my $dbh = $write_obj->dbh;

        my $array_ref = $dbh->selectcol_arrayref($sql_query);

        my $order = 1;
        while ( scalar @{$array_ref} ) {
            my @row = splice @{$array_ref}, 0, 3;
            ($sheet_row) = $write_obj->write_row_direct(
                $sheet,
                {   row       => [ $order, @row ],
                    sheet_row => $sheet_row,
                    sheet_col => $sheet_col,
                }
            );
            $order++;
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

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
        $column_stat->( 'align_length',       'align', 'align_length' );
        $column_stat->( 'indel_length',       'indel', 'indel_length' );
        $column_stat->( 'indel_left_extand',  'indel', 'left_extand' );
        $column_stat->( 'indel_right_extand', 'indel', 'right_extand' );
        $column_stat->(
            'indel_windows', 'isw', 'isw_length', 'WHERE isw_distance <= 0'
        );
        $column_stat->(
            'indel_free_windows', 'isw',
            'isw_length',         'WHERE isw_distance > 0'
        );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- pi_gc_cv
#----------------------------------------------------------#
my $pi_gc_cv = sub {
    {
        my $sheet_name = 'd1_pi_gc_cv';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('common-d1_pi_gc_cv-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    {
        my $sheet_name = 'd2_pi_gc_cv';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('common-d1_pi_gc_cv-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- comb_pi_gc_cv
#----------------------------------------------------------#
my $comb_pi_gc_cv = sub {

    # make combine
    my @combined;
    {
        my $thaw_sql   = $sql_file->retrieve('common-d1_combine-0');
        my $standalone = [ -1, 0 ];
        my %option     = (
            sql_query  => $thaw_sql->as_sql,
            threshold  => $combine,
            standalone => $standalone,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    {
        my $sheet_name = 'd1_comb_pi_gc_cv';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_obj     => $thaw_sql,
                sheet_row   => $sheet_row,
                sheet_col   => $sheet_col,
                combine_col => 'isw.isw_distance',
                combined    => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine_obj( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    # make combine again
    {
        my $thaw_sql = $sql_file->retrieve('common-d1_combine-0');
        $thaw_sql->replace( { distance => 'density' } );
        my $standalone = [ -1, 0 ];
        my %option = (
            sql_query  => $thaw_sql->as_sql,
            threshold  => $combine,
            standalone => $standalone,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    {
        my $sheet_name = 'd2_comb_pi_gc_cv';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_obj     => $thaw_sql,
                sheet_row   => $sheet_row,
                sheet_col   => $sheet_col,
                combine_col => 'isw.isw_density',
                combined    => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine_obj( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- group_distance
#----------------------------------------------------------#
my $group_distance = sub {
    my $sheet_name = 'group_distance';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{AVG_distance AVG_pi COUNT STD_pi SUM_length length_proportion};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # make last portion
    my ( $all_length, $last_portion );
    {
        my $thaw_sql = $sql_file->retrieve('common-d1_combine-0');
        my $portion  = 0.05;
        my %option   = (
            sql_query => $thaw_sql->as_sql,
            portion   => $portion,
        );
        ( $all_length, $last_portion )
            = $write_obj->make_last_portion( \%option );
    }

    my @group_distance = (
        [-1],
        [0],
        [1],
        [2],
        [ 3 .. 5 ],
        [ 6 .. 10 ],
        [ 11 .. 20 ],
        [ 21 .. 50 ],
        [ 51 .. 999 ],
        [ 3 .. 20 ],
        [ 21 .. 999 ],
        [ 101 .. 999 ],
        [0],
        [ 1 .. 5 ],
        [ 6 .. 10 ],
        [ 11 .. 15 ],
        [ 16 .. 20 ],
        [ 21 .. 25 ],
        $last_portion,
    );

    {    # write contents
        my $thaw_sql = $sql_file->retrieve('common-d1_pi_avg-0');
        $thaw_sql->add_select( "SUM(isw_length)", 'SUM_length' );
        $thaw_sql->add_select( "SUM(isw_length) / $all_length * 100",
            'length_proportion' );
        my %option = (
            sql_obj   => $thaw_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group_col => 'isw_distance',
            group     => \@group_distance,
        );
        ($sheet_row) = $write_obj->write_content_group_obj( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- group_density
#----------------------------------------------------------#
my $group_density = sub {
    my $sheet_name = 'group_density';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{AVG_density AVG_pi COUNT STD_pi SUM_length length_proportion};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # make last portion
    my ( $all_length, $last_portion );
    {
        my $thaw_sql = $sql_file->retrieve('common-d1_combine-0');
        $thaw_sql->replace( { distance => 'density' } );
        my $portion = 0.05;
        my %option  = (
            sql_query => $thaw_sql->as_sql,
            portion   => $portion,
        );
        ( $all_length, $last_portion )
            = $write_obj->make_last_portion( \%option );
    }

    my @group_density = (
        [-1],
        [0],
        [1],
        [2],
        [ 3 .. 5 ],
        [ 6 .. 10 ],
        [ 11 .. 20 ],
        [ 21 .. 50 ],
        [ 51 .. 999 ],
        [ 3 .. 20 ],
        [ 21 .. 999 ],
        [ 101 .. 999 ],
        $last_portion,
    );

    {    # write contents
        my $thaw_sql = $sql_file->retrieve('common-d1_pi_avg-0');
        $thaw_sql->add_select( "SUM(isw_length)", 'SUM_length' );
        $thaw_sql->add_select( "SUM(isw_length) / $all_length * 100",
            'length_proportion' );
        $thaw_sql->replace( { distance => 'density' } );
        my %option = (
            sql_obj   => $thaw_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group_col => 'isw_density',
            group     => \@group_density,
        );
        ($sheet_row) = $write_obj->write_content_group_obj( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- comb_coding
#----------------------------------------------------------#
my $comb_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw', 'isw_coding' ) ) {
        return;
    }

    my @type_levels = ( [ 'coding', 1, 1 ], [ 'non_coding', 0, 0 ], );

    my $write_sheet_d1 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        # make combine
        my @combined;
        {
            my $thaw_sql
                = $sql_file->retrieve('common-d1_make_combine_coding-0');
            my $standalone = [ -1, 0 ];
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                threshold  => $combine,
                standalone => $standalone,
                bind_value => [ $feature_1, $feature_2 ],
            );
            @combined = @{ $write_obj->make_combine( \%option ) };
        }

        {
            my $sheet_name = "d1_comb_$name";
            my $sheet;
            my ( $sheet_row, $sheet_col );

            my $thaw_sql = $sql_file->retrieve('common-d1_comb_coding-0');

            {    # write header
                my @headers = $thaw_sql->as_header;
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
                my %option = (
                    sql_obj     => $thaw_sql,
                    sheet_row   => $sheet_row,
                    sheet_col   => $sheet_col,
                    combine_col => 'isw.isw_distance',
                    combined    => \@combined,
                    bind_value  => [ $feature_1, $feature_2 ],
                );
                ($sheet_row)
                    = $write_obj->write_content_combine_obj( $sheet, \%option );
            }

            print "Sheet \"$sheet_name\" has been generated.\n";
        }
    };

    my $write_sheet_d2 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        # make combine
        my @combined;
        {
            my $thaw_sql
                = $sql_file->retrieve('common-d1_make_combine_coding-0');
            $thaw_sql->replace( { distance => 'density' } );
            my $standalone = [ -1, 0 ];
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                threshold  => $combine,
                standalone => $standalone,
                bind_value => [ $feature_1, $feature_2 ],
            );
            @combined = @{ $write_obj->make_combine( \%option ) };
        }

        {
            my $sheet_name = "d2_comb_$name";
            my $sheet;
            my ( $sheet_row, $sheet_col );

            my $thaw_sql = $sql_file->retrieve('common-d1_comb_coding-0');
            $thaw_sql->replace( { distance => 'density' } );

            {    # write header
                my @headers = $thaw_sql->as_header;
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
                my %option = (
                    sql_obj     => $thaw_sql,
                    sheet_row   => $sheet_row,
                    sheet_col   => $sheet_col,
                    combine_col => 'isw.isw_density',
                    combined    => \@combined,
                    bind_value  => [ $feature_1, $feature_2 ],
                );
                ($sheet_row)
                    = $write_obj->write_content_combine_obj( $sheet, \%option );
            }

            print "Sheet \"$sheet_name\" has been generated.\n";
        }
    };

    for (@type_levels) {
        $write_sheet_d1->($_);
    }

    for (@type_levels) {
        $write_sheet_d2->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- comb_slippage
#----------------------------------------------------------#
my $comb_slippage = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel', 'indel_slippage' ) ) {
        return;
    }

    my @type_levels = ( [ 'slippage', 1, 1 ], [ 'non_slippage', 0, 0 ], );

    my $write_sheet_d1 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        # make combine
        my @combined;
        {
            my $thaw_sql
                = $sql_file->retrieve('common-d1_make_combine_slippage-0');
            my $standalone = [ -1, 0 ];
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                threshold  => $combine,
                standalone => $standalone,
                bind_value => [ $feature_1, $feature_2 ],
            );
            @combined = @{ $write_obj->make_combine( \%option ) };
        }

        {
            my $sheet_name = "d1_comb_$name";
            my $sheet;
            my ( $sheet_row, $sheet_col );

            my $thaw_sql = $sql_file->retrieve('common-d1_comb_slippage-0');

            {    # write header
                my @headers = $thaw_sql->as_header;
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
                my %option = (
                    sql_obj     => $thaw_sql,
                    sheet_row   => $sheet_row,
                    sheet_col   => $sheet_col,
                    combine_col => 'isw.isw_distance',
                    combined    => \@combined,
                    bind_value  => [ $feature_1, $feature_2 ],
                );
                ($sheet_row)
                    = $write_obj->write_content_combine_obj( $sheet, \%option );
            }

            print "Sheet \"$sheet_name\" has been generated.\n";
        }
    };

    my $write_sheet_d2 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        # make combine
        my @combined;
        {
            my $thaw_sql
                = $sql_file->retrieve('common-d1_make_combine_slippage-0');
            $thaw_sql->replace( { distance => 'density' } );
            my $standalone = [ -1, 0 ];
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                threshold  => $combine,
                standalone => $standalone,
                bind_value => [ $feature_1, $feature_2 ],
            );
            @combined = @{ $write_obj->make_combine( \%option ) };
        }

        {
            my $sheet_name = "d2_comb_$name";
            my $sheet;
            my ( $sheet_row, $sheet_col );

            my $thaw_sql = $sql_file->retrieve('common-d1_comb_slippage-0');
            $thaw_sql->replace( { distance => 'density' } );

            {    # write header
                my @headers = $thaw_sql->as_header;
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
                my %option = (
                    sql_obj     => $thaw_sql,
                    sheet_row   => $sheet_row,
                    sheet_col   => $sheet_col,
                    combine_col => 'isw.isw_density',
                    combined    => \@combined,
                    bind_value  => [ $feature_1, $feature_2 ],
                );
                ($sheet_row)
                    = $write_obj->write_content_combine_obj( $sheet, \%option );
            }

            print "Sheet \"$sheet_name\" has been generated.\n";
        }
    };

    for (@type_levels) {
        $write_sheet_d1->($_);
    }

    for (@type_levels) {
        $write_sheet_d2->($_);
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
            my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                header     => \@headers,
                query_name => $sheet_name,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        my @dd_density_group = (
            [ 1,  2 ],
            [ 3,  6 ],
            [ 7,  10 ],
            [ 11, 18 ],
            [ 19, 28 ],
            [ 29, 999 ],
        );

        {    # write contents
            my $thaw_sql = $sql_file->retrieve('common-dd_group-4');
            my %option   = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@dd_density_group,
            );
            ($sheet_row) = $write_obj->write_content_dd( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    #----------------------------------------------------------#
    # worksheet -- dd_group_gc
    #----------------------------------------------------------#
    {
        my $sheet_name = 'dd_group_gc';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{isw_distance AVG_gc COUNT STD_gc};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                header     => \@headers,
                query_name => $sheet_name,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        my @dd_density_group
            = ( [ 1, 2 ], [ 3, 10 ], [ 11, 18 ], [ 19, 999 ], );

        {    # write contents
            my $thaw_sql = $sql_file->retrieve('common-dd_group-4');
            $thaw_sql->replace(
                {   AVG_pi => 'AVG_gc',
                    STD_pi => 'STD_gc',
                    isw_pi => 'isw_average_gc',
                }
            );
            my %option = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@dd_density_group,
            );
            ($sheet_row) = $write_obj->write_content_dd( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- indel_size_group
#----------------------------------------------------------#
my $indel_size_group = sub {
    my $sheet_name = 'indel_size_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group = ( [ 1, 5 ], [ 6, 10 ], [ 11, 50 ], [ 51, 300 ], );

    {    # write contents
        my $thaw_sql_R = $sql_file->retrieve('common-indel_size_r-0');
        $thaw_sql_R->add_where(
            'indel.indel_length' => { op => '>=', value => '1' } );
        $thaw_sql_R->add_where(
            'indel.indel_length' => { op => '<=', value => '5' } );

        my $thaw_sql_L = $sql_file->retrieve('common-indel_size_l-0');
        $thaw_sql_L->add_where(
            'indel.indel_length' => { op => '>=', value => '1' } );
        $thaw_sql_L->add_where(
            'indel.indel_length' => { op => '<=', value => '5' } );

        my %option = (
            sql_query_1 => $thaw_sql_R->as_sql,
            sql_query_2 => $thaw_sql_L->as_sql,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_size_asymmetry
#----------------------------------------------------------#
my $indel_size_asymmetry = sub {
    my $sheet_name = 'indel_size_asymmetry';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group = (
        [ 1,  10,  1,  10 ],
        [ 1,  10,  11, 300 ],
        [ 11, 300, 1,  10 ],
        [ 11, 300, 11, 300 ],
    );

    # write contents
    {
        my $sql_query_L = q{
            # indel_size_asymmetry effect for L windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance,
                   AVG(isw_pi) AVG_pi,
                   COUNT(isw_pi) COUNT,
                   STD(isw_pi) STD_pi
            FROM isw, (SELECT i2.indel_id indel_id
                       FROM indel i1, indel i2
                       WHERE i1.indel_id = i2.prev_indel_id
                       AND i1.indel_length BETWEEN ? AND ?
                       AND i2.indel_length BETWEEN ? AND ?
                      ) indel
            WHERE isw.isw_type = 'L'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.indel_id = indel.indel_id
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
        };
        my $sql_query_R = q{
            # indel_size_asymmetry effect for R windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance,
                   AVG(isw_pi) AVG_pi,
                   COUNT(isw_pi) COUNT,
                   STD(isw_pi) STD_pi
            FROM isw, (SELECT i2.indel_id indel_id
                       FROM indel i1, indel i2
                       WHERE i1.indel_id = i2.prev_indel_id
                       AND i1.indel_length BETWEEN ? AND ?
                       AND i2.indel_length BETWEEN ? AND ?
                      ) indel
            WHERE isw.isw_type = 'R'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.indel_id = indel.indel_id
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
            ORDER BY CONCAT(isw.isw_type, isw.isw_distance) DESC
        };
        my %option = (
            sql_query_1 => $sql_query_L,
            sql_query_2 => $sql_query_R,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_extand_group
#----------------------------------------------------------#
my $indel_extand_group = sub {
    my $sheet_name = 'indel_extand_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group
        = ( [ 0, 0 ], [ 1, 2 ], [ 3, 4 ], [ 5, 9 ], [ 10, 19 ], [ 20, 999 ], );

    {    # write contents
        my $thaw_sql_R = $sql_file->retrieve('common-indel_size_r-0');
        $thaw_sql_R->add_where(
            'FLOOR(indel.right_extand / 100)' => { op => '>=', value => '0' } );
        $thaw_sql_R->add_where(
            'FLOOR(indel.right_extand / 100)' => { op => '<=', value => '0' } );

        my $thaw_sql_L = $sql_file->retrieve('common-indel_size_l-0');
        $thaw_sql_L->add_where(
            'FLOOR(indel.left_extand / 100)' => { op => '>=', value => '0' } );
        $thaw_sql_L->add_where(
            'FLOOR(indel.left_extand / 100)' => { op => '<=', value => '0' } );

        my %option = (
            sql_query_1 => $thaw_sql_R->as_sql,
            sql_query_2 => $thaw_sql_L->as_sql,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_extand_asymmetry
#----------------------------------------------------------#
my $indel_extand_asymmetry = sub {
    my $sheet_name = 'indel_extand_asymmetry';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group = (
        [ 0, 4,   0, 4 ],
        [ 0, 4,   5, 999 ],
        [ 5, 999, 0, 4 ],
        [ 5, 999, 5, 999 ],
    );

    # write contents
    {
        my $sql_query_L = q{
            # indel_extand_asymmetry effect for L windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance,
                   AVG(isw_pi) AVG_pi,
                   COUNT(isw_pi) COUNT,
                   STD(isw_pi) STD_pi
            FROM isw, (SELECT i2.indel_id indel_id
                       FROM indel i1, indel i2
                       WHERE i1.indel_id = i2.prev_indel_id
                       AND FLOOR(i1.left_extand / 100) BETWEEN ? AND ?
                       AND FLOOR(i2.right_extand / 100) BETWEEN ? AND ?
                      ) indel
            WHERE isw.isw_type = 'L'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.indel_id = indel.indel_id
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
        };
        my $sql_query_R = q{
            # indel_extand_asymmetry effect for R windows
            SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_distance, AVG(isw_pi) AVG_pi, COUNT(isw_pi) COUNT, STD(isw_pi) STD_pi
            FROM isw, (SELECT i2.indel_id indel_id
                       FROM indel i1, indel i2
                       WHERE i1.indel_id = i2.prev_indel_id
                       AND FLOOR(i1.left_extand / 100) BETWEEN ? AND ?
                       AND FLOOR(i2.right_extand / 100) BETWEEN ? AND ?
                      ) indel
            WHERE isw.isw_type = 'R'
            AND isw.isw_density > 9
            AND isw.isw_distance <= 5
            AND isw.indel_id = indel.indel_id
            GROUP BY CONCAT(isw.isw_type, isw.isw_distance)
            ORDER BY CONCAT(isw.isw_type, isw.isw_distance) DESC
        };
        my %option = (
            sql_query_1 => $sql_query_L,
            sql_query_2 => $sql_query_R,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_position_group
#----------------------------------------------------------#
my $indel_position_group = sub {
    my $sheet_name = 'indel_position_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group
        = ( [ 1, 1, 0, 0 ], [ 1, 1, 1, 1 ], [ 0, 0, 0, 0 ], [ 0, 0, 1, 1 ], );

    {    # write contents
        my $thaw_sql_R = $sql_file->retrieve('common-indel_feature_r-0');
        $thaw_sql_R->add_where(
            'indel.indel_coding' => { op => '>=', value => '0' } );
        $thaw_sql_R->add_where(
            'indel.indel_coding' => { op => '<=', value => '0' } );
        $thaw_sql_R->add_where(
            'indel.indel_repeats' => { op => '>=', value => '0' } );
        $thaw_sql_R->add_where(
            'indel.indel_repeats' => { op => '<=', value => '0' } );

        my $thaw_sql_L = $sql_file->retrieve('common-indel_feature_l-0');
        $thaw_sql_L->add_where(
            'indel.indel_coding' => { op => '>=', value => '0' } );
        $thaw_sql_L->add_where(
            'indel.indel_coding' => { op => '<=', value => '0' } );
        $thaw_sql_L->add_where(
            'indel.indel_repeats' => { op => '>=', value => '0' } );
        $thaw_sql_L->add_where(
            'indel.indel_repeats' => { op => '<=', value => '0' } );

        my %option = (
            sql_query_1 => $thaw_sql_R->as_sql,
            sql_query_2 => $thaw_sql_L->as_sql,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_coding_group
#----------------------------------------------------------#
my $indel_coding_group = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel', 'indel_coding' ) ) {
        return;
    }

    my $sheet_name = 'indel_coding_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group
        = ( [ 0, 0 ], [ 1, 1 ], );

    {    # write contents
        my $thaw_sql_R = $sql_file->retrieve('common-indel_feature_r-0');
        $thaw_sql_R->add_where(
            'indel.indel_coding' => { op => '>=', value => '0' } );
        $thaw_sql_R->add_where(
            'indel.indel_coding' => { op => '<=', value => '0' } );

        my $thaw_sql_L = $sql_file->retrieve('common-indel_feature_l-0');
        $thaw_sql_L->add_where(
            'indel.indel_coding' => { op => '>=', value => '0' } );
        $thaw_sql_L->add_where(
            'indel.indel_coding' => { op => '<=', value => '0' } );

        my %option = (
            sql_query_1 => $thaw_sql_R->as_sql,
            sql_query_2 => $thaw_sql_L->as_sql,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_repeat_group
#----------------------------------------------------------#
my $indel_repeat_group = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel', 'indel_repeats' ) ) {
        return;
    }

    my $sheet_name = 'indel_repeat_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group
        = ( [ 0, 0 ], [ 1, 1 ], );

    {    # write contents
        my $thaw_sql_R = $sql_file->retrieve('common-indel_feature_r-0');
        $thaw_sql_R->add_where(
            'indel.indel_repeats' => { op => '>=', value => '0' } );
        $thaw_sql_R->add_where(
            'indel.indel_repeats' => { op => '<=', value => '0' } );

        my $thaw_sql_L = $sql_file->retrieve('common-indel_feature_l-0');
        $thaw_sql_L->add_where(
            'indel.indel_repeats' => { op => '>=', value => '0' } );
        $thaw_sql_L->add_where(
            'indel.indel_repeats' => { op => '<=', value => '0' } );

        my %option = (
            sql_query_1 => $thaw_sql_R->as_sql,
            sql_query_2 => $thaw_sql_L->as_sql,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_slip_group
#----------------------------------------------------------#
my $indel_slip_group = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel', 'indel_slippage' ) ) {
        return;
    }

    my $sheet_name = 'indel_slip_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group
        = ( [ 0, 0 ], [ 1, 1 ], );

    {    # write contents
        my $thaw_sql_R = $sql_file->retrieve('common-indel_feature_r-0');
        $thaw_sql_R->add_where(
            'indel.indel_slippage' => { op => '>=', value => '0' } );
        $thaw_sql_R->add_where(
            'indel.indel_slippage' => { op => '<=', value => '0' } );

        my $thaw_sql_L = $sql_file->retrieve('common-indel_feature_l-0');
        $thaw_sql_L->add_where(
            'indel.indel_slippage' => { op => '>=', value => '0' } );
        $thaw_sql_L->add_where(
            'indel.indel_slippage' => { op => '<=', value => '0' } );

        my %option = (
            sql_query_1 => $thaw_sql_R->as_sql,
            sql_query_2 => $thaw_sql_L->as_sql,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_gc_group
#----------------------------------------------------------#
my $indel_gc_group = sub {
    my $sheet_name = 'indel_gc_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{isw_distance AVG_pi COUNT STD_pi};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    my @indel_group = ( [ 0, 0.2999 ], [ 0.3, 0.4999 ], [ 0.5, 1 ], );

    {    # write contents
        my $thaw_sql_R = $sql_file->retrieve('common-indel_size_r-0');
        $thaw_sql_R->add_where( 'indel.indel_length' => \'>= 10' );
        $thaw_sql_R->add_where(
            'indel.indel_gc' => { op => '>=', value => '0' } );
        $thaw_sql_R->add_where(
            'indel.indel_gc' => { op => '<=', value => '0' } );

        my $thaw_sql_L = $sql_file->retrieve('common-indel_size_l-0');
        $thaw_sql_L->add_where( 'indel.indel_length' => \'>= 10' );
        $thaw_sql_L->add_where(
            'indel.indel_gc' => { op => '>=', value => '0' } );
        $thaw_sql_L->add_where(
            'indel.indel_gc' => { op => '<=', value => '0' } );

        my %option = (
            sql_query_1 => $thaw_sql_R->as_sql,
            sql_query_2 => $thaw_sql_L->as_sql,
            sheet_row   => $sheet_row,
            sheet_col   => $sheet_col,
            group       => \@indel_group,
        );
        ($sheet_row) = $write_obj->write_content_indel( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_indel_ratio
#----------------------------------------------------------#
my $snp_indel_ratio = sub {
    my $sheet_name = 'snp_indel_ratio';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # create temporary table
    {
        my $sql_query = q{
            DROP TABLE IF EXISTS pi_group
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    {
        my $sql_query = q{
            # create temporary table
            CREATE TABLE pi_group (p_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (p_id))
                ENGINE=MyISAM
                SELECT a.align_pi `pi`,
                       a.align_differences `snp`,
                       i.indel_number `indel`,
                       a.align_gaps `gaps`,
                       a.align_length `align_length`
                FROM (SELECT align_id, COUNT(indel_id) indel_number
                      FROM indel i
                      GROUP BY align_id) i,
                     align a
                WHERE i.align_id = a.align_id
                ORDER BY pi DESC
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    # make group
    my @group_align;
    {
        my $sql_query = q{
            SELECT p_id, align_length
            FROM pi_group
        };
        my %option = (
            sql_query => $sql_query,
            piece     => $piece,
        );
        @group_align = @{ $write_obj->make_combine_piece( \%option ) };
    }

    {    # write header
        my @headers
            = qw{AVG_pi AVG_SNP/Indel COUNT AVG_align_length SUM_align_length AVG_SNP/kb AVG_Indel/kb};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # align query
        my $sql_query = q{
            SELECT AVG(p.pi) `AVG_pi`,
                   AVG(p.snp / p.indel) `AVG_SNP/Indel`,
                   COUNT(*) COUNT,
                   AVG(p.align_length) `AVG_align_length`,
                   SUM(p.align_length) `SUM_align_length`,
                   AVG(p.snp / p.align_length * 1000) `AVG_SNP/kb`,
                   AVG(p.indel / p.align_length * 1000) `AVG_Indel/kb`
            FROM pi_group p
            WHERE p_id IN
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            group     => \@group_align,
        );

        ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
    }

    {    # drop temporary table
        my $sql_query = q{
            DROP TABLE IF EXISTS pi_group
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_length
#----------------------------------------------------------#
my $indel_length = sub {
    my $sheet_name = 'indel_length';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{indel_length indel_number AVG_gc_ratio indel_sum};
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
        my $thaw_sql = $sql_file->retrieve('common-indel_length-0');

        my %option = (
            sql_query => $thaw_sql->as_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_highlight( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_length_100
#----------------------------------------------------------#
my $indel_length_100 = sub {
    my $sheet_name = 'indel_length_100';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{indel_length indel_number AVG_gc_ratio indel_sum};
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
        my $thaw_sql = $sql_file->retrieve('common-indel_length-0');
        $thaw_sql->add_where( 'left_extand'  => \'>= 100' );
        $thaw_sql->add_where( 'right_extand' => \'>= 100' );

        my %option = (
            sql_query => $thaw_sql->as_sql,
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
my $snp_base_change = sub {
    my $sheet_name = 'snp_base_change';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{base_change snp_number};
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $sheet_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        $sheet_row++;
        my $query_name = 'T2Q';
        my $sql_query  = q{
        SELECT s.base_change, s.snp_number / total.total * 100
        FROM (SELECT count(*) total
             FROM snp) total,
             (SELECT CONCAT(target_base, "->", query_base) base_change,
                     COUNT(snp_id) snp_number
             FROM snp
             GROUP BY CONCAT(target_base,  "->",query_base) ) s
        };
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
        my $query_name = 'Q2T';
        my $sql_query  = q{
        SELECT s.base_change, s.snp_number / total.total * 100
        FROM (SELECT count(*) total
             FROM snp) total,
             (SELECT CONCAT(query_base, "->",target_base ) base_change,
                     COUNT(snp_id) snp_number
             FROM snp
             GROUP BY CONCAT(query_base, "->",target_base ) ) s
        };
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
    my @base_pair = qw/A<=>C A<=>G A<=>T C<=>G C<=>T G<=>T/;

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ 'distance', @base_pair ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q{
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND i.isw_distance BETWEEN -1 AND 30
            GROUP BY i.isw_distance
        };
        my $sql_query2 = q{
            # base change
            SELECT i.isw_distance distance, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND CONCAT(target_base, query_base) IN (?, ?)
            AND i.isw_distance BETWEEN -1 AND 30
            GROUP BY i.isw_distance
        };
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- density_snp
#----------------------------------------------------------#
my $density_snp = sub {
    my $sheet_name = 'density_snp';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # Six base pair groups
    my @base_pair = qw/A<=>C A<=>G A<=>T C<=>G C<=>T G<=>T/;

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [ 'density', @base_pair ],
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query1 = q{
            # base change
            SELECT i.isw_density density, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND i.isw_density BETWEEN -1 AND 60
            GROUP BY i.isw_density
        };
        my $sql_query2 = q{
            # base change
            SELECT i.isw_density density, COUNT(s.snp_id) snp_number
            FROM snp s, isw i
            WHERE s.isw_id = i.isw_id
            AND CONCAT(target_base, query_base) IN (?, ?)
            AND i.isw_density BETWEEN -1 AND 60
            GROUP BY i.isw_density
        };
        my %option = (
            sql_query1 => $sql_query1,
            sql_query2 => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            base_pair  => \@base_pair,
        );
        ($sheet_row) = $write_obj->write_content_snp( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- align_coding
#----------------------------------------------------------#
my $align_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'align', 'align_coding' ) ) {
        return;
    }

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT align_coding
            FROM align
            WHERE align_coding IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @coding_levels = (
        [ 1, $quartiles->[0], $quartiles->[1] ],
        [ 2, $quartiles->[1], $quartiles->[2] ],
        [ 3, $quartiles->[2], $quartiles->[3] ],
        [ 4, $quartiles->[3], $quartiles->[4] ],
        [ 9, 0.4,             0.6 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;

        my $sheet_name = "align_coding_$order";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = ( qw{distance AVG_pi COUNT STD_pi}, $low_border,
                $high_border );
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
            my $thaw_sql = $sql_file->retrieve('common-align-0');
            $thaw_sql->add_where(
                'align.align_coding' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'align.align_coding' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $low_border, $high_border ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@coding_levels) {
        &$write_sheet(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- align_repeat
#----------------------------------------------------------#
my $align_repeat = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'align', 'align_repeats' ) ) {
        return;
    }

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT align_repeats
            FROM align
            WHERE align_repeats IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @repeat_levels = (
        [ 1, $quartiles->[0], $quartiles->[1] ],
        [ 2, $quartiles->[1], $quartiles->[2] ],
        [ 3, $quartiles->[2], $quartiles->[3] ],
        [ 4, $quartiles->[3], $quartiles->[4] ],
        [ 9, 0.4,             0.6 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;

        my $sheet_name = "align_repeat_$order";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = ( qw{distance AVG_pi COUNT STD_pi}, $low_border,
                $high_border );
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
            my $thaw_sql = $sql_file->retrieve('common-align-0');
            $thaw_sql->add_where(
                'align.align_repeats' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'align.align_repeats' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $low_border, $high_border ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@repeat_levels) {
        &$write_sheet(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- align_te
#----------------------------------------------------------#
my $align_te = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'align', 'align_te' ) ) {
        return;
    }

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT align_te
            FROM align
            WHERE align_te IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @te_levels = (
        [ 1, $quartiles->[0], $quartiles->[1] ],
        [ 2, $quartiles->[1], $quartiles->[2] ],
        [ 3, $quartiles->[2], $quartiles->[3] ],
        [ 4, $quartiles->[3], $quartiles->[4] ],
        [ 5, 0,               0 ],
        [ 6, 0.0001,          0.1 ],
        [ 7, 0.1,             0.9 ],
        [ 8, 0.9,             1 ],
        [ 9, 0.4,             0.6 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;

        my $sheet_name = "align_te_$order";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = ( qw{distance AVG_pi COUNT STD_pi}, $low_border,
                $high_border );
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
            my $thaw_sql = $sql_file->retrieve('common-align-0');
            $thaw_sql->add_where(
                'align.align_te' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'align.align_te' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $low_border, $high_border ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@te_levels) {
        &$write_sheet(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- align_paralog
#----------------------------------------------------------#
my $align_paralog = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'align', 'align_paralog' ) ) {
        return;
    }

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT align_paralog
            FROM align
            WHERE align_paralog IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @para_levels = (
        [ 1, $quartiles->[0], $quartiles->[1] ],
        [ 2, $quartiles->[1], $quartiles->[2] ],
        [ 3, $quartiles->[2], $quartiles->[3] ],
        [ 4, $quartiles->[3], $quartiles->[4] ],
        [ 5, 0,               0 ],
        [ 6, 0.0001,          0.1 ],
        [ 7, 0.1,             0.9 ],
        [ 8, 0.9,             1 ],
        [ 9, 0.4,             0.6 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;

        my $sheet_name = "align_paralog_$order";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = ( qw{distance AVG_pi COUNT STD_pi}, $low_border,
                $high_border );
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
            my $thaw_sql = $sql_file->retrieve('common-align-0');
            $thaw_sql->add_where(
                'align.align_paralog' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'align.align_paralog' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $low_border, $high_border ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@para_levels) {
        &$write_sheet(@$_);
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$basic; &$process; &$summary; next; }
    if ( $n == 2 )  { &$pi_gc_cv;             next; }
    if ( $n == 3 )  { &$comb_pi_gc_cv;        next; }
    if ( $n == 4 )  { &$group_distance;       &$group_density; next; }
    if ( $n == 5 )  { &$comb_coding;          next; }
    if ( $n == 6 )  { &$comb_slippage;        next; }
    if ( $n == 8 )  { &$dd_group;             next; }
    if ( $n == 9 )  { &$indel_size_group;     &$indel_size_asymmetry; next; }
    if ( $n == 10 ) { &$indel_extand_group;   &$indel_extand_asymmetry; next; }
    if ( $n == 11 ) { &$indel_position_group; next; }
    if ( $n == 12 ) { &$indel_coding_group;   &$indel_repeat_group; next; }
    if ( $n == 13 ) { &$indel_slip_group;     &$indel_gc_group; next; }
    if ( $n == 14 ) { &$snp_indel_ratio;      next; }
    if ( $n == 15 ) { &$indel_length;         &$indel_length_100; next; }
    if ( $n == 16 ) { &$snp_base_change;      next; }
    if ( $n == 17 ) { &$distance_snp;         &$density_snp; next; }

    if ( $n == 51 ) { &$align_coding;  next; }
    if ( $n == 52 ) { &$align_repeat;  next; }
    if ( $n == 53 ) { &$align_te;      next; }
    if ( $n == 54 ) { &$align_paralog; next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    common_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    common_stat_factory.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --output            output filename
        --run               run special analysis

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back
