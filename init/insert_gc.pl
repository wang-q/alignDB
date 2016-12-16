#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;

use AlignDB::GC;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use AlignDB::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
Add GC ralated tables to alignDB

    perl init/insert_gc.pl -d S288cvsRM11_1a --parallel 2

Usage: perl %c [options]
EOF

(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    $description,
    [ 'help|h', 'display this message' ],
    [],
    ['Database init values'],
    [ 'server|s=s',   'MySQL IP/Domain', { default => $conf->{database}{server} }, ],
    [ 'port=i',       'MySQL port',      { default => $conf->{database}{port} }, ],
    [ 'username|u=s', 'username',        { default => $conf->{database}{username} }, ],
    [ 'password|p=s', 'password',        { default => $conf->{database}{password} }, ],
    [ 'db|d=s',       'database name',   { default => $conf->{database}{db} }, ],
    [],
    ['AlignDB::GC init values'],
    [ 'insert_gc=s',      'insert gc or not',      { default => $conf->{gc}{insert_gc} }, ],
    [ 'insert_segment=s', 'insert segment or not', { default => $conf->{gc}{insert_segment} }, ],
    [ 'alt_levelalt_level', 'use alternative segment levels 200..900, 1000..5000', ],
    [ 'one_level',          'use segment level of 100', ],
    [],
    [ 'parallel=i', 'run in parallel mode',       { default => $conf->{generate}{parallel} }, ],
    [ 'batch=i',    '#alignments in one process', { default => $conf->{generate}{batch} }, ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# record config
$stopwatch->record_conf( { opt => $opt, gc => $conf->{gc}, } );

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update GC tables of [$opt->{db}]...");

my @jobs;
{
    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    # empty tables: segment, gsw, extreme
    print "Emptying tables...\n";
    $alignDB->empty_table( 'segment', 'with_window' );
    $alignDB->empty_table( 'gsw',     'with_window' );
    $alignDB->empty_table( 'extreme', 'with_window' );

    @jobs = @{ $alignDB->get_align_ids };
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my @align_ids = @{$chunk_ref};
    my $wid       = MCE->wid;

    $stopwatch->block_message("Process task [$chunk_id] by worker #$wid");

    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );
    my $obj_gc = AlignDB::GC->new(
        wave_window_size => $conf->{gc}{wave_window_size},
        wave_window_step => $conf->{gc}{wave_window_step},
        vicinal_size     => $conf->{gc}{vicinal_size},
        fall_range       => $conf->{gc}{fall_range},
        gsw_size         => $conf->{gc}{gsw_size},
        stat_window_size => $conf->{gc}{stat_window_size},
        stat_window_step => $conf->{gc}{stat_window_step},
    );

    # Database handler
    my DBI $dbh = $alignDB->dbh;

    # alignments' chromosomal location, target_seq and query_seq
    my DBI $sth = $dbh->prepare(
        q{
        SELECT  a.align_comparable_runlist
        FROM    align a
        WHERE   a.align_id = ?
        }
    );

    # for each alignment
    for my $align_id (@align_ids) {
        $alignDB->process_message($align_id);
        $sth->execute($align_id);
        my ($comparable_runlist) = $sth->fetchrow_array;

        # comparable runlist
        my $comparable_set = AlignDB::IntSpan->new($comparable_runlist);

        if ( $opt->{insert_gc} ) {
            insert_extreme( $obj_gc, $alignDB, $align_id, $comparable_set );

            insert_gsw( $obj_gc, $alignDB, $align_id, $comparable_set );
        }

        if ( $opt->{insert_segment} ) {
            my $style = 'normal';
            if ( $opt->{alt_level} ) {
                $style = 'alt_level';
            }
            elsif ( $opt->{one_level} ) {
                $style = 'one_level';
            }
            insert_segment( $obj_gc, $alignDB, $align_id, $comparable_set, $style );
        }
    }

    return;
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $opt->{parallel}, chunk_size => $opt->{batch}, );
$mce->forchunk( \@jobs, $worker, );

$stopwatch->end_message;

# store program's meta info to database
AlignDB::Common->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

exit;

sub insert_extreme {
    my AlignDB::GC $obj_gc              = shift;
    my AlignDB::Common $obj             = shift;
    my $align_id                        = shift;
    my AlignDB::IntSpan $comparable_set = shift;

    my DBI $dbh = $obj->dbh;
    my $seqs_ref = $obj->get_seqs($align_id);
    $seqs_ref = [ $seqs_ref->[0] ];    # only calc target gc

    my @extreme_site;
    my @slidings = $obj_gc->gc_wave( $seqs_ref, $comparable_set );
    for my $s (@slidings) {
        my $flag = $s->{high_low_flag};
        if ( $flag eq 'T' or $flag eq 'C' ) {
            push @extreme_site, $s;
        }
    }

    # get extreme sliding windows' sizes
    my $windows_size = $obj_gc->wave_window_size;
    my $half_length  = int( $windows_size / 2 );

    # left
    my $prev_extreme_middle_right_idx = 1;
    my $prev_extreme_gc               = $slidings[0]->{gc};
    for ( my $i = 0; $i < scalar @extreme_site; $i++ ) {

        # wave_length
        my AlignDB::IntSpan $extreme_set = $extreme_site[$i]->{set};
        my $extreme_middle_left          = $extreme_set->at($half_length);
        my $extreme_middle_left_idx      = $comparable_set->index($extreme_middle_left);
        my $left_wave_length = $extreme_middle_left_idx - $prev_extreme_middle_right_idx + 1;
        $prev_extreme_middle_right_idx = $extreme_middle_left_idx + 1;
        $extreme_site[$i]->{left_wave_length} = $left_wave_length;

        # amplitude
        my $extreme_gc     = $extreme_site[$i]->{gc};
        my $left_amplitude = abs( $extreme_gc - $prev_extreme_gc );
        $extreme_site[$i]->{left_amplitude} = $left_amplitude;
        $prev_extreme_gc = $extreme_gc;
    }

    # right
    my $next_extreme_middle_left_idx = $comparable_set->size;
    my $next_extreme_gc              = $slidings[-1]->{gc};
    for ( my $i = scalar @extreme_site - 1; $i >= 0; $i-- ) {

        # wave_length
        my AlignDB::IntSpan $extreme_set = $extreme_site[$i]->{set};
        my $extreme_middle_right         = $extreme_set->at( $half_length + 1 );
        my $extreme_middle_right_idx     = $comparable_set->index($extreme_middle_right);
        my $right_wave_length = $next_extreme_middle_left_idx - $extreme_middle_right_idx + 1;
        $next_extreme_middle_left_idx = $extreme_middle_right_idx - 1;
        $extreme_site[$i]->{right_wave_length} = $right_wave_length;

        # amplitude
        my $extreme_gc      = $extreme_site[$i]->{gc};
        my $right_amplitude = abs( $extreme_gc - $next_extreme_gc );
        $extreme_site[$i]->{right_amplitude} = $right_amplitude;
        $prev_extreme_gc = $extreme_gc;
    }

    # prepare extreme_insert
    my $extreme_sql = qq{
        INSERT INTO extreme (
            extreme_id, prev_extreme_id, window_id, extreme_type,
            extreme_left_amplitude, extreme_right_amplitude,
            extreme_left_wave_length, extreme_right_wave_length
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?,
            ?, ?            
        )
    };

    my DBI $sth = $dbh->prepare($extreme_sql);

    my $prev_extreme_id = 0;
    for (@extreme_site) {
        my $extreme_set = $_->{set};

        my ($cur_window_id) = $obj->insert_window( $align_id, $extreme_set );

        $sth->execute(
            $prev_extreme_id,     $cur_window_id,        $_->{high_low_flag},
            $_->{left_amplitude}, $_->{right_amplitude}, $_->{left_wave_length},
            $_->{right_wave_length}
        );

        $prev_extreme_id = $obj->last_insert_id;
    }
    $sth->finish;

    return;
}

sub insert_gsw {
    my AlignDB::GC $obj_gc              = shift;
    my AlignDB::Common $obj             = shift;
    my $align_id                        = shift;
    my AlignDB::IntSpan $comparable_set = shift;

    my DBI $dbh = $obj->dbh;

    # get gc sliding windows' sizes
    my $gsw_size = $obj_gc->gsw_size;

    # extreme_id & prev_extreme_id
    my DBI $fetch_ex_id = $dbh->prepare(
        q{
        SELECT e.extreme_id,
                e.prev_extreme_id,
                e.extreme_type
        FROM extreme e, window w
        WHERE e.window_id = w.window_id
        AND w.align_id = ?
        }
    );
    $fetch_ex_id->execute($align_id);

    # extreme_runlist, _amplitude and _gc
    my DBI $fetch_ex_attr = $dbh->prepare(
        q{
        SELECT w.window_runlist,
                w.window_gc,
                e.extreme_left_amplitude,
                e.extreme_left_wave_length
        FROM extreme e, window w
        WHERE e.window_id = w.window_id
        AND e.extreme_id = ?
        }
    );

    # prepare gsw_insert
    my DBI $gsw_insert = $dbh->prepare(
        q{
        INSERT INTO gsw (
            gsw_id, extreme_id, prev_extreme_id, window_id,
            gsw_type, gsw_distance, gsw_distance_crest, gsw_wave_length,
            gsw_amplitude, gsw_trough_gc, gsw_crest_gc, gsw_gradient
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?, ?, ?
        )
        }
    );

    while ( my $ref = $fetch_ex_id->fetchrow_hashref ) {
        my $ex_id      = $ref->{extreme_id};
        my $prev_ex_id = $ref->{prev_extreme_id};
        my $ex_type    = $ref->{extreme_type};

        # bypass the first extreme
        if ( $prev_ex_id == 0 ) {
            next;
        }

        # get attrs of the two extremes
        $fetch_ex_attr->execute($prev_ex_id);
        my ( $prev_ex_runlist, $prev_ex_gc, undef, undef ) = $fetch_ex_attr->fetchrow_array;
        my $prev_ex_set = AlignDB::IntSpan->new($prev_ex_runlist);

        $fetch_ex_attr->execute($ex_id);
        my ( $ex_runlist, $ex_gc, $ex_left_amplitude, $ex_left_wave_length )
            = $fetch_ex_attr->fetchrow_array;
        my $ex_set = AlignDB::IntSpan->new($ex_runlist);

        # determining $gsw_density here, which is different from isw_density
        my $half_length = int( $ex_set->size / 2 );
        my $gsw_density = int( ( $ex_left_wave_length - $half_length ) / $gsw_size );

        # wave length, amplitude, trough_gc and gradient
        my $gsw_wave_length = $ex_left_wave_length;
        my $gsw_amplitude   = $ex_left_amplitude;
        my $gsw_trough_gc   = $ex_type eq 'T' ? $ex_gc : $prev_ex_gc;
        my $gsw_crest_gc    = $ex_type eq 'T' ? $prev_ex_gc : $ex_gc;
        my $gsw_gradient    = int( $gsw_amplitude / $ex_left_wave_length / 0.00001 );

        # determining $gsw_type here, ascend and descent
        my $gsw_type;
        my @gsw_windows;
        if ( $ex_type eq 'T' ) {    # push trough to gsw
            $gsw_type = 'D';        # descend, left of trough, right of crest
            push @gsw_windows,
                {
                type           => 'T',
                set            => $ex_set,
                distance       => 0,
                distance_crest => $gsw_density + 1,
                };
        }
        elsif ( $ex_type eq 'C' ) {    # crest
            $gsw_type = 'A';           # ascend, right of trough, left of crest
            push @gsw_windows,
                {
                type           => 'C',
                set            => $ex_set,
                distance       => $gsw_density + 1,
                distance_crest => 0,
                };
        }
        else {
            warn "extreme_type [$ex_type] error\n";
        }

        {    # More windows will be submitted in the following section

            # distance start counting from trough
            # $sw_start and $sw_end are both index of $comprarable_set
            # $gsw_distance is from 1 to $gsw_density
            # window 0 is trough
            # ..., D2, D1, T0, A1, A2, ...
            my ( $sw_start, $sw_end );
            if ( $gsw_type eq 'A' ) {
                $sw_start = $comparable_set->index( $prev_ex_set->max ) + 1;
                $sw_end   = $sw_start + $gsw_size - 1;
            }
            elsif ( $gsw_type eq 'D' ) {
                $sw_end   = $comparable_set->index( $ex_set->min ) - 1;
                $sw_start = $sw_end - $gsw_size + 1;
            }

            for my $gsw_distance ( 1 .. $gsw_density ) {
                my $gsw_set = $comparable_set->slice( $sw_start, $sw_end );

                push @gsw_windows,
                    {
                    type           => $gsw_type,
                    set            => $gsw_set,
                    distance       => $gsw_distance,
                    distance_crest => $gsw_density - $gsw_distance + 1,
                    };

                if ( $gsw_type eq 'A' ) {
                    $sw_start = $sw_end + 1;
                    $sw_end   = $sw_start + $gsw_size - 1;
                }
                elsif ( $gsw_type eq 'D' ) {
                    $sw_end   = $sw_start - 1;
                    $sw_start = $sw_end - $gsw_size + 1;
                }
            }

            for my $gsw (@gsw_windows) {
                my ($cur_window_id)
                    = $obj->insert_window( $align_id, $gsw->{set} );

                $gsw_insert->execute(
                    $ex_id,           $prev_ex_id,      $cur_window_id,
                    $gsw->{type},     $gsw->{distance}, $gsw->{distance_crest},
                    $gsw_wave_length, $gsw_amplitude,   $gsw_trough_gc,
                    $gsw_crest_gc,    $gsw_gradient,
                );
            }
        }
    }

    return;
}

sub insert_segment {
    my AlignDB::GC $obj_gc              = shift;
    my AlignDB::Common $obj             = shift;
    my $align_id                        = shift;
    my AlignDB::IntSpan $comparable_set = shift;

    my $style = shift || 'normal';

    my DBI $dbh = $obj->dbh;

    my @segment_levels;
    if ( $style eq "alt_level" ) {
        @segment_levels = reverse map { [ $_, $_ * 100, $_ * 100 ] } ( 2 .. 10, 20, 30, 40, 50 );
    }
    elsif ( $style eq "one_level" ) {
        @segment_levels = ( [ 100, 100, 100 ] );
    }
    else {
        @segment_levels
            = ( [ 1, 5000, 5000 ], [ 2, 1000, 1000 ], [ 3, 500, 500 ], );
    }

    for (@segment_levels) {
        my $segment_type = $_->[0];
        my $segment_size = $_->[1];
        my $segment_step = $_->[2];

        my @segment_site = $obj_gc->segment( $comparable_set, $segment_size, $segment_step );

        # prepare segment_insert
        my DBI $sth = $dbh->prepare(
            q{
            INSERT INTO segment (
                segment_id, window_id, segment_type,
                segment_gc_mean, segment_gc_std,
                segment_gc_cv, segment_gc_mdcw
            )
            VALUES (
                NULL, ?, ?,
                ?, ?,
                ?, ?
            )
            }
        );

        for my $segment_set (@segment_site) {

            my ($cur_window_id) = $obj->insert_window( $align_id, $segment_set );

            my $seqs_ref = $obj->get_seqs($align_id);
            $seqs_ref = [ $seqs_ref->[0] ];    # only calc target gc

            if ( $style ne "one_level" ) {
                my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                    = $obj_gc->segment_gc_stat( $seqs_ref, $segment_set );

                $sth->execute( $cur_window_id, $segment_type, $gc_mean, $gc_std, $gc_cv, $gc_mdcw );
            }
            else {
                my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                    = $obj_gc->segment_gc_stat_one( $seqs_ref, $comparable_set, $segment_set );

                $sth->execute( $cur_window_id, $segment_type, $gc_mean, $gc_std, $gc_cv, $gc_mdcw );
            }
        }
        $sth->finish;
    }

    return;
}

__END__
