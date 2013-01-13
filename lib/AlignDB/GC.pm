package AlignDB::GC;
use Moose::Role;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Util qw(:all);

requires 'dbh';

# extreme sliding window size
has 'wave_window_size' => ( is => 'rw', isa => 'Int', default => 100, );

# extreme sliding window step
has 'wave_window_step' => ( is => 'rw', isa => 'Int', default => 50, );

# vicinal size
has 'vicinal_size' => ( is => 'rw', isa => 'Int', default => 500, );

# minimal fall range
has 'fall_range' => ( is => 'rw', isa => 'Num', default => 0.1, );

# gsw window size
has 'gsw_size' => ( is => 'rw', isa => 'Int', default => 100, );

# extreme sliding window size
has 'stat_window_size' => ( is => 'rw', isa => 'Int', default => 100, );

# extreme sliding window step
has 'stat_window_step' => ( is => 'rw', isa => 'Int', default => 100, );

# use 100 .. 900 segment levels
has 'alt_level' => ( is => 'rw', isa => 'Bool', default => 0, );

sub insert_segment {
    my $self           = shift;
    my $align_id       = shift;
    my $comparable_set = shift;

    my $dbh = $self->dbh;

    my @segment_levels = (
        [ 'A', '',   '' ],
        [ 1,   5000, 5000 ],
        [ 2,   1000, 1000 ],
        [ 3,   500,  500 ],
    );

    if ( $self->alt_level ) {
        @segment_levels = reverse map { [ $_, $_ * 100, $_ * 100 ] }
            ( 2 .. 10, 20, 30, 40, 50 );
    }

    for (@segment_levels) {

        my $segment_type = $_->[0];
        my $segment_size = $_->[1];
        my $segment_step = $_->[2];

        my @segment_site
            = $self->segment( $comparable_set, $segment_size, $segment_step );

        # prepare segment_insert
        my $segment_sql = qq{
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
        };

        my $segment_insert = $dbh->prepare($segment_sql);

        foreach my $segment_set (@segment_site) {

            my ($cur_window_id)
                = $self->insert_window( $align_id, $segment_set );

            my $seqs_ref = $self->get_seqs($align_id);
            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $self->segment_gc_stat( $seqs_ref, $segment_set );

            $segment_insert->execute( $cur_window_id, $segment_type, $gc_mean,
                $gc_std, $gc_cv, $gc_mdcw );
        }
        $segment_insert->finish;
    }

    return;
}

sub segment {
    my $self           = shift;
    my $comparable_set = shift;
    my $segment_size   = shift;
    my $segment_step   = shift;

    # if not given $segment_size, do a whole alignment analysis
    my @segment_windows;
    if ($segment_size) {
        @segment_windows
            = $self->sliding_window( $comparable_set, $segment_size,
            $segment_step );
    }
    else {
        push @segment_windows, $comparable_set;
    }

    return @segment_windows;
}

sub sliding_window {
    my $self           = shift;
    my $comparable_set = shift;
    my $window_size    = shift || $self->wave_window_size;
    my $window_step    = shift || $self->wave_window_step;

    my $comparable_number = $comparable_set->cardinality;

    my @sliding_windows;
    my $sliding_start = 1;
    while (1) {
        my $sliding_end = $sliding_start + $window_size - 1;
        last if $sliding_end > $comparable_number;
        my $sliding_window
            = $comparable_set->slice( $sliding_start, $sliding_end );
        $sliding_start += $window_step;

        push @sliding_windows, $sliding_window;
    }

    return @sliding_windows;
}

sub segment_gc_stat {
    my $self        = shift;
    my $seqs_ref    = shift;
    my $segment_set = shift;
    my $window_size = shift || $self->stat_window_size;
    my $window_step = shift || $self->stat_window_step;

    my $segment_number = $segment_set->cardinality;

    # There should be at least 2 sample, i.e., two sliding windows
    #   in the segment
    if ( $segment_number < $window_size + $window_step ) {
        return ( undef, undef, undef, undef );
    }

    my @seqs = @{$seqs_ref};

    my @sliding_windows
        = $self->sliding_window( $segment_set, $window_size, $window_step );
    my @sliding_gcs;

    foreach my $sliding_set (@sliding_windows) {
        my @sliding_seqs = map { $sliding_set->substr_span($_) } @seqs;
        my $sliding_gc = calc_gc_ratio(@sliding_seqs);
        push @sliding_gcs, $sliding_gc;
    }

    my $gc_mean = mean(@sliding_gcs);
    my $gc_std  = stddev(@sliding_gcs);
    my $gc_cv
        = $gc_mean == 0 || $gc_mean == 1 ? undef
        : $gc_mean <= 0.5 ? $gc_std / $gc_mean
        :                   $gc_std / ( 1 - $gc_mean );
    my $gc_mdcw = _mdcw(@sliding_gcs);

    return ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw );
}

sub _mdcw {
    my @array = @_;

    if ( @array <= 1 ) {
        warn "There should be more than 1 elements in the array.\n";
        return;
    }

    my @dcws;
    for ( 1 .. $#array ) {
        my $dcw = $array[$_] - $array[ $_ - 1 ];
        push @dcws, abs $dcw;
    }

    return mean(@dcws);
}

# find local maxima, minima
sub find_extreme_step1 {
    my $self    = shift;
    my $windows = shift;

    my $count = scalar @$windows;

    my $wave_window_size = $self->wave_window_size;
    my $wave_window_step = $self->wave_window_step;
    my $vicinal_size     = $self->vicinal_size;
    my $vicinal_number   = int( $vicinal_size / $wave_window_step );

    for my $i ( 0 .. $count - 1 ) {

        #----------------------------#
        # right vicinal window
        #----------------------------#
        my @right_vicinal_windows;
        my ( $right_low, $right_high, $right_flag );
        foreach ( $i + 1 .. $i + $vicinal_number ) {
            next unless $_ < $count;
            push @right_vicinal_windows, $windows->[$_]->{sw_gc};
        }

        if ( scalar @right_vicinal_windows == 0 ) {
            $right_low  = 0;
            $right_high = 0;
        }
        else {
            if ( $windows->[$i]->{sw_gc} >= max(@right_vicinal_windows) ) {
                $right_low = 1;
            }
            if ( $windows->[$i]->{sw_gc} <= min(@right_vicinal_windows) ) {
                $right_high = 1;
            }
        }

        if ( $right_low and $right_high ) {
            print " " x 4, "Can't determine right vicinity.\n";
            $right_flag = 'N';    # non
        }
        elsif ($right_low) {
            $right_flag = 'L';    # low
        }
        elsif ($right_high) {
            $right_flag = 'H';    # high
        }
        else {
            $right_flag = 'N';
        }

        $windows->[$i]->{right_flag} = $right_flag;

        #----------------------------#
        # left vicinal window
        #----------------------------#
        my @left_vicinal_windows;
        my ( $left_low, $left_high, $left_flag );
        foreach ( $i - $vicinal_number .. $i - 1 ) {
            next if $_ < 0;
            push @left_vicinal_windows, $windows->[$_]->{sw_gc};
        }

        if ( scalar @left_vicinal_windows == 0 ) {
            $left_low  = 0;
            $left_high = 0;
        }
        else {
            if ( $windows->[$i]->{sw_gc} >= max(@left_vicinal_windows) ) {
                $left_low = 1;
            }
            if ( $windows->[$i]->{sw_gc} <= min(@left_vicinal_windows) ) {
                $left_high = 1;
            }
        }

        if ( $left_low and $left_high ) {
            print " " x 4, "Can't determine left vicinity.\n";
            $left_flag = 'N';
        }
        elsif ($left_low) {
            $left_flag = 'L';
        }
        elsif ($left_high) {
            $left_flag = 'H';
        }
        else {
            $left_flag = 'N';
        }

        $windows->[$i]->{left_flag} = $left_flag;

        # high or low window
        my ($high_low_flag);

        if ( $left_flag eq 'H' and $right_flag eq 'H' ) {
            $high_low_flag = 'T';    # trough
        }
        elsif ( $left_flag eq 'L' and $right_flag eq 'L' ) {
            $high_low_flag = 'C';    # crest
        }
        else {
            $high_low_flag = 'N';
        }
        $windows->[$i]->{high_low_flag} = $high_low_flag;
    }

    return;
}

sub find_extreme_step2 {
    my $self    = shift;
    my $windows = shift;

    my $count = scalar @$windows;

    my $fall_range       = $self->fall_range;
    my $vicinal_size     = $self->vicinal_size;
    my $wave_window_step = $self->wave_window_step;
    my $vicinal_number   = int( $vicinal_size / $wave_window_step );

REDO: while (1) {
        my @extreme;
        for my $i ( 0 .. $count - 1 ) {
            my $flag = $windows->[$i]->{high_low_flag};
            if ( $flag eq 'T' or $flag eq 'C' ) {
                push @extreme, $i;
            }
        }

        # Ensure there are no crest--crest or trough--trough
        for my $i ( 0 .. scalar @extreme - 1 ) {
            my $wave = $windows->[ $extreme[$i] ]->{high_low_flag};
            my $gc   = $windows->[ $extreme[$i] ]->{sw_gc};
            my ( $left_wave, $right_wave ) = ( ' ', ' ' );
            my ( $left_gc,   $right_gc )   = ( 0,   0 );

            if ( $i - 1 >= 0 ) {
                $left_gc   = $windows->[ $extreme[ $i - 1 ] ]->{sw_gc};
                $left_wave = $windows->[ $extreme[ $i - 1 ] ]->{high_low_flag};
            }
            if ( $i + 1 < scalar @extreme ) {
                $right_gc   = $windows->[ $extreme[ $i + 1 ] ]->{sw_gc};
                $right_wave = $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag};
            }

            if ( $wave eq 'T' ) {
                if ( $wave eq $left_wave ) {
                    if ( $gc <= $left_gc ) {
                        $windows->[ $extreme[ $i - 1 ] ]->{high_low_flag} = 'N';
                        next REDO;
                    }
                }
                if ( $wave eq $right_wave ) {
                    if ( $gc < $right_gc ) {
                        $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag} = 'N';
                        next REDO;
                    }
                }
            }
            elsif ( $wave eq 'C' ) {
                if ( $wave eq $left_wave ) {
                    if ( $gc >= $left_gc ) {
                        $windows->[ $extreme[ $i - 1 ] ]->{high_low_flag} = 'N';
                        next REDO;
                    }
                }
                if ( $wave eq $right_wave ) {
                    if ( $gc > $right_gc ) {
                        $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag} = 'N';
                        next REDO;
                    }
                }
            }
            else {
                warn "high_low_flag signal errors\n";
            }
        }

        # delete nesting crest--trough
        for my $i ( 0 .. scalar @extreme - 1 ) {
            my $wave = $windows->[ $extreme[$i] ]->{high_low_flag};
            my $gc   = $windows->[ $extreme[$i] ]->{sw_gc};
            my ( $wave1, $wave2, $wave3 ) = ( ' ', ' ', ' ' );
            my ( $gc1,   $gc2,   $gc3 )   = ( 0,   0,   0 );

            if ( $i + 3 < scalar @extreme ) {

                # next 1
                $wave1 = $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag};
                $gc1   = $windows->[ $extreme[ $i + 1 ] ]->{sw_gc};

                # next 2
                $wave2 = $windows->[ $extreme[ $i + 2 ] ]->{high_low_flag};
                $gc2   = $windows->[ $extreme[ $i + 2 ] ]->{sw_gc};

                # next 3
                $wave3 = $windows->[ $extreme[ $i + 3 ] ]->{high_low_flag};
                $gc3   = $windows->[ $extreme[ $i + 3 ] ]->{sw_gc};
            }
            else {
                next;
            }

            if ( abs( $gc1 - $gc2 ) < $fall_range ) {
                if (    max( $gc1, $gc2 ) < max( $gc, $gc3 )
                    and min( $gc1, $gc2 ) > min( $gc, $gc3 ) )
                {
                    $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag} = 'N';
                    $windows->[ $extreme[ $i + 2 ] ]->{high_low_flag} = 'N';
                    next REDO;
                }
            }
        }

        # delete small fall-range crest--trough
        for my $i ( 0 .. scalar @extreme - 1 ) {
            my $wave = $windows->[ $extreme[$i] ]->{high_low_flag};
            my $gc   = $windows->[ $extreme[$i] ]->{sw_gc};
            my ( $wave1, $wave2 ) = ( ' ', ' ' );
            my ( $gc1,   $gc2 )   = ( 0,   0 );

            if ( $i + 2 < scalar @extreme ) {

                # next 1
                $wave1 = $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag};
                $gc1   = $windows->[ $extreme[ $i + 1 ] ]->{sw_gc};

                # next 2
                $wave2 = $windows->[ $extreme[ $i + 2 ] ]->{high_low_flag};
                $gc2   = $windows->[ $extreme[ $i + 2 ] ]->{sw_gc};
            }
            else {
                next;
            }

            if (    abs( $gc1 - $gc ) < $fall_range
                and abs( $gc2 - $gc1 ) < $fall_range )
            {
                $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag} = 'N';
                next REDO;
            }
        }

        # delete vicinal crest--trough
        for my $i ( 0 .. scalar @extreme - 1 ) {
            my $wave = $windows->[ $extreme[$i] ]->{high_low_flag};
            my $gc   = $windows->[ $extreme[$i] ]->{sw_gc};

            if ( $i + 1 < scalar @extreme ) {
                if ( $extreme[ $i + 1 ] - $extreme[$i] <= $vicinal_number ) {
                    $windows->[ $extreme[ $i + 1 ] ]->{high_low_flag} = 'N';
                    next REDO;
                }
            }
        }

        my @extreme2;
        for my $i ( 0 .. $count - 1 ) {
            my $flag = $windows->[$i]->{high_low_flag};
            if ( $flag eq 'T' or $flag eq 'C' ) {
                push @extreme2, $i;
            }
        }

        if ( scalar @extreme == scalar @extreme2 ) {
            last;
        }
    }

    return;
}

sub gc_wave {
    my $self           = shift;
    my $align_id       = shift;
    my $comparable_set = shift;

    my @seqs = @{ $self->get_seqs($align_id) };

    my $comparable_number = $comparable_set->cardinality;

    my $fall_range = $self->fall_range;

    my @sliding_windows = $self->sliding_window($comparable_set);

    my @sliding_attrs;
    for my $i ( 0 .. $#sliding_windows ) {
        my $sliding_window = $sliding_windows[$i];

        my @sw_seqs = map { $sliding_window->substr_span($_) } @seqs;
        my $sw_gc = calc_gc_ratio(@sw_seqs);

        my $sw_length = $sliding_window->cardinality;
        my $sw_span   = scalar $sliding_window->spans;
        my $sw_indel  = $sw_span - 1;

        $sliding_attrs[$i] = {
            sw_set    => $sliding_window,
            sw_length => $sw_length,
            sw_gc     => $sw_gc,
            sw_indel  => $sw_indel,
        };
    }

    $self->find_extreme_step1( \@sliding_attrs );
    $self->find_extreme_step2( \@sliding_attrs );

    my @slidings;
    for (@sliding_attrs) {
        push @slidings,
            {
            set           => $_->{sw_set},
            high_low_flag => $_->{high_low_flag},
            gc            => $_->{sw_gc},
            };
    }

    return @slidings;
}

sub insert_extreme {
    my $self           = shift;
    my $align_id       = shift;
    my $comparable_set = shift;
    my $align_length   = shift;

    my $dbh = $self->dbh;

    # get extreme sliding windows' sizes
    my $windows_size = $self->wave_window_size;
    my $half_length  = int( $windows_size / 2 );

    my @extreme_site;
    my @slidings = $self->gc_wave( $align_id, $comparable_set );

    foreach (@slidings) {
        my $flag = $_->{high_low_flag};
        if ( $flag eq 'T' or $flag eq 'C' ) {
            push @extreme_site, $_;
        }
    }

    my $prev_extreme_middle_right = 1;
    my $prev_extreme_gc           = $slidings[0]->{gc};
    for ( my $i = 0; $i < scalar @extreme_site; $i++ ) {

        # wave_length
        my $extreme_set         = $extreme_site[$i]->{set};
        my $extreme_middle_left = $extreme_set->lookup_index($half_length);
        my $left_wave_length
            = $extreme_middle_left - $prev_extreme_middle_right + 1;
        $prev_extreme_middle_right = $extreme_middle_left + 1;
        $extreme_site[$i]->{left_wave_length} = $left_wave_length;

        # amplitude
        my $extreme_gc     = $extreme_site[$i]->{gc};
        my $left_amplitude = abs( $extreme_gc - $prev_extreme_gc );
        $extreme_site[$i]->{left_amplitude} = $left_amplitude;
        $prev_extreme_gc = $extreme_gc;
    }

    my $next_extreme_middle_left = $align_length;
    my $next_extreme_gc          = $slidings[-1]->{gc};
    for ( my $i = scalar @extreme_site - 1; $i >= 0; $i-- ) {

        # wave_length
        my $extreme_set = $extreme_site[$i]->{set};
        my $extreme_middle_right
            = $extreme_set->lookup_index( $half_length + 1 );
        my $right_wave_length
            = $next_extreme_middle_left - $extreme_middle_right + 1;
        $next_extreme_middle_left = $extreme_middle_right - 1;
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

    my $extreme_insert = $dbh->prepare($extreme_sql);

    my $prev_extreme_id = 0;
    for (@extreme_site) {
        my $extreme_set = $_->{set};

        my ($cur_window_id) = $self->insert_window( $align_id, $extreme_set );

        $extreme_insert->execute(
            $prev_extreme_id,      $cur_window_id,
            $_->{high_low_flag},   $_->{left_amplitude},
            $_->{right_amplitude}, $_->{left_wave_length},
            $_->{right_wave_length}
        );

        $prev_extreme_id = $self->last_insert_id;
    }
    $extreme_insert->finish;

    return;
}

sub insert_gsw {
    my $self           = shift;
    my $align_id       = shift;
    my $comparable_set = shift;

    my $dbh = $self->dbh;

    # get gc sliding windows' sizes
    my $gsw_size  = $self->gsw_size;
    my $gsw0_size = int( $gsw_size / 2 );

    # extreme_id & prev_extreme_id
    my $fetch_ex_id = $dbh->prepare(
        q{
        SELECT e.extreme_id, e.prev_extreme_id, e.extreme_type
        FROM extreme e, window w
        WHERE e.window_id = w.window_id
        AND w.align_id = ?
        }
    );
    $fetch_ex_id->execute($align_id);

    # extreme_runlist, _amplitude and _average_gc
    my $fetch_ex_attr = $dbh->prepare(
        q{
        SELECT w.window_runlist, e.extreme_left_amplitude, w.window_average_gc
        FROM extreme e, window w
        WHERE e.window_id = w.window_id
        AND e.extreme_id = ?
        }
    );

    # prepare gsw_insert
    my $gsw_insert = $dbh->prepare(
        q{
        INSERT INTO gsw (
            gsw_id, extreme_id, prev_extreme_id, window_id,
            gsw_type, gsw_distance, gsw_wave_length,
            gsw_amplitude, gsw_trough_gc, gsw_gradient
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?, ?,
            ?, ?
        )
        }
    );

    while ( my $ref = $fetch_ex_id->fetchrow_hashref ) {
        my $ex_id      = $ref->{extreme_id};
        my $prev_ex_id = $ref->{prev_extreme_id};
        my $ex_type    = $ref->{extreme_type};

        # determining $gsw_type here, ascend and descent
        my $gsw_type;
        if ( $ex_type eq 'C' ) {    # crest
            $gsw_type = 'R';        # ascend, right side of trough
        }
        elsif ( $ex_type eq 'T' ) {    # trough
            $gsw_type = 'L';           # descend, left side of trough
        }
        else {
            warn "extreme_type [$ex_type] error\n";
        }

        # bypass the first extreme
        if ( $prev_ex_id == 0 ) {
            next;
        }

        # get attrs of the two extremes
        $fetch_ex_attr->execute($prev_ex_id);
        my ( $prev_ex_runlist, undef, $prev_ex_gc )
            = $fetch_ex_attr->fetchrow_array;
        my $prev_ex_set = AlignDB::IntSpan->new($prev_ex_runlist);

        $fetch_ex_attr->execute($ex_id);
        my ( $ex_runlist, $ex_left_amplitude, $ex_gc )
            = $fetch_ex_attr->fetchrow_array;
        my $ex_set = AlignDB::IntSpan->new($ex_runlist);

        # find middle points of extremes
        # the interval is between two middle points
        my $half_length        = int( $ex_set->cardinality / 2 );
        my $prev_middle_right  = $prev_ex_set->at( $half_length + 1 );
        my $ex_middle_left     = $ex_set->at($half_length);
        my $interval_start_idx = $comparable_set->index($prev_middle_right) + 1;
        my $interval_end_idx   = $comparable_set->index($ex_middle_left) - 1;

        next if $interval_start_idx > $interval_end_idx;

        my $interval_length = $interval_end_idx - $interval_start_idx + 1;

        # determining $gsw_density here, which is different from isw_density
        # and start at "0", because the first windows is within the trough
        # windows.
        my $gsw_density = int( ( $interval_length - $gsw0_size ) / $gsw_size );

        # wave length, amplitude, trough_gc and gradient
        my $gsw_wave_length = $interval_length;
        my $gsw_amplitude   = $ex_left_amplitude;
        my $gsw_trough_gc   = $ex_type eq 'T' ? $ex_gc : $prev_ex_gc;
        my $gsw_gradient    = $gsw_amplitude / $interval_length;

        {    # More windows will be submitted in the following section

            # distance start counting from trough
            # $sw_start and $sw_end are both index of $comprarable_set
            # $gsw_distance is from 0 to $gsw_density
            # length of window 0 is 50
            # ..., L2, L1, L0, R0, R1, R2, ...
            # Two window 0 have a total length of 100
            my ( $sw_start, $sw_end );
            if ( $gsw_type eq 'R' ) {
                $sw_start = $interval_start_idx;
                $sw_end   = $sw_start + $gsw0_size - 1;
            }
            elsif ( $gsw_type eq 'L' ) {
                $sw_end   = $interval_end_idx;
                $sw_start = $sw_end - $gsw0_size + 1;
            }
            else {
                warn "gsw_type [$gsw_type] error\n";
            }

            for my $gsw_distance ( 0 .. $gsw_density ) {
                my $gsw_set = $comparable_set->slice( $sw_start, $sw_end );

                my ($cur_window_id)
                    = $self->insert_window( $align_id, $gsw_set );

                $gsw_insert->execute(
                    $ex_id,         $prev_ex_id,    $cur_window_id,
                    $gsw_type,      $gsw_distance,  $gsw_wave_length,
                    $gsw_amplitude, $gsw_trough_gc, $gsw_gradient,
                );

                if ( $gsw_type eq 'R' ) {
                    $sw_start = $sw_end + 1;
                    $sw_end   = $sw_start + $gsw_size - 1;
                }
                elsif ( $gsw_type eq 'L' ) {
                    $sw_end   = $sw_start - 1;
                    $sw_start = $sw_end - $gsw_size + 1;
                }
            }
        }
    }

    return;
}

1;
