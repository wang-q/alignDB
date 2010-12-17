package AlignDB::Multi::GC;
use Moose;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Util qw(:all);

extends qw(AlignDB::Multi);

# extreme sliding window size
has 'window_size' => ( is => 'rw', isa => 'Int', default => 100, );

# extreme sliding window step
has 'window_step' => ( is => 'rw', isa => 'Int', default => 100, );

sub insert_segment {
    my $self           = shift;
    my $align_id       = shift;
    my $comparable_set = shift;

    my $dbh = $self->dbh;

    my @segment_levels = (
        [ 'A', '',    '' ],
        [ 0,   10000, 10000 ],
        [ 1,   5000,  5000 ],
        [ 2,   1000,  1000 ],
        [ 3,   500,   500 ],
    );

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
    my $window_size    = shift || $self->window_size;
    my $window_step    = shift || $self->window_step;

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
    my $window_size = shift || $self->window_size;
    my $window_step = shift || $self->window_step;

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
        = $gc_mean <= 0.5
        ? $gc_std / $gc_mean
        : $gc_std / ( 1 - $gc_mean );
    my $gc_mdcw = _mdcw(@sliding_gcs);

    return ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw );
}

sub _mdcw {
    my @array = @_;

    if ( @array <= 1 ) {
        confess "There should be more than 1 elements in the array.\n";
        return;
    }

    my @dcws;
    for ( 1 .. $#array ) {
        my $dcw = $array[$_] - $array[ $_ - 1 ];
        push @dcws, abs $dcw;
    }

    return mean(@dcws);
}

1;
