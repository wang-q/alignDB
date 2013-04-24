package AlignDB::Ofg;
use Moose::Role;
use Moose::Util::TypeConstraints;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use AlignDB::DeltaG;

requires 'dbh';
requires 'get_target_info';
requires 'get_seqs';
requires 'process_message';
requires 'insert_window';

use FindBin;
use lib "$FindBin::Bin/..";
use AlignDB::Position;

# sliding windows' size, default is 100
has 'sw_size' => ( is => 'rw', isa => 'Int', default => sub {100}, );

# max outside window distance
has 'max_out_distance' => (
    is      => 'rw',
    isa     => 'Int',
    default => sub {30},
);

# max inside window distance
has 'max_in_distance' => (
    is      => 'rw',
    isa     => 'Int',
    default => sub {30},
);

# edge
# ------+------------------------------------+--------
#    2 1 -1 -2     -89 -90  -90 -89     -2 -1 1 2
#
# edge_only
# ------+-----------------------+--------
#    2 1 -1 -2             -2 -1 1 2
#
# center
# ------+-----------------+--------
#  7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8
#
# center_intact
# ------+-----------------+--------
#  3 2 1        0          1 2 3 4
enum 'Styles', [qw(edge edge_only center center_intact)];
has 'style' => ( is => 'rw', isa => 'Styles', default => 'edge', );

# dG calculator
has 'insert_dG' => ( is => 'rw', isa => 'Bool', default => 0, );
has 'dG_calc' =>
    ( is => 'ro', isa => 'Object', default => sub { AlignDB::DeltaG->new }, );

sub insert_ofg {
    my $self        = shift;
    my $align_ids   = shift;
    my $all_ofgs    = shift;
    my $chr_ofg_set = shift;

    my $dbh = $self->dbh;
    my $pos_finder = AlignDB::Position->new( dbh => $dbh );

    # insert into ofg
    my $ofg_insert_sth = $dbh->prepare(
        q{
        INSERT INTO ofg (
            ofg_id, window_id, ofg_tag, ofg_type
        )
        VALUES (
            NULL, ?, ?, ?
        )
        }
    );

    # for each alignment
    for my $align_id (@$align_ids) {
        my $target_info    = $self->get_target_info($align_id);
        my $chr_name       = $target_info->{chr_name};
        my $chr_start      = $target_info->{chr_start};
        my $chr_end        = $target_info->{chr_end};
        my $align_length   = $target_info->{align_length};
        my $target_runlist = $target_info->{seq_runlist};
        my ( $target_seq, $query_seq ) = @{ $self->get_seqs($align_id) };

        next if $chr_name =~ /rand|un|contig|hap|scaf/i;

        $self->process_message($align_id);

        $chr_name =~ s/chr0?//i;
        my $chr_set = AlignDB::IntSpan->new("$chr_start-$chr_end");

        # chr_ofg_set has intersect with chr_set
        #   ? there are ofgs in this alignmet
        #   : there is no ofg
        next unless exists $chr_ofg_set->{$chr_name};
        next if $chr_ofg_set->{$chr_name}->intersect($chr_set)->is_empty;
        my @align_ofg;
        for (@$all_ofgs) {
            next if $_->{chr} ne $chr_name;
            my $iset = $_->{set}->intersect($chr_set);
            if ( $iset->is_not_empty ) {
                print ' ' x 4, "Find ofg: $chr_name $iset\n";
                push @align_ofg,
                    { set => $iset, tag => $_->{tag}, type => $_->{type} };
            }
        }
        if ( @align_ofg == 0 ) {
            warn "Match wrong ofg positions\n";
            next;
        }

        # target runlist
        my $target_set = AlignDB::IntSpan->new($target_runlist);

        # insert internal indels, that are, indels in target_set
        # indels in query_set is equal to spans of target_set minus one
        my $internal_indel_flag = 1;

        #----------------------------#
        # INSERT INTO ofg
        #----------------------------#
        for my $ofg (@align_ofg) {
            my $ofg_set  = $ofg->{set};
            my $ofg_tag  = $ofg->{tag};
            my $ofg_type = $ofg->{type};

            # ofg position set
            my $ofg_start = $pos_finder->at_align( $align_id, $ofg_set->min );
            my $ofg_end   = $pos_finder->at_align( $align_id, $ofg_set->max );
            next if $ofg_start > $ofg_end;
            $ofg_set   = AlignDB::IntSpan->new("$ofg_start-$ofg_end");
            $ofg_set   = $ofg_set->intersect($target_set);
            $ofg_start = $ofg_set->min;
            $ofg_end   = $ofg_set->max;

            # window
            my ($cur_window_id)
                = $self->insert_window( $align_id, $ofg_set,
                $internal_indel_flag );

            # insert to table
            $ofg_insert_sth->execute( $cur_window_id, $ofg_tag, $ofg_type );
        }

        #----------------------------#
        # INSERT INTO ofgsw
        #----------------------------#
        $self->insert_ofgsw( $align_id, $self->style );
    }

    return;
}

sub insert_ofgsw {
    my $self     = shift;
    my $align_id = shift;
    my $style    = shift;

    my $dbh          = $self->dbh;
    my $window_maker = AlignDB::Window->new(
        sw_size          => $self->sw_size,
        max_out_distance => $self->max_out_distance,
        max_in_distance  => $self->max_in_distance,
    );

    my $target_info    = $self->get_target_info($align_id);
    my $chr_name       = $target_info->{chr_name};
    my $chr_start      = $target_info->{chr_start};
    my $chr_end        = $target_info->{chr_end};
    my $align_length   = $target_info->{align_length};
    my $target_runlist = $target_info->{seq_runlist};
    my ( $target_seq, $query_seq ) = @{ $self->get_seqs($align_id) };

    # target runlist
    my $target_set = AlignDB::IntSpan->new($target_runlist);

    # insert internal indels, that are, indels in target_set
    # indels in query_set is equal to spans of target_set minus one
    my $internal_indel_flag = 1;

    # ofg_id
    my $fetch_ofg_id = $dbh->prepare(
        q{
        SELECT ofg_id
          FROM ofg o, window w
         WHERE o.window_id = w.window_id
           AND w.align_id = ?
        }
    );
    $fetch_ofg_id->execute($align_id);

    # ofg_info
    my $fetch_ofg_info = $dbh->prepare(
        q{
        SELECT w.window_runlist
          FROM ofg o, window w
         WHERE o.window_id = w.window_id
           AND o.ofg_id = ?
        }
    );

    # prepare ofgsw_insert
    my $ofgsw_insert = $dbh->prepare(
        q{
        INSERT INTO ofgsw (
            ofgsw_id, window_id, ofg_id,
            ofgsw_type, ofgsw_distance, ofgsw_dG
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?
        )
        }
    );

    my $ofgsw_size = $window_maker->sw_size;

    while ( my ($ofg_id) = $fetch_ofg_id->fetchrow_array ) {

        $fetch_ofg_info->execute($ofg_id);
        my ($ofg_runlist) = $fetch_ofg_info->fetchrow_array;
        my $ofg_set = AlignDB::IntSpan->new($ofg_runlist);

        if ( $style =~ /^edge/ ) {

            # outside rsw
            my @out_rsw
                = $window_maker->outside_window( $target_set, $ofg_set->min,
                $ofg_set->max );

            for my $outside_rsw (@out_rsw) {
                my ($cur_window_id)
                    = $self->insert_window( $align_id, $outside_rsw->{set},
                    $internal_indel_flag );

                my $deltaG
                    = $self->insert_dG
                    ? $self->_calc_deltaG( $align_id, $outside_rsw->{set} )
                    : undef;
                $ofgsw_insert->execute( $cur_window_id, $ofg_id,
                    $outside_rsw->{type}, $outside_rsw->{distance}, $deltaG, );
            }

            # inside rsw
            my @in_rsw
                = $window_maker->inside_window( $target_set, $ofg_set->min,
                $ofg_set->max );

            for my $inside_rsw (@in_rsw) {
                my ($cur_window_id)
                    = $self->insert_window( $align_id, $inside_rsw->{set},
                    $internal_indel_flag );

                my $deltaG
                    = $self->insert_dG
                    ? $self->_calc_deltaG( $align_id, $inside_rsw->{set} )
                    : undef;
                $ofgsw_insert->execute( $cur_window_id, $ofg_id,
                    $inside_rsw->{type}, $inside_rsw->{distance}, $deltaG, );
            }

            if ( $style ne 'edge_only' ) {

                # inside rsw 2
                # rsw2 start from -90, so there will be no conflicts with rsw
                my @in_rsw2
                    = $window_maker->inside_window2( $target_set, $ofg_set->min,
                    $ofg_set->max );

                for my $inside_rsw (@in_rsw2) {
                    my ($cur_window_id)
                        = $self->insert_window( $align_id, $inside_rsw->{set},
                        $internal_indel_flag );

                    my $deltaG
                        = $self->insert_dG
                        ? $self->_calc_deltaG( $align_id, $inside_rsw->{set} )
                        : undef;
                    $ofgsw_insert->execute( $cur_window_id, $ofg_id,
                        $inside_rsw->{type}, $inside_rsw->{distance}, $deltaG,
                    );
                }
            }
        }
        elsif ( $style eq 'center' ) {
            my @center_rsw
                = $window_maker->center_window( $target_set, $ofg_set->min,
                $ofg_set->max );

            for my $rsw (@center_rsw) {
                my ($cur_window_id)
                    = $self->insert_window( $align_id, $rsw->{set},
                    $internal_indel_flag );

                my $deltaG
                    = $self->insert_dG
                    ? $self->_calc_deltaG( $align_id, $rsw->{set} )
                    : undef;
                $ofgsw_insert->execute( $cur_window_id, $ofg_id, $rsw->{type},
                    $rsw->{distance}, $deltaG, );
            }
        }
        elsif ( $style eq 'center_intact' ) {
            my @center_rsw = $window_maker->center_intact_window( $target_set,
                $ofg_set->min, $ofg_set->max );

            for my $rsw (@center_rsw) {
                my ($cur_window_id)
                    = $self->insert_window( $align_id, $rsw->{set},
                    $internal_indel_flag );

                my $deltaG
                    = $self->insert_dG
                    ? $self->_calc_deltaG( $align_id, $rsw->{set} )
                    : undef;
                $ofgsw_insert->execute( $cur_window_id, $ofg_id, $rsw->{type},
                    $rsw->{distance}, $deltaG, );
            }
        }
    }

    return;
}

sub get_tags {
    my $self = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT DISTINCT o.ofg_tag
        FROM ofg o
        ORDER BY o.ofg_tag
    };

    my $ary_ref = $dbh->selectcol_arrayref($query);

    return $ary_ref;
}

sub get_types {
    my $self = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT DISTINCT o.ofg_type
        FROM ofg o
        ORDER BY o.ofg_type
    };

    my $ary_ref = $dbh->selectcol_arrayref($query);

    return $ary_ref;
}

sub get_tts {
    my $self = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT DISTINCT CONCAT(o.ofg_tag, "_", o.ofg_type)
        FROM ofg o
        ORDER BY CONCAT(o.ofg_tag, "_", o.ofg_type)
    };

    my $ary_ref = $dbh->selectcol_arrayref($query);

    return $ary_ref;
}

sub _calc_deltaG {
    my $self     = shift;
    my $align_id = shift;
    my $set      = shift;

    my $seqs_ref   = $self->get_seqs($align_id);
    my @seq_slices = map { $set->substr_span($_) } @{$seqs_ref};
    my @seq_dG     = map { $self->dG_calc->polymer_deltaG($_) } @seq_slices;

    return mean(@seq_dG);
}

1;

__END__
