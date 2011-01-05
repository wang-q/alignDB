package AlignDB::Multi::Position;
use Moose;
use Carp;

use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;

#----------------------------------------------------------#
# Object attributes
#----------------------------------------------------------#
has 'dbh' => ( is => 'ro', isa => 'Object' );

has 'target_align_id' => ( is => 'ro', isa => 'Int' );
has 'target_info'     => ( is => 'ro', isa => 'HashRef' );

has 'target_chr_pos'   => ( is => 'ro', isa => 'ArrayRef' );
has 'target_align_pos' => ( is => 'ro', isa => 'HashRef' );

#----------------------------------------------------------#
# Methods
#----------------------------------------------------------#
# cache target info
sub update_target_info {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $stored_align_id = $self->target_align_id;
    if ( defined $stored_align_id and $align_id == $stored_align_id ) {
        return;
    }

    my $align_query = q{
        SELECT a.align_length,
               s.chr_start,
               s.chr_end,
               s.seq_runlist
        FROM align a
        inner join sequence s on a.align_id = s.align_id
        inner join target t on s.seq_id = t.seq_id
        WHERE 1 = 1
        AND a.align_id = ?
    };

    my $align_sth = $dbh->prepare($align_query);
    $align_sth->execute($align_id);

    my ($align_length,   $target_chr_start,
        $target_chr_end, $target_runlist,
    ) = $align_sth->fetchrow_array;

    my $target_set = AlignDB::IntSpan->new($target_runlist);

    $self->{target_info} = {
        align_length => $align_length,
        chr_start    => $target_chr_start,
        chr_end      => $target_chr_end,
        set          => $target_set,
    };

    return;
}

# Give a target chr position, return an align position
# start at '1'
# when position is beyond alignment range, force return
#  ________CCATG--------CAGA--ATTTACCGA______________
#  |      ||   |        |  |  |       |             |
#  |      |1   5       14 17 20      28             |
# -7      0                                        42
# aa:  M  Q  N  L  P
# dna: ATGCAGAATTTACCG
# align runlist:  1-5,14-17,20-28
# coding runlist: 3-5,14-17,20-27
# Now calculate distance to translation start point
# 3:  index(3) - index(3) + 1 = 3 - 3 + 1 = 1
# 5:  index(5) - index(3) + 1 = 3
# 21: index(21) - index(3) + 1 = 11 - 3 + 1= 9
# -7: (min() - (-7)) + (3 - 1 + 1) = 8 + 3 = 11
# 42: (42 - max()) + (cardinality() - 3 + 1) = 14 + 16 = 30
sub at_align {
    my $self     = shift;
    my $align_id = shift;
    my $chr_pos  = shift;

    $self->update_target_info($align_id);
    my $target_info = $self->target_info;

    my $align_length     = $target_info->{align_length};
    my $target_chr_start = $target_info->{chr_start};
    my $target_chr_end   = $target_info->{chr_end};
    my $target_set       = $target_info->{set};

    my $align_pos;
    if ( $chr_pos < $target_chr_start ) {    # 0 or negtive interger
        $align_pos = $chr_pos - $target_chr_start + 1;
    }
    elsif ( $chr_pos > $target_chr_end ) {    # larger than align_length
        $align_pos = $chr_pos - $target_chr_end + $align_length;
    }
    else {                                    # normal
        my $target_pos = $chr_pos - $target_chr_start + 1;
        $align_pos = $target_set->at($target_pos);
    }

    return $align_pos;
}

# Give an align position, return a target chr position
# start at '1'
sub at_target_chr {
    my $self      = shift;
    my $align_id  = shift;
    my $align_pos = shift;

    $self->update_target_info($align_id);
    my $target_info = $self->target_info;

    my $align_length     = $target_info->{align_length};
    my $target_chr_start = $target_info->{chr_start};
    my $target_chr_end   = $target_info->{chr_end};
    my $target_set       = $target_info->{set};

    local $SIG{__WARN__} = sub {
        my ($sig) = @_;
        print "$sig\n";
        print Dump {
            align_id         => $align_id,
            align_pos        => $align_pos,
            align_length     => $align_length,
            target_chr_start => $target_chr_start,
            target_chr_end   => $target_chr_end,
        };
    };

    if ( $align_pos > $align_length or $align_length < 1 ) {
        return;
    }

    my $target_pos     = $target_set->index($align_pos);
    my $target_chr_pos = $target_chr_start + $target_pos - 1;

    return $target_chr_pos;
}

# example:
## chr_position and align_position transforming
#$obj->target_chr_align_transform($align_id);
#$obj->query_chr_align_transform($align_id);
#my $target_chr_pos   = $obj->target_chr_pos;    # align to chr array
#my $target_align_pos = $obj->target_align_pos;  # chr to align hash
sub target_chr_align_transform {
    my $self     = shift;
    my $align_id = shift;

    $self->update_target_info($align_id);
    my $target_info = $self->target_info;

    my $align_length     = $target_info->{align_length};
    my $target_chr_start = $target_info->{chr_start};
    my $target_chr_end   = $target_info->{chr_end};
    my $target_set       = $target_info->{set};

    # align_position to chr_position transforming array
    # the first element [0] will be ignored
    my @target_chr_pos;
    {
        my $target_cursor = 0;    # target
        for ( my $i = 1; $i <= $align_length; $i++ ) {
            if ( $target_set->contain($i) ) {
                $target_cursor++;
            }
            $target_chr_pos[$i] = $target_chr_start + $target_cursor - 1;
        }
    }

    # chr_position to align_position transforming hash
    my %target_align_pos;
    {
        for ( my $i = $target_chr_start; $i <= $target_chr_end; $i++ ) {
            my $j = $i - $target_chr_start + 1;
            $target_align_pos{$i} = $target_set->at($j);
        }
    }

    $self->{target_chr_pos}   = \@target_chr_pos;
    $self->{target_align_pos} = \%target_align_pos;

    return;
}

sub positioning_align {
    my $self  = shift;
    my $chr   = shift;
    my $start = shift;
    my $end   = shift || $start;

    my $dbh = $self->dbh;

    if ( $chr !~ /chr/i ) {
        $chr = 'chr' . $chr;
    }
    my $align_id_query = q{
        # align contain a chr position
        SELECT s.align_id
        FROM sequence s
        inner join target t on t.seq_id = s.seq_id
        inner join chromosome c on c.chr_id = s.chr_id
        WHERE 1 = 1
        AND c.chr_name like ?
        AND s.chr_start <= ?
        AND s.chr_end >= ?
    };
    my $align_id_sth = $dbh->prepare($align_id_query);
    $align_id_sth->execute( $chr, $start, $end );
    my @align_ids;
    while ( my @row = $align_id_sth->fetchrow_array ) {
        push @align_ids, $row[0];
    }

    return \@align_ids;
}

sub positioning_align_chr_id {
    my $self   = shift;
    my $chr_id = shift;
    my $start  = shift;
    my $end    = shift || $start;

    my $dbh = $self->dbh;

    my $align_id_query = q{
        # align contain a chr position
        SELECT s.align_id
        FROM sequence s
        inner join target t on s.seq_id = t.seq_id
        WHERE 1 = 1
        AND s.chr_id = ?
        AND s.chr_start <= ?
        AND s.chr_end >= ?
    };
    my $align_id_sth = $dbh->prepare($align_id_query);
    $align_id_sth->execute( $chr_id, $start, $end );
    my @align_ids;
    while ( my @row = $align_id_sth->fetchrow_array ) {
        push @align_ids, $row[0];
    }

    return \@align_ids;
}

sub positioning_align_chr_name {
    my $self   = shift;
    my $chr_name = shift;
    my $start  = shift;
    my $end    = shift || $start;

    my $dbh = $self->dbh;

    my $align_id_query = q{
        # align contain a chr position
        SELECT s.align_id
        FROM sequence s
        inner join target t on s.seq_id = t.seq_id
        inner join chromosome c on c.chr_id = s.chr_id
        WHERE 1 = 1
        AND c.chr_name = ?
        AND s.chr_start <= ?
        AND s.chr_end >= ?
    };
    my $align_id_sth = $dbh->prepare($align_id_query);
    $align_id_sth->execute( $chr_name, $start, $end );
    my @align_ids;
    while ( my @row = $align_id_sth->fetchrow_array ) {
        push @align_ids, $row[0];
    }

    return \@align_ids;
}

1;
