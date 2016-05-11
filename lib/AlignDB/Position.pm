package AlignDB::Position;
use Moose;

use YAML::Syck;

use AlignDB::IntSpan;

#----------------------------------------------------------#
# Object attributes
#----------------------------------------------------------#
has 'dbh' => ( is => 'ro', isa => 'Object' );

has 'target_align_id' => ( is => 'ro', isa => 'Int' );
has 'target_info'     => ( is => 'ro', isa => 'HashRef' );

#----------------------------------------------------------#
# Methods
#----------------------------------------------------------#
# cache target info
sub update_target_info {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;

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
        WHERE s.seq_role = "T"
        AND a.align_id = ?
    };

    my DBI $sth = $dbh->prepare($align_query);
    $sth->execute($align_id);

    my ( $align_length, $target_chr_start, $target_chr_end, $target_runlist, )
        = $sth->fetchrow_array;

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
# 42: (42 - max()) + (size() - 3 + 1) = 14 + 16 = 30
sub at_align {
    my $self     = shift;
    my $align_id = shift;
    my $chr_pos  = shift;

    $self->update_target_info($align_id);
    my $target_info = $self->target_info;

    my $align_length                = $target_info->{align_length};
    my $target_chr_start            = $target_info->{chr_start};
    my $target_chr_end              = $target_info->{chr_end};
    my AlignDB::IntSpan $target_set = $target_info->{set};

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
# if the position in target is located in a gap, then return the left base's
# position. (Just like GATK's indel left align)
sub at_target_chr {
    my $self      = shift;
    my $align_id  = shift;
    my $align_pos = shift;

    $self->update_target_info($align_id);
    my $target_info = $self->target_info;

    my $align_length                = $target_info->{align_length};
    my $target_chr_start            = $target_info->{chr_start};
    my $target_chr_end              = $target_info->{chr_end};
    my AlignDB::IntSpan $target_set = $target_info->{set};

    local $SIG{__WARN__} = sub {
        my ($sig) = @_;
        print "$sig\n";
        print YAML::Syck::Dump {
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

    my $target_pos;
    if ( $target_set->contains($align_pos) ) {
        $target_pos = $target_set->index($align_pos);
    }
    else {
        my @sets = $target_set->sets;
        for my $i ( 0 .. $#sets ) {
            if ( $sets[$i]->max < $align_pos ) {
                next;
            }
            elsif ( $i > 0 ) {
                $target_pos = $sets[ $i - 1 ]->max;
                last;
            }
            else {

                # still wrong, we will get messages
            }
        }
    }
    my $target_chr_pos = $target_chr_start + $target_pos - 1;

    return $target_chr_pos;
}

1;
