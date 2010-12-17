package AlignDB::Join;
use Moose;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all);
use Math::Combinatorics;
use Statistics::Descriptive;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#use AlignDB;
#use AlignDB::Position;

has 'goal'   => ( is => 'rw', isa => 'Str' );
has 'ref'    => ( is => 'rw', isa => 'Str' );
has 'target' => ( is => 'rw', isa => 'Str' );
has 'queries' =>
    ( is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] } );
has 'dbs'   => ( is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] } );
has 'names' => ( is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] } );
has 'db_info_of' => ( is => 'ro', isa => 'HashRef', default => sub { {} } );
has 'info_of'    => ( is => 'ro', isa => 'HashRef', default => sub { {} } );
has 'target_obj' => ( is => 'rw', isa => 'Object' );
has 'ref_obj'    => ( is => 'rw', isa => 'Object' );

has 'inter_chr' => ( is => 'ro', isa => 'HashRef', default => sub { {} } );

has 'server' => ( is => 'ro', isa => 'Str' );    # e.g. '202.119.43.5'
has 'user'   => ( is => 'ro', isa => 'Str' );    # database username
has 'passwd' => ( is => 'ro', isa => 'Str' );    # database password

has reduce_end       => ( is => 'rw', isa => 'Int', default => 0 );
has length_threshold => ( is => 'rw', isa => 'Int', default => 1000 );
has discard_paralog  => ( is => 'rw', isa => 'Int', default => 0 );
has discard_distant  => ( is => 'rw', isa => 'Int', default => 0 );
has percentile_discard => ( is => 'rw', isa => 'Num' );

has indel_expand => ( is => 'rw', isa => 'Int', default => 25 );
has indel_join   => ( is => 'rw', isa => 'Int', default => 25 );

has raw_fasta => ( is => 'rw', isa => 'Int', default => 0 );

sub BUILD {
    my $self = shift;
    my @names = ( $self->ref, $self->target, @{ $self->queries } );
    $self->{names} = \@names;

    $self->{target} =~ /^(\d+)(.+)/;
    $self->{target_db} = $self->{dbs}->[$1];
    $self->{ref} =~ /^(\d+)(.+)/;
    $self->{ref_db} = $self->{dbs}->[$1];

    if ( $self->discard_distant ) {
        $self->_calc_discard;
    }

    $self->_db_info;
    $self->_info_of;
    $self->_inter_chr;
    $self->_segments;

}

sub _calc_discard {
    my $self = shift;

    my $ref        = $self->ref;
    my @all_dbs    = @{ $self->dbs };
    my %db_info_of = %{ $self->db_info_of };

    $ref =~ /^(\d+)(.+)/;
    my $db_name_idx = $1;
    my $db_name     = $all_dbs[$db_name_idx];

    my $per_idn_query = qq{
        SELECT  a.identities / a.align_length
        FROM align a
    };

    my $dbh = $db_info_of{$db_name}->{dbh};
    my $sth = $dbh->prepare($per_idn_query);

    my @data;
    $sth->execute;
    while ( my @row = $sth->fetchrow_array ) {
        push @data, $row[0];
    }

    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@data);
    $self->{percentile_discard} = $stat->percentile( $self->discard_distant );
}

sub _db_info {
    my $self = shift;

    my $server = $self->server;
    my $user   = $self->user;
    my $passwd = $self->passwd;

    my %db_info_of;
    for ( @{ $self->dbs } ) {
        my $cur_obj = AlignDB->new(
            mysql  => "$_:$server",
            user   => $user,
            passwd => $passwd,
        );
        my $cur_dbh = $cur_obj->dbh;
        my $cur_pos_obj = AlignDB::Position->new( dbh => $cur_dbh );
        $db_info_of{$_} = {
            target => {
                taxon_id => '',
                name     => '',
            },
            query => {
                taxon_id => '',
                name     => '',
            },
            obj     => $cur_obj,
            pos_obj => $cur_pos_obj,
        };

        (   $db_info_of{$_}->{target}{taxon_id},
            $db_info_of{$_}->{query}{taxon_id},
        ) = $cur_obj->get_taxon_ids;

        ( $db_info_of{$_}->{target}{name}, $db_info_of{$_}->{query}{name}, )
            = $cur_obj->get_names;

        my $chr_ref = $cur_obj->get_chrs('target');

        for my $ref ( @{$chr_ref} ) {
            my ( $chr_id, $chr_name, $chr_length ) = @{$ref};
            my $chr_set = $self->_build_chr_set( $cur_dbh, $chr_id );
            $db_info_of{$_}->{chrs}{$chr_name} = $chr_set;
        }
    }

    $self->{db_info_of} = \%db_info_of;

    $self->{target} =~ /^(\d+)(.+)/;
    my $target_db = $self->{dbs}->[$1];
    $self->{target_obj} = $db_info_of{$target_db}->{obj};
    $self->{ref} =~ /^(\d+)(.+)/;
    my $ref_db = $self->{dbs}->[$1];
    $self->{ref_obj} = $db_info_of{$ref_db}->{obj};
}

sub _info_of {
    my $self = shift;

    my @all_names  = @{ $self->names };
    my @all_dbs    = @{ $self->dbs };
    my %db_info_of = %{ $self->db_info_of };

    #----------------------------#
    # build %info_of all_names hash
    #----------------------------#
    my %info_of;
    for my $name (@all_names) {
        $name =~ /^(\d+)(.+)/;
        my $db_name_idx = $1;
        my $torq        = $2;
        if ( not( $torq =~ /^t/i or $torq =~ /^q/i ) ) {
            die "$torq is not equal to target or query\n";
        }
        my $db_name = $all_dbs[$db_name_idx];
        $info_of{$name} = $db_info_of{$db_name}->{$torq};
    }

    $self->{info_of} = \%info_of;

}

sub _inter_chr {
    my $self = shift;

    my %inter_chr;

    my $target_db  = $self->target_obj->db;
    my $db_info_of = $self->db_info_of;

    my @chrs = sort keys %{ $db_info_of->{$target_db}{chrs} };

    for my $chr (@chrs) {
        my $inter_chr_set = AlignDB::IntSpan->new;
        for my $db ( @{ $self->dbs } ) {
            my $cur_chr_set = $db_info_of->{$target_db}{chrs}{$chr};
            if ( $inter_chr_set->is_empty ) {
                $inter_chr_set = $cur_chr_set;
            }
            else {
                $inter_chr_set = $inter_chr_set->intersect($cur_chr_set);
            }
        }
        $inter_chr{$chr} = $inter_chr_set;
    }

    $self->{inter_chr} = \%inter_chr;
}

sub _build_chr_set {
    my $self   = shift;
    my $dbh    = shift;
    my $chr_id = shift;

    my $reduce_end = $self->reduce_end;

    my $chr_set = AlignDB::IntSpan->new;

    my $chr_query = qq{
        SELECT  s.chr_start + $reduce_end,
                s.chr_end - $reduce_end
        FROM sequence s, chromosome c
        WHERE c.chr_id = ?
        AND s.chr_id = c.chr_id
    };

    # build $chr_set
    my $chr_sth = $dbh->prepare($chr_query);
    $chr_sth->execute($chr_id);
    while ( my @row = $chr_sth->fetchrow_array ) {
        my ( $chr_start, $chr_end ) = @row;
        next if $chr_start > $chr_end;
        $chr_set->add_range( $chr_start, $chr_end );
    }

    return $chr_set;
}

sub _segments {
    my $self      = shift;
    my $chr_name  = shift;
    my $seg_start = shift;
    my $seg_end   = shift;

    my $seg_length = $seg_end - $seg_start + 1;
    return $seg_length <= $self->length_threshold;

    print "\n$seg_length,$seg_start-$seg_end\n";

    my @all_dbs    = @{ $self->dbs };
    my @all_names  = @{ $self->names };
    my %db_info_of = %{ $self->db_info_of };
    my %info_of    = %{ $self->info_of };

    # target chr position
    for my $db_name (@all_dbs) {
        $db_info_of{$db_name}->{target}{chr_start} = $seg_start;
        $db_info_of{$db_name}->{target}{chr_end}   = $seg_end;
    }

    for my $db_name (@all_dbs) {
        my $pos_obj = $db_info_of{$db_name}->{pos_obj};
        my ( $align_id, $dummy ) = @{
            $pos_obj->positioning_align_chr_name( $chr_name, $seg_start,
                $seg_end )
            };

        if ( !defined $align_id ) {
            warn " " x 4, "Find no align in $db_name, jump to next\n";
            return;
        }
        elsif ( defined $dummy ) {
            warn " " x 4, "Overlapped alignment in $db_name!\n";
            if ( $self->discard_paralog ) {
                return;
            }
        }
        $db_info_of{$db_name}->{align_id} = $align_id;
    }

    #----------------------------#
    # get seq, use align coordinates
    #----------------------------#
    for my $db_name (@all_dbs) {
        print " " x 4, "build $db_name seq\n";
        my $align_id = $db_info_of{$db_name}->{align_id};

        my $error
            = $self->build_seq( $db_info_of{$db_name}, $seg_start, $seg_end );
        if ($error) {
            warn $error . " in $db_name $align_id\n";
            return;
        }
    }

    #----------------------------#
    # discard alignments which have low percentage identity to outgroup
    #----------------------------#
    if ( $self->discard_distant ) {
        my $ref = $self->ref;
        $ref =~ /^(\d+)(.+)/;
        my $db_name_idx = $1;
        my $db_name     = $all_dbs[$db_name_idx];
        my $result      = pair_seq_stat(
            $db_info_of{$db_name}->{target}{seqs},
            $db_info_of{$db_name}->{query}{seqs},
        );
        my $seq_legnth = $result->[0];
        my $identities = $result->[2];
        my $per_idn    = $identities / $seq_legnth;

        if ( $per_idn >= $self->percentile_discard ) {
            warn " " x 4 . "Low percentage identity with outgroup\n";
            return;
        }
    }

    #----------------------------#
    # start peusdo-alignment, according to common sequences
    #----------------------------#
    print " " x 4, "start peusdo-alignment\n";
    my $pos_count = 0;
    while (1) {
        $pos_count++;
        my $max_length = 0;
        for my $db_name (@all_dbs) {
            $max_length = max( $max_length,
                length $db_info_of{$db_name}->{target}{seqs} );
        }
        if ( $pos_count >= $max_length ) {
            last;
        }

        my @target_bases;
        for my $db_name (@all_dbs) {
            push @target_bases,
                substr( $db_info_of{$db_name}->{target}{seqs},
                $pos_count - 1, 1 );
        }

        if ( all { $_ eq $target_bases[0] } @target_bases ) {
            next;
        }
        elsif ( all { $_ ne '-' } @target_bases ) {
            warn " " x 8 . "align error in $pos_count, [@target_bases]\n";
            return;
        }

        # insert a '-' in current position
        for ( 0 .. @all_dbs - 1 ) {
            my $db_name = $all_dbs[$_];
            if ( $target_bases[$_] eq '-' ) {
                next;
            }
            else {
                substr(
                    $db_info_of{$db_name}->{target}{seqs},
                    $pos_count - 1,
                    0, '-'
                );
                substr(
                    $db_info_of{$db_name}->{query}{seqs},
                    $pos_count - 1,
                    0, '-'
                );
            }
        }
    }

    #----------------------------#
    # clustalw realign indel_flank region
    #----------------------------#
    {
        print " " x 4, "start finding realign region\n";
        $self->realign;
    }

    #----------------------------#
    # output a raw fasta alignment for further use
    #----------------------------#
    if ( $self->raw_fasta ) {
        my $goal_db = $self->goal_db;
        unless ( -e $goal_db ) {
            mkdir $goal_db, 0777
                or die "Cannot create \"$goal_db\" directory: $!";
        }
        my $first_taxon_id = $info_of{ $all_names[1] }->{taxon_id};
        my $outfile
            = "$goal_db/" . "raw_"
            . "id$first_taxon_id" . "_"
            . $chr_name . "_"
            . $seg_start . "_"
            . $seg_end . ".fas";
        print " " x 4, "$outfile\n";
        open my $out_fh, '>', $outfile
            or die("Cannot open output file $outfile");
        for my $name (@all_names) {
            my $seq = $info_of{$name}->{seqs};
            print {$out_fh} ">", $info_of{$name}->{name}, "\n";
            print {$out_fh} $seq, "\n";
        }
        close $out_fh;
    }

    return;
}

# get seq, use align coordinates
sub _build_seq {
    my $self      = shift;
    my $db_info   = shift;
    my $seg_start = shift;
    my $seg_end   = shift;

    my $target_seq_query = q{
        SELECT  t.target_seq,
                s.chr_id,
                s.chr_strand
        FROM target t, sequence s
        WHERE t.seq_id = s.seq_id
        AND t.align_id = ?
    };
    my $query_seq_query = q{
        SELECT  q.query_seq,
                q.query_strand,
                s.chr_id,
                s.chr_strand
        FROM query q, sequence s
        WHERE q.seq_id = s.seq_id
        AND q.align_id = ?
    };

    my $dbh      = $db_info->{dbh};
    my $pos_obj  = $db_info->{pos_obj};
    my $align_id = $db_info->{align_id};

    my $target_sth = $dbh->prepare($target_seq_query);
    $target_sth->execute($align_id);
    (   $db_info->{target}{full_seqs},
        $db_info->{target}{chr_id},
        $db_info->{target}{chr_strand},
    ) = $target_sth->fetchrow_array;

    my $query_sth = $dbh->prepare($query_seq_query);
    $query_sth->execute($align_id);
    (   $db_info->{query}{full_seqs}, $db_info->{query}{query_strand},
        $db_info->{query}{chr_id},    $db_info->{query}{chr_strand},
    ) = $query_sth->fetchrow_array;

    my $align_start = $pos_obj->at_align( $align_id, $seg_start );
    my $align_end   = $pos_obj->at_align( $align_id, $seg_end );

    # align_start and align_end should must be available
    unless ( $align_start and $align_end ) {
        return " " x 8 . "align_start or align_end error";
    }

    my $align_length = $align_end - $align_start + 1;

    $db_info->{target}{seqs} = substr(
        $db_info->{target}{full_seqs},
        $align_start - 1,
        $align_length
    );
    $db_info->{query}{seqs} = substr(
        $db_info->{query}{full_seqs},
        $align_start - 1,
        $align_length
    );

    unless (length $db_info->{target}{seqs} == length $db_info->{query}{seqs}
        and length $db_info->{target}{seqs} > 0 )
    {
        return " " x 8 . "seq-length error";
    }

    delete $db_info->{target}{full_seqs};
    delete $db_info->{query}{full_seqs};

    return;
}

sub realign {
    my $self = shift;

    my %info_of   = %{ $self->info_of };
    my @all_names = @{ $self->all_names };

    # use AlignDB::IntSpan to find nearby indels
    #   expand indel by a range of $indel_expand

    my %indel_sets;
    for (@all_names) {
        $indel_sets{$_}
            = find_indel_set( $info_of{$_}->{seqs}, $self->indel_expand );
    }

    my $realign_region = AlignDB::IntSpan->new;
    my $combinat       = Math::Combinatorics->new(
        count => 2,
        data  => \@all_names,
    );
    while ( my @combo = $combinat->next_combination ) {
        print " " x 8, "pairwise correction @combo\n";
        my $intersect_set = AlignDB::IntSpan->new;
        my $union_set     = AlignDB::IntSpan->new;
        $intersect_set
            = $indel_sets{ $combo[0] }->intersect( $indel_sets{ $combo[1] } );
        $union_set
            = $indel_sets{ $combo[0] }->union( $indel_sets{ $combo[1] } );

        for my $span ( $union_set->runlists ) {
            my $flag_set = $intersect_set->intersect($span);
            if ( $flag_set->is_not_empty ) {
                $realign_region->add($span);
            }
        }
    }

    # join adjacent realign regions
    $realign_region = $realign_region->join_span( $self->indel_join );

    # realign all segments in realign_region
    my @realign_region_spans = $realign_region->spans;
    for ( reverse @realign_region_spans ) {
        my $seg_start = $_->[0];
        my $seg_end   = $_->[1];
        my @segments;
        for (@all_names) {
            my $seg = substr(
                $info_of{$_}->{seqs},
                $seg_start - 1,
                $seg_end - $seg_start + 1
            );
            push @segments, $seg;
        }

        my $realign_segments = clustal_align( \@segments );

        for (@all_names) {
            my $seg = shift @$realign_segments;
            substr(
                $info_of{$_}->{seqs},
                $seg_start - 1,
                $seg_end - $seg_start + 1, $seg
            );
        }
    }

    return;
}

1;

