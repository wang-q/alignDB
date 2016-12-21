package AlignDB::Outgroup;
use Moose;
use autodie;

extends qw(AlignDB::Common);

sub add_align {
    my $self      = shift;
    my $info_refs = shift;

    for my $info ( @{$info_refs} ) {
        $info->{seq} = uc $info->{seq};
    }
    my $seq_refs = [ map { $_->{seq} } @{$info_refs} ];

    my $target_idx = 0;

    # check align length
    my $align_length = length $seq_refs->[$target_idx];
    for ( @{$seq_refs} ) {
        if ( ( length $_ ) != $align_length ) {
            Carp::confess "Sequences should have the same length!\n";
            return;
        }
    }

    # check seq number
    my $seq_count = scalar @{$seq_refs};
    if ( $seq_count < 3 ) {
        Carp::confess "Too few sequences [$seq_count]\n";
        return;
    }

    # appoint outgroup
    my $outgroup_idx = scalar( @{$info_refs} ) - 1;
    my $outgroup_seq = $seq_refs->[$outgroup_idx];

    #----------------------------#
    # INSERT INTO align
    #----------------------------#
    my $align_id = $self->_insert_align( $seq_refs, "outgroup" );
    printf "Prosess align [%s] at %s.%s(%s):%s-%s\n", $align_id,
        $info_refs->[$target_idx]{name},
        $info_refs->[$target_idx]{chr},
        $info_refs->[$target_idx]{strand},
        $info_refs->[$target_idx]{start},
        $info_refs->[$target_idx]{end};

    #----------------------------#
    # UPDATE align, INSERT sequence, target, queries and outgroup
    #----------------------------#
    $self->_insert_set_and_sequence( $align_id, $info_refs, $seq_refs, "outgroup" );

    #----------------------------#
    # INSERT INTO indel
    #----------------------------#
    $self->_insert_indel($align_id, $outgroup_seq);

    #----------------------------#
    # INSERT INTO snp
    #----------------------------#
    $self->_insert_snp($align_id, $outgroup_seq);

    return;
}

sub update_D_values {
    my $self     = shift;
    my $align_id = shift;

    # Get database handle
    my DBI $dbh = $self->dbh;

    my $isw_id_ref = $dbh->selectcol_arrayref(
        q{
        SELECT isw_id
        FROM isw w, indel i
        where 1 = 1
        AND w.isw_indel_id = i.indel_id
        AND i.align_id = ?
        },
        {},
        $align_id
    );

    my DBI $read_sth = $dbh->prepare(
        q{
        SELECT s.snp_outgroup_base, s.snp_all_bases, w.isw_length, i.indel_occured
        FROM snp s, isw w, indel i
        WHERE 1 = 1
        AND s.isw_id = w.isw_id
        AND w.isw_indel_id = i.indel_id
        AND i.indel_occured <> 'unknown'
        AND w.isw_id = ?
        }
    );

    my DBI $update_sth = $dbh->prepare(
        q{
        UPDATE isw
        SET isw_d_indel = ?,
            isw_d_noindel = ?,
            isw_d_bii = ?,
            isw_d_bnn = ?,
            isw_d_complex = ?,
            isw_d_indel2 = ?,
            isw_d_noindel2 = ?,
            isw_d_bii2 = ?,
            isw_d_bnn2 = ?,
            isw_d_complex2 = ?,
            isw_d_indel3 = ?,
            isw_d_noindel3 = ?,
            isw_d_bii3 = ?,
            isw_d_bnn3 = ?,
            isw_d_complex3 = ?
        WHERE isw_id = ?
        }
    );

    for my $isw_id ( @{$isw_id_ref} ) {
        my $window_length;
        my ( $d_indel,  $d_noindel,  $d_bii,  $d_bnn,  $d_complex )  = (0) x 5;
        my ( $d_indel2, $d_noindel2, $d_bii2, $d_bnn2, $d_complex2 ) = (0) x 5;
        my ( $d_indel3, $d_noindel3, $d_bii3, $d_bnn3, $d_complex3 ) = (0) x 5;

        my $group_i = AlignDB::IntSpan->new;
        my $group_n = AlignDB::IntSpan->new;
        my $ref_seq;
        my @sequences;

        # removes all mutations on the deepest indel branches
        my $ref_seq2;
        my @sequences2;

        # removes all mutations on the deepest indel branches without recombinations
        my $ref_seq3;
        my @sequences3;

        $read_sth->execute($isw_id);
        while ( my @row = $read_sth->fetchrow_array ) {
            $window_length = $row[2];
            my $snp_outgroup_base = $row[0];

            $ref_seq .= $snp_outgroup_base;
            my @snp_base    = split '', $row[1];
            my @indel_occur = split '', $row[3];
            my $align_cnt   = scalar @indel_occur;
            for my $i ( 0 .. $align_cnt - 1 ) {
                if ( $indel_occur[$i] eq '1' ) {
                    $group_i->add($i);
                }
                elsif ( $indel_occur[$i] eq '0' ) {
                    $group_n->add($i);
                }
                else {
                    Carp::confess "$indel_occur[$i]\n";
                }
                $sequences[$i] .= $snp_base[$i];
            }

            # find mutations on the deepest branches
            my @group_i = $group_i->elements;
            my @group_n = $group_n->elements;
            my @i_snp   = App::Fasops::Common::uniq( @snp_base[@group_i] );
            my @n_snp   = App::Fasops::Common::uniq( @snp_base[@group_n] );
            if ( @i_snp == 1 and $i_snp[0] ne $snp_outgroup_base ) {

                # removes all mutations on the deepest indel branches
                # add bases with recombination events
                if ( App::Fasops::Common::any { $_ eq $i_snp[0] } @n_snp ) {
                    $ref_seq3 .= $snp_outgroup_base;
                    for my $i ( 0 .. $align_cnt - 1 ) {
                        $sequences3[$i] .= $snp_base[$i];
                    }
                }
            }
            elsif ( @n_snp == 1 and $n_snp[0] ne $snp_outgroup_base ) {

                # removes all mutations on the deepest noindel branches
                # add bases with recombination events
                if ( App::Fasops::Common::any { $_ eq $n_snp[0] } @i_snp ) {
                    $ref_seq3 .= $snp_outgroup_base;
                    for my $i ( 0 .. $align_cnt - 1 ) {
                        $sequences3[$i] .= $snp_base[$i];
                    }
                }
            }
            else {
                $ref_seq2 .= $snp_outgroup_base;
                for my $i ( 0 .. $align_cnt - 1 ) {
                    $sequences2[$i] .= $snp_base[$i];
                }
                $ref_seq3 .= $snp_outgroup_base;
                for my $i ( 0 .. $align_cnt - 1 ) {
                    $sequences3[$i] .= $snp_base[$i];
                }
            }
        }

        if ( !( $group_i->is_empty and $group_n->is_empty ) ) {
            ( $d_indel, $d_noindel, $d_bii, $d_bnn, $d_complex )
                = _two_group_D( $group_i, $group_n, $ref_seq, \@sequences, $window_length );

            if ( @sequences2 > 0 and length $sequences2[0] > 0 ) {
                ( $d_indel2, $d_noindel2, $d_bii2, $d_bnn2, $d_complex2 )
                    = _two_group_D( $group_i, $group_n, $ref_seq2, \@sequences2, $window_length );
            }

            if ( @sequences3 > 0 and length $sequences3[0] > 0 ) {
                ( $d_indel3, $d_noindel3, $d_bii3, $d_bnn3, $d_complex3 )
                    = _two_group_D( $group_i, $group_n, $ref_seq3, \@sequences3, $window_length );
            }
        }
        $update_sth->execute(
            $d_indel,    $d_noindel, $d_bii,      $d_bnn,      $d_complex, $d_indel2,
            $d_noindel2, $d_bii2,    $d_bnn2,     $d_complex2, $d_indel3,  $d_noindel3,
            $d_bii3,     $d_bnn3,    $d_complex3, $isw_id
        );
    }

    return;
}

#----------------------------------------------------------#
# Internal Subroutines
#----------------------------------------------------------#
sub _D_indels {
    my $seq_refs = shift;

    my $length = length $seq_refs->[0];
    my ( $d1, $d2, $dc ) = App::Fasops::Common::ref_pair_D($seq_refs);
    for ( $d1, $d2, $dc ) {
        $_ *= $length;
    }

    return ( $d1, $d2, $dc );
}

sub _two_group_D {
    my AlignDB::IntSpan $group1 = shift;
    my AlignDB::IntSpan $group2 = shift;
    my $outgroup_seq            = shift;    # string
    my $sequences               = shift;    # array_ref of strings
    my $window_length = shift || length $outgroup_seq;

    my ( $d_1, $d_2, $d_b11, $d_b22, $d_complex ) = (0) x 5;
    my ( @d1, @d2, @db11, @db22, @dc );

    my ( @g1_seqs, @g2_seqs );
    for ( $group1->elements ) {
        push @g1_seqs, $sequences->[$_];
    }
    for ( $group2->elements ) {
        push @g2_seqs, $sequences->[$_];
    }
    for my $g1_side_seq (@g1_seqs) {
        for my $g2_side_seq (@g2_seqs) {
            my ( $di, $dn, $dc )
                = _D_indels( [ $g1_side_seq, $g2_side_seq, $outgroup_seq, ] );
            push @d1, $di;
            push @d2, $dn;
            push @dc, $dc;
        }
    }
    $d_1 = App::Fasops::Common::mean(@d1) / $window_length;
    $d_2 = App::Fasops::Common::mean(@d2) / $window_length;

    my $i = 0;
    while ( $g1_seqs[ $i + 1 ] ) {
        my $j = $i + 1;
        while ( $g1_seqs[$j] ) {
            my ( $d1, $d2, $dc )
                = _D_indels( [ $g1_seqs[$i], $g1_seqs[$j], $outgroup_seq, ] );
            push @db11, ( $d1 + $d2 );
            push @dc, $dc;
            $j++;
        }
        $i++;
    }
    $i = 0;
    while ( $g2_seqs[ $i + 1 ] ) {
        my $j = $i + 1;
        while ( $g2_seqs[$j] ) {
            my ( $d1, $d2, $dc )
                = _D_indels( [ $g2_seqs[$i], $g2_seqs[$j], $outgroup_seq, ] );
            push @db22, ( $d1 + $d2 );
            push @dc, $dc;
            $j++;
        }
        $i++;
    }
    $d_b11     = App::Fasops::Common::mean(@db11) / $window_length;
    $d_b22     = App::Fasops::Common::mean(@db22) / $window_length;
    $d_complex = App::Fasops::Common::mean(@dc) / $window_length;

    return ( $d_1, $d_2, $d_b11, $d_b22, $d_complex );
}

1;

__END__
