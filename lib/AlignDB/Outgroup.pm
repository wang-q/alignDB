package AlignDB::Outgroup;
use Moose;
use autodie;

use IO::Zlib;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any uniq);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../";
extends qw(AlignDB);

sub _insert_ref_sequences {
    my $self     = shift;
    my $align_id = shift;
    my $info_of  = shift;
    my $ref_name = shift;
    my $ref_seq  = shift;

    my $dbh = $self->dbh;

    my $align_length = length $ref_seq;
    my $align_set    = AlignDB::IntSpan->new("1-$align_length");

    my $ref_index = -1;
    {
        $info_of->{$ref_name}{align_id} = $align_id;
        $info_of->{$ref_name}{seq}      = $ref_seq;
        $info_of->{$ref_name}{gc}       = calc_gc_ratio($ref_seq);
        my $seq_indel_set = find_indel_set($ref_seq);
        my $seq_set       = $align_set->diff($seq_indel_set);
        $info_of->{$ref_name}{runlist} = $seq_set->runlist;
        $info_of->{$ref_name}{length}  = $seq_set->cardinality;

        # for polarize indels
        $info_of->{$ref_name}{indel_set} = $seq_indel_set;
    }

    my $insert = $dbh->prepare(
        q{
        INSERT INTO reference (
            ref_id, seq_id, ref_raw_seq, ref_complex_indel
        )
        VALUES (
            NULL, ?, ?, ?
        )
        }
    );
    my $seq_id = $self->_insert_seq( $info_of->{$ref_name} );
    $insert->execute( $seq_id, $ref_seq, '-' );
    $insert->finish;

    return;
}

sub _polarize_indel {
    my $self          = shift;
    my $align_id      = shift;
    my $ref_indel_set = shift;

    my $dbh = $self->dbh;

    # Complex indels determined without outgroup are still complex
    # polarize clear ones
    my $indel_info_sth = $dbh->prepare(
        q{
        SELECT indel_id, indel_start, indel_end, indel_type, indel_occured
        FROM indel
        WHERE 1 = 1
        AND indel_type <> 'C'
        AND indel_occured <> 'unknown'
        AND align_id = ?
        }
    );
    my $update_indel_sth = $dbh->prepare(
        q{
        UPDATE indel
        SET indel_type = ?,
            indel_occured = ?,
            indel_freq = ?
        WHERE indel_id = ?
        }
    );

    $indel_info_sth->execute($align_id);
    while ( my @row = $indel_info_sth->fetchrow_array ) {
        my ( $indel_id, $indel_start, $indel_end, $old_indel_type,
            $old_indel_occured )
            = @row;

        my $indel_set = AlignDB::IntSpan->new("$indel_start-$indel_end");

        my ( $indel_type, $indel_occured, $indel_freq );

        if ( $ref_indel_set->intersect($indel_set)->is_empty ) {

            # ref NNNN
            #     NNNN
            #     N--N
            $indel_type = 'D';
        }
        elsif ( $ref_indel_set->superset($indel_set) ) {
            my $island = $ref_indel_set->find_islands($indel_start);
            if ( $island->equal($indel_set) ) {

                # ref N--N
                #     NNNN
                #     N--N
                $indel_type = 'I';
            }
            else {
                # ref N--N
                #     NNNN
                #     N-NN
                $indel_type = 'C';
            }
        }
        else {
            # ref N-NN
            #     NNNN
            #     N--N
            $indel_type = 'C';
        }

        if ( $indel_type eq 'C' ) {
            $indel_occured = 'unknown';
            $indel_freq    = -1;
        }
        else {
            #       |    I    |  D  new
            # ------+---------+--------
            #       |     -   |     N
            #  I    | N x N o | N x N x
            #       | - o - x | - o - o
            # ------+---------+--------
            #       |     -   |     N
            #  D    | - x - x | - x - o
            #  old  | N o N o | N o N x
            $indel_occured = $old_indel_occured;
            if ( $old_indel_type eq $indel_type ) {
                $indel_occured =~ tr/ox/xo/;
            }
            $indel_freq = $indel_occured =~ tr/o/o/;
        }

        $update_indel_sth->execute( $indel_type, $indel_occured, $indel_freq,
            $align_id );
    }

    $update_indel_sth->finish;
    $indel_info_sth->finish;

    return;
}

sub _polarize_snp {
    my $self     = shift;
    my $align_id = shift;
    my $ref_seq  = shift;

    my $dbh = $self->dbh;

    my $snp_info_sth = $dbh->prepare(
        q{
        SELECT snp_id, snp_pos, all_bases
        FROM snp
        WHERE 1 = 1
        AND align_id = ?
        }
    );
    my $update_snp_sth = $dbh->prepare(
        q{
        UPDATE snp
        SET ref_base = ?,
            mutant_to = ?,
            snp_freq = ?,
            snp_occured = ?
        WHERE snp_id = ?
        }
    );

    $snp_info_sth->execute($align_id);
    while ( my @row = $snp_info_sth->fetchrow_array ) {
        my ( $snp_id, $snp_pos, $all_bases ) = @row;

        my $ref_base = substr $ref_seq, $snp_pos - 1, 1;

        my @nts = split '', $all_bases;
        my @class;
        for my $nt (@nts) {
            my $class_bool = 0;
            for (@class) {
                if ( $_ eq $nt ) { $class_bool = 1; }
            }
            unless ($class_bool) {
                push @class, $nt;
            }
        }

        my ( $mutant_to, $snp_freq, $snp_occured );

        if ( scalar @class < 2 ) {
            confess "Not a real SNP\n";
        }
        elsif ( scalar @class == 2 ) {
            for my $nt (@nts) {
                if ( $nt eq $ref_base ) {
                    $snp_occured .= 'x';
                }
                else {
                    $snp_occured .= 'o';
                    $snp_freq++;
                    $mutant_to = "$ref_base->$nt";
                }
            }
        }
        else {
            $snp_freq    = -1;
            $mutant_to   = 'Complex';
            $snp_occured = 'unknown';
        }

        # ref_base is not equal to any nts
        if ( $snp_occured eq ( 'o' x ( length $snp_occured ) ) ) {
            $snp_freq    = -1;
            $mutant_to   = 'Complex';
            $snp_occured = 'unknown';
        }

        $update_snp_sth->execute( $ref_base, $mutant_to, $snp_freq,
            $snp_occured, $snp_id );
    }

    return;
}

sub add_align {
    my $self     = shift;
    my $info_of  = shift;
    my $names    = shift;
    my $seq_refs = shift;

    my $dbh = $self->dbh;

    # check align length
    my $align_length = length $seq_refs->[0];
    for ( @{$seq_refs} ) {
        if ( ( length $_ ) != $align_length ) {
            confess "Sequences should have the same length!\n";
        }
    }

    if ( scalar @{$names} != scalar @{$seq_refs} ) {
        warn Dump $names;
        confess "Names and Sequences should have the same number!\n";
    }
    my $seq_count = scalar @{$seq_refs};
    if ( $seq_count < 3 ) {
        confess "Too few sequences [$seq_count]\n";
    }

    # appoint reference/outgroup
    my $ref_name = $names->[-1];
    my $ref_seq  = $seq_refs->[-1];

    # exclude outgroup
    my $ingroup_names = [ @{$names}[ 0 .. $seq_count - 2 ] ];
    my $ingroup_seqs  = [ @{$seq_refs}[ 0 .. $seq_count - 2 ] ];

    #----------------------------#
    # INSERT INTO align
    #----------------------------#
    my $align_id = $self->_insert_align( @{$ingroup_seqs} );
    printf "Prosess align %s in %s %s - %s\n", $align_id,
        $info_of->{ $names->[0] }{chr_name},
        $info_of->{ $names->[0] }{chr_start},
        $info_of->{ $names->[0] }{chr_end};

    #----------------------------#
    # UPDATE align, INSERT INTO sequence, target, queries
    #----------------------------#
    $self->_insert_set_and_sequence( $align_id, $info_of, $ingroup_names,
        $ingroup_seqs );

    #----------------------------#
    # INSERT INTO ref
    #----------------------------#
    $self->_insert_ref_sequences( $align_id, $info_of, $ref_name, $ref_seq );

    #----------------------------#
    # INSERT INTO indel
    #----------------------------#
    $self->_insert_indel($align_id);
    $self->_polarize_indel( $align_id, $info_of->{$ref_name}{indel_set} );

    #----------------------------#
    # INSERT INTO snp
    #----------------------------#
    $self->_insert_snp($align_id);
    $self->_polarize_snp( $align_id, $ref_seq );

    return;
}

sub update_D_values {
    my $self     = shift;
    my $align_id = shift;

    # Get database handle
    my $dbh = $self->dbh;

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

    my $read_sql = $dbh->prepare(
        q{
        SELECT s.ref_base, s.all_bases, w.isw_length, i.indel_occured
        FROM snp s, isw w, indel i
        WHERE 1 = 1
        AND s.isw_id = w.isw_id
        AND w.isw_indel_id = i.indel_id
        AND i.indel_occured <> 'unknown'
        AND w.isw_id = ?
        }
    );

    my $update_sql = $dbh->prepare(
        'UPDATE isw
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
        WHERE isw_id = ?'
    );

ISW: for my $isw_id (@{$isw_id_ref}) {
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

        # removes all mutations on the deepest indel branches
        # without recombinatio
        my $ref_seq3;
        my @sequences3;

        $read_sql->execute($isw_id);
        while ( my @row = $read_sql->fetchrow_array ) {
            $window_length = $row[2];
            my $ref_base = $row[0];

            $ref_seq .= $ref_base;
            my @snp_base    = split '', $row[1];
            my @indel_occur = split '', $row[3];
            my $align_cnt   = scalar @indel_occur;
            for my $i ( 0 .. $align_cnt - 1 ) {
                if ( $indel_occur[$i] eq 'o' ) {
                    $group_i->add($i);
                }
                elsif ( $indel_occur[$i] eq 'x' ) {
                    $group_n->add($i);
                }
                else {
                    die "$indel_occur[$i]\n";
                }
                $sequences[$i] .= $snp_base[$i];
            }

            # find mutations on the deepest branches
            my @group_i = $group_i->elements;
            my @group_n = $group_n->elements;
            my @i_snp   = uniq( @snp_base[@group_i] );
            my @n_snp   = uniq( @snp_base[@group_n] );
            if ( @i_snp == 1 and $i_snp[0] ne $ref_base ) {

                # removes all mutations on the deepest indel branches
                # add bases with recombination events
                if ( any { $_ eq $i_snp[0] } @n_snp ) {
                    $ref_seq3 .= $ref_base;
                    for my $i ( 0 .. $align_cnt - 1 ) {
                        $sequences3[$i] .= $snp_base[$i];
                    }
                }
            }
            elsif ( @n_snp == 1 and $n_snp[0] ne $ref_base ) {

                # removes all mutations on the deepest noindel branches
                # add bases with recombination events
                if ( any { $_ eq $n_snp[0] } @i_snp ) {
                    $ref_seq3 .= $ref_base;
                    for my $i ( 0 .. $align_cnt - 1 ) {
                        $sequences3[$i] .= $snp_base[$i];
                    }
                }
            }
            else {
                $ref_seq2 .= $ref_base;
                for my $i ( 0 .. $align_cnt - 1 ) {
                    $sequences2[$i] .= $snp_base[$i];
                }
                $ref_seq3 .= $ref_base;
                for my $i ( 0 .. $align_cnt - 1 ) {
                    $sequences3[$i] .= $snp_base[$i];
                }
            }
        }

        if ( !( $group_i->empty and $group_n->empty ) ) {
            ( $d_indel, $d_noindel, $d_bii, $d_bnn, $d_complex )
                = _two_group_D( $group_i, $group_n, $ref_seq, \@sequences,
                $window_length );

            if ( @sequences2 > 0 and length $sequences2[0] > 0 ) {
                ( $d_indel2, $d_noindel2, $d_bii2, $d_bnn2, $d_complex2 )
                    = _two_group_D( $group_i, $group_n, $ref_seq2, \@sequences2,
                    $window_length );
            }

            if ( @sequences3 > 0 and length $sequences3[0] > 0 ) {
                ( $d_indel3, $d_noindel3, $d_bii3, $d_bnn3, $d_complex3 )
                    = _two_group_D( $group_i, $group_n, $ref_seq3, \@sequences3,
                    $window_length );
            }
        }
        $update_sql->execute(
            $d_indel,   $d_noindel,  $d_bii,      $d_bnn,
            $d_complex, $d_indel2,   $d_noindel2, $d_bii2,
            $d_bnn2,    $d_complex2, $d_indel3,   $d_noindel3,
            $d_bii3,    $d_bnn3,     $d_complex3, $isw_id
        );
    }

    return;
}

#----------------------------------------------------------#
# Internal Subroutines
#----------------------------------------------------------#
sub _D_indels {
    my ( $ref_seq, $first_seq, $second_seq ) = @_;

    my $length = length $ref_seq;
    my ( $d1, $d2, $dc ) = ref_pair_D( $ref_seq, $first_seq, $second_seq );
    for ( $d1, $d2, $dc ) {
        $_ *= $length;
    }

    return ( $d1, $d2, $dc );
}

sub _two_group_D {
    my $group1        = shift;                      # AlignDB::IntSpan object
    my $group2        = shift;                      # AlignDB::IntSpan object
    my $ref_seq       = shift;                      # string
    my $sequences     = shift;                      # array_ref of strings
    my $window_length = shift || length $ref_seq;

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
                = _D_indels( $ref_seq, $g1_side_seq, $g2_side_seq );
            push @d1, $di;
            push @d2, $dn;
            push @dc, $dc;
        }
    }
    $d_1 = average(@d1) / $window_length;
    $d_2 = average(@d2) / $window_length;

    my $i = 0;
    while ( $g1_seqs[ $i + 1 ] ) {
        my $j = $i + 1;
        while ( $g1_seqs[$j] ) {
            my ( $d1, $d2, $dc )
                = _D_indels( $ref_seq, $g1_seqs[$i], $g1_seqs[$j] );
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
                = _D_indels( $ref_seq, $g2_seqs[$i], $g2_seqs[$j] );
            push @db22, ( $d1 + $d2 );
            push @dc, $dc;
            $j++;
        }
        $i++;
    }
    $d_b11     = average(@db11) / $window_length;
    $d_b22     = average(@db22) / $window_length;
    $d_complex = average(@dc) / $window_length;

    return ( $d_1, $d_2, $d_b11, $d_b22, $d_complex );
}

1;

__END__
