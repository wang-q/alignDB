package AlignDB::Multi;
use Moose;
use Carp;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any uniq);
use Statistics::Lite qw(mean);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Util qw(:all);

extends qw(AlignDB);

sub add_align {
    my $self        = shift;
    my $ref_info    = shift;
    my $target_info = shift;
    my @query_info  = @_;

}

sub parse_fasta_file {
    my $self       = shift;
    my $fasta_file = shift;

    my $fasta_qr = qr{
        id(\d+)         # target taxon id
        [\_\-]          # spacer
        ((?:chr)?\w+)   # target chr name
        [\_\-]          # spacer
        (\d+)           # target chr start
        [\_\-]          # spacer
        (\d+)           # target chr end
        \.fa            # suffix
    }xi;
    $fasta_file =~ $fasta_qr;

    my $target_taxon_id  = $1;
    my $target_chr_name  = $2;
    my $target_chr_start = $3;
    my $target_chr_end   = $4;

    my $target_chr_id
        = $self->get_chr_id_hash($target_taxon_id)->{$target_chr_name};

    # Get database handle
    my $dbh = $self->dbh;

    my $align_insert_sql = $dbh->prepare(
        q{
        INSERT INTO align (
            align_id, 
            align_length, align_comparables, align_identities,
            align_differences, align_gaps, align_ns,
            align_error, align_pi
        )
        VALUES (
            NULL,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?
        )
        }
    );

    my $align_update_sql = $dbh->prepare(
        q{
        UPDATE align
        SET align_indels = ?,
            align_target_gc = ?,
            align_average_gc = ?,
            align_comparable_runlist = ?,
            align_indel_runlist = ?
        WHERE align_id = ?
        }
    );
    my ( $align_indels, $align_target_gc, $align_average_gc,
        $align_comparable_runlist, $align_indel_runlist );

    # read in fasta file
    open my $fasta_fh, '<', $fasta_file or die $!;
    my $fasta_content = do { local $/; <$fasta_fh> };
    close $fasta_fh;

    my @fasta_heads = ( $fasta_content =~ />(\w+)\r?\n/g );
    my @fasta_lines = ( $fasta_content =~ />\w+\r?\n([\-\w]+)\r?\n/g );
    my $seq_count   = scalar @fasta_heads;

    my $align_length = length $fasta_lines[0];
    for (@fasta_lines) {
        if ( ( length $_ ) ne $align_length ) {
            croak "sequences should have the same length!\n";
        }
    }

    $align_target_gc  = calc_gc_ratio( $fasta_lines[1] );
    $align_average_gc = calc_gc_ratio( @fasta_lines[ 1 .. $#fasta_lines ] );

    my $align_stat = multi_seq_stat(@fasta_lines);
    $align_insert_sql->execute(
        $align_stat->[0], $align_stat->[1], $align_stat->[2],
        $align_stat->[3], $align_stat->[4], $align_stat->[5],
        $align_stat->[6], $align_stat->[7]
    );
    $align_insert_sql->finish;
    my $align_id = $self->last_insert_id;

    #----------------------------------------------------------#
    # INSERT INTO sequence, target, query
    #----------------------------------------------------------#
    {
        my @seq_infos;    # store info of each seqs
        for my $i ( 0 .. $seq_count - 1 ) {
            my $seq_info = {};
            $seq_info->{align_id} = $align_id;
            $seq_info->{name}     = $fasta_heads[$i];
            $seq_info->{seq}      = $fasta_lines[$i];
            $seq_info->{length}   = $align_length;
            $seq_info->{gc}       = calc_gc_ratio( $fasta_lines[$i] );

            my $seq_set   = AlignDB::IntSpan->new("1-$align_length");
            my $indel_set = find_indel_set( $fasta_lines[$i] );
            $seq_set = $seq_set->diff($indel_set);
            $seq_info->{runlist} = $seq_set->runlist;

            push @seq_infos, $seq_info;
        }

        # ref
        my $ref_insert = $dbh->prepare(
            q{
            INSERT INTO reference (
                ref_id, seq_id, ref_raw_seq, ref_complex_indel
            )
            VALUES (
                NULL, ?, ?, ?
            )
            }
        );
        my $ref_seq_id = $self->_insert_seq( $seq_infos[0] );
        $ref_insert->execute( $ref_seq_id, $seq_infos[0]->{seq}, '-' );
        $ref_insert->finish;

        # target
        my $target_insert = $dbh->prepare(
            q{
            INSERT INTO target ( target_id, seq_id )
            VALUES ( NULL, ? )
            }
        );
        $seq_infos[1]->{chr_id}     = $target_chr_id;
        $seq_infos[1]->{chr_start}  = $target_chr_start;
        $seq_infos[1]->{chr_end}    = $target_chr_end;
        $seq_infos[1]->{chr_strand} = '+';
        my $target_seq_id = $self->_insert_seq( $seq_infos[1] );
        $target_insert->execute($target_seq_id);
        $target_insert->finish;

        # and queries
        my $query_insert = $dbh->prepare(
            q{
            INSERT INTO query ( query_id, seq_id, query_strand )
            VALUES ( NULL, ?, ? )
            }
        );
        for my $i ( 2 .. $seq_count - 1 ) {
            my $query_seq_id = $self->_insert_seq( $seq_infos[$i] );
            $query_insert->execute( $query_seq_id, '+', );
        }
        $query_insert->finish;
    }

    # find indel
    my $indel_set = AlignDB::IntSpan->new;
    my $reference = shift @fasta_lines;
    for my $cur_line (@fasta_lines) {
        my $cur_indel_set = find_indel_set($cur_line);
        $indel_set->merge($cur_indel_set);
    }
    unshift @fasta_lines, $reference;

    my $align_set = AlignDB::IntSpan->new("1-$align_length");

    my @indel_spans = $indel_set->spans;
    $align_indels = scalar @indel_spans;
    my $comparable_set = $align_set->diff($indel_set);
    $align_indel_runlist      = $indel_set->runlist;
    $align_comparable_runlist = $comparable_set->runlist;

    $align_update_sql->execute( $align_indels, $align_target_gc,
        $align_average_gc, $align_comparable_runlist, $align_indel_runlist,
        $align_id );

    # insert indels
    if ( $align_indels > 0 ) {
        my $indel_insert_sql = $dbh->prepare(
            q{
            INSERT INTO indel (
                indel_id, align_id, prev_indel_id,
                indel_start, indel_end, indel_length,
                indel_seq,  indel_gc, indel_freq,
                indel_occured, indel_type
            )
            VALUES (
                NULL, ?, ?,
                ?, ?, ?,
                ?, ?, ?,
                ?, ?
            )
            }
        );

        my $prev_indel_id = 0;
        for my $cur_indel (@indel_spans) {
            my ( $indel_start, $indel_end ) = @{$cur_indel};
            my $indel_length = $indel_end - $indel_start + 1;

            my @indel_seqs;
            for my $cur_line (@fasta_lines) {
                push @indel_seqs,
                    ( substr $cur_line, $indel_start - 1, $indel_length );
            }

            my $indel_seq = '';
            my $indel_type;
            my @indel_class;
            for my $seq (@indel_seqs) {
                unless ( $seq =~ /-/ || $seq =~ /N/i ) {
                    if ( $indel_seq =~ /-/ ) { die "aaaa$seq\n"; }
                    $indel_seq = $seq;
                }
                my $class_bool = 0;
                for (@indel_class) {
                    if ( $_ eq $seq ) { $class_bool = 1; }
                }
                unless ($class_bool) {
                    push @indel_class, $seq;
                }
            }

            if ( scalar @indel_class < 2 ) {
                die "no indel!\n";
            }
            elsif ( scalar @indel_class > 2 ) {
                $indel_type = 'C';
            }
            my $ref_seq = shift @indel_seqs;
            unless ($indel_type) {
                if ( $ref_seq eq ( '-' x ( length $ref_seq ) ) ) {
                    $indel_type = 'I';
                }
                elsif ( !( $ref_seq =~ /-/ ) ) {
                    $indel_type = 'D';
                }
                else {
                    die $cur_indel;
                }
            }
            my $indel_frequency = 0;
            my $indel_occured;
            if ( $indel_type eq 'C' ) {
                $indel_frequency = -1;
                $indel_occured   = 'unknown';
            }
            else {
                for (@indel_seqs) {
                    if ( $ref_seq ne $_ ) {
                        $indel_frequency++;
                        $indel_occured .= 'o';
                    }
                    else {
                        $indel_occured .= 'x';
                    }
                }
            }
            if ( $indel_occured eq ( 'o' x ( length $indel_occured ) ) ) {
                my $drop_set
                    = AlignDB::IntSpan->new("$indel_start-$indel_end");
                $indel_set = $indel_set->diff($drop_set);
                next;
            }
            my $indel_gc = calc_gc_ratio($indel_seq);

            $indel_insert_sql->execute(
                $align_id,  $prev_indel_id,   $indel_start,
                $indel_end, $indel_length,    $indel_seq,
                $indel_gc,  $indel_frequency, $indel_occured,
                $indel_type,
            );
            ($prev_indel_id) = $self->last_insert_id;
        }
    }

    #intevals
    if ( $align_indels > 1 ) {
        my $fetch_indel_id_sql = $dbh->prepare(
            q{
            SELECT indel_id, prev_indel_id
            FROM indel
            WHERE align_id = ?
              AND indel_start = ?
            }
        );
        my $isw_insert_sql = $dbh->prepare(
            q{
            INSERT INTO isw (
                isw_id, indel_id, prev_indel_id, isw_indel_id,
                isw_start, isw_end, isw_length, isw_type,
                isw_distance, isw_density, isw_differences,
                isw_pi, isw_target_gc, isw_average_gc
            )
            VALUES (
                NULL, ?, ?, ?,
                ?, ?, ?, ?,
                ?, ?, ?,
                ?, ?, ?
            )
            }
        );
        my $snp_insert_sql = $dbh->prepare(
            q{
            INSERT INTO snp (
                snp_id, isw_id, align_id, snp_pos,
                target_base, ref_base, all_bases
            )
            VALUES (
                NULL, ?, ?, ?,
                ?, ?, ?
            )
            }
        );

        my $window_maker = $self->window_maker;

        my @interval_spans = $comparable_set->spans;
        shift @interval_spans;
        pop @interval_spans;

        for my $cur_inteval (@interval_spans) {
            my ( $interval_start, $interval_end ) = @{$cur_inteval};
            my $indel_start = $interval_end + 1;

            $fetch_indel_id_sql->execute( $align_id, $indel_start );
            my ( $indel_id, $prev_indel_id, )
                = $fetch_indel_id_sql->fetchrow_array;

            my @isws
                = $window_maker->interval_window( $align_set, $interval_start,
                $interval_end );

            for my $isw (@isws) {
                my $isw_set    = $isw->{set};
                my $isw_start  = $isw_set->min;
                my $isw_end    = $isw_set->max;
                my $isw_length = $isw_end - $isw_start + 1;

                my $isw_distance = $isw->{distance};
                my $isw_density  = $isw->{density};
                my $isw_type     = $isw->{type};
                my @isw_seq;

                for my $cur_line (@fasta_lines) {
                    push @isw_seq,
                        ( substr $cur_line, $isw_start - 1, $isw_length );
                }
                my $isw_ref_seq    = shift @isw_seq;
                my $isw_average_gc = &calc_gc_ratio(@isw_seq);
                my $isw_target_gc  = &calc_gc_ratio( $isw_seq[0] );
                my $isw_stat       = multi_seq_stat(@isw_seq);
                my $isw_difference = $isw_stat->[3];
                my $isw_pi         = $isw_stat->[7];

                my $isw_indel_id;
                if ( $isw_type eq 'L' ) {
                    $isw_indel_id = $prev_indel_id;
                }
                elsif ( $isw_type eq 'R' ) {
                    $isw_indel_id = $indel_id;
                }
                elsif ( $isw_type eq 'S' ) {
                    $isw_indel_id = $prev_indel_id;
                }

                $isw_insert_sql->execute(
                    $indel_id,       $prev_indel_id, $isw_indel_id,
                    $isw_start,      $isw_end,       $isw_length,
                    $isw_type,       $isw_distance,  $isw_density,
                    $isw_difference, $isw_pi,        $isw_target_gc,
                    $isw_average_gc
                );
                my $isw_id = $self->last_insert_id;

                for my $snp_pos ( 0 .. $isw_length - 1 ) {
                    my @bases;
                    for my $cur_isw_seq (@isw_seq) {
                        my $cur_nuc = substr $cur_isw_seq, $snp_pos, 1;
                        push @bases, $cur_nuc;
                    }
                    next if any { $_ !~ /[agct]/i } @bases;
                    my $class = scalar uniq(@bases);
                    if ( $class > 1 ) {
                        my $snp_position = $snp_pos + $isw_start;
                        my $snp_ref_base = substr $isw_ref_seq, $snp_pos, 1;
                        my $snp_target_base = substr $isw_seq[0], $snp_pos, 1;
                        my $snp_other_base = join '', @bases;

                        $snp_insert_sql->execute(
                            $isw_id,       $align_id,
                            $snp_position, $snp_target_base,
                            $snp_ref_base, $snp_other_base
                        );
                    }
                }
            }
        }
    }

    return;
}

sub _insert_seq {
    my $self     = shift;
    my $seq_info = shift;

    croak "Pass a seq to this method!\n" if !defined $seq_info->{seq};

    for my $key (qw{chr_id chr_start chr_end chr_strand length gc runlist}) {
        if ( !defined $seq_info->{$key} ) {
            $seq_info->{$key} = undef;
        }
    }

    # Get database handle
    my $dbh = $self->dbh;

    my $seq_insert = $dbh->prepare(
        q{
        INSERT INTO sequence (
            seq_id, chr_id, align_id, chr_start, chr_end,
            chr_strand, seq_length, seq_seq, seq_gc, seq_runlist
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?, ?, ?
        )
        }
    );

    $seq_insert->execute(
        $seq_info->{chr_id},     $seq_info->{align_id},
        $seq_info->{chr_start},  $seq_info->{chr_end},
        $seq_info->{chr_strand}, $seq_info->{length},
        $seq_info->{seq},        $seq_info->{gc},
        $seq_info->{runlist},
    );

    my $seq_id = $self->last_insert_id;

    return $seq_id;
}

sub update_misc {
    my $self = shift;

    $self->update_snp;
    $self->update_indel;
    $self->update_D_values;

    return;
}

sub update_snp {
    my $self = shift;

    # Get database handle
    my $dbh = $self->dbh;

    my $snp_info_sth = $dbh->prepare(
        'SELECT snp_id, ref_base, all_bases
        FROM snp'
    );
    my $update_snp_sth = $dbh->prepare(
        'UPDATE snp
        SET mutant_to = ?, snp_freq = ?, snp_occured = ?
        WHERE snp_id = ?'
    );

    $snp_info_sth->execute;
    while ( my @row = $snp_info_sth->fetchrow_array ) {
        my ( $snp_id, $ref_base, $all_bases ) = @row;

        my @all_nucs = split '', $all_bases;
        my @class;
        for my $nuc (@all_nucs) {
            my $class_null = 0;
            for (@class) {
                if ( $_ eq $nuc ) { $class_null = 1; }
            }
            unless ($class_null) {
                push @class, $nuc;
            }
        }

        my ( $mutant_to, $snp_freq, $snp_occured );

        if ( scalar @class < 2 ) {
            croak "Not a real SNP\n";
        }
        elsif ( scalar @class == 2 ) {
            for my $nuc (@all_nucs) {
                if ( $nuc eq $ref_base ) {
                    $snp_occured .= 'x';
                }
                else {
                    $snp_occured .= 'o';
                    $snp_freq++;
                    $mutant_to = "$ref_base->$nuc";
                }
            }
        }
        else {
            $snp_freq    = -1;
            $mutant_to   = 'C';
            $snp_occured = 'unknown';
        }
        if ( $snp_occured eq ( 'o' x ( length $snp_occured ) ) ) {
            $snp_freq    = -1;
            $mutant_to   = 'C';
            $snp_occured = 'unknown';
        }

        $update_snp_sth->execute( $mutant_to, $snp_freq, $snp_occured,
            $snp_id );
    }

    return;
}

sub update_indel {
    my $self = shift;

    # Get database handle
    my $dbh = $self->dbh;

    my $align_info_sth = $dbh->prepare(
        'SELECT a.align_id, a.align_length, s.seq_seq
        FROM align a, target t, sequence s
        WHERE a.align_id = s.align_id
        AND s.seq_id = t.seq_id'
    );
    my $indel_info_sth = $dbh->prepare(
        'SELECT indel_id, indel_start, indel_end,
                indel_length, indel_seq
        FROM indel
        WHERE align_id = ?'
    );
    my $update_indel_sth = $dbh->prepare(
        'UPDATE indel
        SET left_extand = ?, right_extand = ?
        WHERE indel_id = ?'
    );
    my $update_slippage_sth = $dbh->prepare(
        'UPDATE indel
        SET indel_slippage = ?
        WHERE indel_id = ?'
    );

    #motif-repeat parameters
    my $min_reps = {
        1 => 4,    # mononucl. with >= 4 repeats
    };

    $align_info_sth->execute;
    while ( my @row = $align_info_sth->fetchrow_array ) {
        my ( $align_id, $align_length, $target_seq ) = @row;

        my $align_start = 1;
        my $align_end   = $align_length;

        my $indel_info;
        my @indel_ids;
        $indel_info_sth->execute($align_id);
        while ( my @row = $indel_info_sth->fetchrow_array ) {
            my ( $id, $start, $end, $length, $seq ) = @row;
            $indel_info->{$id}->{start}  = $start;
            $indel_info->{$id}->{end}    = $end;
            $indel_info->{$id}->{length} = $length;
            $indel_info->{$id}->{seq}    = $seq;
            push @indel_ids, $id;
        }

        #----------------------------#
        # left_ and right_extand
        #----------------------------#
        for my $indel_id (@indel_ids) {
            my $former_end;
            my $latter_start;
            if ( $indel_info->{ $indel_id - 1 }->{end} ) {
                $former_end = $indel_info->{ $indel_id - 1 }->{end} + 1;
            }
            else {
                $former_end = $align_start;
            }
            if ( $indel_info->{ $indel_id + 1 }->{start} ) {
                $latter_start = $indel_info->{ $indel_id + 1 }->{start} - 1;
            }
            else {
                $latter_start = $align_end;
            }
            my $left_extand = $indel_info->{$indel_id}->{start} - $former_end;
            if ( $left_extand < 0 ) {
                my $start = $indel_info->{$indel_id}->{start};
                my $end   = $indel_info->{ $indel_id - 1 }->{end};
                croak "$indel_id,$end,$start\n";
            }
            my $right_extand
                = $latter_start - $indel_info->{$indel_id}->{end};

            $indel_info->{$indel_id}{left_extand}  = $left_extand;
            $indel_info->{$indel_id}{right_extand} = $right_extand;
            $update_indel_sth->execute( $left_extand, $right_extand,
                $indel_id );
        }

        #----------------------------#
        # indel slippage
        #----------------------------#
        for my $indel_id (@indel_ids) {
            my $indel_start    = $indel_info->{$indel_id}{length};
            my $indel_end      = $indel_info->{$indel_id}{end};
            my $indel_length   = $indel_info->{$indel_id}{length};
            my $indel_seq      = $indel_info->{$indel_id}{seq};
            my $left_extand    = $indel_info->{$indel_id}{left_extand};
            my $right_extand   = $indel_info->{$indel_id}{right_extand};
            my $indel_slippate = 0;

            next unless $indel_seq;

            if ( exists $min_reps->{$indel_length} ) {
                my $reps         = $min_reps->{$indel_length};
                my $fland_length = $indel_length * $reps;

                my $left_flank = " ";    # avoid warning from $flank
                if ( $fland_length <= $left_extand ) {
                    $left_flank
                        = substr( $target_seq,
                        $indel_start - $fland_length - 1,
                        $fland_length );
                }

                my $right_flank = " ";
                if ( $fland_length <= $right_extand ) {
                    $right_flank
                        = substr( $target_seq, $indel_end, $fland_length );
                }

                my $flank = $left_flank . $indel_seq . $right_flank;
                my $regex = $indel_seq . "{$reps,}";
                $regex = quotemeta $regex;

                if ( $flank =~ /$regex/ ) {
                    $indel_slippate = 1;
                }
            }
            else {

                # indel 23-28, length 6: substr 17-22
                # seq start at 1 and string start at 0, so minus 1
                # substr(..., 16, 6)
                my $left_flank;
                if ( $indel_length <= $left_extand ) {
                    $left_flank
                        = substr( $target_seq,
                        $indel_start - $indel_length - 1,
                        $indel_length );
                }

                # indel 23-28, length 6: substr 29-34
                # substr(..., 28, 6)
                my $right_flank;
                if ( $indel_length <= $right_extand ) {
                    $right_flank
                        = substr( $target_seq, $indel_end, $indel_length );
                }

                if ( $left_flank and $indel_seq eq $left_flank ) {
                    $indel_slippate = 1;
                }
                elsif ( $right_flank and $indel_seq eq $right_flank ) {
                    $indel_slippate = 1;
                }
            }
            $update_slippage_sth->execute( $indel_slippate, $indel_id );
        }
    }

    return;
}

sub update_D_values {
    my $self = shift;

    # Get database handle
    my $dbh = $self->dbh;

    my $isw_id_ref = $dbh->selectcol_arrayref('SELECT isw_id FROM isw');
    my @isw_ids    = @{$isw_id_ref};

    my $read_sql = $dbh->prepare(
        'SELECT s.ref_base, s.all_bases, w.isw_length, i.indel_occured
        FROM snp s, isw w, indel i
        WHERE s.isw_id = w.isw_id
        AND w.isw_indel_id = i.indel_id
        AND w.isw_id = ?'
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

ISW: for my $isw_id (@isw_ids) {
        my $length;
        my $window_length;
        my ( $d_indel, $d_noindel, $d_bii, $d_bnn, $d_complex ) = (0) x 5;
        my ( $d_indel2, $d_noindel2, $d_bii2, $d_bnn2, $d_complex2 )
            = (0) x 5;
        my ( $d_indel3, $d_noindel3, $d_bii3, $d_bnn3, $d_complex3 )
            = (0) x 5;

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
            $length        = $row[1];
            next ISW if $row[3] eq 'unknown';
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
                    = _two_group_D( $group_i, $group_n, $ref_seq2,
                    \@sequences2, $window_length );
            }

            if ( @sequences3 > 0 and length $sequences3[0] > 0 ) {
                ( $d_indel3, $d_noindel3, $d_bii3, $d_bnn3, $d_complex3 )
                    = _two_group_D( $group_i, $group_n, $ref_seq3,
                    \@sequences3, $window_length );
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

sub get_seqs {
    my $self     = shift;
    my $align_id = shift;
    my $flag;
    if ( defined $self->caching_id ) {
        if ( $self->caching_id == $align_id ) {
            $flag = 0;
        }
        else {
            $flag = 1;
        }
    }
    else {
        $flag = 1;
    }

    if ($flag) {

        # Get database handle
        my $dbh = $self->dbh;

        my $target_sth = $dbh->prepare(
            q{
            SELECT s.seq_seq
            FROM target t, sequence s
            WHERE s.align_id = ?
            AND s.seq_id = t.seq_id
            }
        );

        my $query_sth = $dbh->prepare(
            q{
            SELECT s.seq_seq
            FROM query q, sequence s
            WHERE s.align_id = ?
            AND s.seq_id = q.seq_id
            }
        );

        $target_sth->execute($align_id);
        my ($target_seq) = $target_sth->fetchrow_array;
        $target_sth->finish;

        $query_sth->execute($align_id);
        my (@query_seqs) = $query_sth->fetchrow_array;
        $query_sth->finish;

        $self->{caching_id} = $align_id;
        $self->{caching_seqs} = [ $target_seq, @query_seqs ];
    }

    return $self->caching_seqs;
}

sub get_slice_stat {
    my $self     = shift;
    my $align_id = shift;
    my $set      = shift;

    my $seqs_ref = $self->get_seqs($align_id);

    my @seq_slices = map { $set->substr_span($_) } @$seqs_ref;

    my $result = multi_seq_stat(@seq_slices);

    return $result;
}

=method insert_window

      Usage : $self->insert_window(
            :     $align_id, $window_set, $internal_indel,
            : );
    Purpose : Add a window into malignDB4
    Returns : Int, window_id
            : each member is a hash_ref
 Parameters : $align_id       : align_id
            : $window_set     : AlignDB::IntSpan object
            : $internal_indel : count all indels in this set (flag)
     Throws : no exceptions
   Comments : none 
   See Also : n/a

=cut

sub insert_window {
    my $self           = shift;
    my $align_id       = shift;
    my $window_set     = shift;
    my $internal_indel = shift;

    # Get database handle
    my $dbh = $self->dbh;

    my $window_sql = qq{
        INSERT INTO window (
            window_id, align_id, window_start, window_end, window_length,
            window_runlist, window_comparables, window_identities,
            window_differences, window_indel, window_pi,
            window_target_gc, window_average_gc,
            window_coding, window_repeats
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?,
            NULL, NULL
        )
    };
    my $window_insert = $dbh->prepare($window_sql);

    my $window_start   = $window_set->min;
    my $window_end     = $window_set->max;
    my $window_length  = $window_set->cardinality;
    my $window_runlist = $window_set->runlist;

    # do or do not count internal indels within window_set
    # $set_indel is equal to $window_span - 1
    my $window_indel;
    if ($internal_indel) {
        my ( $set_indel, $real_indel )
            = $self->get_slice_indel( $align_id, $window_set );
        $window_indel = $set_indel + $real_indel;
    }
    else {
        $window_indel = $window_set->span_size - 1;
    }

    my $window_stat = $self->get_slice_stat( $align_id, $window_set );

    $window_insert->execute(
        $align_id,         $window_start,     $window_end,
        $window_length,    $window_runlist,   $window_stat->[1],
        $window_stat->[2], $window_stat->[3], $window_indel,
        $window_stat->[7], $window_stat->[8], $window_stat->[9],
    );

    my $cur_window_id = $self->last_insert_id;

    return $cur_window_id;
}

sub get_align_ids {
    my $self = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT a.align_id
        FROM align a
    };

    my $align_ids = $dbh->selectcol_arrayref($query);

    return $align_ids;
}

sub get_align_ids_of_chr {
    my $self   = shift;
    my $chr_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT a.align_id
        FROM align a, target t, sequence s
        WHERE a.align_id = s.align_id
        AND t.seq_id = s.seq_id
        AND s.chr_id = ?
        ORDER BY a.align_id
    };

    my @align_ids;

    my $sth = $dbh->prepare($query);
    $sth->execute($chr_id);
    while ( my @row = $sth->fetchrow_array ) {
        my ($align_id) = @row;
        push @align_ids, $align_id;
    }

    return \@align_ids;
}

sub get_target_info {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT c.chr_name,
               s.chr_start,
               s.chr_end,
               s.seq_runlist
        FROM sequence s
        inner join target t on s.seq_id = t.seq_id
        inner join chromosome c on s.chr_id = c.chr_id
        WHERE 1 = 1
        AND s.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my ( $chr_name, $chr_start, $chr_end, $target_runlist )
        = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_name, $chr_start, $chr_end, $target_runlist );
}

#----------------------------------------------------------#
# Internal Subroutines
#----------------------------------------------------------#
sub _D_indels {
    my ( $ref_seq, $first_seq, $second_seq ) = @_;

    my $length = length $ref_seq;
    my ( $d1, $d2, $dc ) = &ref_pair_D( $ref_seq, $first_seq, $second_seq );
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
                = &_D_indels( $ref_seq, $g1_side_seq, $g2_side_seq );
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
                = &_D_indels( $ref_seq, $g1_seqs[$i], $g1_seqs[$j] );
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
                = &_D_indels( $ref_seq, $g2_seqs[$i], $g2_seqs[$j] );
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
