package AlignDB;
use Moose;
use DBI;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::DeltaG;
use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Util qw(:all);

has 'mysql'  => ( is => 'ro', isa => 'Str' );    # e.g. 'alignDB:202.119.43.5'
has 'server' => ( is => 'ro', isa => 'Str' );    # e.g. '202.119.43.5'
has 'db'     => ( is => 'ro', isa => 'Str' );    # e.g. 'alignDB'
has 'user'   => ( is => 'ro', isa => 'Str' );    # database username
has 'passwd' => ( is => 'ro', isa => 'Str' );    # database password
has 'dbh'    => ( is => 'ro', isa => 'Object' ); # store database handle here
has 'dG_calc'      => ( is => 'ro', isa => 'Object' ); # dG calculator
has 'window_maker' => ( is => 'ro', isa => 'Object' ); # sliding windows maker
has 'insert_dG' => ( is => 'ro', isa => 'Bool', default => 0 );
has 'insert_ssw' => ( is => 'ro', isa => 'Bool', default => 0 );
has 'threshold'  => ( is => 'ro', isa => 'Int',  default => 10_000 );

has 'caching_id' => ( is => 'ro', isa => 'Int' );      # caching seqs
has 'caching_seqs' =>
    ( is => 'ro', isa => 'ArrayRef[Str]', default => sub { [] } );
has 'caching_refs' => ( is => 'ro', isa => 'Ref' );

sub BUILD {
    my $self = shift;

    # Connect to the alignDB database
    if ( $self->mysql ) {
        my ( $server, $db ) = split ':', $self->mysql;
        $self->{server} ||= $server;
        $self->{db}     ||= $db;
    }
    elsif ( $self->server and $self->db ) {
        $self->{mysql} = $self->db . ':' . $self->server;
    }
    else {
        confess "You should provide either mysql or db-server\n";
    }

    my $mysql  = $self->mysql;
    my $user   = $self->user;
    my $passwd = $self->passwd;
    my $server = $self->server;
    my $db     = $self->db;

    my $dbh;
    $dbh = DBI->connect( "dbi:mysql:$mysql", $user, $passwd )
        or confess "Cannot connect to MySQL database at $mysql";
    $self->{dbh} = $dbh;

    my $dG_calc = AlignDB::DeltaG->new;
    $self->{dG_calc} = $dG_calc;

    my $window_maker = AlignDB::Window->new;
    $self->{window_maker} = $window_maker;

    return;
}

sub _tvsq_id {
    my ( $self, $target_taxon_id, $target_name, $query_taxon_id, $query_name,
        $ref_taxon_id, $ref_name )
        = @_;

    my $dbh = $self->dbh;

    my $tvsq_id;
    my $tvsq = $dbh->prepare(
        'SELECT tvsq_id FROM tvsq
        WHERE target_taxon_id = ?
        AND query_taxon_id = ?'
    );
    $tvsq->execute( $target_taxon_id, $query_taxon_id );
    while ( my @row = $tvsq->fetchrow_array ) {
        ($tvsq_id) = @row;
    }
    $tvsq->finish;

    unless ( defined $tvsq_id ) {
        my $tvsq_insert = $dbh->prepare(
            'INSERT INTO tvsq (
                tvsq_id, target_taxon_id, target_name,
                query_taxon_id, query_name,
                ref_taxon_id, ref_name
            )
            VALUES (
                NULL, ?, ?,
                ?, ?,
                ?, ?
            )'
        );
        $tvsq_insert->execute(
            $target_taxon_id, $target_name,  $query_taxon_id,
            $query_name,      $ref_taxon_id, $ref_name
        );
        ($tvsq_id) = $self->last_insert_id;
    }

    return $tvsq_id;
}

sub _insert_align {
    my ( $self, $tvsq_id, $target_seq_ref, $query_seq_ref ) = @_;

    my $dbh = $self->dbh;

    my $align_insert = $dbh->prepare(
        'INSERT INTO align (
            align_id, tvsq_id, align_length,
            comparable_bases, identities, differences,
            gaps, ns, align_error,
            pi, align_target_gc, align_average_gc,
            comparable_runlist, indel_runlist
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            NULL, NULL
        )'
    );

    my $result = pair_seq_stat( $$target_seq_ref, $$query_seq_ref );

    $align_insert->execute(
        $tvsq_id,     $result->[0], $result->[1], $result->[2],
        $result->[3], $result->[4], $result->[5], $result->[6],
        $result->[7], $result->[8], $result->[9],
    );
    $align_insert->finish;

    my $align_id = $self->last_insert_id;

    return $align_id;
}

sub _insert_isw {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my ( $align_set, $comparable_set, $indel_set )
        = $self->get_sets($align_id);

    # indel_id & foregoing_indel_id
    my $fetch_indel_id_isw = $dbh->prepare(
        'SELECT indel_id, foregoing_indel_id
        FROM indel
        WHERE align_id = ?'
    );
    $fetch_indel_id_isw->execute($align_id);

    # indel_end
    my $fetch_foregoing_indel_end = $dbh->prepare(
        'SELECT indel_end
        FROM indel
        WHERE indel_id = ?'
    );

    # foregoing_indel_start
    my $fetch_indel_start = $dbh->prepare(
        'SELECT indel_start
        FROM indel
        WHERE indel_id = ?'
    );

    # prepare isw_insert
    my $isw_insert = $dbh->prepare(
        'INSERT INTO isw (
            isw_id, indel_id, foregoing_indel_id,
            isw_start, isw_end, isw_length, 
            isw_type, isw_distance, isw_density,
            isw_pi, isw_target_gc, isw_average_gc,
            isw_target_dG, isw_query_dG,
            isw_d_indel, isw_d_noindel, isw_d_complex
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?,
            NULL, NULL, NULL
        )'
    );

    while ( my $ref = $fetch_indel_id_isw->fetchrow_hashref ) {
        my $indel_id           = $ref->{indel_id};
        my $foregoing_indel_id = $ref->{foregoing_indel_id};

        # bypass the first indel
        if ( $foregoing_indel_id == 0 ) {
            next;
        }

        my ( $interval_start, $interval_end, $interval_length ) = ( 0, 0, 0 );

        $fetch_foregoing_indel_end->execute($foregoing_indel_id);
        ($interval_start) = $fetch_foregoing_indel_end->fetchrow_array;
        $interval_start++;

        $fetch_indel_start->execute($indel_id);
        ($interval_end) = $fetch_indel_start->fetchrow_array;
        $interval_end--;

        if ( $interval_start > $interval_end ) {
            print Dump(
                {   interval_start  => $interval_start,
                    interval_end    => $interval_end,
                    interval_length => $interval_length,
                }
            );
            print "start $interval_start > end $interval_end.\n";
            next;
        }

        my $window_maker = $self->window_maker;

        my @isws
            = $window_maker->interval_window( $align_set, $interval_start,
            $interval_end );

        foreach my $isw (@isws) {
            my $isw_set    = $isw->{set};
            my $isw_start  = $isw_set->min;
            my $isw_end    = $isw_set->max;
            my $isw_length = $isw_end - $isw_start + 1;

            $isw_insert->execute(
                $indel_id,
                $foregoing_indel_id,
                $isw_start,
                $isw_end,
                $isw_length,
                $isw->{type},
                $isw->{distance},
                $isw->{density},
                ( $self->get_slice_pi_gc_dG( $align_id, $isw_set ) )
            );
        }
    }

    $isw_insert->finish;
    $fetch_foregoing_indel_end->finish;
    $fetch_indel_start->finish;
    $fetch_indel_id_isw->finish;

    return;

}

sub _insert_snp {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my ( $target_seq, $query_seq ) = @{ $self->get_seqs($align_id) };
    my ( $align_set, $comparable_set, $indel_set )
        = $self->get_sets($align_id);
    my $ref_info = $self->caching_refs;

    my %snp_site = %{ pair_snp_sites( $target_seq, $query_seq ) };
    my $snp_insert = $dbh->prepare(
        'INSERT INTO snp (
            snp_id, isw_id, align_id, snp_pos,
            target_base, query_base, ref_base, snp_occured
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?, ?, ?
        )'
    );

    # isw_id
    my $fetch_isw_id = $dbh->prepare(
        'SELECT isw_id
        FROM isw, indel
        WHERE isw.indel_id = indel.indel_id
        AND indel.align_id = ?
        AND isw.isw_start <= ?
        AND isw.isw_end >= ?'
    );

    # %snp_site keys are snp positions,
    #   where is generated at align section
    foreach ( sort { $a <=> $b } keys %snp_site ) {
        $fetch_isw_id->execute( $align_id, $_, $_ );
        my ($isw_id)    = $fetch_isw_id->fetchrow_array;
        my $ref_base    = undef;
        my $snp_occured = undef;
        if ( defined $ref_info->{seqs} ) {
            $ref_base = substr( $ref_info->{seqs}, $_ - 1, 1 );
            if ( $ref_base eq $snp_site{$_}{target_base} ) {
                $snp_occured = "Q";
            }
            elsif ( $ref_base eq $snp_site{$_}{query_base} ) {
                $snp_occured = "T";
            }
            else {
                $snp_occured = "N";
            }
        }
        $snp_insert->execute(
            $isw_id, $align_id, $_,
            $snp_site{$_}{target_base},
            $snp_site{$_}{query_base},
            $ref_base, $snp_occured,
        );
    }
    $snp_insert->finish;

    return;
}

sub _modify_isw {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    # get sliding windows' sizes
    my $window_maker = $self->window_maker;
    my $windows_size = $window_maker->sw_size;

    # indel_id & foregoing_indel_id
    my $fetch_indel_id = $dbh->prepare(
        'SELECT i1.indel_id,
                i1.indel_occured,
                i2.indel_occured
        FROM indel i1, indel i2
        WHERE i1.foregoing_indel_id = i2.indel_id
        AND i1.align_id = ?
        AND i1.left_extand >= ?'
    );

    # isw_id
    my $fetch_isw_id = $dbh->prepare(
        'SELECT isw_id, isw_type
        FROM isw
        WHERE indel_id = ?'
    );

    # snp
    my $fetch_snp = $dbh->prepare(
        'SELECT snp_occured, COUNT(*)
        FROM snp
        WHERE isw_id = ?
        GROUP BY snp_occured'
    );

    # update isw
    my $update_isw = $dbh->prepare(
        'UPDATE isw
        SET isw_d_indel = ? / isw_length,
            isw_d_noindel = ? /isw_length,
            isw_d_complex = ? /isw_length
        WHERE isw_id = ?'
    );

    $fetch_indel_id->execute( $align_id, $windows_size );
    while ( my @row = $fetch_indel_id->fetchrow_array ) {
        my ( $indel_id, $indel_occured, $foregoing_indel_occured ) = @row;
        $fetch_isw_id->execute($indel_id);
        while ( my @row = $fetch_isw_id->fetchrow_array ) {
            my ( $isw_id, $isw_type ) = @row;
            my %occured;
            my ( $d_indel, $d_noindel, $d_complex );
            $fetch_snp->execute($isw_id);
            while ( my @row = $fetch_snp->fetchrow_array ) {
                my ( $snp_occured, $number ) = @row;
                $occured{$snp_occured} = $number;
            }

            # When there is not a snp,
            #   $occured{$snp_occured} will be undef.
            foreach (qw{T Q N}) {
                $occured{$_} ||= 0;
            }
            if ( $indel_occured eq "T" and $isw_type eq "R" ) {
                $d_indel   = $occured{T};
                $d_noindel = $occured{Q};
                $d_complex = $occured{N};
            }
            elsif ( $indel_occured eq "Q" and $isw_type eq "R" ) {
                $d_indel   = $occured{Q};
                $d_noindel = $occured{T};
                $d_complex = $occured{N};
            }
            elsif ( $foregoing_indel_occured eq "T"
                and $isw_type eq "L" )
            {
                $d_indel   = $occured{T};
                $d_noindel = $occured{Q};
                $d_complex = $occured{N};
            }
            elsif ( $foregoing_indel_occured eq "Q"
                and $isw_type eq "L" )
            {
                $d_indel   = $occured{Q};
                $d_noindel = $occured{T};
                $d_complex = $occured{N};
            }
            else {
                next;
            }

            $update_isw->execute( $d_indel, $d_noindel, $d_complex, $isw_id );
        }
    }

    return;
}

sub _insert_ssw {
    my $self     = shift;
    my $align_id = shift;

    my ( $align_set, $comparable_set, $indel_set )
        = $self->get_sets($align_id);

    my $dbh = $self->dbh;

    # get sliding windows' sizes
    my $window_maker = $self->window_maker;
    my $ssw_size     = $window_maker->sw_size;

    my $ssw_max_distance = 10;
    my $ssw_size_window0 = int( $ssw_size / 2 );

    # extreme_id & foregoing_extreme_id
    my $fetch_snp_id = $dbh->prepare(
        'SELECT s.snp_id, s.snp_pos, s.snp_occured
        FROM snp s
        WHERE s.align_id = ?'
    );
    $fetch_snp_id->execute($align_id);

    # prepare ssw_insert
    my $ssw_insert = $dbh->prepare(
        'INSERT INTO ssw (
            ssw_id, snp_id, window_id,
            ssw_type, ssw_distance, 
            ssw_d_snp, ssw_d_nosnp, ssw_d_complex
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, 
            ?, ?, ?
        )'
    );

    # store snp_info, use snp_index as key
    my $snp_info = {};
    while ( my @row = $fetch_snp_id->fetchrow_array ) {
        my ( $snp_id, $snp_pos, $snp_occured ) = @row;

        # index of snp in the $comparable_set
        my $snp_index = lookup_member( $comparable_set, $snp_pos );

        $snp_info->{$snp_index}{snp_id}      = $snp_id;
        $snp_info->{$snp_index}{snp_pos}     = $snp_pos;
        $snp_info->{$snp_index}{snp_occured} = $snp_occured;
    }

    foreach my $snp_index ( sort { $a <=> $b } keys %$snp_info ) {
        my $snp_id      = $snp_info->{$snp_index}{snp_id};
        my $snp_pos     = $snp_info->{$snp_index}{snp_pos};
        my $snp_occured = $snp_info->{$snp_index}{snp_occured};

        #print "snp_id: $snp_id | snp_occured: $snp_occured\n";

        unless ( $snp_occured eq 'T' or $snp_occured eq 'Q' ) {

            #print " " x 4, "jump to next\n";
            next;
        }

        # ssw has two types: L & R
        # More windows will be submitted in the following for-loop
        #   and if section
        foreach my $ssw_type (qw/L R/) {

            # $ssw_start and $ssw_end are both index of $comprarable_set
            my ( $ssw_start, $ssw_end );
            if ( $ssw_type eq 'R' ) {
                $ssw_start = $snp_index + 1;
                $ssw_end   = $ssw_start + $ssw_size_window0 - 1;
            }
            elsif ( $ssw_type eq 'L' ) {
                $ssw_end   = $snp_index - 1;
                $ssw_start = $ssw_end - $ssw_size_window0 + 1;
            }
            else {
                print "ssw_type \"$ssw_type\" error\n";
            }

            # $gsw_distance is from 0 to $gsw_density
        SSW: foreach my $i ( 0 .. $ssw_max_distance ) {
                my $ssw_set
                    = subset_span( $comparable_set, $ssw_start, $ssw_end );
                my $ssw_set_member_number = $ssw_set->cardinality;
                unless ($ssw_set_member_number) {
                    last SSW;
                }
                my $ssw_distance = $i;

                my ( $d_snp, $d_nosnp, $d_complex ) = ( 0, 0, 0 );

                foreach ( $ssw_start .. $ssw_end ) {
                    if ( exists $snp_info->{$_} ) {
                        my $occured = $snp_info->{$_}{snp_occured};
                        if ( $occured eq 'N' ) {
                            $d_complex++;
                        }
                        elsif ( $occured eq $snp_occured ) {
                            $d_snp++;
                        }
                        elsif ( $occured ne $snp_occured ) {
                            $d_nosnp++;
                        }
                    }
                }
                $d_snp     /= $ssw_set_member_number;
                $d_nosnp   /= $ssw_set_member_number;
                $d_complex /= $ssw_set_member_number;

                my ($cur_window_id)
                    = $self->insert_window( $align_id, $ssw_set );

                $ssw_insert->execute( $snp_id, $cur_window_id, $ssw_type,
                    $ssw_distance, $d_snp, $d_nosnp, $d_complex );

                if ( $ssw_type eq 'R' ) {
                    $ssw_start = $ssw_end + 1;
                    $ssw_end   = $ssw_start + $ssw_size - 1;
                }
                elsif ( $ssw_type eq 'L' ) {
                    $ssw_end   = $ssw_start - 1;
                    $ssw_start = $ssw_end - $ssw_size + 1;
                }
            }
        }
    }

    return;
}

##################################################
# Usage      : $self->add_align(
#            :     $target_info, $query_info,
#            :     $ref_info, $all_indel,
#            : );
# Purpose    : Add a new alignment and all its underlings
#            : to the database
# Returns    : $align_id
# Parameters : $target_info:    target info hash_ref
#            : $query_info:     query info hash_ref
#            : $ref_info:       reference info hash_ref
#            : $all_indel:      pre-defined indels
# Throws     : no exceptions
# Comments   : This method is the most important one in this module.
#            : All generating operations are performed here.
# See Also   : n/a
sub add_align {
    my ( $self, $target_info, $query_info, $ref_info, $all_indel ) = @_;

    $self->{caching_refs} = $ref_info;

    my $target_seq = $target_info->{seqs};
    my $query_seq  = $query_info->{seqs};

    my ( $ref_seq, $ref_raw_seq );
    if ( defined $ref_info->{seqs} ) {
        $ref_seq     = $ref_info->{seqs};
        $ref_raw_seq = $ref_info->{raw_seqs};
    }

    my $dbh = $self->dbh;

    # deltaG calculator
    my $dG_calc = $self->dG_calc;

    #----------------------------------------------------------#
    # GET seq info
    #----------------------------------------------------------#
    # taxon_id
    my $target_taxon_id = $target_info->{taxon_id} || 0;
    my $query_taxon_id  = $query_info->{taxon_id}  || 0;
    my $ref_taxon_id    = $ref_info->{taxon_id}    || 0;

    # name
    my $target_name = $target_info->{name} || $target_taxon_id;
    my $query_name  = $query_info->{name}  || $query_taxon_id;
    my $ref_name    = $ref_info->{name}    || $ref_taxon_id;

    # chr_id
    my $target_chr_id = $target_info->{chr_id} || undef;
    my $query_chr_id  = $query_info->{chr_id}  || undef;

    # chr_name
    my $target_chr_name = $target_info->{chr_name} || undef;
    my $query_chr_name  = $query_info->{chr_name}  || undef;

    # chr_start
    my $target_chr_start = $target_info->{chr_start} || undef;
    my $query_chr_start  = $query_info->{chr_start}  || undef;

    # chr_end
    my $target_chr_end = $target_info->{chr_end} || undef;
    my $query_chr_end  = $query_info->{chr_end}  || undef;

    # chr_strand
    my $target_chr_strand = $target_info->{chr_strand} || undef;
    my $query_chr_strand  = $query_info->{chr_strand}  || undef;

    # query_strand
    my $query_strand = $query_info->{query_strand} || "+";

    # Get tvsq_id from table tvsq or insert a new one
    my $tvsq_id = $self->_tvsq_id(
        $target_taxon_id, $target_name,  $query_taxon_id,
        $query_name,      $ref_taxon_id, $ref_name
    );

    # align length
    my $align_length = length $target_seq;
    if ( length $query_seq != $align_length ) {
        print "The two sequences has not equal length.\n";
        return;
    }

    #----------------------------#
    # INSERT INTO align
    #----------------------------#
    my $align_id
        = $self->_insert_align( $tvsq_id, \$target_seq, \$query_seq );
    print
        "Prosess align $align_id ",
        "in $target_chr_name $target_chr_start - $target_chr_end\n";

    #----------------------------#
    # INSERT INTO sequence, target, query
    #----------------------------#
    {
        my $sequence_insert = $dbh->prepare(
            'INSERT INTO sequence (
                seq_id, chr_id, chr_start, chr_end,
                chr_strand, seq_length
            )
            VALUES (
                NULL, ?, ?, ?,
                ?, ?
            )'
        );
        my $first_sequence = $target_seq;
        $first_sequence =~ s/\-//g;
        my $first_sequence_length = length $first_sequence;
        $sequence_insert->execute(
            $target_chr_id,  $target_chr_start,
            $target_chr_end, $target_chr_strand,
            $first_sequence_length,
        );
        my $target_seq_id = $self->last_insert_id;
        my $target_insert = $dbh->prepare(
            'INSERT INTO target (
                target_id, align_id, seq_id,
                target_seq, target_runlist
            )
            VALUES (
                NULL, ?, ?,
                ?, NULL
            )'
        );
        $target_insert->execute( $align_id, $target_seq_id, $target_seq, );
        $target_insert->finish;

        my $second_sequence = $query_seq;
        $second_sequence =~ s/\-//g;
        my $second_sequence_length = length $second_sequence;
        $sequence_insert->execute( $query_chr_id, $query_chr_start,
            $query_chr_end, $query_chr_strand, $second_sequence_length, );
        my $query_seq_id = $self->last_insert_id;
        my $query_insert = $dbh->prepare(
            'INSERT INTO query (
                query_id, seq_id, align_id,
                query_seq, query_strand, query_runlist
            )
            VALUES (
                NULL, ?, ?,
                ?, ?, NULL
            )'
        );
        $query_insert->execute( $query_seq_id, $align_id, $query_seq,
            $query_strand, );
        $query_insert->finish;
        $sequence_insert->finish;
    }

    #----------------------------#
    # INSERT INTO reference
    #----------------------------#
    if ( defined $ref_info->{seqs} ) {
        my $ref_complex_indel = $ref_info->{complex};
        my $ref_insert        = $dbh->prepare(
            'INSERT INTO reference (
                ref_id, align_id,
                ref_seq, ref_raw_seq, ref_complex_indel
            )
            VALUES (
                NULL, ?,
                ?, ?, ?
            )'
        );
        $ref_insert->execute( $align_id, $ref_seq, $ref_raw_seq,
            $ref_complex_indel );
        $ref_insert->finish;
    }

    #----------------------------#
    # INSERT INTO indel
    #----------------------------#
    {
        my $indel_insert = $dbh->prepare(
            'INSERT INTO indel (
                indel_id, foregoing_indel_id, align_id,
                indel_start, indel_end, indel_length,
                indel_seq, indel_insert, left_extand, right_extand,
                indel_gc_ratio, indel_dG, indel_occured, indel_type
            )
            VALUES (
                NULL, ?, ?,
                ?, ?, ?,
                ?, ?, ?, ?,
                ?, ?, ?, ?
            )'
        );
        my $indel_regex = qr{^\-+$};
        my $base_regex  = qr{^\w+$};

        my @indel_site
            = @{ pair_indel_sites( $target_seq, $query_seq, $all_indel ) };

        my $complex_set = AlignDB::IntSpan->new;
        if ( defined $ref_info->{seqs} ) {
            $complex_set->add( $ref_info->{complex} );
        }

        my $foregoing_indel_id = 0;
        foreach (@indel_site) {

            # $indel_occured:
            #   'C': complex indel
            #   'T': occured in target seq
            #   'Q': occured in query seq
            #   'N': occured in other place, noindel
            # $indel_type:
            #   'C': complex indel
            #   'D': deletion
            #   'I': insertion
            #   'N': noindel
            my $indel_occured = undef;
            my $indel_type    = undef;
            if ( defined $ref_info->{seqs} ) {
                my $start  = $_->{start};
                my $end    = $_->{end};
                my $length = $_->{length};

                if ( $complex_set->superset("$start-$end") ) {
                    $indel_occured = "C";
                    $indel_type    = "C";
                }
                else {
                    my $r_str = substr( $ref_seq,    $start - 1, $length );
                    my $t_str = substr( $target_seq, $start - 1, $length );
                    my $q_str = substr( $query_seq,  $start - 1, $length );
                    ( $indel_occured, $indel_type )
                        = ref_indel_type( $r_str, $t_str, $q_str );
                }
            }
            $indel_insert->execute(
                $foregoing_indel_id, $align_id,         $_->{start},
                $_->{end},           $_->{length},      $_->{seq},
                $_->{insert},        $_->{left_extand}, $_->{right_extand},
                $_->{gc_ratio},      $_->{dG},          $indel_occured,
                $indel_type,
            );
            ($foregoing_indel_id) = $self->last_insert_id;
        }
        $indel_insert->finish;
    }

    #----------------------------#
    # UPDATE target, query with runlist
    #----------------------------#
    {
        my $align_set      = AlignDB::IntSpan->new("1-$align_length");
        my $comparable_set = AlignDB::IntSpan->new;
        my $indel_set      = AlignDB::IntSpan->new;

        my $deletion_set  = AlignDB::IntSpan->new;
        my $insertion_set = AlignDB::IntSpan->new;

        # indel
        my $fetch_indel = $dbh->prepare(
            'SELECT indel_start, indel_end
            FROM indel
            WHERE align_id = ?
            AND indel_insert = ?'
        );

        # for target insert, i.e. deletion in query
        $fetch_indel->execute( $align_id, 'D' );
        while ( my @row = $fetch_indel->fetchrow_array ) {
            my ( $indel_start, $indel_end ) = @row;
            $deletion_set->add("$indel_start-$indel_end");
        }

        # for query insert, i.e. insertion in query
        $fetch_indel->execute( $align_id, 'I' );
        while ( my @row = $fetch_indel->fetchrow_array ) {
            my ( $indel_start, $indel_end ) = @row;
            $insertion_set->add("$indel_start-$indel_end");
        }

        $indel_set = $insertion_set->union($deletion_set);
        my $target_set = $align_set->diff($insertion_set);
        my $query_set  = $align_set->diff($deletion_set);
        $comparable_set = $align_set->diff($indel_set);

        my $target_update = $dbh->prepare(
            'UPDATE target
            SET target_runlist = ?
            WHERE align_id = ?'
        );
        $target_update->execute( $target_set->runlist, $align_id );

        my $query_update = $dbh->prepare(
            'UPDATE query
            SET query_runlist = ?
            WHERE align_id = ?'
        );
        $query_update->execute( $query_set->runlist, $align_id );

        my $align_update = $dbh->prepare(
            'UPDATE align
            SET comparable_runlist = ?,
                indel_runlist = ?
            WHERE align_id = ?'
        );
        $align_update->execute( $comparable_set->runlist, $indel_set->runlist,
            $align_id );
    }

    #----------------------------#
    # INSERT INTO isw
    #----------------------------#
    $self->_insert_isw($align_id);

    #----------------------------#
    # INSERT INTO snp
    #----------------------------#
    $self->_insert_snp($align_id);

    #----------------------------#
    # MODIFY isw
    #----------------------------#
    if ( defined $ref_info->{seqs} ) {
        $self->_modify_isw($align_id);
    }

    #----------------------------#
    # INSERT INTO ssw
    #----------------------------#
    if ( defined $ref_info->{seqs} and $self->insert_ssw ) {
        $self->_insert_ssw($align_id);
    }

    return $align_id;
}

##################################################
# Usage      : $self->parse_axt_file( $opt );
# Purpose    : read in alignments and chromosome position
#            :   info from .axt file
#            : then pass them to add_align method
# Returns    : none
# Parameters : $opt:        option hashref
# Throws     : no exceptions
# Comments   : This method is the most important one in this module.
#            : All generating operations are performed here.
# See Also   : n/a
sub parse_axt_file {
    my $self = shift;
    my $opt  = shift;

    my $axt_file        = $opt->{axt_file};
    my $target_taxon_id = $opt->{target_taxon_id};
    my $target_name     = $opt->{target_name};
    my $query_taxon_id  = $opt->{query_taxon_id};
    my $query_name      = $opt->{query_name};
    my $threshold       = $opt->{threshold};
    my $not_db          = $opt->{not_db};

    # minimal length
    $threshold ||= $self->threshold;

    open my $axt_fh, '<', $axt_file
        or confess "Cannot open align file $axt_file";

    my $target_chr_id = $self->get_chr_id_hash($target_taxon_id);
    my $query_chr_id  = $self->get_chr_id_hash($query_taxon_id);

    my %all_info_of;

    while (1) {
        my $summary_line = <$axt_fh>;
        unless ($summary_line) {
            last;
        }
        if ( $summary_line =~ /^#/ ) {
            next;
        }
        chomp $summary_line;
        chomp( my $first_line = <$axt_fh> );
        $first_line = uc $first_line;
        chomp( my $second_line = <$axt_fh> );
        $second_line = uc $second_line;
        my $dummy = <$axt_fh>;

        unless ( length $first_line > $threshold ) {
            next;
        }

        my ($align_serial, $first_chr,    $first_start,
            $first_end,    $second_chr,   $second_start,
            $second_end,   $query_strand, $align_score,
        ) = split /\s+/, $summary_line;

        #print "$align_serial\n";

        if ( $first_chr =~ /scaff|contig|super|bac|ran/i ) {
            $first_chr = 'chrUn';
        }
        if ( $second_chr =~ /scaff|contig|super|bac|ran/i ) {
            $second_chr = 'chrUn';
        }

        my $first_strand;
        if ( $first_end > $first_start ) {
            $first_strand = "+";
        }
        else {
            $first_strand = "-";
        }
        my $second_strand;
        if ( $second_end > $second_start ) {
            $second_strand = "+";
        }
        else {
            $second_strand = "-";
        }

        my $target_info = {
            taxon_id   => $target_taxon_id,
            name       => $target_name,
            chr_id     => $target_chr_id->{$first_chr},
            chr_name   => $first_chr,
            chr_start  => $first_start,
            chr_end    => $first_end,
            chr_strand => $first_strand,
            seqs       => $first_line,
        };
        my $query_info = {
            taxon_id     => $query_taxon_id,
            name         => $query_name,
            chr_id       => $query_chr_id->{$second_chr},
            chr_name     => $second_chr,
            chr_start    => $second_start,
            chr_end      => $second_end,
            chr_strand   => $second_strand,
            query_strand => $query_strand,
            seqs         => $second_line,
        };

        if ($not_db) {
            $all_info_of{$align_serial}           = {};
            $all_info_of{$align_serial}->{target} = $target_info;
            $all_info_of{$align_serial}->{query}  = $query_info;
        }
        else {
            $self->add_align( $target_info, $query_info );
        }
    }

    close $axt_fh;

    return \%all_info_of;
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

        my $dbh = $self->dbh;

        my $target_sth = $dbh->prepare(
            'SELECT t.target_seq
            FROM target t
            WHERE t.align_id = ?'
        );

        my $query_sth = $dbh->prepare(
            'SELECT q.query_seq
            FROM query q
            WHERE q.align_id = ?'
        );

        $target_sth->execute($align_id);
        my ($target_seq) = $target_sth->fetchrow_array;
        $target_sth->finish;

        $query_sth->execute($align_id);
        my ($query_seqs) = $query_sth->fetchrow_array;
        $query_sth->finish;

        $self->{caching_id} = $align_id;
        $self->{caching_seqs} = [ $target_seq, $query_seqs ];
    }

    return $self->caching_seqs;
}

sub get_sets {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $sth = $dbh->prepare(
        'SELECT align_length, comparable_runlist, indel_runlist
        FROM align 
        WHERE align_id = ?'
    );
    $sth->execute($align_id);
    my ( $align_length, $comparable_runlist, $indel_runlist )
        = $sth->fetchrow_array;
    $sth->finish;

    my $align_set      = AlignDB::IntSpan->new("1-$align_length");
    my $comparable_set = AlignDB::IntSpan->new($comparable_runlist);
    my $indel_set      = AlignDB::IntSpan->new($indel_runlist);

    return ( $align_set, $comparable_set, $indel_set );
}

sub get_chr_id_hash {
    my ( $self, $taxon_id ) = @_;

    my %chr_id = ();
    my $dbh    = $self->dbh;
    my $chromosome
        = $dbh->prepare(q{SELECT * FROM chromosome WHERE taxon_id = ?});
    $chromosome->execute($taxon_id);
    while ( my $ref = $chromosome->fetchrow_hashref ) {
        $chr_id{ $ref->{chr_name} } = $ref->{chr_id};
    }

    $chromosome->finish;

    return \%chr_id;
}

sub get_slice_stat {
    my $self     = shift;
    my $align_id = shift;
    my $set      = shift;

    my $dG_calc = $self->dG_calc;

    my $seqs_ref = $self->get_seqs($align_id);

    my @seq_slices = map { $set->substr_span($_) } @$seqs_ref;

    my $result = pair_seq_stat(@seq_slices);

    if ( $self->insert_dG ) {
        my @seq_dG = map { $dG_calc->polymer_deltaG($_) } @seq_slices;

        push @$result, (@seq_dG);
    }
    else {
        push @$result, ( undef, undef );
    }

    return $result;
}

sub get_slice_indel {
    my $self     = shift;
    my $align_id = shift;
    my $set      = shift;

    # slice indel; every gap is treated as an indel
    my $set_indel = $set->span_size - 1;

    # real indels in this alignment slice
    my $seqs_ref       = $self->get_seqs($align_id);
    my @seq_slices     = map { $set->substr_span($_) } @$seqs_ref;
    my @seq_indel_sets = map { find_indel_set($_) } @seq_slices;

    my $indel_set = AlignDB::IntSpan->new;
    for (@seq_indel_sets) {
        $indel_set->add($_);
    }
    my $real_indel = $indel_set->span_size - 1;

    # avoid empty sets which have span_size of 0
    $set_indel  = 0 if $set_indel < 0;
    $real_indel = 0 if $real_indel < 0;

    return ( $set_indel, $real_indel );
}

sub get_slice_pi_gc_dG {
    my $self     = shift;
    my $align_id = shift;
    my $set      = shift;

    my $slice_stat = $self->get_slice_stat( $align_id, $set );
    my $pi         = $slice_stat->[7];
    my $target_gc  = $slice_stat->[8];
    my $average_gc = $slice_stat->[9];
    my $target_dG  = $slice_stat->[10];
    my $query_dG   = $slice_stat->[11];

    return ( $pi, $target_gc, $average_gc, $target_dG, $query_dG );
}

##################################################
# Usage      : internal method
#            : $self->insert_window(
#            :     $align_id, $window_set, $internal_indel,
#            : );
# Purpose    : Add a window into alignDB,
# Returns    : a window_id
# Parameters : $align_id       : align_id
#            : $window_set     : AlignDB::IntSpan object
#            : $internal_indel : count all indels in this set (flag)
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub insert_window {
    my $self           = shift;
    my $align_id       = shift;
    my $window_set     = shift;
    my $internal_indel = shift;

    my $dbh = $self->dbh;

    my $window_sql = qq{
        INSERT INTO window (
            window_id, align_id, window_start, window_end, window_length,
            window_runlist, window_comparables, window_identities,
            window_differences, window_indel, window_pi,
            window_target_gc, window_average_gc,
            window_target_dG, window_query_dG,
            window_feature1, window_feature2, window_feature3
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?,
            ?, ?,
            NULL, NULL, NULL
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

        # avoid empty sets which have span_size of 0
        $window_indel = 0 if $window_indel < 0;
    }

    my $window_stat = $self->get_slice_stat( $align_id, $window_set );
    $window_insert->execute(
        $align_id,          $window_start,     $window_end,
        $window_length,     $window_runlist,   $window_stat->[1],
        $window_stat->[2],  $window_stat->[3], $window_indel,
        $window_stat->[7],  $window_stat->[8], $window_stat->[9],
        $window_stat->[10], $window_stat->[11],
    );

    return $self->last_insert_id;
}

sub last_insert_id {
    my $self = shift;

    my $dbh         = $self->dbh;
    my $last_insert = $dbh->prepare(q{SELECT LAST_INSERT_ID()});
    $last_insert->execute;
    my ($last_insert_id) = $last_insert->fetchrow_array;

    return $last_insert_id;
}

sub execute_sql {
    my ( $self, $sql_query, $bind_value ) = @_;

    # init
    my $dbh = $self->dbh;
    my $sth = $dbh->prepare($sql_query);

    # bind value
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    $sth->execute(@$bind_value);
}

##################################################
# Usage      : $self->empty_table(
#                  $table,
#                  $with_window
#              );
# Purpose    : Clean a table for new insertion
# Returns    : none
# Parameters : $table       : table name
#            : $with_window : delete corresponding window rows
# Throws     : Exception if $table does not exist
# Comments   : none
# See Also   : n/a
sub empty_table {
    my ( $self, $table, $with_window ) = @_;

    my $dbh = $self->dbh;

    # check table existing
    my @table_names = $dbh->tables( '', '', '' );

    # returned table names are quoted by `` (back-quotes)
    unless ( any { $_ =~ qr{`$table`} } @table_names ) {
        print "Table $table does not exist\n";
        return;
    }

    if ( !$with_window ) {
        $dbh->do(qq{DELETE FROM $table});
    }
    else {

        # In MySQL 4.1, you must use the alias (if one was given)
        #   when referring to a table name
        my $sql = qq{
            DELETE t,
                   w
            FROM $table t,
                 window w
            WHERE t.window_id = w.window_id
        };
        $dbh->do($sql);
    }

    # set AUTO_INCREMENT to 1 for this table
    # MyISAM table will memory last increment even this table is emptied
    $dbh->do("ALTER TABLE $table AUTO_INCREMENT = 1");

    return;
}

##################################################
# Usage      : $self->create_column(
#            :     $table,
#            :     $column,
#            :     $column_definition
#            : );
# Purpose    : Add $column to $table
#            :   with $column_definition as properties
# Returns    : none
# Parameters : $table             : table name
#            : $column            : cojumn name
#            : $column_definition : column properties
# Throws     : Exception if $table does not exist
# Comments   : none
# See Also   : n/a
sub create_column {
    my ( $self, $table, $column, $column_definition ) = @_;

    $column_definition ||= "DOUBLE";

    my $dbh = $self->dbh;

    # check table existing
    {
        my @table_names = $dbh->tables( '', '', '' );

        # table names are quoted by ` (back-quotes) which is the
        #   quote_identifier
        my $table_name = "`$table`";
        unless ( any { $_ =~ /$table_name/i } @table_names ) {
            print "Table $table does not exist\n";
            return;
        }
    }

    # check column existing
    # then create column
    {
        my $sql_query = qq{
            SHOW FIELDS
            FROM $table
            LIKE "$column"
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute;
        my ($field) = $sth->fetchrow_array;

        if ($field) {
            $dbh->do(qq{ALTER TABLE $table DROP COLUMN $column});
        }

        $dbh->do(
            qq{ALTER TABLE $table ADD COLUMN $column $column_definition} );
    }

    return;
}

sub check_column {
    my ( $self, $table, $column ) = @_;

    # init
    my $dbh = $self->dbh;

    # check table existing
    {
        my @table_names = $dbh->tables( '', '', '' );

        # table names are quoted by ` (back-quotes) which is the
        #   quote_identifier
        my $table_name = "`$table`";
        unless ( any { $_ =~ /$table_name/i } @table_names ) {
            print " " x 4, "Table $table does not exist\n";
            return;
        }
    }

    # check column existing
    {
        my $sql_query = qq{
            SHOW FIELDS
            FROM $table
            LIKE "$column"
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute;
        my ($field) = $sth->fetchrow_array;

        if ( not $field ) {
            print " " x 4, "Column $column does not exist\n";
            return;
        }
    }

    # check values in column
    # return row number (maybe 0 row)
    {
        my $sql_query = qq{
            SELECT COUNT($column)
            FROM $table
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute;
        my ($count) = $sth->fetchrow_array;

        if ( not $count ) {
            print " " x 4, "Column $column has no records\n";
        }

        return $count;
    }
}

##################################################
# Usage      : $self->get_names;
# Purpose    : get target_name, query_name & ref_name
# Returns    : ( $target_name, $query_name, $ref_name )
# Parameters : none
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub get_names {
    my $self = shift;
    my $align_id = shift || 1;

    my $dbh = $self->dbh;

    my $query = qq{
        SELECT  t.target_name,
                t.query_name,
                t.ref_name
        FROM tvsq t, align a
        WHERE t.tvsq_id = a.tvsq_id
        AND a.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my ( $target_name, $query_name, $ref_name ) = $sth->fetchrow_array;
    $sth->finish;

    return ( $target_name, $query_name, $ref_name );
}

##################################################
# Usage      : $self->get_taxon_ids;
# Purpose    : get target_taxon_id, query_taxon_id & ref_taxon_id
# Returns    : ( $target_taxon_id, $query_taxon_id, $ref_taxon_id )
# Parameters : none
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub get_taxon_ids {
    my $self = shift;
    my $align_id = shift || 1;

    my $dbh = $self->dbh;

    my $query = qq{
        SELECT  t.target_taxon_id,
                t.query_taxon_id,
                t.ref_taxon_id
        FROM tvsq t, align a
        WHERE t.tvsq_id = a.tvsq_id
        AND a.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my ( $target_taxon_id, $query_taxon_id, $ref_taxon_id )
        = $sth->fetchrow_array;
    $sth->finish;

    return ( $target_taxon_id, $query_taxon_id, $ref_taxon_id );
}

##################################################
# Usage      : $self->get_chr_info($chr_id);
# Purpose    : get the chr_name and chr_length of a given chr_id
# Returns    : ($chr_name, $chr_length)
# Parameters : $chr_id
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub get_chr_info {
    my $self   = shift;
    my $chr_id = shift;

    my $dbh = $self->dbh;

    my $query = qq{
        SELECT c.chr_name, c.chr_length
        FROM chromosome c
        WHERE c.chr_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($chr_id);
    my ( $chr_name, $chr_length ) = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_name, $chr_length );
}

##################################################
# Usage      : $self->get_align_ids;
# Purpose    : get an array of align_ids
# Returns    : $tbl_ary_ref ( \@align_ids)
# Parameters : none
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
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

##################################################
# Usage      : $self->get_align_ids_of_chr($chr_id);
# Purpose    : get an array of align_ids of a chr
# Returns    : $tbl_ary_ref ( \@align_ids)
# Parameters : none
# Throws     : no exceptions
# Comments   : none
# See Also   : n/a
sub get_align_ids_of_chr {
    my $self   = shift;
    my $chr_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT a.align_id
        FROM align a, target t, sequence s
        WHERE a.align_id = t.align_id
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
               t.target_runlist
        FROM target t, sequence s, chromosome c
        WHERE t.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        AND t.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my ( $chr_name, $chr_start, $chr_end, $target_runlist )
        = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_name, $chr_start, $chr_end, $target_runlist );
}

sub get_query_info {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT c.chr_name,
               s.chr_start,
               s.chr_end,
               q.query_runlist,
               q.query_strand
        FROM query q, sequence s, chromosome c
        WHERE q.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        AND q.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my ( $chr_name, $chr_start, $chr_end, $query_runlist, $query_strand )
        = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_name, $chr_start, $chr_end, $query_runlist, $query_strand );
}

##################################################
# Usage      : $self->get_target_chr_info($align_id);
# Purpose    : get target chr_id, chr_name, chr_length
# Returns    : ( $chr_id, $chr_name, $chr_length )
# Parameters : $align_id
# Throws     : no exceptions
# Comments   : none
# See Also   : get_query_chr_info
sub get_target_chr_info {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT c.chr_id, c.chr_name, c.chr_length
        FROM target t, sequence s, chromosome c
        WHERE t.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        AND t.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my ( $chr_id, $chr_name, $chr_length ) = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_id, $chr_name, $chr_length );
}

##################################################
# Usage      : $self->get_query_chr_info($align_id);
# Purpose    : get query chr_id, chr_name, chr_length
# Returns    : ( $chr_id, $chr_name, $chr_length )
# Parameters : $align_id
# Throws     : no exceptions
# Comments   : none
# See Also   : get_target_chr_info
sub get_query_chr_info {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT c.chr_id, c.chr_name, c.chr_length
        FROM query q, sequence s, chromosome c
        WHERE q.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        AND q.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my ( $chr_id, $chr_name, $chr_length ) = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_id, $chr_name, $chr_length );
}

##################################################
# Usage      : $self->get_chrs($goal);
# Purpose    : get all chrs of $goal (target or query)
# Returns    : arrayref of [ $chr_id, $chr_name, $chr_length ]
# Parameters : none
# Throws     : no exceptions
# Comments   : none
# See Also   : get_chr_info
sub get_chrs {
    my $self = shift;
    my $goal = shift;

    my $dbh = $self->dbh;

    # select all _GOAL_ chromosomes in this database
    my $chr_query = qq{
        # Processed chr
        SELECT c.chr_id, c.chr_name, c.chr_length
        FROM _GOAL_ G, sequence s, chromosome c
        WHERE G.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        GROUP BY c.chr_id
        ORDER BY c.chr_id
    };
    $chr_query =~ s/_GOAL_/$goal/;
    my $sth = $dbh->prepare($chr_query);
    $sth->execute;

    my $array_ref = $sth->fetchall_arrayref;
    $sth->finish;

    return $array_ref;
}

sub add_meta {
    my $self      = shift;
    my $meta_hash = shift;

    my $dbh = $self->dbh;

    my $query = q{
        INSERT INTO meta ( meta_id, meta_key, meta_value )
        VALUES ( NULL, ?, ? )
    };

    my $sth = $dbh->prepare($query);

    for my $key ( sort keys %$meta_hash ) {
        $sth->execute( $key, $meta_hash->{$key} );
    }

    $sth->finish;

    return;
}

sub add_meta_stopwatch {
    my $self      = shift;
    my $stopwatch = shift;

    my $dbh     = $self->dbh;
    my $ary_ref = $dbh->selectcol_arrayref(
        q{
        SELECT meta_value
        from meta
        where meta_key = 'uuid'
        }
    );

    my $uuid = $stopwatch->uuid;
    if ( !any { $_ eq $uuid } @$ary_ref ) {
        $self->add_meta(
            {   '------'      => '------',
                a_operation   => $stopwatch->operation,
                b_start_time  => $stopwatch->start_time2,
                c_store_time  => $stopwatch->time_now,
                d_duration    => $stopwatch->duration_now,
                e_cmd_line    => $stopwatch->cmd_line,
                f_init_config => $stopwatch->init_config,
                uuid          => $uuid,
            }
        );
    }

    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
__END__

=head1 NAME

    AlignDB - convert alignment filea to an indel-concentrated RDBMS
              (Default format is blastZ F<.axt>)

=head1 SYNOPSIS

    use AlignDB;
    my $obj = AlignDB->new(
        mysql          => "$db:$server",
        user           => $username,
        passwd         => $password,
        insert_dG      => $insert_dG,
    );
    $obj->parse_axt_file(
        {   axt_file        => $infile,
            target_taxon_id => $target_taxon_id,
            target_name     => $target_name,
            query_taxon_id  => $query_taxon_id,
            query_name      => $query_name,
            threshold       => $axt_threshold,
        }
    );

=head1 DESCRIPTION

C<AlignDB> is a simple class to convert alignment files to an indel-concentrated
RDBMS.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Authors' emails to help us keep track of the bugs
and their resolution. 

=head1 AUTHOR

B<Wang Qiang>

Email: wangqiang1997{at}nju{dot}edu{dot}cn

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=head1 METHODS

=cut

=head2 new

 Title   :  new
 Usage   :  my $obj = AlignDB->new;
 Function:  Builds a new AlignDB object 
 Returns :  AlignDB initialized with the correct format
 Args    :  mysql  => dsn, e.g. 'alignDB:202.119.43.248'
            user   => user
            passwd => passwd
            dbh    => database handle from DBI->connect
            insert_dG      => dG, default 0
            insert_ssw     => table: ssw, default 0

The constructor C<new> creates and returns an C<AlignDB> object.

=cut

