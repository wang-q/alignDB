package AlignDB;
use Moose;
use autodie;
use DBI;

use IO::Zlib;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all uniq);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Util qw(:all);

has 'mysql'  => ( is => 'ro', isa => 'Str' );    # e.g. 'alignDB:202.119.43.5'
has 'server' => ( is => 'ro', isa => 'Str' );    # e.g. '202.119.43.5'
has 'db'     => ( is => 'ro', isa => 'Str' );    # e.g. 'alignDB'
has 'user'   => ( is => 'ro', isa => 'Str' );    # database username
has 'passwd' => ( is => 'ro', isa => 'Str' );    # database password
has 'dbh'    => ( is => 'ro', isa => 'Ref' );    # store database handle here
has 'window_maker' => ( is => 'ro', isa => 'Object' );   # sliding windows maker
has 'threshold' => ( is => 'ro', isa => 'Int', default => sub {5_000} );

has 'caching_id' => ( is => 'ro', isa => 'Int' );        # caching seqs
has 'caching_seqs' =>
    ( is => 'ro', isa => 'ArrayRef[Str]', default => sub { [] } );

# target info
has 'caching_info' => ( is => 'ro', isa => 'HashRef', default => sub { {} } );

# don't connect mysql
has 'mocking' => ( is => 'ro', isa => 'Bool', default => 0 );

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
    elsif ( $self->mocking ) {

        # do nothing
    }
    else {
        confess "You should provide either mysql or db-server\n";
    }

    my $mysql  = $self->mysql;
    my $user   = $self->user;
    my $passwd = $self->passwd;
    my $server = $self->server;
    my $db     = $self->db;

    my $dbh = {};
    if ( !$self->mocking ) {
        $dbh = DBI->connect( "dbi:mysql:$mysql", $user, $passwd )
            or confess "Cannot connect to MySQL database at $mysql";
    }
    $self->{dbh} = $dbh;

    my $window_maker = AlignDB::Window->new;
    $self->{window_maker} = $window_maker;

    return;
}

sub _insert_align {
    my $self = shift;
    my @seqs = @_;

    my $dbh = $self->dbh;

    my $align_insert = $dbh->prepare(
        q{
        INSERT INTO align (
            align_id, align_length,
            align_comparables, align_identities, align_differences,
            align_gaps, align_ns, align_error,
            align_pi, align_target_gc, align_average_gc
        )
        VALUES (
            NULL, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?
        )
        }
    );

    my $result = multi_seq_stat(@seqs);

    $align_insert->execute(
        $result->[0], $result->[1], $result->[2], $result->[3], $result->[4],
        $result->[5], $result->[6], $result->[7], $result->[8], $result->[9],
    );
    $align_insert->finish;

    my $align_id = $self->last_insert_id;

    return $align_id;
}

sub _insert_seq {
    my $self     = shift;
    my $seq_info = shift;

    confess "Pass a seq to this method!\n" if !defined $seq_info->{seq};

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

sub _insert_set_and_sequence {
    my $self      = shift;
    my $align_id  = shift;
    my $info_refs = shift;
    my $seq_refs  = shift;

    my $dbh = $self->dbh;

    my $seq_number   = scalar @{$seq_refs};
    my $align_length = length $seq_refs->[0];

    my $align_set      = AlignDB::IntSpan->new("1-$align_length");
    my $indel_set      = AlignDB::IntSpan->new;
    my $comparable_set = AlignDB::IntSpan->new;

    for my $i ( 0 .. $seq_number - 1 ) {
        my $name = $info_refs->[$i]{name};
        $info_refs->[$i]{align_id} = $align_id;
        $info_refs->[$i]{seq}      = $seq_refs->[$i];
        $info_refs->[$i]{gc}       = calc_gc_ratio( $seq_refs->[$i] );
        my $seq_indel_set = find_indel_set( $seq_refs->[$i] );
        my $seq_set       = $align_set->diff($seq_indel_set);
        $info_refs->[$i]{runlist} = $seq_set->runlist;
        $info_refs->[$i]{length}  = $seq_set->cardinality;

        $indel_set->merge($seq_indel_set);
    }

    $comparable_set = $align_set->diff($indel_set);

    {    # sets
        my $align_update = $dbh->prepare(
            q{
            UPDATE align
            SET align_indels = ?,
                align_comparable_runlist = ?,
                align_indel_runlist = ?
            WHERE align_id = ?
            }
        );
        $align_update->execute( scalar $indel_set->spans,
            $comparable_set->runlist, $indel_set->runlist, $align_id );

    }

    {    # target
        my $insert = $dbh->prepare(
            q{
            INSERT INTO target ( target_id, seq_id )
            VALUES ( NULL, ? )
            }
        );
        my $seq_id = $self->_insert_seq( $info_refs->[0] );
        $insert->execute($seq_id);
        $insert->finish;
    }

    {    # and queries
        my $insert = $dbh->prepare(
            q{
            INSERT INTO query (
                query_id, seq_id, query_strand, query_position
            )
            VALUES ( NULL, ?, ?, ? )
            }
        );
        for my $i ( 1 .. $seq_number - 1 ) {
            my $seq_id = $self->_insert_seq( $info_refs->[$i] );
            $insert->execute( $seq_id, $info_refs->[$i]{chr_strand}, $i - 1 );
        }
        $insert->finish;
    }

    $self->{caching_id}   = $align_id;
    $self->{caching_seqs} = $seq_refs;

    return;
}

sub _insert_indel {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $indel_insert = $dbh->prepare(
        q{
        INSERT INTO indel (
            indel_id, prev_indel_id, align_id,
            indel_start, indel_end, indel_length,
            indel_seq, indel_all_seqs, left_extand, right_extand,
            indel_gc, indel_freq, indel_occured, indel_type
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?, ?, ?
        )
        }
    );

    my ( $align_set, undef, $indel_set ) = @{ $self->get_sets($align_id) };

    my $seq_refs  = $self->get_seqs($align_id);
    my $seq_count = scalar @{$seq_refs};

    my @indel_sites;
    for my $cur_indel ( $indel_set->spans ) {
        my ( $indel_start, $indel_end ) = @{$cur_indel};
        my $indel_length = $indel_end - $indel_start + 1;

        my @indel_seqs;
        for my $seq ( @{$seq_refs} ) {
            push @indel_seqs, ( substr $seq, $indel_start - 1, $indel_length );
        }
        my $indel_all_seqs = join "|", @indel_seqs;

        my $indel_type;
        my @uniq_indel_seqs = uniq(@indel_seqs);

        # seqs with least '-' char wins
        my ($indel_seq) = map { $_->[0] }
            sort { $a->[1] <=> $b->[1] }
            map { [ $_, tr/-/-/ ] } @uniq_indel_seqs;

        if ( scalar @uniq_indel_seqs < 2 ) {
            confess "no indel!\n";
        }
        elsif ( scalar @uniq_indel_seqs > 2 ) {
            $indel_type = 'C';
        }
        elsif ( $indel_seq =~ /-/ ) {
            $indel_type = 'C';
        }
        else {
            #   'D': means deletion relative to target/first seq
            #        target is ----
            #   'I': means insertion relative to target/first seq
            #        target is AAAA
            if ( $indel_seqs[0] eq $indel_seq ) {
                $indel_type = 'I';
            }
            else {
                $indel_type = 'D';
            }
        }

        my $indel_freq = 0;
        my $indel_occured;
        if ( $indel_type eq 'C' ) {
            $indel_freq    = -1;
            $indel_occured = 'unknown';
        }
        else {
            for (@indel_seqs) {

                # same as target 'x', not 'o'
                if ( $indel_seqs[0] eq $_ ) {
                    $indel_freq++;
                    $indel_occured .= 'o';
                }
                else {
                    $indel_occured .= 'x';
                }
            }
        }

        # here freq is the minor allele freq
        $indel_freq = min( $indel_freq, $seq_count - $indel_freq );

        my $indel_gc = calc_gc_ratio($indel_seq);

        push @indel_sites,
            {
            start    => $indel_start,
            end      => $indel_end,
            length   => $indel_length,
            seq      => $indel_seq,
            all_seqs => $indel_all_seqs,
            gc       => $indel_gc,
            freq     => $indel_freq,
            occured  => $indel_occured,
            type     => $indel_type,
            };
    }

    my $anterior_indel_end = 0;
    for my $i ( 0 .. scalar @indel_sites - 1 ) {
        my $cur_indel_start = $indel_sites[$i]->{start};
        my $cur_left_extand = $cur_indel_start - 1 - $anterior_indel_end;
        my $cur_indel_end   = $indel_sites[$i]->{end};
        $anterior_indel_end = $cur_indel_end;
        $indel_sites[$i]->{left_extand} = $cur_left_extand;
    }

    my $posterior_indel_start = $align_set->max + 1;
    for my $i ( reverse( 0 .. scalar @indel_sites - 1 ) ) {
        my $cur_indel_end    = $indel_sites[$i]->{end};
        my $cur_right_extand = $posterior_indel_start - $cur_indel_end - 1;
        my $cur_indel_start  = $indel_sites[$i]->{start};
        $posterior_indel_start = $cur_indel_start;
        $indel_sites[$i]->{right_extand} = $cur_right_extand;
    }

    my $prev_indel_id = 0;
    for (@indel_sites) {
        $indel_insert->execute(
            $prev_indel_id, $align_id,         $_->{start},
            $_->{end},      $_->{length},      $_->{seq},
            $_->{all_seqs}, $_->{left_extand}, $_->{right_extand},
            $_->{gc},       $_->{freq},        $_->{occured},
            $_->{type},
        );
        ($prev_indel_id) = $self->last_insert_id;
    }

    $indel_insert->finish;

    return;
}

sub _insert_snp {
    my $self     = shift;
    my $align_id = shift;

    my $dbh        = $self->dbh;
    my $snp_insert = $dbh->prepare(
        q{
        INSERT INTO snp (
            snp_id, align_id, snp_pos,
            target_base, query_base, all_bases,
            mutant_to, snp_freq, snp_occured
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?,
            ?, ?, ?
        )
        }
    );

    my ( $align_set, $comparable_set, $indel_set )
        = @{ $self->get_sets($align_id) };

    my $seq_refs  = $self->get_seqs($align_id);
    my $seq_count = scalar @{$seq_refs};

    my $snp_site = multi_snp_site( @{$seq_refs} );

    # %{$snp_site} keys are snp positions
    for my $pos ( sort { $a <=> $b } keys %{$snp_site} ) {

        my @bases = @{ $snp_site->{$pos} };

        my $target_base = $bases[0];
        my $all_bases = join '', @bases;

        my $query_base;
        my $mutant_to;
        my $snp_freq = 0;
        my $snp_occured;
        my @class = uniq(@bases);
        if ( scalar @class < 2 ) {
            confess "no snp!\n";
        }
        elsif ( scalar @class > 2 ) {
            $snp_freq    = -1;
            $snp_occured = 'unknown';
        }
        else {
            for (@bases) {
                if ( $target_base ne $_ ) {
                    $snp_freq++;
                    $snp_occured .= 'o';
                }
                else {
                    $snp_occured .= 'x';
                }
            }
            ($query_base) = grep { $_ ne $target_base } @bases;
            $mutant_to = $target_base . '<->' . $query_base;
        }

        # here freq is the minor allele freq
        $snp_freq = min( $snp_freq, $seq_count - $snp_freq );

        $snp_insert->execute( $align_id, $pos, $target_base, $query_base,
            $all_bases, $mutant_to, $snp_freq, $snp_occured, );
    }
    $snp_insert->finish;

    return;
}

sub insert_isw {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my ( $align_set, $comparable_set, $indel_set )
        = @{ $self->get_sets($align_id) };

    # indel_id & prev_indel_id
    my $fetch_indel_id_isw = $dbh->prepare(
        q{
        SELECT indel_id, prev_indel_id
        FROM indel
        WHERE align_id = ?
        }
    );
    $fetch_indel_id_isw->execute($align_id);

    # indel_end
    my $fetch_prev_indel_end = $dbh->prepare(
        q{
        SELECT indel_end
        FROM indel
        WHERE indel_id = ?
        }
    );

    # prev_indel_start
    my $fetch_indel_start = $dbh->prepare(
        q{
        SELECT indel_start
        FROM indel
        WHERE indel_id = ?
        }
    );

    # prepare isw_insert
    my $isw_insert = $dbh->prepare(
        q{
        INSERT INTO isw (
            isw_id, indel_id, prev_indel_id, isw_indel_id,
            isw_start, isw_end, isw_length, 
            isw_type, isw_distance, isw_density,
            isw_differences, isw_pi,
            isw_target_gc, isw_average_gc
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?,
            ?, ?
        )
        }
    );

    while ( my $ref = $fetch_indel_id_isw->fetchrow_hashref ) {
        my $indel_id      = $ref->{indel_id};
        my $prev_indel_id = $ref->{prev_indel_id};

        # bypass the first indel
        if ( $prev_indel_id == 0 ) {
            next;
        }

        my ( $interval_start, $interval_end, $interval_length ) = ( 0, 0, 0 );

        $fetch_prev_indel_end->execute($prev_indel_id);
        ($interval_start) = $fetch_prev_indel_end->fetchrow_array;
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

        my @isws = $window_maker->interval_window( $align_set, $interval_start,
            $interval_end );

        for my $isw (@isws) {
            my $isw_set    = $isw->{set};
            my $isw_start  = $isw_set->min;
            my $isw_end    = $isw_set->max;
            my $isw_length = $isw_end - $isw_start + 1;

            my $isw_indel_id;
            if ( $isw->{type} eq 'L' ) {
                $isw_indel_id = $prev_indel_id;
            }
            elsif ( $isw->{type} eq 'R' ) {
                $isw_indel_id = $indel_id;
            }
            elsif ( $isw->{type} eq 'S' ) {
                $isw_indel_id = $prev_indel_id;
            }
            my $isw_stat = $self->get_slice_stat( $align_id, $isw_set );
            $isw_insert->execute(
                $indel_id,      $prev_indel_id,   $isw_indel_id,
                $isw_start,     $isw_end,         $isw_length,
                $isw->{type},   $isw->{distance}, $isw->{density},
                $isw_stat->[3], $isw_stat->[7],   $isw_stat->[8],
                $isw_stat->[9],
            );
        }
    }

    $isw_insert->finish;
    $fetch_prev_indel_end->finish;
    $fetch_indel_start->finish;
    $fetch_indel_id_isw->finish;

    return;
}

sub isw_snp_fk {
    my $self     = shift;
    my $align_id = shift;

    my $dbh          = $self->dbh;
    my $fetch_snp_id = $dbh->prepare(
        q{
        SELECT s.snp_id, s.snp_pos
        FROM snp s
        WHERE s.align_id = ?
        }
    );

    my $fetch_isw_id = $dbh->prepare(
        q{
        SELECT isw_id
        FROM isw, indel
        WHERE isw.indel_id = indel.indel_id
        AND indel.align_id = ?
        AND isw.isw_start <= ?
        AND isw.isw_end >= ?
        }
    );

    my $snp_update = $dbh->prepare(
        q{
        UPDATE snp
        SET isw_id = ?
        WHERE snp_id = ?
        }
    );

    # process each snp
    $fetch_snp_id->execute($align_id);
    while ( my ( $snp_id, $snp_pos ) = $fetch_snp_id->fetchrow_array ) {
        $fetch_isw_id->execute( $align_id, $snp_pos, $snp_pos );
        my ($isw_id) = $fetch_isw_id->fetchrow_array;

        if ($isw_id) {
            $snp_update->execute( $isw_id, $snp_id );
        }
    }

    $snp_update->finish;
    $fetch_isw_id->finish;
    $fetch_snp_id->finish;

    return;
}

sub _modify_isw {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    # get sliding windows' sizes
    my $window_maker = $self->window_maker;
    my $windows_size = $window_maker->sw_size;

    # indel_id & prev_indel_id
    my $fetch_indel_id = $dbh->prepare(
        q{
        SELECT i1.indel_id,
                i1.indel_occured,
                i2.indel_occured
        FROM indel i1, indel i2
        WHERE i1.prev_indel_id = i2.indel_id
        AND i1.align_id = ?
        AND i1.left_extand >= ?
        }
    );

    # isw_id
    my $fetch_isw_id = $dbh->prepare(
        q{
        SELECT isw_id, isw_type
        FROM isw
        WHERE indel_id = ?
        }
    );

    # snp
    my $fetch_snp = $dbh->prepare(
        q{
        SELECT snp_occured, COUNT(*)
        FROM snp
        WHERE isw_id = ?
        GROUP BY snp_occured
        }
    );

    # update isw
    my $update_isw = $dbh->prepare(
        q{
        UPDATE isw
        SET isw_d_indel = ? / isw_length,
            isw_d_noindel = ? /isw_length,
            isw_d_complex = ? /isw_length
        WHERE isw_id = ?
        }
    );

    $fetch_indel_id->execute( $align_id, $windows_size );
    while ( my @row = $fetch_indel_id->fetchrow_array ) {
        my ( $indel_id, $indel_occured, $prev_indel_occured ) = @row;
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
            elsif ( $prev_indel_occured eq "T"
                and $isw_type eq "L" )
            {
                $d_indel   = $occured{T};
                $d_noindel = $occured{Q};
                $d_complex = $occured{N};
            }
            elsif ( $prev_indel_occured eq "Q"
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

sub insert_ssw {
    my $self     = shift;
    my $align_id = shift;

    my ( $align_set, $comparable_set, $indel_set )
        = @{ $self->get_sets($align_id) };

    my $dbh = $self->dbh;

    # get sliding windows' sizes
    my $window_maker = $self->window_maker;
    my $ssw_size     = $window_maker->sw_size;

    my $ssw_max_distance = 10;
    my $ssw_size_window0 = int( $ssw_size / 2 );

    my $fetch_snp_id = $dbh->prepare(
        q{
        SELECT s.snp_id, s.snp_pos, s.snp_occured
        FROM snp s
        WHERE s.align_id = ?
        }
    );
    $fetch_snp_id->execute($align_id);

    # prepare ssw_insert
    my $ssw_insert = $dbh->prepare(
        q{
        INSERT INTO ssw (
            ssw_id, snp_id, window_id,
            ssw_type, ssw_distance, 
            ssw_d_snp, ssw_d_nosnp, ssw_d_complex
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, 
            ?, ?, ?
        )
        }
    );

    # store snp_info, use snp_index as key
    my $snp_info = {};
    while ( my @row = $fetch_snp_id->fetchrow_array ) {
        my ( $snp_id, $snp_pos, $snp_occured ) = @row;

        # index of snp in the $comparable_set
        my $snp_index = $comparable_set->index($snp_pos);

        $snp_info->{$snp_index}{snp_id}      = $snp_id;
        $snp_info->{$snp_index}{snp_pos}     = $snp_pos;
        $snp_info->{$snp_index}{snp_occured} = $snp_occured;
    }

    for my $snp_index ( sort { $a <=> $b } keys %$snp_info ) {
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
        for my $ssw_type (qw/L R/) {

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

            # distance is from 0 to density
        SSW: for my $i ( 0 .. $ssw_max_distance ) {
                my $ssw_set = $comparable_set->slice( $ssw_start, $ssw_end );
                my $ssw_set_member_number = $ssw_set->cardinality;
                unless ($ssw_set_member_number) {
                    last SSW;
                }
                my $ssw_distance = $i;

                my ( $d_snp, $d_nosnp, $d_complex ) = ( 0, 0, 0 );

                for ( $ssw_start .. $ssw_end ) {
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

# Add a new alignment to the database
# This method is the most important one in this module.
# All generating operations are performed here.
sub add_align {
    my $self      = shift;
    my $info_refs = shift;
    my $seq_refs  = shift;

    my $dbh = $self->dbh;

    my $target_idx = 0;

    # check align length
    my $align_length = length $seq_refs->[$target_idx];
    for ( @{$seq_refs} ) {
        if ( ( length $_ ) != $align_length ) {
            confess "Sequences should have the same length!\n";
        }
    }

    # check seq number
    my $seq_number = scalar @{$seq_refs};
    if ( $seq_number < 2 ) {
        confess "Too few sequences [$seq_number]\n";
    }

    # check info and seq numbers
    if ( $seq_number != scalar @{$info_refs} ) {
        confess "Number of infos is not equal to seqs!\n";
    }

    #----------------------------#
    # INSERT INTO align
    #----------------------------#
    my $align_id = $self->_insert_align( @{$seq_refs} );
    printf "Prosess align [%s] at %s.%s(%s):%s-%s\n", $align_id,
        $info_refs->[$target_idx]{name},
        $info_refs->[$target_idx]{chr_name},
        $info_refs->[$target_idx]{chr_strand},
        $info_refs->[$target_idx]{chr_start},
        $info_refs->[$target_idx]{chr_end};

    #----------------------------#
    # UPDATE align, INSERT INTO sequence, target, queries
    #----------------------------#
    $self->_insert_set_and_sequence( $align_id, $info_refs, $seq_refs );

    #----------------------------#
    # INSERT INTO indel
    #----------------------------#
    $self->_insert_indel($align_id);

    #----------------------------#
    # INSERT INTO snp
    #----------------------------#
    $self->_insert_snp($align_id);

    return $align_id;
}

# read in alignments and chromosome position info from .axt file then pass them
#   to add_align method
sub parse_axt_file {
    my $self   = shift;
    my $infile = shift;
    my $opt    = shift;

    my $target_taxon_id = $opt->{target_taxon_id};
    my $query_taxon_id  = $opt->{query_taxon_id};
    my $threshold       = $opt->{threshold};
    my $gzip            = $opt->{gzip};

    my $target_name      = $self->get_name_of($target_taxon_id);
    my $query_name       = $self->get_name_of($query_taxon_id);
    my $target_chr_id_of = $self->get_chr_id_hash($target_taxon_id);
    my $query_chr_id_of  = $self->get_chr_id_hash($query_taxon_id);

    # minimal length
    $threshold ||= $self->threshold;

    my $in_fh;
    if ( !$gzip ) {
        open $in_fh, '<', $infile;
    }
    else {
        $in_fh = IO::Zlib->new( $infile, "rb" );
    }

    while (1) {
        my $summary_line = <$in_fh>;
        last unless $summary_line;
        next if $summary_line =~ /^#/;

        chomp $summary_line;
        chomp( my $first_line = <$in_fh> );
        $first_line = uc $first_line;
        chomp( my $second_line = <$in_fh> );
        $second_line = uc $second_line;
        my $dummy = <$in_fh>;

        next if length $first_line < $threshold;

        my ($align_serial, $first_chr,    $first_start,
            $first_end,    $second_chr,   $second_start,
            $second_end,   $query_strand, $align_score,
        ) = split /\s+/, $summary_line;

        my $info_refs = [
            {   taxon_id   => $target_taxon_id,
                name       => $target_name,
                chr_name   => $first_chr,
                chr_id     => $target_chr_id_of->{$first_chr},
                chr_start  => $first_start,
                chr_end    => $first_end,
                chr_strand => '+',
            },
            {   taxon_id   => $query_taxon_id,
                name       => $query_name,
                chr_name   => $second_chr,
                chr_id     => $query_chr_id_of->{$second_chr},
                chr_start  => $second_start,
                chr_end    => $second_end,
                chr_strand => $query_strand,
            },
        ];

        $self->add_align( $info_refs, [ $first_line, $second_line ], );
    }

    if ( !$gzip ) {
        close $in_fh;
    }
    else {
        $in_fh->close;
    }

    return;
}

# blocked fasta format
sub parse_block_fasta_file {
    my $self   = shift;
    my $in_file = shift;
    my $opt    = shift;

    my $id_of     = $opt->{id_of};
    my $threshold = $opt->{threshold};
    my $gzip      = $opt->{gzip};

    my $in_fh;
    if ( !$gzip ) {
        open $in_fh, '<', $in_file;
    }
    else {
        $in_fh = IO::Zlib->new( $in_file, "rb" );
    }

    my $content = '';
    while ( my $line = <$in_fh> ) {
        if ( $line =~ /^\s+$/ and $content =~ /\S/ ) {
            my @lines = grep {/\S/} split /\n/, $content;
            $content = '';
            die "headers not equal to seqs\n" if @lines % 2;
            die "Two few lines in block\n" if @lines < 4;

            my ( @headers, @seqs );
            while (@lines) {
                my $header = shift @lines;
                $header =~ s/^\>//;
                chomp $header;
                my $seq = shift @lines;
                chomp $seq;
                $seq = uc $seq;
                push @headers, $header;
                push @seqs,    $seq;
            }

            next if length $seqs[0] < $threshold;

            #S288C:
            #  chr_end: 667886
            #  chr_id: 265
            #  chr_name: chrIV
            #  chr_start: 652404
            #  chr_strand: +
            #  name: S288C
            #  taxon_id: 4932
            my $info_refs = [];
            for my $header (@headers) {
                my $info_ref = decode_header($header);
                $info_ref->{taxon_id} = $id_of->{ $info_ref->{name} };
                $info_ref->{chr_id}
                    = $self->get_chr_id_hash( $info_ref->{taxon_id} )
                    ->{ $info_ref->{chr_name} };

                push @{$info_refs}, $info_ref;
            }

            $self->add_align( $info_refs, \@seqs );
        }
        else {
            $content .= $line;
        }
    }

    if ( !$gzip ) {
        close $in_fh;
    }
    else {
        $in_fh->close;
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

        my $dbh = $self->dbh;

        my $target_sth = $dbh->prepare(
            q{
            SELECT s.seq_seq
            FROM sequence s
            INNER JOIN target t on s.seq_id = t.seq_id
            WHERE s.align_id = ?
            }
        );

        my $query_sth = $dbh->prepare(
            q{
            SELECT s.seq_seq
            FROM sequence s
            INNER JOIN query q on s.seq_id = q.seq_id
            WHERE s.align_id = ?
            ORDER BY q.query_position
            }
        );

        $target_sth->execute($align_id);
        my ($target_seq) = $target_sth->fetchrow_array;
        $target_sth->finish;

        $query_sth->execute($align_id);
        my @query_seqs;
        while ( my ($query_seq) = $query_sth->fetchrow_array ) {
            push @query_seqs, $query_seq;
        }
        $query_sth->finish;

        $self->{caching_id} = $align_id;
        $self->{caching_seqs} = [ $target_seq, @query_seqs ];
    }

    return $self->caching_seqs;
}

sub get_seq_ref {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $sth = $dbh->prepare(
        q{
        SELECT s.seq_seq
        FROM sequence s
        INNER JOIN reference r on s.seq_id = r.seq_id
        WHERE s.align_id = ?
        }
    );
    $sth->execute($align_id);
    my ($seq) = $sth->fetchrow_array;
    $sth->finish;

    return $seq;
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

    my @names;
    for my $table (qw{target query reference}) {
        my $sql = qq{
            SELECT 
                c.common_name
            FROM
                sequence s
                    INNER JOIN
                _TABLE_ ON s.seq_id = _TABLE_.seq_id
                    INNER JOIN
                (SELECT 
                    c.taxon_id, t.common_name, c.chr_id
                FROM
                    chromosome c
                INNER JOIN taxon t ON c.taxon_id = t.taxon_id) c ON c.chr_id = s.chr_id
            WHERE
                s.align_id = ?
        };

        $sql =~ s/_TABLE_/$table/g;

        $sql .= "ORDER BY query.query_position" if $table eq 'query';

        my $sth = $dbh->prepare($sql);
        $sth->execute($align_id);
        while ( my ($name) = $sth->fetchrow_array ) {
            push @names, $name;
        }
        $sth->finish;
    }

    return (@names);
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

    my @ids;
    for my $table (qw{target query reference}) {
        my $query = qq{
            SELECT 
                c.taxon_id
            FROM
                sequence s
                    INNER JOIN
                _TABLE_ ON s.seq_id = _TABLE_.seq_id
                    INNER JOIN
                (SELECT 
                    c.taxon_id, t.common_name, c.chr_id
                FROM
                    chromosome c
                INNER JOIN taxon t ON c.taxon_id = t.taxon_id) c ON c.chr_id = s.chr_id
            WHERE
                s.align_id = ?
        };

        $query =~ s/_TABLE_/$table/g;

        my $sth = $dbh->prepare($query);
        $sth->execute($align_id);
        while ( my ($id) = $sth->fetchrow_array ) {
            push @ids, $id;
        }
        $sth->finish;
    }

    return (@ids);
}

sub update_names {
    my $self    = shift;
    my $name_of = shift;

    my $dbh = $self->dbh;

    my $query = q{
        UPDATE taxon
        SET common_name = ?
        WHERE taxon_id = ?
    };

    my $sth = $dbh->prepare($query);
    for my $taxon_id ( keys %{$name_of} ) {
        $sth->execute( $name_of->{$taxon_id}, $taxon_id );
    }
    $sth->finish;

    return;
}

sub get_name_of {
    my $self     = shift;
    my $taxon_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        select common_name
        from taxon
        where taxon_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($taxon_id);
    my ($name) = $sth->fetchrow_array;
    $sth->finish;

    return $name;
}

sub get_sets {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $sth = $dbh->prepare(
        q{
        SELECT align_length, align_comparable_runlist, align_indel_runlist
        FROM align 
        WHERE align_id = ?
        }
    );
    $sth->execute($align_id);
    my ( $align_length, $comparable_runlist, $indel_runlist )
        = $sth->fetchrow_array;
    $sth->finish;

    my $align_set      = AlignDB::IntSpan->new("1-$align_length");
    my $comparable_set = AlignDB::IntSpan->new($comparable_runlist);
    my $indel_set      = AlignDB::IntSpan->new($indel_runlist);

    return [ $align_set, $comparable_set, $indel_set ];
}

sub get_chr_id_hash {
    my $self     = shift;
    my $taxon_id = shift;

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

    my $seqs_ref   = $self->get_seqs($align_id);
    my @seq_slices = map { $set->substr_span($_) } @$seqs_ref;
    my $result     = multi_seq_stat(@seq_slices);

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

    my $window_sql = q{
        INSERT INTO window (
            window_id, align_id, window_start, window_end, window_length,
            window_runlist, window_comparables, window_identities,
            window_differences, window_indel, window_pi,
            window_target_gc, window_average_gc
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?
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
        $align_id,         $window_start,     $window_end,
        $window_length,    $window_runlist,   $window_stat->[1],
        $window_stat->[2], $window_stat->[3], $window_indel,
        $window_stat->[7], $window_stat->[8], $window_stat->[9],
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

sub index_isw_indel_id {
    my $self = shift;

    my $query = "CREATE INDEX indel_isw_id_FK ON isw ( isw_indel_id );";
    $self->execute_sql($query);
    return;
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

        $dbh->do(qq{ALTER TABLE $table ADD COLUMN $column $column_definition});
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
        FROM sequence s
        INNER JOIN target t ON s.seq_id = t.seq_id
        INNER JOIN align a ON s.align_id = a.align_id
        WHERE s.chr_id = ?
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

sub get_align_ids_of_chr_name {
    my $self     = shift;
    my $chr_name = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT a.align_id
        FROM sequence s
        INNER JOIN target t ON s.seq_id = t.seq_id
        INNER JOIN align a ON s.align_id = a.align_id
        INNER JOIN chromosome c on s.chr_id = c.chr_id
        WHERE c.chr_name = ?
        ORDER BY a.align_id
    };

    my @align_ids;

    my $sth = $dbh->prepare($query);
    $sth->execute($chr_name);
    while ( my @row = $sth->fetchrow_array ) {
        my ($align_id) = @row;
        push @align_ids, $align_id;
    }

    return \@align_ids;
}

sub get_target_info {
    my $self     = shift;
    my $align_id = shift;

    my $caching_info = $self->caching_info;

    if ( !exists $caching_info->{$align_id} ) {
        my $dbh = $self->dbh;

        my $query = q{
            SELECT c.taxon_id,
                   c.chr_id,
                   c.chr_name,
                   c.chr_length,
                   s.chr_start,
                   s.chr_end,
                   s.chr_strand,
                   s.seq_length,
                   s.seq_gc,
                   s.seq_runlist,
                   a.align_length
            FROM sequence s 
            INNER JOIN target t ON s.seq_id = t.seq_id
            LEFT JOIN chromosome c ON s.chr_id = c.chr_id
            INNER JOIN align a ON s.align_id = a.align_id
            WHERE s.align_id = ?
        };

        my $sth = $dbh->prepare($query);
        $sth->execute($align_id);
        my $hash_ref = $sth->fetchrow_hashref;
        $sth->finish;
        $caching_info->{$align_id} = $hash_ref;
    }

    return $caching_info->{$align_id};
}

# general purpose
sub get_queries_info {
    my $self     = shift;
    my $align_id = shift;

    my $dbh = $self->dbh;

    my $query = q{
        SELECT c.taxon_id,
               c.chr_id,
               c.chr_name,
               c.chr_length,
               s.chr_start,
               s.chr_end,
               s.chr_strand,
               s.seq_length,
               s.seq_gc,
               s.seq_runlist,
               a.align_length,
               q.query_strand,
               q.query_position
        FROM sequence s 
        INNER JOIN query q ON s.seq_id = q.seq_id
        LEFT JOIN chromosome c ON s.chr_id = c.chr_id
        INNER JOIN align a ON s.align_id = a.align_id
        WHERE s.align_id = ?
    };

    my $sth = $dbh->prepare($query);
    $sth->execute($align_id);
    my @array;
    while ( my $hash_ref = $sth->fetchrow_hashref ) {
        push @array, $hash_ref;
    }
    $sth->finish;

    return @array;
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
        FROM sequence s 
        INNER JOIN target t ON s.seq_id = t.seq_id
        INNER JOIN chromosome c ON s.chr_id = c.chr_id
        WHERE s.align_id = ?
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
        FROM sequence s 
        INNER JOIN query q ON s.seq_id = q.seq_id
        INNER JOIN chromosome c ON s.chr_id = c.chr_id
        WHERE s.align_id = ?
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
    my $goal = shift || 'target';

    my $dbh = $self->dbh;

    # select all _GOAL_ chromosomes in this database
    my $chr_query = q{
        # Processed chr
        SELECT c.chr_id, c.chr_name, c.chr_length
        FROM sequence s 
        INNER JOIN _GOAL_ G ON s.seq_id = G.seq_id
        INNER JOIN chromosome c ON s.chr_id = c.chr_id
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

sub get_freq {
    my $self = shift;

    my $dbh = $self->dbh;

    my $sql_query = q{
        SELECT DISTINCT COUNT(q.query_id) + 1
        FROM  query q, sequence s
        WHERE q.seq_id = s.seq_id
        GROUP BY s.align_id
    };
    my $sth = $dbh->prepare($sql_query);

    my @counts;
    $sth->execute;
    while ( my ($count) = $sth->fetchrow_array ) {
        push @counts, $count;
    }
    if ( scalar @counts > 1 ) {
        warn "Database corrupts, freqs are not consistent\n";
    }

    return $counts[0];
}

sub process_message {
    my $self     = shift;
    my $align_id = shift;

    my $info = $self->get_target_info($align_id);

    printf "Process align [%s] at %s(%s):%s-%s\n", $align_id, $info->{chr_name},
        $info->{chr_strand}, $info->{chr_start}, $info->{chr_end};

    return;
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
    );

=head1 DESCRIPTION

C<AlignDB> is a simple class to convert alignment files to an indel-concentrated
RDBMS.

=head1 AUTHOR

B<Wang Qiang>

Email: wangq{at}nju{dot}edu{dot}cn

=head1 APPENDIX

Internal methods are usually preceded with a _


