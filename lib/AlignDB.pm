package AlignDB;
use Moose;
use autodie;

use Carp;
use DBI;
use IO::Zlib;
use List::Util;
use List::MoreUtils::PP;
use YAML::Syck;

use AlignDB::IntSpan;
use AlignDB::Window;

use App::RL::Common;
use App::Fasops::Common;

has 'mysql'  => ( is => 'ro', isa => 'Str' );    # e.g. 'alignDB:202.119.43.5'
has 'server' => ( is => 'ro', isa => 'Str' );    # e.g. '202.119.43.5'
has 'db'     => ( is => 'ro', isa => 'Str' );    # e.g. 'alignDB'
has 'user'   => ( is => 'ro', isa => 'Str' );    # database username
has 'passwd' => ( is => 'ro', isa => 'Str' );    # database password

has 'dbh'          => ( is => 'ro', isa => 'Object' );    # store database handle here
has 'window_maker' => ( is => 'ro', isa => 'Object' );    # sliding windows maker
has 'threshold' => ( is => 'ro', isa => 'Int', default => sub {5_000} );

has 'caching_id' => ( is => 'ro', isa => 'Int' );                                     # caching seqs
has 'caching_seqs' => ( is => 'ro', isa => 'ArrayRef[Str]', default => sub { [] } );

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
        Carp::confess "You should provide either mysql or db:server\n";
    }

    my $mysql  = $self->mysql;
    my $user   = $self->user;
    my $passwd = $self->passwd;

    my DBI $dbh = {};
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
    my $self     = shift;
    my $seq_refs = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        INSERT INTO align (
            align_id, align_length,
            align_comparables, align_identities, align_differences,
            align_gaps, align_ns, align_error,
            align_pi, align_gc, align_average_gc
        )
        VALUES (
            NULL, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?
        )
        }
    );

    my $result     = App::Fasops::Common::multi_seq_stat($seq_refs);
    my $target_gc  = App::Fasops::Common::calc_gc_ratio( [ $seq_refs->[0] ] );
    my $average_gc = App::Fasops::Common::calc_gc_ratio($seq_refs);

    $sth->execute(
        $result->[0], $result->[1], $result->[2], $result->[3], $result->[4],
        $result->[5], $result->[6], $result->[7], $target_gc,   $average_gc,
    );
    $sth->finish;

    my $align_id = $self->last_insert_id;

    return $align_id;
}

sub get_chr_id {
    my $self        = shift;
    my $common_name = shift;
    my $chr_name    = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT c.chr_id
        FROM chromosome c
        WHERE 1=1
        and c.common_name = ?
        and c.chr_name = ?
        }
    );
    $sth->execute( $common_name, $chr_name );
    my ($chr_id) = $sth->fetchrow_array;
    $sth->finish;

    return $chr_id;
}

sub _insert_seq {
    my $self     = shift;
    my $align_id = shift;
    my $seq_info = shift;

    Carp::confess "Need align_id\n" if !defined $align_id;
    Carp::confess "Need seq_info\n" if !defined $seq_info || ref $seq_info ne "HASH";

    Carp::confess "Need common_name\n"             if !defined $seq_info->{name};
    Carp::confess "Need chr_name\n"                if !defined $seq_info->{chr_name};
    Carp::confess "T, Q, or O\n"                   if !defined $seq_info->{seq_role};
    Carp::confess "Need positions in alignments\n" if !defined $seq_info->{seq_position};
    Carp::confess "Pass a seq to this method\n"    if !defined $seq_info->{seq};

    for my $key (qw{chr_start chr_end chr_strand length gc runlist}) {
        if ( !defined $seq_info->{$key} ) {
            $seq_info->{$key} = undef;
        }
    }

    my $chr_id = $self->get_chr_id( $seq_info->{name}, $seq_info->{chr_name} );

    # Get database handle
    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        INSERT INTO sequence (
            seq_id, align_id, chr_id,
            seq_role, seq_position, common_name,
            chr_name, chr_start, chr_end, chr_strand,
            seq_length, seq_gc, seq_runlist, seq_seq
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?, ?, ?
        )
        }
    );

    $sth->execute(
        $align_id,                 $chr_id,              $seq_info->{seq_role},
        $seq_info->{seq_position}, $seq_info->{name},    $seq_info->{chr_name},
        $seq_info->{chr_start},    $seq_info->{chr_end}, $seq_info->{chr_strand},
        $seq_info->{length},       $seq_info->{gc},      $seq_info->{runlist},
        $seq_info->{seq},
    );
}

sub _insert_set_and_sequence {
    my $self      = shift;
    my $align_id  = shift;
    my $info_refs = shift;
    my $seq_refs  = shift;

    my DBI $dbh = $self->dbh;

    my $seq_number   = scalar @{$seq_refs};
    my $align_length = length $seq_refs->[0];

    {    # sets
        my $align_set      = AlignDB::IntSpan->new("1-$align_length");
        my $indel_set      = AlignDB::IntSpan->new;
        my $comparable_set = AlignDB::IntSpan->new;
        for my $i ( 0 .. $seq_number - 1 ) {
            $info_refs->[$i]{seq} = $seq_refs->[$i];
            $info_refs->[$i]{gc} = App::Fasops::Common::calc_gc_ratio( [ $seq_refs->[$i] ] );
            my $seq_indel_set = App::Fasops::Common::indel_intspan( $seq_refs->[$i] );
            my $seq_set       = $align_set->diff($seq_indel_set);
            $info_refs->[$i]{runlist} = $seq_set->runlist;
            $info_refs->[$i]{length}  = $seq_set->size;

            $indel_set->merge($seq_indel_set);
        }
        $comparable_set = $align_set->diff($indel_set);

        my DBI $sth = $dbh->prepare(
            q{
            UPDATE align
            SET align_indels = ?,
                align_comparable_runlist = ?,
                align_indel_runlist = ?
            WHERE align_id = ?
            }
        );
        $sth->execute( scalar $indel_set->span_size,
            $comparable_set->runlist, $indel_set->runlist, $align_id );
    }

    {    # target
        $info_refs->[0]{seq_role}     = "T";
        $info_refs->[0]{seq_position} = 0;
        $self->_insert_seq( $align_id, $info_refs->[0] );
    }

    {    # and queries
        for my $i ( 1 .. $seq_number - 1 ) {
            $info_refs->[$i]{seq_role}     = "Q";
            $info_refs->[$i]{seq_position} = $i;
            $self->_insert_seq( $align_id, $info_refs->[$i] );
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

    my $seq_refs = [ map { $_->{seq} } @{$info_refs} ];

    my $target_idx = 0;

    # check align length
    my $align_length = length $seq_refs->[$target_idx];
    for ( @{$seq_refs} ) {
        if ( ( length $_ ) != $align_length ) {
            Carp::confess "Sequences should have the same length!\n";
        }
    }

    # check seq number
    my $seq_count = scalar @{$seq_refs};
    if ( $seq_count < 2 ) {
        Carp::confess "Too few sequences [$seq_count]\n";
    }

    #----------------------------#
    # INSERT INTO align
    #----------------------------#
    my $align_id = $self->_insert_align($seq_refs);
    printf "Prosess align [%s] at %s.%s(%s):%s-%s\n", $align_id,
        $info_refs->[$target_idx]{name},
        $info_refs->[$target_idx]{chr_name},
        $info_refs->[$target_idx]{chr_strand},
        $info_refs->[$target_idx]{chr_start},
        $info_refs->[$target_idx]{chr_end};
    $self->{caching_id}   = $align_id;
    $self->{caching_seqs} = $seq_refs;

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

sub parse_axt_file {
    my $self   = shift;
    my $infile = shift;
    my $opt    = shift;

    # minimal length
    my $threshold = $opt->{threshold};
    $threshold ||= $self->threshold;

    my $query_length_of = $self->get_chr_legnth_of( $opt->{qname} );

    my $in_fh = IO::Zlib->new( $infile, "rb" );
    my $content = '';    # content of one block
    while (1) {
        last if $in_fh->eof and $content eq '';
        my $line = '';
        if ( !$in_fh->eof ) {
            $line = $in_fh->getline;
        }

        if ( substr( $line, 0, 1 ) eq "#" ) {
            next;
        }
        elsif ( ( $line eq '' or $line =~ /^\s+$/ ) and $content ne '' ) {
            my $info_refs = App::Fasops::Common::parse_axt_block( $content, $query_length_of );
            $content = '';

            next if length( $info_refs->[0]{seq} ) < $threshold;

            $info_refs->[0]{name} = $opt->{tname};
            $info_refs->[1]{name} = $opt->{qname};

            $self->add_align($info_refs);
        }
        else {
            $content .= $line;
        }
    }
    $in_fh->close;

    return;
}

# blocked fasta format
sub parse_fas_file {
    my $self   = shift;
    my $infile = shift;
    my $opt    = shift;

    # minimal length
    my $threshold = $opt->{threshold};
    $threshold ||= $self->threshold;

    my $in_fh = IO::Zlib->new( $infile, "rb" );

    my $content = '';    # content of one block
    while (1) {
        last if $in_fh->eof and $content eq '';
        my $line = '';
        if ( !$in_fh->eof ) {
            $line = $in_fh->getline;
        }
        if ( ( $line eq '' or $line =~ /^\s+$/ ) and $content ne '' ) {
            my $info_of = App::Fasops::Common::parse_block($content);
            $content = '';

            my $info_refs = [ map { $info_of->{$_} } keys %{$info_of} ];
            next if length( $info_refs->[0]{seq} ) < $threshold;

            $self->add_align($info_refs);
        }
        else {
            $content .= $line;
        }
    }

    $in_fh->close;

    return;
}

sub _insert_indel {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
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

    my AlignDB::IntSpan( $align_set, undef, $indel_set ) = @{ $self->get_sets($align_id) };

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
        my @uniq_indel_seqs = List::MoreUtils::PP::uniq(@indel_seqs);

        # seqs with least '-' char wins
        my ($indel_seq) = map { $_->[0] }
            sort { $a->[1] <=> $b->[1] }
            map { [ $_, tr/-/-/ ] } @uniq_indel_seqs;

        if ( scalar @uniq_indel_seqs < 2 ) {
            warn "no indel!\n";
            next;
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
        $indel_freq = List::Util::min( $indel_freq, $seq_count - $indel_freq );

        my $indel_gc = App::Fasops::Common::calc_gc_ratio( [$indel_seq] );

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
        $sth->execute(
            $prev_indel_id,     $align_id, $_->{start},    $_->{end},
            $_->{length},       $_->{seq}, $_->{all_seqs}, $_->{left_extand},
            $_->{right_extand}, $_->{gc},  $_->{freq},     $_->{occured},
            $_->{type},
        );
        ($prev_indel_id) = $self->last_insert_id;
    }

    $sth->finish;

    return;
}

sub _insert_snp {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
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

    my $seq_refs  = $self->get_seqs($align_id);
    my $seq_count = scalar @{$seq_refs};
    my $length    = length $seq_refs->[0];

    my $snp_site = {};
    for my $pos ( 1 .. $length ) {
        my @bases;
        for my $i ( 0 .. $seq_count - 1 ) {
            my $base = substr( $seq_refs->[$i], $pos - 1, 1 );
            push @bases, $base;
        }

        if ( List::MoreUtils::PP::all { $_ =~ /[agct]/i } @bases ) {
            if ( List::MoreUtils::PP::any { $_ ne $bases[0] } @bases ) {
                $snp_site->{$pos} = \@bases;
            }
        }
    }

    # %{$snp_site} keys are snp positions
    for my $pos ( sort { $a <=> $b } keys %{$snp_site} ) {

        my @bases = @{ $snp_site->{$pos} };

        my $target_base = $bases[0];
        my $all_bases = join '', @bases;

        my $query_base;
        my $mutant_to;
        my $snp_freq = 0;
        my $snp_occured;
        my @class = List::MoreUtils::PP::uniq(@bases);
        if ( scalar @class < 2 ) {
            Carp::confess "no snp\n";
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
        $snp_freq = List::Util::min( $snp_freq, $seq_count - $snp_freq );

        $sth->execute( $align_id, $pos, $target_base, $query_base,
            $all_bases, $mutant_to, $snp_freq, $snp_occured, );
    }
    $sth->finish;

    return;
}

sub insert_isw {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;

    my ( $align_set, undef, undef )
        = @{ $self->get_sets($align_id) };

    # indel_id & prev_indel_id
    my DBI $fetch_indel_id_isw = $dbh->prepare(
        q{
        SELECT indel_id, prev_indel_id
        FROM indel
        WHERE align_id = ?
        }
    );
    $fetch_indel_id_isw->execute($align_id);

    # indel_end
    my DBI $fetch_prev_indel_end = $dbh->prepare(
        q{
        SELECT indel_end
        FROM indel
        WHERE indel_id = ?
        }
    );

    # prev_indel_start
    my DBI $fetch_indel_start = $dbh->prepare(
        q{
        SELECT indel_start
        FROM indel
        WHERE indel_id = ?
        }
    );

    # prepare insert isw
    my DBI $sth = $dbh->prepare(
        q{
        INSERT INTO isw (
            isw_id, indel_id, prev_indel_id, align_id, isw_indel_id,
            isw_start, isw_end, isw_length, isw_type,
            isw_distance, isw_density, isw_differences, isw_pi,
            isw_target_gc, isw_average_gc
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?, ?, ?,
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
            print YAML::Syck::Dump(
                {   interval_start  => $interval_start,
                    interval_end    => $interval_end,
                    interval_length => $interval_length,
                }
            );
            print "start $interval_start > end $interval_end.\n";
            next;
        }

        my AlignDB::Window $window_maker = $self->window_maker;

        my @isws = $window_maker->interval_window( $align_set, $interval_start, $interval_end );

        for my $isw (@isws) {
            my AlignDB::IntSpan $isw_set = $isw->{set};
            my $isw_start                = $isw_set->min;
            my $isw_end                  = $isw_set->max;
            my $isw_length               = $isw_end - $isw_start + 1;

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
            $sth->execute(
                $indel_id,        $prev_indel_id,  $align_id,      $isw_indel_id,
                $isw_start,       $isw_end,        $isw_length,    $isw->{type},
                $isw->{distance}, $isw->{density}, $isw_stat->[3], $isw_stat->[7],
                $isw_stat->[8],   $isw_stat->[9],
            );
        }
    }

    $sth->finish;
    $fetch_prev_indel_end->finish;
    $fetch_indel_start->finish;
    $fetch_indel_id_isw->finish;

    return;
}

sub isw_snp_fk {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh          = $self->dbh;
    my DBI $fetch_snp_id = $dbh->prepare(
        q{
        SELECT s.snp_id, s.snp_pos
        FROM snp s
        WHERE s.align_id = ?
        }
    );

    my DBI $fetch_isw_id = $dbh->prepare(
        q{
        SELECT isw_id
        FROM isw
        WHERE 1 = 1
        AND align_id = ?
        AND isw_start <= ?
        AND isw_end >= ?
        }
    );

    my DBI $sth = $dbh->prepare(
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
            $sth->execute( $isw_id, $snp_id );
        }
    }

    $sth->finish;
    $fetch_isw_id->finish;
    $fetch_snp_id->finish;

    return;
}

sub insert_ssw {
    my $self     = shift;
    my $align_id = shift;

    my AlignDB::IntSpan( undef, $comparable_set, undef )
        = @{ $self->get_sets($align_id) };

    my DBI $dbh = $self->dbh;

    # get sliding windows' sizes
    my AlignDB::Window $window_maker = $self->window_maker;
    my $ssw_size = $window_maker->sw_size;

    my $ssw_max_distance = 10;
    my $ssw_size_window0 = int( $ssw_size / 2 );

    my DBI $fetch_snp_id = $dbh->prepare(
        q{
        SELECT s.snp_id, s.snp_pos, s.snp_occured
        FROM snp s
        WHERE s.align_id = ?
        }
    );
    $fetch_snp_id->execute($align_id);

    # prepare ssw_insert
    my DBI $sth = $dbh->prepare(
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
        my $snp_id = $snp_info->{$snp_index}{snp_id};

        #        my $snp_pos     = $snp_info->{$snp_index}{snp_pos};
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
                my $ssw_set_member_number = $ssw_set->size;
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

                my ($cur_window_id) = $self->insert_window( $align_id, $ssw_set );

                $sth->execute( $snp_id, $cur_window_id, $ssw_type,
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

        my DBI $dbh = $self->dbh;

        my DBI $target_sth = $dbh->prepare(
            q{
            SELECT s.seq_seq
            FROM sequence s
            INNER JOIN target t on s.seq_id = t.seq_id
            WHERE s.align_id = ?
            }
        );

        my DBI $query_sth = $dbh->prepare(
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

sub get_ref_seq {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
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

# ( $target_name, $query_name, $ref_name )
sub get_names {
    my $self = shift;
    my $align_id = shift || 1;

    my DBI $dbh = $self->dbh;

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
                chromosome c ON c.chr_id = s.chr_id
            WHERE
                s.align_id = ?
        };

        $sql =~ s/_TABLE_/$table/g;

        $sql .= "ORDER BY query.query_position" if $table eq 'query';

        my DBI $sth = $dbh->prepare($sql);
        $sth->execute($align_id);
        while ( my ($name) = $sth->fetchrow_array ) {
            push @names, $name;
        }
        $sth->finish;
    }

    return @names;
}

#@returns AlignDB::IntSpan
sub get_sets {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;

    my DBI $sth = $dbh->prepare(
        q{
        SELECT align_length, align_comparable_runlist, align_indel_runlist
        FROM align 
        WHERE align_id = ?
        }
    );
    $sth->execute($align_id);
    my ( $align_length, $comparable_runlist, $indel_runlist ) = $sth->fetchrow_array;
    $sth->finish;

    my $align_set      = AlignDB::IntSpan->new("1-$align_length");
    my $comparable_set = AlignDB::IntSpan->new($comparable_runlist);
    my $indel_set      = AlignDB::IntSpan->new($indel_runlist);

    return [ $align_set, $comparable_set, $indel_set ];
}

sub get_chr_legnth_of {
    my $self        = shift;
    my $common_name = shift;

    my %length_of = ();
    my DBI $dbh   = $self->dbh;
    my DBI $std   = $dbh->prepare(
        q{
        SELECT * FROM chromosome WHERE common_name = ?
        }
    );
    $std->execute($common_name);
    while ( my $ref = $std->fetchrow_hashref ) {
        $length_of{ $ref->{chr_name} } = $ref->{chr_length};
    }
    $std->finish;

    return \%length_of;
}

sub get_slice_stat {
    my $self                 = shift;
    my $align_id             = shift;
    my AlignDB::IntSpan $set = shift;

    my $seqs_ref   = $self->get_seqs($align_id);
    my @seq_slices = map { $set->substr_span($_) } @$seqs_ref;
    my $result     = App::Fasops::Common::multi_seq_stat( \@seq_slices );
    my $target_gc  = App::Fasops::Common::calc_gc_ratio( [ $seq_slices[0] ] );

    push @{$result}, $target_gc;
    return $result;
}

sub get_slice_indel {
    my $self                 = shift;
    my $align_id             = shift;
    my AlignDB::IntSpan $set = shift;

    # slice indel; every gap is treated as an indel
    my $set_indel = $set->span_size - 1;

    # real indels in this alignment slice
    my $seqs_ref       = $self->get_seqs($align_id);
    my @seq_slices     = map { $set->substr_span($_) } @$seqs_ref;
    my @seq_indel_sets = map { App::Fasops::Common::indel_intspan($_) } @seq_slices;

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

# internal method
# my $window_id = $self->insert_window(
#     $align_id, $window_set, $internal_indel,
# );
sub insert_window {
    my $self                        = shift;
    my $align_id                    = shift;
    my AlignDB::IntSpan $window_set = shift;
    my $internal_indel              = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        INSERT INTO window (
            window_id, align_id, window_start, window_end, window_length,
            window_runlist, window_comparables, window_identities,
            window_differences, window_indel, window_pi,
            window_gc
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?
        )
        }
    );

    my $window_start   = $window_set->min;
    my $window_end     = $window_set->max;
    my $window_length  = $window_set->size;
    my $window_runlist = $window_set->runlist;

    # do or do not count internal indels within window_set
    # $set_indel is equal to $window_span - 1
    my $window_indel;
    if ($internal_indel) {
        my ( $set_indel, $real_indel ) = $self->get_slice_indel( $align_id, $window_set );
        $window_indel = $set_indel + $real_indel;
    }
    else {
        $window_indel = $window_set->span_size - 1;

        # avoid empty sets which have span_size of 0
        $window_indel = 0 if $window_indel < 0;
    }

    my $window_stat = $self->get_slice_stat( $align_id, $window_set );
    $sth->execute(
        $align_id,       $window_start,     $window_end,       $window_length,
        $window_runlist, $window_stat->[1], $window_stat->[2], $window_stat->[3],
        $window_indel,   $window_stat->[7], $window_stat->[8],
    );

    return $self->last_insert_id;
}

sub last_insert_id {
    my $self = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT LAST_INSERT_ID()
        }
    );
    $sth->execute;
    my ($last_insert_id) = $sth->fetchrow_array;

    return $last_insert_id;
}

sub execute_sql {
    my ( $self, $sql_query, $bind_value ) = @_;

    # init
    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare($sql_query);

    # bind value
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    $sth->execute(@$bind_value);
}

# Clean a table for new insertion
sub empty_table {
    my ( $self, $table, $with_window ) = @_;

    my DBI $dbh = $self->dbh;

    # check table existing
    my @table_names = $dbh->tables( '', '', '' );

    # returned table names are quoted by `` (back-quotes)
    unless ( List::MoreUtils::PP::any { $_ =~ qr{`$table`} } @table_names ) {
        print "Table $table does not exist\n";
        return;
    }

    if ( !$with_window ) {
        $dbh->do(
            qq{
            DELETE FROM $table
            }
        );
    }
    else {

        # In MySQL 4.1, you must use the alias (if one was given) when referring to a table name
        $dbh->do(
            qq{
            DELETE t,
                   w
            FROM $table t,
                 window w
            WHERE t.window_id = w.window_id
            }
        );
    }

    # set AUTO_INCREMENT to 1 for this table
    # MyISAM table will memory last increment even this table is emptied
    $dbh->do("ALTER TABLE $table AUTO_INCREMENT = 1");

    return;
}

# $self->create_column(
#     $table,
#     $column,
#     $column_definition
# );
# Add $column to $table with $column_definition as properties
sub create_column {
    my ( $self, $table, $column, $column_definition ) = @_;

    $column_definition ||= "DOUBLE";

    my DBI $dbh = $self->dbh;

    # check table existing
    {
        my @table_names = $dbh->tables( '', '', '' );

        # table names are quoted by ` (back-quotes) which is the
        #   quote_identifier
        my $table_name = "`$table`";
        unless ( List::MoreUtils::PP::any { $_ =~ /$table_name/i } @table_names ) {
            print "Table $table does not exist\n";
            return;
        }
    }

    # check column existing
    # then create column
    {
        my DBI $sth = $dbh->prepare(
            qq{
            SHOW FIELDS
            FROM $table
            LIKE "$column"
            }
        );
        $sth->execute;
        my ($field) = $sth->fetchrow_array;

        if ($field) {
            $dbh->do(
                qq{
                ALTER TABLE $table DROP COLUMN $column
                }
            );
        }

        $dbh->do(
            qq{
            ALTER TABLE $table ADD COLUMN $column $column_definition
            }
        );
    }

    return;
}

sub check_column {
    my ( $self, $table, $column ) = @_;

    # init
    my DBI $dbh = $self->dbh;

    # check table existing
    {
        my @table_names = $dbh->tables( '', '', '' );

        # table names are quoted by ` (back-quotes) which is the
        #   quote_identifier
        my $table_name = "`$table`";
        unless ( List::MoreUtils::PP::any { $_ =~ /$table_name/i } @table_names ) {
            print " " x 4, "Table $table does not exist\n";
            return;
        }
    }

    # check column existing
    {
        my DBI $sth = $dbh->prepare(
            qq{
            SHOW FIELDS
            FROM $table
            LIKE "$column"
            }
        );
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
        my DBI $sth = $dbh->prepare(
            qq{
            SELECT COUNT($column)
            FROM $table
            }
        );
        $sth->execute;
        my ($count) = $sth->fetchrow_array;

        if ( not $count ) {
            print " " x 4, "Column $column has no records\n";
        }

        return $count;
    }
}

# get the chr_name and chr_length of a given chr_id
sub get_chr_info {
    my $self   = shift;
    my $chr_id = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        qq{
        SELECT c.chr_name, c.chr_length
        FROM chromosome c
        WHERE c.chr_id = ?
        }
    );
    $sth->execute($chr_id);
    my ( $chr_name, $chr_length ) = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_name, $chr_length );
}

# get an array_ref of align_ids
sub get_align_ids {
    my $self = shift;

    my DBI $dbh = $self->dbh;
    my $align_ids = $dbh->selectcol_arrayref(
        q{
        SELECT a.align_id
        FROM align a
        }
    );

    return $align_ids;
}

sub get_align_ids_of_chr_name {
    my $self     = shift;
    my $chr_name = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT a.align_id
        FROM sequence s
        INNER JOIN target t ON s.seq_id = t.seq_id
        INNER JOIN align a ON s.align_id = a.align_id
        INNER JOIN chromosome c on s.chr_id = c.chr_id
        WHERE c.chr_name = ?
        ORDER BY a.align_id
        }
    );
    $sth->execute($chr_name);

    my @align_ids;
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
        my DBI $dbh = $self->dbh;
        my DBI $sth = $dbh->prepare(
            q{
            SELECT c.common_name,
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
            }
        );
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

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT c.common_name,
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
        }
    );
    $sth->execute($align_id);

    my @array;
    while ( my $hash_ref = $sth->fetchrow_hashref ) {
        push @array, $hash_ref;
    }
    $sth->finish;

    return @array;
}

sub get_target_chr_info {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT c.chr_id, c.chr_name, c.chr_length
        FROM sequence s
        INNER JOIN target t ON s.seq_id = t.seq_id
        INNER JOIN chromosome c ON s.chr_id = c.chr_id
        WHERE s.align_id = ?
        }
    );
    $sth->execute($align_id);
    my ( $chr_id, $chr_name, $chr_length ) = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_id, $chr_name, $chr_length );
}

sub get_query_chr_info {
    my $self     = shift;
    my $align_id = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT c.chr_id, c.chr_name, c.chr_length
        FROM sequence s
        INNER JOIN query q ON s.seq_id = q.seq_id
        INNER JOIN chromosome c ON s.chr_id = c.chr_id
        WHERE s.align_id = ?
        }
    );
    $sth->execute($align_id);

    my ( $chr_id, $chr_name, $chr_length ) = $sth->fetchrow_array;
    $sth->finish;

    return ( $chr_id, $chr_name, $chr_length );
}

sub get_chrs {
    my $self = shift;
    my $goal = shift || 'target';

    my DBI $dbh = $self->dbh;

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
    my DBI $sth = $dbh->prepare($chr_query);
    $sth->execute;

    my $array_ref = $sth->fetchall_arrayref;
    $sth->finish;

    return $array_ref;
}

sub get_freq {
    my $self = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT DISTINCT COUNT(q.query_id) + 1
        FROM  query q, sequence s
        WHERE q.seq_id = s.seq_id
        GROUP BY s.align_id
        }
    );
    $sth->execute;

    my @counts;
    while ( my ($count) = $sth->fetchrow_array ) {
        push @counts, $count;
    }
    if ( scalar @counts > 1 ) {
        warn "Database corrupts, freqs are not consistent\n";
    }

    return $counts[0];
}

sub find_align {
    my $self     = shift;
    my $chr_name = shift;
    my $start    = shift;
    my $end      = shift || $start;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT s.align_id
        FROM sequence s
        INNER JOIN target t ON s.seq_id = t.seq_id
        WHERE 1 = 1
        AND s.chr_name = ?
        AND s.chr_start <= ?
        AND s.chr_end >= ?
        }
    );
    $sth->execute( $chr_name, $start, $end );

    my @align_ids;
    while ( my @row = $sth->fetchrow_array ) {
        push @align_ids, $row[0];
    }

    return \@align_ids;
}

sub find_align_chr_id {
    my $self   = shift;
    my $chr_id = shift;
    my $start  = shift;
    my $end    = shift || $start;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        SELECT s.align_id
        FROM sequence s
        INNER JOIN target t ON s.seq_id = t.seq_id
        WHERE 1 = 1
        AND s.chr_id = ?
        AND s.chr_start <= ?
        AND s.chr_end >= ?
        }
    );
    $sth->execute( $chr_id, $start, $end );

    my @align_ids;
    while ( my @row = $sth->fetchrow_array ) {
        push @align_ids, $row[0];
    }

    return \@align_ids;
}

sub process_message {
    my $self     = shift;
    my $align_id = shift;

    my $info = $self->get_target_info($align_id);

    printf "Process align [%s] at %s.%s(%s):%s-%s\n", $align_id, $info->{common_name},
        $info->{chr_name},
        $info->{chr_strand}, $info->{chr_start}, $info->{chr_end};

    return;
}

sub add_meta {
    my $self      = shift;
    my $meta_hash = shift;

    my DBI $dbh = $self->dbh;
    my DBI $sth = $dbh->prepare(
        q{
        INSERT INTO meta ( meta_id, meta_key, meta_value )
        VALUES ( NULL, ?, ? )
        }
    );

    for my $key ( sort keys %$meta_hash ) {
        $sth->execute( $key, $meta_hash->{$key} );
    }

    $sth->finish;

    return;
}

sub add_meta_stopwatch {
    my $self = shift;
    my AlignDB::Stopwatch $stopwatch = shift;

    my DBI $dbh = $self->dbh;
    my $ary_ref = $dbh->selectcol_arrayref(
        q{
        SELECT meta_value
        from meta
        where meta_key = 'uuid'
        }
    );

    my $uuid = $stopwatch->uuid;
    if ( !List::MoreUtils::PP::any { $_ eq $uuid } @$ary_ref ) {
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

=head1 SYNOPSIS

    use AlignDB;
    my $obj = AlignDB->new(
        mysql          => "$db:$server",
        user           => $username,
        passwd         => $password,
    );

=head1 AUTHOR

B<Qiang Wang>

=cut
