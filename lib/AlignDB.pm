package AlignDB;
use Moose;
use Carp;
use autodie;
use DBI;

use IO::Zlib;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any uniq);
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Window;
use AlignDB::Util qw(:all);

has 'mysql'  => ( is => 'ro', isa => 'Str' );      # e.g. 'alignDB:202.119.43.5'
has 'server' => ( is => 'ro', isa => 'Str' );      # e.g. '202.119.43.5'
has 'db'     => ( is => 'ro', isa => 'Str' );      # e.g. 'alignDB'
has 'user'   => ( is => 'ro', isa => 'Str' );      # database username
has 'passwd' => ( is => 'ro', isa => 'Str' );      # database password
has 'dbh'    => ( is => 'ro', isa => 'Object' );   # store database handle here
has 'window_maker' => ( is => 'ro', isa => 'Object' );   # sliding windows maker
has 'threshold' => ( is => 'ro', isa => 'Int', default => 10_000 );

has 'caching_id' => ( is => 'ro', isa => 'Int' );        # caching seqs
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

    my $window_maker = AlignDB::Window->new;
    $self->{window_maker} = $window_maker;

    return;
}

sub _insert_align {
    my ( $self, $target_seq_ref, $query_seq_ref ) = @_;

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

    my $result = pair_seq_stat( $$target_seq_ref, $$query_seq_ref );

    $align_insert->execute(
        $result->[0], $result->[1], $result->[2], $result->[3], $result->[4],
        $result->[5], $result->[6], $result->[7], $result->[8], $result->[9],
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
            isw_pi, isw_target_gc, isw_average_gc,
            isw_d_indel, isw_d_noindel, isw_d_complex
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
            NULL, NULL, NULL
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
                $isw_stat->[7], $isw_stat->[8],   $isw_stat->[9],
            );
        }
    }

    $isw_insert->finish;
    $fetch_prev_indel_end->finish;
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
        = @{ $self->get_sets($align_id) };
    my $ref_info = $self->caching_refs;

    my %snp_site = %{ pair_snp_sites( $target_seq, $query_seq ) };
    my $snp_insert = $dbh->prepare(
        q{
        INSERT INTO snp (
            snp_id, isw_id, align_id, snp_pos,
            target_base, query_base, ref_base, snp_occured
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?, ?, ?
        )
        }
    );

    # isw_id
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

    # %snp_site keys are snp positions,
    #   where is generated at align section
    for ( sort { $a <=> $b } keys %snp_site ) {
        $fetch_isw_id->execute( $align_id, $_, $_ );
        my ($isw_id)    = $fetch_isw_id->fetchrow_array;
        my $ref_base    = undef;
        my $snp_occured = undef;
        if ( defined $ref_info->{seq} ) {
            $ref_base = substr( $ref_info->{seq}, $_ - 1, 1 );
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

    # extreme_id & prev_extreme_id
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
        my $snp_index = $comparable_set->lookup_member($snp_pos);

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

            # $gsw_distance is from 0 to $gsw_density
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
    my $self   = shift;
    my $infile = shift;
    my $opt    = shift;

    my $target_taxon_id = $opt->{target_taxon_id};
    my $query_taxon_id  = $opt->{query_taxon_id};
    my $threshold       = $opt->{threshold};
    my $gzip            = $opt->{gzip};

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
        unless ($summary_line) {
            last;
        }
        if ( $summary_line =~ /^#/ ) {
            next;
        }
        chomp $summary_line;
        chomp( my $first_line = <$in_fh> );
        $first_line = uc $first_line;
        chomp( my $second_line = <$in_fh> );
        $second_line = uc $second_line;
        my $dummy = <$in_fh>;

        unless ( length $first_line > $threshold ) {
            next;
        }

        my ($align_serial, $first_chr,    $first_start,
            $first_end,    $second_chr,   $second_start,
            $second_end,   $query_strand, $align_score,
        ) = split /\s+/, $summary_line;

        my $target_info = {
            taxon_id   => $target_taxon_id,
            chr_name   => $first_chr,
            chr_start  => $first_start,
            chr_end    => $first_end,
            chr_strand => '+',
            seq        => $first_line,
        };
        my $query_info = {
            taxon_id     => $query_taxon_id,
            chr_name     => $second_chr,
            chr_start    => $second_start,
            chr_end      => $second_end,
            chr_strand   => $query_strand,
            query_strand => $query_strand,
            seq          => $second_line,
        };

        $self->add_align( $target_info, $query_info );
    }

    if ( !$gzip ) {
        close $in_fh;
    }
    else {
        $in_fh->close;
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

    my $target_seq = $target_info->{seq};
    my $query_seq  = $query_info->{seq};

    my ( $ref_seq, $ref_raw_seq );
    if ( defined $ref_info->{seq} ) {
        $ref_seq     = $ref_info->{seq};
        $ref_raw_seq = $ref_info->{raw_seq};
    }

    my $dbh = $self->dbh;

    #----------------------------#
    # fill empty key
    #----------------------------#
    for my $key (qw{chr_name chr_start chr_end chr_strand}) {
        if ( !defined $query_info->{$key} ) {
            $query_info->{$key} = undef;
        }
        if ($ref_seq) {
            if ( !defined $ref_info->{$key} ) {
                $ref_info->{$key} = undef;
            }
        }
    }
    $query_info->{query_strand} ||= "+";

    {
        my $target_chr_id_of
            = $self->get_chr_id_hash( $target_info->{taxon_id} );
        if ( exists $target_chr_id_of->{ $target_info->{chr_name} } ) {
            $target_info->{chr_id}
                = $target_chr_id_of->{ $target_info->{chr_name} };
        }
        else {
            $target_info->{chr_id} = $target_chr_id_of->{'chrUn'};
        }

        my $query_chr_id_of = $self->get_chr_id_hash( $query_info->{taxon_id} );
        if ( defined $query_info->{chr_name}
            and exists $query_chr_id_of->{ $query_info->{chr_name} } )
        {
            $query_info->{chr_id}
                = $query_chr_id_of->{ $query_info->{chr_name} };
        }
        else {
            $query_info->{chr_id} = $query_chr_id_of->{'chrUn'};
        }

        if ($ref_seq) {
            my $ref_chr_id_of = $self->get_chr_id_hash( $ref_info->{taxon_id} );
            if ( defined $ref_info->{chr_name}
                and exists $ref_chr_id_of->{ $ref_info->{chr_name} } )
            {
                $ref_info->{chr_id} = $ref_chr_id_of->{ $ref_info->{chr_name} };
            }
            else {
                $ref_info->{chr_id} = $ref_chr_id_of->{'chrUn'};
            }
        }
    }

    # check align length
    my $align_length = length $target_seq;
    if ( length $query_seq != $align_length ) {
        print "The two sequences has not equal length.\n";
        return;
    }

    #----------------------------#
    # INSERT INTO align
    #----------------------------#
    my $align_id = $self->_insert_align( \$target_seq, \$query_seq );
    printf "Prosess align %s in %s %s - %s\n", $align_id,
        $target_info->{chr_name}, $target_info->{chr_start},
        $target_info->{chr_end};

    #----------------------------#
    # INSERT INTO target, query
    #----------------------------#
    my $align_set = AlignDB::IntSpan->new("1-$align_length");
    {    # target
        my $target_insert = $dbh->prepare(
            q{
            INSERT INTO target ( target_id, seq_id )
            VALUES ( NULL, ? )
            }
        );
        $target_info->{align_id} = $align_id;
        $target_info->{gc}       = calc_gc_ratio($target_seq);
        my $seq_set = $align_set->diff( find_indel_set($target_seq) );
        $target_info->{runlist} = $seq_set->runlist;
        $target_info->{length}  = $seq_set->cardinality;
        my $target_seq_id = $self->_insert_seq($target_info);
        $target_insert->execute($target_seq_id);
        $target_insert->finish;
    }

    {    # and queries
        my $query_insert = $dbh->prepare(
            q{
            INSERT INTO query (
                query_id, seq_id, query_strand, query_position
            )
            VALUES ( NULL, ?, ?, 1 )
            }
        );
        $query_info->{align_id} = $align_id;
        $query_info->{gc}       = calc_gc_ratio($query_seq);
        my $seq_set = $align_set->diff( find_indel_set($query_seq) );
        $query_info->{runlist} = $seq_set->runlist;
        $query_info->{length}  = $seq_set->cardinality;
        my $query_seq_id = $self->_insert_seq($query_info);
        $query_insert->execute( $query_seq_id, $query_info->{query_strand} );
        $query_insert->finish;
    }

    #----------------------------#
    # INSERT INTO reference
    #----------------------------#
    if ($ref_seq) {
        my $ref_complex_indel = $ref_info->{complex};
        my $ref_insert        = $dbh->prepare(
            q{
            INSERT INTO reference (
                ref_id, seq_id, ref_raw_seq, ref_complex_indel
            )
            VALUES (
                NULL, ?, ?, ?
            )
            }
        );
        $ref_info->{align_id} = $align_id;
        $ref_info->{gc}       = calc_gc_ratio($ref_seq);
        my $seq_set = $align_set->diff( find_indel_set($ref_seq) );
        $ref_info->{runlist} = $seq_set->runlist;
        $ref_info->{length}  = $seq_set->cardinality;
        my $ref_seq_id = $self->_insert_seq($ref_info);
        $ref_insert->execute( $ref_seq_id, $ref_raw_seq, $ref_complex_indel );
        $ref_insert->finish;
    }

    #----------------------------#
    # UPDATE align with runlist
    #----------------------------#
    {
        my $align_set      = AlignDB::IntSpan->new("1-$align_length");
        my $target_gap_set = $align_set->diff( $target_info->{runlist} );
        my $query_gap_set  = $align_set->diff( $query_info->{runlist} );
        my $indel_set      = $target_gap_set->union($query_gap_set);
        my $comparable_set = $align_set->diff($indel_set);

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

    #----------------------------#
    # INSERT INTO indel
    #----------------------------#
    {
        my $indel_insert = $dbh->prepare(
            'INSERT INTO indel (
                indel_id, prev_indel_id, align_id,
                indel_start, indel_end, indel_length,
                indel_seq, left_extand, right_extand,
                indel_gc, indel_freq, indel_occured, indel_type
            )
            VALUES (
                NULL, ?, ?,
                ?, ?, ?,
                ?, ?, ?,
                ?, ?, ?, ?
            )'
        );

        my @indel_site
            = @{ pair_indel_sites( $target_seq, $query_seq, $all_indel ) };

        my $complex_set = AlignDB::IntSpan->new;
        if ( defined $ref_info->{seq} ) {
            $complex_set->add( $ref_info->{complex} );
        }

        my $prev_indel_id = 0;
        for (@indel_site) {

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
            my $indel_occured   = undef;
            my $indel_type      = undef;
            my $indel_frequency = undef;
            if ( defined $ref_info->{seq} ) {
                my $start  = $_->{start};
                my $end    = $_->{end};
                my $length = $_->{length};

                if ( $complex_set->superset("$start-$end") ) {
                    $indel_occured   = "C";
                    $indel_type      = "C";
                    $indel_frequency = -1;
                }
                else {
                    my $r_str = substr( $ref_seq,    $start - 1, $length );
                    my $t_str = substr( $target_seq, $start - 1, $length );
                    my $q_str = substr( $query_seq,  $start - 1, $length );
                    ( $indel_occured, $indel_type )
                        = ref_indel_type( $r_str, $t_str, $q_str );
                    $indel_frequency = 1;
                }
            }
            $indel_insert->execute(
                $prev_indel_id,    $align_id,          $_->{start},
                $_->{end},         $_->{length},       $_->{seq},
                $_->{left_extand}, $_->{right_extand}, $_->{gc},
                $indel_frequency,  $indel_occured,     $indel_type,
            );
            ($prev_indel_id) = $self->last_insert_id;
        }
        $indel_insert->finish;
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
    if ( defined $ref_info->{seq} ) {
        $self->_modify_isw($align_id);
    }

    return $align_id;
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
        my $query = qq{
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

        $query =~ s/_TABLE_/$table/g;

        my $sth = $dbh->prepare($query);
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
    my $result     = pair_seq_stat(@seq_slices);

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
            window_target_gc, window_average_gc,
            window_coding, window_repeats, window_ns_indel
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?,
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

    return $hash_ref;
}

sub get_query_info {
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
    my $hash_ref = $sth->fetchrow_hashref;
    $sth->finish;

    return $hash_ref;
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
        die "Database corrupts, freqs are not consistent\n";
    }

    return $counts[0];
}

sub process_message {
    my $self     = shift;
    my $align_id = shift;

    my $target_info = $self->get_target_info($align_id);

    printf "Prosess align %s in %s %s - %s\n", $align_id,
        $target_info->{chr_name}, $target_info->{chr_start},
        $target_info->{chr_end};

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


