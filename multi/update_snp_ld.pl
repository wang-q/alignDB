#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all);
use Scalar::Util qw(looks_like_number);

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Multi;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update LD of $db...");

#----------------------------#
# Add columns
#----------------------------#
my $all_freq;
my @jobs;
{
    my $obj = AlignDB::Multi->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    my $dbh = $obj->dbh;

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
    $all_freq = $counts[0];
    print "all_freq is $all_freq\n";

    # r and D' with nearest indel
    $obj->create_column( "snp", "snp_r",      "DOUBLE" );
    $obj->create_column( "snp", "snp_dprime", "DOUBLE" );

    #  r and D' with nearest snp
    ##$obj->create_column( "snp", "snp_nearest_snp_id", "int" );
    #$obj->create_column( "snp", "snp_r_s",             "DOUBLE" );
    #$obj->create_column( "snp", "snp_dprime_s",        "DOUBLE" );

    # r2 and |D'| with near snps
    # include self
    $obj->create_column( "snp", "snp_near_snp_number", "int" );
    $obj->create_column( "snp", "snp_r2_s",            "DOUBLE" );
    $obj->create_column( "snp", "snp_dprime_abs_s",    "DOUBLE" );

    ## indel group
    #$obj->create_column( "snp", "snp_r_i",      "DOUBLE" );
    #$obj->create_column( "snp", "snp_dprime_i", "DOUBLE" );
    #
    ## nonindel group
    #$obj->create_column( "snp", "snp_r_ni",      "DOUBLE" );
    #$obj->create_column( "snp", "snp_dprime_ni", "DOUBLE" );

    print "Table snp altered\n";

    my @align_ids = @{ $obj->get_align_ids };

    while ( scalar @align_ids ) {
        my @batching = splice @align_ids, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $worker = sub {
    my $job       = shift;
    my @align_ids = @$job;

    my $obj = AlignDB::Multi->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # select all indels in this alignment
    my $indel_query = q{
        SELECT i.indel_id, i.indel_freq, i.indel_occured
        FROM indel i
        WHERE 1=1
        AND i.indel_type != 'C'
        AND i.align_id = ?
    };
    my $indel_sth = $dbh->prepare($indel_query);

    # select all isw in this indel
    my $snp_query = q{
        SELECT s.snp_id, s.snp_freq, s.snp_occured, s.snp_pos
        FROM isw w, snp s
        WHERE 1 = 1
        AND w.isw_id = s.isw_id
        AND s.snp_occured != 'unknown'
        AND w.isw_indel_id = ?
    };
    my $snp_sth = $dbh->prepare($snp_query);

    # update snp table in the new feature columns
    #my $update_query = q{
    #    UPDATE snp s
    #    SET s.snp_r         = ?,
    #        s.snp_dprime    = ?,
    #        s.snp_r_i       = ?,
    #        s.snp_dprime_i  = ?,
    #        s.snp_r_ni      = ?,
    #        s.snp_dprime_ni = ?
    #    WHERE s.snp_id = ?
    #};
    #my $update_query = q{
    #    UPDATE
    #        snp s
    #    SET
    #        s.snp_r         = ?,
    #        s.snp_dprime    = ?,
    #        s.snp_nearest_snp_id = ?,
    #        s.snp_r_s       = ?,
    #        s.snp_dprime_s  = ?
    #    WHERE s.snp_id = ?
    #};
    my $update_query = q{
        UPDATE
            snp s
        SET
            s.snp_r         = ?,
            s.snp_dprime    = ?,
            s.snp_near_snp_number = ?,
            s.snp_r2_s       = ?,
            s.snp_dprime_abs_s  = ?
        WHERE s.snp_id = ?
    };
    my $update_sth = $dbh->prepare($update_query);

    # for each align
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);

        my $all_snps = $dbh->selectall_arrayref(
            q{
            SELECT s.snp_id, s.snp_freq, s.snp_occured, s.snp_pos
            FROM snp s
            WHERE 1 = 1
            AND s.snp_occured != 'unknown'
            AND s.align_id = ?
            }, {}, $align_id
        );

        #print Dump $all_snps;
        #last;

        $indel_sth->execute($align_id);
        while ( my @row = $indel_sth->fetchrow_array ) {
            my ( $indel_id, $indel_freq, $indel_occured ) = @row;

            #my $group_i  = AlignDB::IntSpan->new;
            #my $group_ni = AlignDB::IntSpan->new;
            #for my $i ( 0 .. $all_freq - 1 ) {
            #    my $ichar = substr $indel_occured, $i, 1;
            #    if ( $ichar eq 'o' ) {
            #        $group_i->add( $i + 1 );
            #    }
            #    elsif ( $ichar eq 'x' ) {
            #        $group_ni->add( $i + 1 );
            #    }
            #    else {
            #        die "$i $ichar\n";
            #    }
            #}

            $snp_sth->execute($indel_id);
            while ( my @row = $snp_sth->fetchrow_array ) {
                my ( $snp_id, $snp_freq, $snp_occured, $snp_pos ) = @row;

                # LD with nearest indel
                my ( $r, $dprime );
                ( $r, $dprime ) = calc_ld( $indel_occured, $snp_occured );

                #$update_sth->execute( $r, $dprime, $snp_id );

                # LD with nearest snp
                #my ( $nearest_snp_id, $r_s, $dprime_s );
                #if ( scalar @{$all_snps} > 1 ) {
                #    my $nearest_snp = find_nearest_snp( $all_snps, $snp_pos );
                #    $nearest_snp_id = $nearest_snp->[0];
                #    ( $r_s, $dprime_s )
                #        = calc_ld( $snp_occured, $nearest_snp->[2] );
                #}
                #$update_sth->execute( $r, $dprime, $nearest_snp_id, $r_s,
                #    $dprime_s, $snp_id );

                # average LD with near snps
                my ( $near_snp_number, $r2_s, $dprime_abs_s );
                if ( scalar @{$all_snps} > 1 ) {
                    my $near_snps = find_near_snps( $all_snps, $snp_pos );
                    $near_snp_number = scalar @{$near_snps};

                    if ( $near_snp_number > 1 ) {
                        my ( @r2_s, @dprime_abs_s );
                        my $combinat = Math::Combinatorics->new(
                            count => 2,
                            data  => $near_snps,
                        );
                        while ( my @combo = $combinat->next_combination ) {
                            my ( $pair_r_s, $pair_dprime_s )
                                = calc_ld( $combo[0]->[2], $combo[1]->[2] );
                            push @r2_s,         $pair_r_s**2;
                            push @dprime_abs_s, abs($pair_dprime_s);
                        }
                        $r2_s         = average(@r2_s);
                        $dprime_abs_s = average(@dprime_abs_s);
                    }
                }
                $update_sth->execute( $r, $dprime, $near_snp_number, $r2_s,
                    $dprime_abs_s, $snp_id );

                #my ( $r_i, $dprime_i, $r_ni, $dprime_ni, );

                #
                #if ( $group_i->size > 1 and $group_i->size < $all_freq - 1 ) {
                #
                #    # LD inside indel group
                #    ( $r_i, $dprime_i ) = calc_ld(
                #        $group_i->substr_span($indel_occured),
                #        $group_i->substr_span($snp_occured)
                #    );
                #
                #    # LD inside nonindel group
                #    ( $r_ni, $dprime_ni ) = calc_ld(
                #        $group_ni->substr_span($indel_occured),
                #        $group_ni->substr_span($snp_occured)
                #    );
                #}
                #$update_sth->execute( $r, $dprime, $r_i, $dprime_i, $r_ni,
                #    $dprime_ni, $snp_id );
            }
        }
    }
    $update_sth->finish;
    $snp_sth->finish;
    $indel_sth->finish;
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB::Multi->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}
exit;

sub calc_ld {
    my $strA = shift;
    my $strB = shift;

    if ( length $strA != length $strB ) {
        warn "length not equal for $strA and $strA\n";
        return;
    }

    for ( $strA, $strB ) {
        if (/[^ox]/) {
            warn "$_ contains illegal chars\n";
            return;
        }
    }

    my $size = length $strA;

    my $A_count = $strA =~ tr/o/o/;
    my $fA      = $A_count / $size;
    my $fa      = 1 - $fA;

    my $B_count = $strB =~ tr/o/o/;
    my $fB      = $B_count / $size;
    my $fb      = 1 - $fB;

    if ( any { $_ == 0 } ( $fA, $fa, $fB, $fb ) ) {
        return ( undef, undef );
    }

    # o in strA and o in strB as fAB
    my ( $AB_count, $fAB ) = ( 0, 0 );
    for my $i ( 1 .. $size ) {
        my $ichar = substr $strA, $i - 1, 1;
        my $schar = substr $strB, $i - 1, 1;
        if ( $ichar eq 'o' and $schar eq 'o' ) {
            $AB_count++;
        }
    }
    $fAB = $AB_count / $size;

    my $DAB = $fAB - $fA * $fB;

    my ( $r, $dprime );
    $r = $DAB / sqrt( $fA * $fa * $fB * $fb );

    if ( $DAB < 0 ) {
        $dprime = $DAB / min( $fA * $fB, $fa * $fb );
    }
    else {
        $dprime = $DAB / min( $fA * $fb, $fa * $fB );
    }

    return ( $r, $dprime );
}

sub find_nearest_snp {
    my $all_ary = shift;
    my $pos     = shift;

    my @sorted = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        map { [ $_, abs( $_->[3] - $pos ) ] }
        grep { $_->[3] != $pos } @{$all_ary};

    return $sorted[0];
}

sub find_near_snps {
    my $all_ary = shift;
    my $pos     = shift;

    my @sorted = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        grep { $_->[1] <= 50 }
        map { [ $_, abs( $_->[3] - $pos ) ] } @{$all_ary};

    return \@sorted;
}

__END__

=head1 NAME

    update_snp_ld.pl -  Add snp LD to the nearest indel
    
=head1 SYNOPSIS

    update_snp_ld.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password

=cut

