#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all);
use Scalar::Util qw(looks_like_number);

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

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

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'        => \( my $server         = $Config->{database}{server} ),
    'port|P=i'          => \( my $port           = $Config->{database}{port} ),
    'db|d=s'            => \( my $db             = $Config->{database}{db} ),
    'username|u=s'      => \( my $username       = $Config->{database}{username} ),
    'password|p=s'      => \( my $password       = $Config->{database}{password} ),
    'insert_segment=s'  => \( my $insert_segment = 1 ),
    'near|near_range=i' => \( my $near_range     = 100 ),
    'parallel=i'        => \( my $parallel       = $Config->{generate}{parallel} ),
    'batch=i'           => \( my $batch_number   = $Config->{generate}{batch} ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update LD of $db...");

#----------------------------#
# Add columns
#----------------------------#
my $seq_count;
my @jobs;
{
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    $seq_count = $obj->get_seq_count;

    # r and D' with nearest indel
    $obj->create_column( "snp", "snp_r",      "DOUBLE" );
    $obj->create_column( "snp", "snp_dprime", "DOUBLE" );

    # r2 and |D'| with near snps
    # include self
    $obj->create_column( "snp", "snp_near_snp_number", "int" );
    $obj->create_column( "snp", "snp_r2_s",            "DOUBLE" );
    $obj->create_column( "snp", "snp_dprime_abs_s",    "DOUBLE" );

    # indel group
    $obj->create_column( "snp", "snp_r2_i",         "DOUBLE" );
    $obj->create_column( "snp", "snp_dprime_abs_i", "DOUBLE" );

    # nonindel group
    $obj->create_column( "snp", "snp_r2_ni",         "DOUBLE" );
    $obj->create_column( "snp", "snp_dprime_abs_ni", "DOUBLE" );

    if ($insert_segment) {
        $obj->create_column( "segment", "segment_r2_s",         "DOUBLE" );
        $obj->create_column( "segment", "segment_dprime_abs_s", "DOUBLE" );
    }

    #  r and D' with nearest snp
    ##$obj->create_column( "snp", "snp_nearest_snp_id", "int" );
    #$obj->create_column( "snp", "snp_r_s",             "DOUBLE" );
    #$obj->create_column( "snp", "snp_dprime_s",        "DOUBLE" );

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

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my DBI $dbh = $obj->dbh;

    # select all indels in this alignment
    my DBI $indel_sth = $dbh->prepare(
        q{
        SELECT i.indel_id, i.indel_freq, i.indel_occured
        FROM indel i
        WHERE 1=1
        AND i.indel_occured != 'unknown'
        AND i.align_id = ?
        }
    );

    # select all isw in this indel
    my DBI $snp_sth = $dbh->prepare(
        q{
        SELECT s.snp_id, s.snp_freq, s.snp_occured, s.snp_pos
        FROM isw w, snp s
        WHERE 1 = 1
        AND w.isw_id = s.isw_id
        AND s.snp_occured != 'unknown'
        AND w.isw_indel_id = ?
        }
    );

    my DBI $update_indel_snp_sth = $dbh->prepare(
        q{
        UPDATE
            snp s
        SET
            s.snp_r         = ?,
            s.snp_dprime    = ?
        WHERE s.snp_id = ?
        }
    );

    my DBI $update_snps_sth = $dbh->prepare(
        q{
        UPDATE
            snp s
        SET
            s.snp_near_snp_number = ?,
            s.snp_r2_s            = ?,
            s.snp_dprime_abs_s    = ?
        WHERE s.snp_id = ?
        }
    );

    my DBI $update_sub_snps_sth = $dbh->prepare(
        q{
        UPDATE
            snp s
        SET
            s.snp_r2_i          = ?,
            s.snp_dprime_abs_i  = ?,
            s.snp_r2_ni         = ?,
            s.snp_dprime_abs_ni = ?
        WHERE s.snp_id = ?
        }
    );

    my DBI $update_segment_snps_sth = $dbh->prepare(
        q{
        UPDATE
            segment s
        SET
            s.segment_r2_s         = ?,
            s.segment_dprime_abs_s = ?
        WHERE s.segment_id = ?
        }
    );

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

        # average LD with near snps
        if ( scalar @{$all_snps} > 1 ) {
            for my $cur_snp ( @{$all_snps} ) {
                my $snp_id          = $cur_snp->[0];
                my $snp_pos         = $cur_snp->[3];
                my $near_snps       = find_near_snps( $all_snps, $snp_pos, $near_range );
                my $near_snp_number = scalar @{$near_snps};
                my ( $r2_s, $dprime_abs_s ) = combo_ld($near_snps);

                $update_snps_sth->execute( $near_snp_number, $r2_s, $dprime_abs_s, $snp_id );
            }
        }

        if ($insert_segment) {

            # snps in segment
            my $all_segments = $dbh->selectall_arrayref(
                q{
                SELECT s.segment_id, w.window_start, w.window_end
                FROM segment s
                    INNER JOIN window w ON s.window_id = w.window_id
                WHERE w.window_length = 500
                AND w.align_id = ?
                }, {}, $align_id
            );

            for my $cur_segment ( @{$all_segments} ) {
                my $segment_id = $cur_segment->[0];
                my $segment_set
                    = AlignDB::IntSpan->new( $cur_segment->[1] . '-' . $cur_segment->[2] );
                my $set_snps = find_set_snps( $all_snps, $segment_set );
                my ( $r2_seg, $dprime_abs_seg ) = combo_ld($set_snps);

                $update_segment_snps_sth->execute( $r2_seg, $dprime_abs_seg, $segment_id );
            }
        }

        # LD with nearest indel
        $indel_sth->execute($align_id);
        while ( my @row = $indel_sth->fetchrow_array ) {
            my ( $indel_id, $indel_freq, $indel_occured ) = @row;

            my $group_i  = AlignDB::IntSpan->new;
            my $group_ni = AlignDB::IntSpan->new;
            for my $i ( 0 .. $seq_count - 1 ) {
                my $ichar = substr $indel_occured, $i, 1;
                if ( $ichar eq 'o' ) {
                    $group_i->add( $i + 1 );
                }
                elsif ( $ichar eq 'x' ) {
                    $group_ni->add( $i + 1 );
                }
                else {
                    die "indel occured error [$i] [$ichar]\n";
                }
            }

            $snp_sth->execute($indel_id);
            while ( my @row = $snp_sth->fetchrow_array ) {
                my ( $snp_id, $snp_freq, $snp_occured, $snp_pos ) = @row;

                # LD with nearest indel
                my ( $r, $dprime );
                ( $r, $dprime ) = calc_ld( $indel_occured, $snp_occured );
                $update_indel_snp_sth->execute( $r, $dprime, $snp_id );

                my ( $r2_i, $dprime_abs_i, $r2_ni, $dprime_abs_ni, );
                my $near_snps = find_near_snps( $all_snps, $snp_pos, $near_range );
                my $near_snp_number = scalar @{$near_snps};
                if ( $near_snp_number > 1 ) {
                    my @near_snps_i
                        = map { [ $_->[0], $_->[1], $group_i->substr_span( $_->[2] ), $_->[3], ] }
                        @{$near_snps};
                    my @near_snps_ni
                        = map { [ $_->[0], $_->[1], $group_ni->substr_span( $_->[2] ), $_->[3], ] }
                        @{$near_snps};

                    ( $r2_i,  $dprime_abs_i )  = combo_ld( \@near_snps_i );
                    ( $r2_ni, $dprime_abs_ni ) = combo_ld( \@near_snps_ni );
                }
                $update_sub_snps_sth->execute( $r2_i, $dprime_abs_i, $r2_ni,
                    $dprime_abs_ni, $snp_id );
            }
        }
    }
    $update_segment_snps_sth->finish;
    $update_sub_snps_sth->finish;
    $update_snps_sth->finish;
    $update_indel_snp_sth->finish;
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
    AlignDB->new(
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
    my $range   = shift;

    my @sorted = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        grep { $_->[1] <= $range }
        map { [ $_, abs( $_->[3] - $pos ) ] } @{$all_ary};

    return \@sorted;
}

sub find_set_snps {
    my $all_ary = shift;
    my $set     = shift;    # AlignDB::IntSpan object

    my @sorted = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        grep { $set->contain( $_->[1] ) }
        map { [ $_, $_->[3] ] } @{$all_ary};

    return \@sorted;
}

sub combo_ld {
    my $snps = shift;

    my ( $r2, $dprime_abs );
    if ( scalar @{$snps} > 1 ) {
        my ( @r2, @dprime_abs );
        my $combinat = Math::Combinatorics->new(
            count => 2,
            data  => $snps,
        );
        while ( my @combo = $combinat->next_combination ) {
            my ( $pair_r, $pair_dprime )
                = calc_ld( $combo[0]->[2], $combo[1]->[2] );
            next unless defined $pair_r;
            push @r2,         $pair_r**2;
            push @dprime_abs, abs($pair_dprime);
        }
        if ( scalar @r2 > 0 ) {
            $r2         = average(@r2);
            $dprime_abs = average(@dprime_abs);
        }
    }

    return ( $r2, $dprime_abs );
}

__END__
