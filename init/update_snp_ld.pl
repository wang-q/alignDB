#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use List::Util;
use List::MoreUtils::PP;
use MCE;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use App::Fasops::Common;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use AlignDB::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
LD values of SNPs to the nearest indel and to nearby SNPs

    perl init/update_snp_ld.pl -d S288cvsRM11_1a --parallel 2

Usage: perl %c [options]
EOF

(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    $description,
    [ 'help|h', 'display this message' ],
    [],
    ['Database init values'],
    [ 'server|s=s',   'MySQL IP/Domain', { default => $conf->{database}{server} }, ],
    [ 'port=i',       'MySQL port',      { default => $conf->{database}{port} }, ],
    [ 'username|u=s', 'username',        { default => $conf->{database}{username} }, ],
    [ 'password|p=s', 'password',        { default => $conf->{database}{password} }, ],
    [ 'db|d=s',       'database name',   { default => $conf->{database}{db} }, ],
    [],
    [ 'insert_segment=s', 'LD in segments',       { default => 1 }, ],
    [ 'insert_ofg=s',     'LD in ofgs',           { default => 1 }, ],
    [ 'near=i',           'near range',           { default => 100 }, ],
    [ 'parallel=i',       'run in parallel mode', { default => $conf->{generate}{parallel} }, ],
    [ 'batch=i', '#alignments in one process', { default => $conf->{generate}{batch} }, ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# record config
$stopwatch->record_conf($opt);

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update LD of [$opt->{db}]...");

my $seq_count;
my @jobs;
{
    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    $seq_count = $alignDB->get_seq_count;

    print "Add columns to tables\n";

    # r and D' with nearest indel
    $alignDB->create_column( "snp", "snp_r",      "DOUBLE" );
    $alignDB->create_column( "snp", "snp_dprime", "DOUBLE" );

    # r2 and |D'| with near snps
    # include self
    $alignDB->create_column( "snp", "snp_near_snp_number", "int" );
    $alignDB->create_column( "snp", "snp_r2_s",            "DOUBLE" );
    $alignDB->create_column( "snp", "snp_dprime_abs_s",    "DOUBLE" );

    # indel group
    $alignDB->create_column( "snp", "snp_r2_i",         "DOUBLE" );
    $alignDB->create_column( "snp", "snp_dprime_abs_i", "DOUBLE" );

    # nonindel group
    $alignDB->create_column( "snp", "snp_r2_ni",         "DOUBLE" );
    $alignDB->create_column( "snp", "snp_dprime_abs_ni", "DOUBLE" );

    if ( $opt->{insert_segment} ) {
        $alignDB->create_column( "segment", "segment_r2_s",         "DOUBLE" );
        $alignDB->create_column( "segment", "segment_dprime_abs_s", "DOUBLE" );
    }

    if ( $opt->{insert_ofg} ) {
        $alignDB->create_column( "ofgsw", "ofgsw_r2_s",         "DOUBLE" );
        $alignDB->create_column( "ofgsw", "ofgsw_dprime_abs_s", "DOUBLE" );
    }

    #  r and D' with nearest snp
    ##$obj->create_column( "snp", "snp_nearest_snp_id", "int" );
    #$obj->create_column( "snp", "snp_r_s",             "DOUBLE" );
    #$obj->create_column( "snp", "snp_dprime_s",        "DOUBLE" );

    print "Table snp altered\n";

    @jobs = @{ $alignDB->get_align_ids };
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my @align_ids = @{$chunk_ref};
    my $wid       = MCE->wid;

    $stopwatch->block_message("Process task [$chunk_id] by worker #$wid");

    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    # Database handler
    my DBI $dbh = $alignDB->dbh;

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

    # select all SNPs belongs to this indel
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

    my DBI $update_ofgsw_snps_sth = $dbh->prepare(
        q{
        UPDATE
            ofgsw sw
        SET
            sw.ofgsw_r2_s         = ?,
            sw.ofgsw_dprime_abs_s = ?
        WHERE sw.ofgsw_id = ?
        }
    );

    # for each align
    for my $align_id (@align_ids) {
        $alignDB->process_message($align_id);

        # all SNPs in this alignment
        my $all_snps = $dbh->selectall_arrayref(
            q{
            SELECT s.snp_id, s.snp_freq, s.snp_occured, s.snp_pos
            FROM snp s
            WHERE 1 = 1
            AND s.snp_occured != 'unknown'
            AND s.align_id = ?
            }, {}, $align_id
        );

        # average LD with nearby snps
        if ( scalar @{$all_snps} > 1 ) {
            for my $cur_snp ( @{$all_snps} ) {
                my $snp_id          = $cur_snp->[0];
                my $snp_pos         = $cur_snp->[3];
                my $near_snps       = find_near_snps( $all_snps, $snp_pos, $opt->{near} );
                my $near_snp_number = scalar @{$near_snps};
                my ( $r2_s, $dprime_abs_s ) = combo_ld($near_snps);

                $update_snps_sth->execute( $near_snp_number, $r2_s, $dprime_abs_s, $snp_id );
            }
        }

        if ( $opt->{insert_segment} ) {

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

        if ( $opt->{insert_ofg} ) {

            # snps in ofgsw
            my $all_ofgsws = $dbh->selectall_arrayref(
                q{
                SELECT s.ofgsw_id, w.window_start, w.window_end
                FROM ofgsw s
                    INNER JOIN window w ON s.window_id = w.window_id
                WHERE w.align_id = ?
                }, {}, $align_id
            );

            for my $cur_ofgsw ( @{$all_ofgsws} ) {
                my $ofgsw_id = $cur_ofgsw->[0];
                my $ofgsw_set
                    = AlignDB::IntSpan->new( $cur_ofgsw->[1] . '-' . $cur_ofgsw->[2] );
                my $set_snps = find_set_snps( $all_snps, $ofgsw_set );
                my ( $r2_ofgsw, $dprime_abs_ofgsw ) = combo_ld($set_snps);

                $update_ofgsw_snps_sth->execute( $r2_ofgsw, $dprime_abs_ofgsw, $ofgsw_id );
            }
        }

        # LD with nearest indel
        $indel_sth->execute($align_id);
        while ( my @row = $indel_sth->fetchrow_array ) {
            my ( $indel_id, undef, $indel_occured ) = @row;

            my $group_i  = AlignDB::IntSpan->new;
            my $group_ni = AlignDB::IntSpan->new;
            for my $i ( 0 .. $seq_count - 1 ) {
                my $ichar = substr $indel_occured, $i, 1;
                if ( $ichar eq '1' ) {
                    $group_i->add( $i + 1 );
                }
                elsif ( $ichar eq '0' ) {
                    $group_ni->add( $i + 1 );
                }
                else {
                    die "indel occured error [$i] [$ichar]\n";
                }
            }

            $snp_sth->execute($indel_id);
            while ( my @row2 = $snp_sth->fetchrow_array ) {
                my ( $snp_id, undef, $snp_occured, $snp_pos ) = @row2;

                # LD with nearest indel
                my ( $r, $dprime );
                ( $r, $dprime ) = App::Fasops::Common::calc_ld( $indel_occured, $snp_occured );
                $update_indel_snp_sth->execute( $r, $dprime, $snp_id );

                my ( $r2_i, $dprime_abs_i, $r2_ni, $dprime_abs_ni, );
                my $near_snps = find_near_snps( $all_snps, $snp_pos, $opt->{near} );
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

    # finish all sths
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
my $mce = MCE->new( max_workers => $opt->{parallel}, chunk_size => $opt->{batch}, );
$mce->forchunk( \@jobs, $worker, );

$stopwatch->end_message;

# store program's meta info to database
AlignDB::Common->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
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
    my AlignDB::IntSpan $set = shift;

    my @sorted = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        grep { $set->contains( $_->[1] ) }
        map { [ $_, $_->[3] ] } @{$all_ary};

    return \@sorted;
}

sub combo_ld {
    my $snps = shift;

    my $count = scalar @{$snps};

    my ( $r2, $dprime_abs );
    if ( $count >= 2 ) {
        my ( @r2, @dprime_abs );
        for my $i ( 0 .. $count - 1 ) {
            for my $j ( $i + 1 .. $count - 1 ) {
                my ( $pair_r, $pair_dprime )
                    = App::Fasops::Common::calc_ld( $snps->[$i][2], $snps->[$j][2] );
                next unless defined $pair_r;
                push @r2,         $pair_r**2;
                push @dprime_abs, abs($pair_dprime);
            }
        }
        if ( scalar @r2 > 0 ) {
            $r2         = App::Fasops::Common::mean(@r2);
            $dprime_abs = App::Fasops::Common::mean(@dprime_abs);
        }
    }

    return ( $r2, $dprime_abs );
}

__END__
