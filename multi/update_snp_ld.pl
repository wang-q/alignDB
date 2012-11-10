#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Scalar::Util qw(looks_like_number);

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

    # r
    $obj->create_column( "snp", "snp_r", "DOUBLE" );

    # D'
    $obj->create_column( "snp", "snp_dprime", "DOUBLE" );

    print "Table snp altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{

    my $obj = AlignDB::Multi->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    my @align_ids = @{ $obj->get_align_ids };

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
        SELECT s.snp_id, s.snp_freq, s.snp_occured
        FROM isw w, snp s
        WHERE 1 = 1
        AND w.isw_id = s.isw_id
        AND s.snp_occured != 'unknown'
        AND w.isw_indel_id = ?
    };
    my $snp_sth = $dbh->prepare($snp_query);

    # update snp table in the new feature columns
    my $update_query = q{
        UPDATE snp s
        SET s.snp_r         = ?,
            s.snp_dprime    = ?
        WHERE s.snp_id = ?
    };
    my $update_sth = $dbh->prepare($update_query);

    # for each align
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);

        $indel_sth->execute($align_id);
        while ( my @row = $indel_sth->fetchrow_array ) {
            my ( $indel_id, $indel_freq, $indel_occured ) = @row;

            $snp_sth->execute($indel_id);
            while ( my @row = $snp_sth->fetchrow_array ) {
                my ( $snp_id, $snp_freq, $snp_occured ) = @row;

                my ( $r, $dprime ) = calc_ld( $indel_occured, $snp_occured );
                $update_sth->execute( $r, $dprime, $snp_id );
            }
        }
    }
    $update_sth->finish;
    $snp_sth->finish;
    $indel_sth->finish;
}

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

    # indel and snp as fAB
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

