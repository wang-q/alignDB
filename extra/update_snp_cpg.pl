#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Update CpG info of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# get reference names via AlignDB methods
my ( undef, undef, $ref_name ) = $obj->get_names;
if ( !$ref_name or $ref_name eq 'NULL' ) {
    die "$db is not a three-way alignDB\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    my @align_ids = @{ $obj->get_align_ids };

    # select all snps in this alignment
    my $snp_query = q{
        SELECT  snp_id, snp_pos, target_base, query_base,
                ref_base, snp_occured,
                IF(snp_occured = "T", target_base, query_base) `mutated_base`
        FROM snp
        WHERE align_id = ?
        AND ref_base IN ("C", "G")
        AND snp_occured IN ("T", "Q")
    };
    my $snp_sth = $dbh->prepare($snp_query);

    # update snp table in the new feature column
    my $snp_update = q{
        UPDATE snp
        SET snp_cpg = ?
        WHERE snp_id = ?
    };
    my $snp_update_sth = $dbh->prepare($snp_update);

    for my $align_id (@align_ids) {
        print "Processing align_id $align_id\n";

        my ( $target_seq, $query_seq ) = @{ $obj->get_seqs($align_id) };

        $snp_sth->execute($align_id);
        while ( my @row = $snp_sth->fetchrow_array ) {
            my ($snp_id,   $snp_pos,     $target_base, $query_base,
                $ref_base, $snp_occured, $mutated_base
            ) = @row;

            my ( %left_base, %right_base );
            $left_base{T}  = substr( $target_seq, $snp_pos - 2, 1 );
            $right_base{T} = substr( $target_seq, $snp_pos,     1 );
            $left_base{Q}  = substr( $query_seq,  $snp_pos - 2, 1 );
            $right_base{Q} = substr( $query_seq,  $snp_pos,     1 );

            my $snp_cpg = 0;

            # CpG to TpG, C to T transition
            # On the reverse strand, is CpG to CpA
            if ( $ref_base eq "C" ) {    # original base is C
                if ( $mutated_base eq "T" ) {
                    if ( $right_base{T} eq "G" and $right_base{Q} eq "G" ) {
                        $snp_cpg = 1;
                    }
                }
            }
            elsif ( $ref_base eq "G" ) {    # original base is G
                if ( $mutated_base eq "A" ) {
                    if ( $left_base{T} eq "C" and $left_base{Q} eq "C" ) {
                        $snp_cpg = 1;
                    }
                }
            }

            next unless $snp_cpg;

            $snp_update_sth->execute( $snp_cpg, $snp_id );
        }

    }

    # update NULL value of snp_cpg to 0
    my $snp_null = q{
        UPDATE snp
        SET snp_cpg = 0
        WHERE snp_cpg IS NULL
    };
    my $snp_null_sth = $dbh->prepare($snp_null);
    $snp_null_sth->execute;

    $snp_null_sth->finish;
    $snp_update_sth->finish;
    $snp_sth->finish;

}

{
    print "Processing isw_cpg\n";

    # select all snps in this alignment
    my $isw_query = q{
        SELECT  i.isw_id id,
                COUNT(*) /i.isw_length * 1.0 cpg
        FROM isw i, snp s
        WHERE i.isw_id = s.isw_id
        AND s.snp_cpg = 1
        GROUP BY i.isw_id
    };
    my $isw_sth = $dbh->prepare($isw_query);

    # update isw table in the new feature column
    my $isw_update = q{
        UPDATE isw
        SET isw_cpg_pi = ?
        WHERE isw_id = ?
    };
    my $isw_update_sth = $dbh->prepare($isw_update);

    # for isw
    $isw_sth->execute;
    while ( my @row = $isw_sth->fetchrow_array ) {
        my ( $isw_id, $cpg ) = @row;
        $isw_update_sth->execute( $cpg, $isw_id );
    }

    # update NULL value of isw_cpg to 0
    my $isw_null = q{
        UPDATE isw
        SET isw_cpg_pi = 0
        WHERE isw_cpg_pi IS NULL
    };
    my $isw_null_sth = $dbh->prepare($isw_null);
    $isw_null_sth->execute;

}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    update_snp_cpg.pl - Add additional CpG info to alignDB
                        1 for CpG and 0 for non

=head1 SYNOPSIS

    update_snp_cpg.pl [options]
        Options:
            --help            brief help message
            --man             full documentation
            --server          MySQL server IP/Domain name
            --db              database name
            --username        username
            --password        password

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<update_snp_cpg.pl> will Add additional CpG info to alignDB,
1 for CpG and 0 for non.

=cut

