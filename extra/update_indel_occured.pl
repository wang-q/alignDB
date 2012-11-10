#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

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

# occured parameters
my $aim_db = $Config->{occured}{aim_db};
my $ref_db = $Config->{occured}{ref_db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'username=s' => \$username,
    'password=s' => \$password,
    'aim_db=s'   => \$aim_db,
    'ref_db=s'   => \$ref_db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $aim_obj = AlignDB->new(
    mysql  => "$aim_db:$server",
    user   => $username,
    passwd => $password,
);

my $ref_obj = AlignDB->new(
    mysql  => "$ref_db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $aim_dbh = $aim_obj->dbh;
my $ref_dbh = $ref_obj->dbh;

#----------------------------#
# Add columns
#----------------------------#
{
    $aim_obj->create_column( "indel", "indel_other_occured", "CHAR(1)" );
    print "Table indel altered\n";

    $aim_obj->create_column( "snp", "snp_other_occured", "CHAR(1)" );
    print "Table snp altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{

    # alignments
    my $aim_align_query = q{
        SELECT a.align_id, chr_id, chr_start, chr_end
        FROM align a, target t, sequence s
        WHERE a.align_id = t.align_id
        AND t.seq_id = s.seq_id
    };
    my $aim_align_sth = $aim_dbh->prepare($aim_align_query);

    my $ref_align_query = q{
        SELECT a.align_id, chr_id, chr_start, chr_end
        FROM align a, target t, sequence s
        WHERE a.align_id = t.align_id
        AND t.seq_id = s.seq_id
        AND a.align_id = ?
    };
    my $ref_align_sth = $ref_dbh->prepare($ref_align_query);

    # select all indels in this alignment
    my $indel_query = q{
        SELECT indel_id, indel_start, indel_end, indel_occured, indel_seq
        FROM indel
        WHERE align_id = ?
    };
    my $aim_indel_sth = $aim_dbh->prepare($indel_query);
    my $ref_indel_sth = $ref_dbh->prepare($indel_query);

    # update indel table in the new feature column
    my $indel_update = q{
        UPDATE indel
        SET indel.indel_other_occured = ?
        WHERE 1 = 1
        AND indel.indel_id = ?
    };
    my $indel_update_sth = $aim_dbh->prepare($indel_update);

    # select all snps in this alignment
    my $snp_query = q{
        SELECT snp_id, snp_pos, snp_occured, target_base, query_base
        FROM snp
        WHERE align_id = ?
    };
    my $aim_snp_sth = $aim_dbh->prepare($snp_query);
    my $ref_snp_sth = $ref_dbh->prepare($snp_query);

    # update snp table in the new feature column
    my $snp_update = q{
        UPDATE snp
        SET snp.snp_other_occured = ?
        WHERE 1 = 1
        AND snp.snp_id = ?
    };
    my $snp_update_sth = $aim_dbh->prepare($snp_update);

    $aim_align_sth->execute;

    # foreach alignment
ALN: while ( my $row_hashref = $aim_align_sth->fetchrow_hashref ) {
        my $aim_hashref = $row_hashref;
        my $align_id    = $aim_hashref->{align_id};
        print "Processing align_id $align_id\n";

        $ref_align_sth->execute($align_id);
        my $ref_hashref = $ref_align_sth->fetchrow_hashref;

        # if any align info in aim and ref is not identical,
        #   abandon this alignment
        for ( sort keys %$aim_hashref ) {
            if ( $aim_hashref->{$_} ne $ref_hashref->{$_} ) {
                print " " x 4, "$_ does not match; ";
                print "align info error, jump to next\n";
                next ALN;
            }
            else {
                print " " x 4, "$_ matched\n";
            }
        }

        #----------------------------#
        # indel
        #----------------------------#
        # store indel info into a hash.
        # because indel_start and indel_end are unique, use them as keys
        # To be more accurate, indel_seq should be used as a constraint
        $aim_indel_sth->execute($align_id);
        my %aim_indel_info;
        while ( my @row = $aim_indel_sth->fetchrow_array ) {
            my ($indel_id,      $indel_start, $indel_end,
                $indel_occured, $indel_seq
            ) = @row;
            $aim_indel_info{"$indel_start-$indel_end-$indel_seq"} = {
                indel_id      => $indel_id,
                indel_occured => $indel_occured,
            };
        }

        $ref_indel_sth->execute($align_id);
        my %ref_indel_info;
        while ( my @row = $ref_indel_sth->fetchrow_array ) {
            my ($indel_id,      $indel_start, $indel_end,
                $indel_occured, $indel_seq
            ) = @row;
            $ref_indel_info{"$indel_start-$indel_end-$indel_seq"} = {
                indel_id      => $indel_id,
                indel_occured => $indel_occured,
            };
        }

        # find corresponding indels
        foreach my $indel_key ( keys %aim_indel_info ) {
            if ( exists $ref_indel_info{$indel_key} ) {
                my $aim_indel_id = $aim_indel_info{$indel_key}->{indel_id};
                my $indel_other_occured
                    = $ref_indel_info{$indel_key}->{indel_occured};
                $indel_update_sth->execute( $indel_other_occured,
                    $aim_indel_id );
            }
        }

        #----------------------------#
        # snp
        #----------------------------#
        # store snp info into a hash.
        # use snp_pos, target_base and query_base as keys
        $aim_snp_sth->execute($align_id);
        my %aim_snp_info;
        while ( my @row = $aim_snp_sth->fetchrow_array ) {
            my ( $snp_id, $snp_pos, $snp_occured, $target_base, $query_base )
                = @row;
            $aim_snp_info{"$snp_pos-$target_base-$query_base"} = {
                snp_id      => $snp_id,
                snp_occured => $snp_occured,
            };
        }

        $ref_snp_sth->execute($align_id);
        my %ref_snp_info;
        while ( my @row = $ref_snp_sth->fetchrow_array ) {
            my ( $snp_id, $snp_pos, $snp_occured, $target_base, $query_base )
                = @row;
            $ref_snp_info{"$snp_pos-$target_base-$query_base"} = {
                snp_id      => $snp_id,
                snp_occured => $snp_occured,
            };
        }

        # find corresponding snps
        foreach my $snp_key ( keys %aim_snp_info ) {
            if ( exists $ref_snp_info{$snp_key} ) {
                my $aim_snp_id        = $aim_snp_info{$snp_key}->{snp_id};
                my $snp_other_occured = $ref_snp_info{$snp_key}->{snp_occured};
                $snp_update_sth->execute( $snp_other_occured, $aim_snp_id );
            }
        }
    }

    $snp_update_sth->finish;
    $ref_snp_sth->finish;
    $aim_snp_sth->finish;

    $indel_update_sth->finish;
    $ref_indel_sth->finish;
    $aim_indel_sth->finish;

    $ref_align_sth->finish;
    $aim_align_sth->finish;
}

__END__

=head1 NAME

    update_indel_occured.pl - Find indel occured in different lineage
                                when compared with different outgroup

=head1 SYNOPSIS

    update_indel_slippage.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --username        username
       --password        password
       --aim_db          aim database name
       --ref_db          ref database name
       

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut

