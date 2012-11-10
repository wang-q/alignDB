#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

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

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    my @names = $aim_obj->get_names;
    my ( undef, undef, $ref_ref_name ) = $ref_obj->get_names;
    push @names, $ref_ref_name;

    #----------------------------#
    # SQL queries
    #----------------------------#
    my @aim_align_ids = @{ $aim_obj->get_align_ids };

    # sequences
    my $aim_seq_query = q{
        SELECT t.target_seq, q.query_seq, r.ref_seq
        FROM align a, target t, query q, reference r
        WHERE a.align_id = ?
        AND a.align_id = t.align_id
        AND a.align_id = q.align_id
        AND a.align_id = r.align_id
    };
    my $aim_seq_sth = $aim_dbh->prepare($aim_seq_query);

    my $ref_seq_query = q{
        SELECT r.ref_seq
        FROM reference r
        where r.align_id = ?
    };
    my $ref_seq_sth = $ref_dbh->prepare($ref_seq_query);

    # find all clearly occured indels
    my $occured_indel = q{
        SELECT i.indel_id, indel_start, indel_end
        FROM   indel i
        WHERE  1 = 1
        AND i.align_id = ?
        AND i.indel_occured = ?
        And i.indel_other_occured = ?
    };
    my $occured_indel_sth = $aim_dbh->prepare($occured_indel);

    # find clearly occured indels' vicinity regions
    my $indel_vicinity = q{
        SELECT isw.isw_id, isw.isw_start, isw.isw_end
        FROM   isw,
               (SELECT isw_id
                FROM   isw
                WHERE isw_type = 'R'
                AND isw.indel_id = ?
                UNION 
                SELECT isw_id
                FROM   isw
                WHERE  isw_type = 'L'
                AND isw.prev_indel_id = ?) i
        WHERE  isw_distance >= 0
        AND isw_distance <= 5
        AND isw.isw_id = i.isw_id
    };
    my $indel_vicinity_sth = $aim_dbh->prepare($indel_vicinity);

    my @occured_groups
        = ( [ 'T', 'T' ], [ 'Q', 'Q' ], [ 'T', 'Q' ], [ 'Q', 'T' ], );

    #----------------------------#
    # write indel vicinity sequences to file
    #----------------------------#
    # vicinity region is composed by isws
    for my $g (@occured_groups) {
        my $group_name = $g->[0] . $g->[1];
        my @group_segments;
        print "For group $group_name\n";

        # foreach alignment
        for my $align_id (@aim_align_ids) {
            print " " x 4, "Processing align_id $align_id\n";

            $aim_seq_sth->execute($align_id);
            my ( $aim_target_seq, $aim_query_seq, $aim_ref_seq )
                = $aim_seq_sth->fetchrow_array;

            $ref_seq_sth->execute($align_id);
            my ($ref_ref_seq) = $ref_seq_sth->fetchrow_array;

            $occured_indel_sth->execute( $align_id, $g->[0], $g->[1] );
            while ( my @row = $occured_indel_sth->fetchrow_array ) {
                my ( $indel_id, $indel_start, $indel_end ) = @row;

                # indel vicinity region including indel itself
                my $vicinity_region = AlignDB::IntSpan->new;
                $vicinity_region->add("$indel_start-$indel_end");

                # merge all isws of this indel to vicinity region
                $indel_vicinity_sth->execute( $indel_id, $indel_id );
                while ( my @row2 = $indel_vicinity_sth->fetchrow_array ) {
                    my ( $isw_id, $isw_start, $isw_end ) = @row2;
                    $vicinity_region->add("$isw_start-$isw_end");
                }
                next if $vicinity_region->cardinality < 100;

                if ( $vicinity_region->spans > 1 ) {
                    warn "Vicinity of indel [$indel_id] is interrupted\n";
                }

                # merge this vicinity region to occured groups
                my @segments;
                for (
                    $aim_target_seq, $aim_query_seq,
                    $aim_ref_seq,    $ref_ref_seq
                    )
                {
                    my $seg = $vicinity_region->substr_span($_);
                    push @segments, $seg;
                }

                #@segments = @{ &clustal_align( \@segments ) };

                for ( 0 .. 3 ) {
                    $group_segments[$_] .= $segments[$_];
                }
            }
        }

        $ref_seq_sth->finish;
        $aim_seq_sth->finish;

        my $outfile = "$aim_db-$ref_ref_name" . "-$group_name.fas";
        open my $outfh, ">", $outfile;
        for ( 0 .. 3 ) {
            print {$outfh} ">", $names[$_], "\n";
            print {$outfh} $group_segments[$_], "\n";
        }
        close $outfh;
    }

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

