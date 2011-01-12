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

# reference database name
my $aim_db = $Config->{gene}{aim_db};
my $ref_db = $Config->{gene}{ref_db};

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
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Init $aim_db...");

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

# add a column gene_feature5 to gene
{
    $aim_obj->create_column( "gene", "gene_feature5", "DOUBLE" );
    print "Table gene altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{

    # gene divergence
    # count every full exons' difference and comparable base number
    # use ref_db
    my $gene_ref_query = q{
        SELECT  ge.gene_stable_id,
                SUM(e.diff),
                SUM(e.comp)
        FROM    (SELECT e.exon_stable_id, COUNT(*) count,
                        SUM(w.window_differences) diff,
                        SUM(w.window_comparables) comp
                FROM exon e, window w
                WHERE e.exon_is_full = 1
                AND e.window_id = w.window_id
                GROUP BY e.exon_stable_id) e,
                (SELECT g.gene_stable_id, e.exon_stable_id
                FROM gene g, exon e
                WHERE g.gene_id = e.gene_id
                GROUP BY g.gene_stable_id, e.exon_stable_id) ge
        WHERE ge.exon_stable_id = e.exon_stable_id
        GROUP BY ge.gene_stable_id
    };
    my $gene_ref_query_sth = $ref_dbh->prepare($gene_ref_query);

    print "Get gene divergences in ref_db: $ref_db\n";
    $gene_ref_query_sth->execute;
    my %gene_divergence;
    while ( my @row = $gene_ref_query_sth->fetchrow_array ) {
        my ( $gene_stable_id, $difference, $comparable ) = @row;
        next if $comparable <= 0;
        $gene_divergence{$gene_stable_id} = $difference / $comparable;
    }
    $gene_ref_query_sth->finish;

    #DumpFile("gene_D.yaml", \%gene_divergence);

    # update handler in db
    my $gene_update_query = q{
        UPDATE gene
        SET gene_feature5 = ?
        WHERE gene_stable_id = ?
    };
    my $gene_update_sth = $aim_dbh->prepare($gene_update_query);

    # update each gene
    print "Update gene_feature5 in aim_db: $aim_db\n";
    foreach my $gene_stable_id ( sort keys %gene_divergence ) {
        my $divergence = $gene_divergence{$gene_stable_id};
        $gene_update_sth->execute( $divergence, $gene_stable_id, );
    }
    $gene_update_sth->finish;
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    update_gene_D.pl - Add Divergence of genes to a refernce, gene_feature5 

=head1 SYNOPSIS

    update_gene_D.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --aim_db            aim database name
        --ref_db            ref database name

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

