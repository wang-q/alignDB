#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Statistics::Lite qw(median);

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

my $quan_file = "$FindBin::Bin/nbt.1551-S2.csv";

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
$stopwatch->start_message("Update $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# add a column gene_feature10 to gene, Quantification 
{
    $obj->create_column( "gene", "gene_feature10", "DOUBLE" );
    print "Table gene altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

{    # update 
    my @name_quan = @{ read_in_quan($quan_file) };

    my $gene_query = q{
        UPDATE  gene
        SET gene_feature10 = ?
        WHERE gene_stable_id = ?
    };
    my $gene_sth = $dbh->prepare($gene_query);

    foreach my $item (@name_quan) {
        $gene_sth->execute($item->[1], $item->[0]);
    }

    $gene_sth->finish;
}

$stopwatch->end_message;

exit;

sub read_in_quan {
    my $infile = shift;

    my @name_quan;
    open my $infh, '<', $infile;
    while (<$infh>) {
        next unless /^\d+/;
        my ($name, $avg) = ( split /,/ )[1,9];
        push @name_quan, [$name, $avg];
    }
    close $infh;

    return \@name_quan;
}

__END__

=head1 NAME

    update_segment.pl - Add additional slippage-like info to alignDB
                        1 for slippage-like and 0 for non

=head1 SYNOPSIS

    update_segment.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password

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

