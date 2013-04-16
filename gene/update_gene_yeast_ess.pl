#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Statistics::Lite qw(median);
use Text::CSV_XS;

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

my $ess_file  = "$FindBin::Bin/Essential_ORFs.txt";
my $rich_file = "$FindBin::Bin/Deutschbauer-2005.csv";
my $chem_file = "$FindBin::Bin/Hillenmeyer-2008.csv";
my $rec_file  = "$FindBin::Bin/forWebORFs.txt";
my $interact_file  = "$FindBin::Bin/sgadata_costanzo2009_stringentCutoff_101120.txt";

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

# add a column gene_feature4 to gene, Essential_ORFs
# add a column gene_feature6 to gene, forWebORFs
# add a column gene_feature6 to gene, interact number
{
    $obj->create_column( "gene", "gene_feature4", "char(8)" );
    $obj->create_column( "gene", "gene_feature6", "DOUBLE" );
    $obj->create_column( "gene", "gene_feature7", "DOUBLE" );
    print "Table gene altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

{    # update essential

    my $update = sub {
        my $ary_ref = shift;
        my $value = shift;
        
        my $sql = qq{
            UPDATE  gene
            SET gene_feature4 = ?
            WHERE gene_stable_id = ?
            AND gene_feature4 IS NULL
        };
        my $sth = $dbh->prepare($sql);

        for my $name (@{$ary_ref}) {
            $sth->execute($value, $name);
        }

        $sth->finish;
    };
    
    my $ess_names = read_in_ess($ess_file);
    printf "Ess genes %d\n", scalar @{$ess_names};
    $update->($ess_names, 'ESS');

    my $rich_names = read_in_rich($rich_file);
    printf "Rich genes %d\n", scalar @{$rich_names};
    $update->($rich_names, 'RICH');

    my $chem_names = read_in_chem($chem_file);
    printf "Chem genes %d\n", scalar @{$chem_names};
    $update->($chem_names, 'CHEM');

    $dbh->do(
        q{
        # let NULL to be zero
        UPDATE  gene
        SET gene_feature4 = 'OTHER'
        WHERE gene_feature4 IS NULL
        }
    );
}

{    # update interact  number
    my %interact_of = %{ read_in_interact($interact_file) };
    printf "Interact %d\n", scalar keys %interact_of;

    my $gene_query = q{
        UPDATE  gene
        SET gene_feature7 = ?
        WHERE gene_stable_id = ?
    };
    my $gene_sth = $dbh->prepare($gene_query);

    foreach my $gene ( keys %interact_of ) {
        $gene_sth->execute( $interact_of{$gene}, $gene );
    }

    $gene_sth->finish;
}

{    # update recombination rate
    my %rec_of = %{ read_in_rec($rec_file) };

    my $gene_query = q{
        UPDATE  gene
        SET gene_feature6 = ?
        WHERE gene_stable_id = ?
    };
    my $gene_sth = $dbh->prepare($gene_query);

    foreach my $gene ( keys %rec_of ) {
        $gene_sth->execute( $rec_of{$gene}, $gene );
    }

    $gene_sth->finish;
}

$stopwatch->end_message;

exit;

sub read_in_ess {
    my $infile = shift;

    my @names;
    open my $infh, '<', $infile;
    while (<$infh>) {
        next unless /^\d+/;
        push @names, ( split /\s+/ )[1];
    }
    close $infh;

    return \@names;
}

sub read_in_rich {
    my $infile = shift;

    my @names;
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $infile;
    $csv->getline($csv_fh);    # bypass title line
    while ( my $row = $csv->getline($csv_fh) ) {
        push @names, $row->[0];
    }
    close $csv_fh;

    return \@names;
}

sub read_in_chem {
    my $infile = shift;

    my @names;
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $infile;
    $csv->getline($csv_fh);    # bypass title line
    while ( my $row = $csv->getline($csv_fh) ) {
        next if $row->[2] > 5;
        push @names, $row->[0];
    }
    close $csv_fh;

    return \@names;
}

sub read_in_interact {
    my $infile = shift;

    my %names;
    open my $infh, '<', $infile;
    while (<$infh>) {
        next if /^[^Y]/;    # gene stable id start with a 'Y'
        my @fields = split /\s+/;
        my $stable_id = shift @fields;
        $names{$stable_id}++;
    }
    close $infh;

    return \%names;
}

sub read_in_rec {
    my $infile = shift;

    my %rec_of;
    open my $infh, '<', $infile;
    while (<$infh>) {
        next if /^[^Y]/;    # gene stable id start with a 'Y'
        my @fields = split /\s+/;
        next if @fields < 2;
        my $stable_id = shift @fields;
        my $rec_rate  = median(@fields);
        $rec_of{$stable_id} = $rec_rate;
    }
    close $infh;

    return \%rec_of;
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

