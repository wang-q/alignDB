#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
$MongoDB::BSON::utf8_flag_on      = 0;
use Text::CSV_XS;

use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new;

# Database init values
my $server = "localhost";
my $port   = 27017;
my $dbname = "alignDB";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'   => \$help,
    'man'      => \$man,
    'server=s' => \$server,
    'port=i'   => \$port,
    'db=s'     => \$dbname,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init object
#----------------------------------------------------------#
$stopwatch->start_message("Init $dbname...");

my $mongo = MongoDB::MongoClient->new(
    host          => $server,
    port          => $port,
    query_timeout => -1,
);
my $db = $mongo->get_database($dbname);
$db->drop;

#----------------------------------------------------------#
# taxon
#----------------------------------------------------------#
{
    my $file = "$FindBin::Bin/../data/taxon.csv";
    print "Use $file to Init table taxon\n";

    my $coll = $db->get_collection('taxon');

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $file;
    my @fields = @{ $csv->getline($csv_fh) };    # title line
    while ( my $row = $csv->getline($csv_fh) ) {
        my $data = {};
        for my $i ( 0 .. $#fields ) {
            $data->{ $fields[$i] } = $row->[$i];
        }
        $coll->insert( $data, { safe => 1 }, );
    }
    close $csv_fh;
    $coll->ensure_index( { 'taxon_id' => 1 }, { 'unique' => 1 } );

    print "There are @{[$coll->count]} documents in collection taxon\n";

    # a test
    my $result = $coll->find_one( { 'common_name' => 'human' } );
    if ($result) {
        printf
            "Found a row: \n\ttaxon_id = %s, genus = %s, species = %s, common_name = %s\n\n",
            $result->{taxon_id}, $result->{genus}, $result->{species},
            $result->{common_name};
    }
}

#----------------------------------------------------------#
# chromosome
#----------------------------------------------------------#
{
    my $file = "$FindBin::Bin/../data/chr_length.csv";
    print "Use $file to Init table chromosome\n";

    my $coll_chr   = $db->get_collection('chromosome');
    my $coll_taxon = $db->get_collection('taxon');

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $file;
    my @fields = qw{taxon_id chr_name chr_length};
    $csv->getline($csv_fh);    # bypass title line

    my @chrs;
    while ( my $row = $csv->getline($csv_fh) ) {
        my $data = {};
        for my $i ( 0 .. $#fields ) {
            $data->{ $fields[$i] } = $row->[$i];
        }
        my $taxon
            = $coll_taxon->find_one( { 'taxon_id' => $data->{taxon_id} } );
        if ( !$taxon ) {
            printf "Taxon %s doesn't exist\n", $data->{taxon_id};
        }
        $data->{taxon} = $taxon;

        push @chrs, $data;
    }
    close $csv_fh;

    $coll_chr->batch_insert( \@chrs, { safe => 1 }, );

    $coll_chr->ensure_index( { 'taxon_id'  => 1 } );
    $coll_chr->ensure_index( { 'taxon._id' => 1 } );

    print
        "There are @{[$coll_chr->count]} documents in collection chromosome\n";

    # a test
    my $result = $coll_chr->find_one( { 'taxon.taxon_id' => 31033 } );
    if ($result) {
        printf
            "Found a row: \n\ttaxon_id = %s, chr_name = %s, chr_length = %s\n\n",
            $result->{taxon}{taxon_id}, $result->{chr_name},
            $result->{chr_length};
    }
}

$stopwatch->end_message;

__END__


=head1 NAME

    init_mg.pl - Initiate alignDB

=head1 SYNOPSIS

    init_alignDB.pl [options]
      Options:
        --help          brief help message
        --man           full documentation
        --server        MySQL server IP/Domain name
        --port          MySQL server port
        --db            database name
        --username      username
        --password      password
        --init_sql      init sql filename
      
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
