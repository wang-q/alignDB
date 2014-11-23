#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
$MongoDB::BSON::utf8_flag_on      = 0;
use MongoDB::OID;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::WriteExcel;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new;

# Database init values
my $server = "localhost";
my $port   = 27017;
my $dbname = "alignDB";

my $outfile;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    's|server=s' => \$server,
    'P|port=i'   => \$port,
    'd|db=s'     => \$dbname,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$dbname.mg.xlsx" unless $outfile;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Do stat for $dbname...");

my $mongo = MongoDB::MongoClient->new(
    host          => $server,
    port          => $port,
    query_timeout => -1,
);
my $db = $mongo->get_database($dbname);

#----------------------------------------------------------#
# worksheet -- distance_to_trough
#----------------------------------------------------------#
#select
#    gsw_distance distance,
#    gsw_wave_length wave_length,
#    gsw_amplitude amplitude,
#    gsw_trough_gc trough_gc,
#    gsw_gradient gradient,
#    window_average_gc window_gc,
#    gsw_cv window_cv,
#    gsw_intra_cv intra_cv,
#    window_indel window_indel,
#    if(window_indel > 0, 1, 0) has_indel
#from
#    gsw g
#        inner join
#    window w ON g.window_id = w.window_id
#where
#    gsw_amplitude >= 0.1;
my $output = "$dbname.gsw.csv";
open my $fh, '>', $output;
my @headers = qw{
    distance_to_trough distance_to_crest window_gc window_cv trough_gc crest_gc
    gradient value flag
};

print {$fh} join( ",", @headers ), "\n";

my $coll = $db->get_collection('gsw');
my $cursor = $coll->find;
#my $cursor = $coll->find( { chr_name => 'chr1' } );
while ( my $object = $cursor->next ) {
    print {$fh} join( ",",
        $object->{distance},  $object->{distance_crest},
        $object->{gc},        $object->{gc_cv},
        $object->{trough_gc}, $object->{crest_gc},
        $object->{gradient},  $object->{bed_count},
        $object->{bed_count} > 0 ? 1 : 0, ),
        "\n";
}

$stopwatch->end_message;
exit;
