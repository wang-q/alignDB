#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;

use MongoDB;
$MongoDB::BSON::looks_like_number = 1;
$MongoDB::BSON::utf8_flag_on      = 0;
use MongoDB::OID;

use MCE;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

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

# dir
my $dir;

my $target;    # target sequence

# mongodb has write lock, so we should make perl busy
my $truncated_length = 500_000;

my $fill       = 50;
my $min_length = 5000;

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    's|server=s' => \$server,
    'P|port=i'   => \$port,
    'd|db=s'     => \$dbname,
    'dir=s'      => \$dir,
    'target=s'   => \$target,
    'length=i'   => \$truncated_length,
    'fill=i'     => \$fill,
    'min=i'      => \$min_length,
    'parallel=i' => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Search for all files and push their paths to @axt_files
#----------------------------------------------------------#
my @files
    = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )->in($dir);
printf "\n----Total .fa Files: %4s----\n\n", scalar @files;

{    # update names
    my ( $target_taxon_id, $target_name ) = split ",", $target;
    $target_name = $target_taxon_id unless $target_name;

    my $mongo = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    );
    my $db = $mongo->get_database($dbname);

    my $coll_taxon = $db->get_collection('taxon');
    $coll_taxon->update( { 'taxon_id' => $target_taxon_id },
        { '$set' => { 'common_name' => $target_name } } );

    # chromosomes' taxon info is denormalized
    my $coll_chr = $db->get_collection('chromosome');
    $coll_chr->update(
        { 'taxon_id' => $target_taxon_id },
        { '$set'     => { 'taxon.common_name' => $target_name } },
        { 'multiple' => 1 },
    );
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $infile = $chunk_ref->[0];
    
    my $wid = MCE->wid;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message(
        "* Process task [$chunk_id] by worker #$wid\n* File [$infile]...");

    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my ( $target_taxon_id, $target_name ) = split ",", $target;

    die "target_taxon_id not defined\n" unless $target_taxon_id;
    $target_name = $target_taxon_id unless $target_name;
    my $chr_name = basename($infile);
    $chr_name =~ s/\..+?$//;

    my ( $seq_of, $seq_names ) = read_fasta($infile);
    my $chr_seq    = $seq_of->{ $seq_names->[0] };
    my $chr_length = length $chr_seq;

    # find chromosome OID
    my $coll_chr = $db->get_collection('chromosome');
    my $chr_id
        = $coll_chr->find_one(
        { 'taxon.taxon_id' => $target_taxon_id, 'chr_name' => $chr_name } );
    return unless $chr_id;
    $chr_id = $chr_id->{_id};

    my $ambiguous_set = AlignDB::IntSpan->new;
    for ( my $pos = 0; $pos < $chr_length; $pos++ ) {
        my $base = substr $chr_seq, $pos, 1;
        if ( $base =~ /[^ACGT-]/i ) {
            $ambiguous_set->add( $pos + 1 );
        }
    }

    print "Ambiguous chromosome region for $chr_name:\n    "
        . $ambiguous_set->runlist . "\n";

    my $valid_set = AlignDB::IntSpan->new("1-$chr_length");
    $valid_set->subtract($ambiguous_set);
    $valid_set = $valid_set->fill( $fill - 1 );   # fill gaps smaller than $fill

    print "Valid chromosome region for $chr_name:\n    "
        . $valid_set->runlist . "\n";

    my @regions;    # ([start, end], [start, end], ...)
    for my $set ( $valid_set->sets ) {
        my $size = $set->size;
        next if $size < $min_length;

        my @set_regions;
        my $pos = $set->min;
        my $max = $set->max;
        while ( $max - $pos + 1 > $truncated_length ) {
            push @set_regions, [ $pos, $pos + $truncated_length - 1 ];
            $pos += $truncated_length;
        }
        if ( scalar @set_regions > 0 ) {
            $set_regions[-1]->[1] = $max;
        }
        else {
            @set_regions = ( [ $pos, $max ] );
        }
        push @regions, @set_regions;
    }

    # collections
    my $coll_seq   = $db->get_collection('sequence');
    my $coll_align = $db->get_collection('align');

    for my $region (@regions) {
        my ( $start, $end ) = @{$region};
        my $length  = $end - $start + 1;
        my $runlist = AlignDB::IntSpan->new("$start-$end")->runlist;
        my $seq     = substr $chr_seq, $start - 1, $length;

        my $data_seq = {
            taxon_id   => $target_taxon_id,
            name       => $target_name,
            chr_id     => $chr_id,
            chr_name   => $chr_name,
            chr_start  => $start,
            chr_end    => $end,
            chr_strand => '+',
            seq        => $seq,
            length     => $length,
            runlist    => $runlist,
            level      => 1,                  # top level
        };
        my $seq_id = $coll_seq->insert($data_seq);

        my $data_align = {
            chr_name   => $chr_name,
            chr_start  => $start,
            chr_end    => $end,
            chr_strand => '+',
            length     => $length,
            runlist    => $runlist,
            seq_id     => $seq_id,
        };
        $coll_align->insert($data_align);
    }

    $inner_watch->block_message( "$infile has been processed.", "duration" );

    return;
};

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $parallel, );
$mce->foreach( \@files, $worker, ); # foreach implies chunk_size => 1.

# indexing
{
    my $db = MongoDB::MongoClient->new(
        host          => $server,
        port          => $port,
        query_timeout => -1,
    )->get_database($dbname);

    my $coll_seq = $db->get_collection('sequence');
    $coll_seq->ensure_index( { 'chr_name'  => 1 } );
    $coll_seq->ensure_index( { 'chr_start' => 1 } );
    $coll_seq->ensure_index( { 'chr_end'   => 1 } );
    $coll_seq->ensure_index( { 'level'     => 1 } );

    my $coll_align = $db->get_collection('align');
    $coll_align->ensure_index( { 'chr_name'  => 1 } );
    $coll_align->ensure_index( { 'chr_start' => 1 } );
    $coll_align->ensure_index( { 'chr_end'   => 1 } );
    $coll_align->ensure_index( { 'seq_id'    => 1 } );
}

$stopwatch->end_message( "All files have been processed.", "duration" );

exit;

__END__

=head1 NAME

    gen_alignDB_genome.pl - Generate alignDB from fasta files

=head1 SYNOPSIS

    gen_alignDB_genome.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --username          username
        --password          password
        --db                database name
        --dir               genome files' directory
        --target            "target_taxon_id,target_name"
        --parallel          run in parallel mode

=cut

perl init/init_alignDB.pl -d Athvsself
perl init/gen_alignDB_genome.pl -d Athvsself -t "3702,Ath" --dir /home/wangq/data/alignment/arabidopsis19/ath_65  --parallel 4

>perl init_alignDB.pl -d nipvsself
>perl gen_alignDB_genome.pl -d nipvsself -t "39947,Nip" --dir e:\data\alignment\rice\nip_58\  --parallel 4

>perl init_alignDB.pl -d 9311vsself
>perl gen_alignDB_genome.pl -d 9311vsself -t "39946,9311" --dir e:\data\alignment\rice\9311_58\  --parallel 4

perl init/init_alignDB.pl -d S288Cvsself
perl init/gen_alignDB_genome.pl -d S288Cvsself -t "4932,S288C" --dir /home/wangq/data/alignment/yeast65/S288C/  --parallel 4
perl init/insert_gc.pl -d S288Cvsself --parallel 4

perl mg/init_mg.pl -d alignDB
perl mg/gen_mg.pl -d alignDB -t "4932,S288C" --dir d:\data\alignment\self_alignment\S288C\  --parallel 1
