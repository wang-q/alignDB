#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use Path::Tiny;

use AlignDB::Run;
use AlignDB::Stopwatch;

use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Outgroup;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

gen_alignDB_fas.pl - Generate alignDB from fas files

=head1 SYNOPSIS

    perl gen_alignDB_fas.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --dir_align -da STR     .fas files' directory
        --outgroup  -o          alignments have an outgroup
        --block                 input is blocked fasta
        --length    -l  INT     threshold of alignment length
        --parallel      INT     run in parallel mode
        --batch         INT     number of alignments in one child process
        --gzip                  open .axt.gz files


    perl ~/Scripts/alignDB/multi/fasta_malignDB.pl -d S288CvsRM11Spar --block --id ~/data/alignment/yeast_combine/id2name.csv --dir ~/data/alignment/yeast_combine/S288CvsRM11Spar_mafft --length 5000 --paralle 1
    
    perl d:/wq/Scripts/alignDB/init/init_alignDB.pl -d Acetobacter_pasteurianus
    perl d:/wq/Scripts/alignDB/init/gen_alignDB_fas.pl -d Acetobacter_pasteurianus --block --id d:\data\alignment\bac_new\Acetobacter_pasteurianus\round2\id2name.csv --dir d:\data\alignment\bac_new\Acetobacter_pasteurianus\round2\Acetobacter_pasteurianus_mft\ --length 5000 --paralle 1
    
=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'         => \( my $server           = $Config->{database}{server} ),
    'port|P=i'           => \( my $port             = $Config->{database}{port} ),
    'db|d=s'             => \( my $db               = $Config->{database}{db} ),
    'username|u=s'       => \( my $username         = $Config->{database}{username} ),
    'password|p=s'       => \( my $password         = $Config->{database}{password} ),
    'dir_align|dir|da=s' => \( my $dir_align        = '' ),
    'outgroup|o'         => \my $outgroup,
    'block'              => \my $block,
    'length|lt|l=i'      => \( my $length_threshold = 5000 ),
    'parallel=i'         => \( my $parallel         = $Config->{generate}{parallel} ),
    'batch=i'            => \( my $batch_number     = $Config->{generate}{batch} ),
    'gzip'               => \my $gzip,
) or HelpMessage(1);

#----------------------------------------------------------#
# Search for all files and push their paths to @files
#----------------------------------------------------------#
my @files;
if ( !$gzip ) {
    @files = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )->in($dir_align);
    printf "\n----Total .fas Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files
        = sort File::Find::Rule->file->name( '*.fa.gz', '*.fas.gz', '*.fasta.gz' )->in($dir_align);
    printf "\n----Total .fas.gz Files: %4s----\n\n", scalar @files;
    $gzip++;
}

my @jobs;
if ( !$block ) {
    while ( scalar @files ) {
        my @batching = splice @files, 0, $batch_number;
        push @jobs, [@batching];
    }
}
else {
    @jobs = map { [$_] } @files;
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job = shift;
    my $opt = shift;

    my @infiles = @$job;

    my $obj;
    if ( !$outgroup ) {
        $obj = AlignDB->new(
            mysql  => "$db:$server",
            user   => $username,
            passwd => $password,
        );
    }
    else {
        $obj = AlignDB::Outgroup->new(
            mysql  => "$db:$server",
            user   => $username,
            passwd => $password,
        );
    }

    for my $infile (@infiles) {
        print "process " . path($infile)->basename . "\n";
        $obj->parse_block_fasta_file( $infile, $opt );
        print "Done.\n\n";
    }

    return;
};

#----------------------------------------------------------#
# start insert
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
    opt      => {
        threshold => $length_threshold,
    },
);
$run->run;

$stopwatch->end_message( "All files have been processed.", "duration" );

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}

exit;

__END__
