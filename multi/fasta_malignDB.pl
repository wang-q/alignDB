#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;

use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Multi;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# program parameters
my $dir_fa           = '';
my $length_thredhold = 5000;

my $block;         # input is galaxy style blocked fasta
my $file_id_of;    # taxon_id-name mapping file

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $gzip;          # open .fas.gz

my $help = 0;
my $man  = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    's|server=s'    => \$server,
    'P|port=i'      => \$port,
    'd|db=s'        => \$db,
    'u|username=s'  => \$username,
    'p|password=s'  => \$password,
    'l|lt|length=i' => \$length_thredhold,
    'dir=s'         => \$dir_fa,
    'id|id_of=s'    => \$file_id_of,
    'block'         => \$block,
    'parallel=i'    => \$parallel,
    'batch=i'       => \$batch_number,
    'gzip'          => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# call init_alignDB.pl
#----------------------------------------------------------#
system "perl $FindBin::Bin/../init/init_alignDB.pl"
    . " -d=$db -i=$FindBin::Bin/../minit.sql";

my $id_of = {};
{
    my $name_of = {};
    open my $fh, '<', $file_id_of;
    while (<$fh>) {
        chomp;
        my ( $id, $name ) = split /,/;
        $id_of->{$name} = $id;
        $name_of->{$id} = $name;
    }
    close $fh;

    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->update_names($name_of);
}

#----------------------------------------------------------#
# Search for all files and push their paths to @files
#----------------------------------------------------------#
my @files;
if ( !$gzip ) {
    @files = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )
        ->in($dir_fa);
    printf "\n----Total .fas Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files = sort File::Find::Rule->file->name( '*.fa.gz', '*.fas.gz',
        '*.fasta.gz' )->in($dir_fa);
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

    my $obj = AlignDB::Multi->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    for my $infile (@infiles) {
        print "process " . basename($infile) . "\n";
        $obj->parse_block_fasta_file( $infile, $opt );
        print "Done.\n\n";
    }
};

#----------------------------------------------------------#
# start insert
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
    opt      => {
        id_of     => $id_of,
        threshold => $length_thredhold,
        gzip      => $gzip,
    },
);
$run->run;

$stopwatch->block_message( "All files have been processed.", "duration" );

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    $stopwatch->block_message("Update $db...");
    my $obj = AlignDB::Multi->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    $obj->update_misc;
}

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}

__END__

=head1 NAME

    fasta_malignDB.pl - Initiate, generate and update malignDB

=head1 SYNOPSIS

    fasta_malignDB.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --db                database name
        --username          username
        --password          password
        --dir               .fas files' directory
        --length            threshold of alignment length
        --parallel          run in parallel mode

=cut

perl ~/Scripts/alignDB/multi/fasta_malignDB.pl -d S288CvsRM11Spar --block --id ~/data/alignment/yeast_combine/id2name.csv --dir ~/data/alignment/yeast_combine/S288CvsRM11Spar_mafft --length 5000 --paralle 1
