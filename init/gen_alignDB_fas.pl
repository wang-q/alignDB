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
use AlignDB::Outgroup;

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

# dir of alignments
my $dir_align = '';

# alignments have an outgroup
my $outgroup;

# program parameters
my $length_threshold = 5000;

my $block;         # input is galaxy style blocked fasta
my $file_id_of;    # taxon_id-name mapping file

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $gzip;          # open .gz

my $help = 0;
my $man  = 0;

GetOptions(
    'help|?'             => \$help,
    'man'                => \$man,
    's|server=s'         => \$server,
    'P|port=i'           => \$port,
    'u|username=s'       => \$username,
    'p|password=s'       => \$password,
    'd|db=s'             => \$db,
    'da|dir|dir_align=s' => \$dir_align,
    'o|outgroup'         => \$outgroup,
    'l|lt|length=i'      => \$length_threshold,
    'id|id_of=s'         => \$file_id_of,
    'block'              => \$block,
    'parallel=i'         => \$parallel,
    'batch=i'            => \$batch_number,
    'gzip'               => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# update names
#----------------------------------------------------------#
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
        ->in($dir_align);
    printf "\n----Total .fas Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files = sort File::Find::Rule->file->name( '*.fa.gz', '*.fas.gz',
        '*.fasta.gz' )->in($dir_align);
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
        print "process " . basename($infile) . "\n";
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
        id_of     => $id_of,
        threshold => $length_threshold,
        gzip      => $gzip,
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

=head1 NAME

    gen_alignDB_fas.pl - Generate alignDB from fas files

=head1 SYNOPSIS

    gen_alignDB_fas.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --username          username
        --password          password
        --db                database name
        --dir_align         .axt files' directory
        --length            threshold of alignment length
        --parallel          run in parallel mode

=cut

perl ~/Scripts/alignDB/multi/fasta_malignDB.pl -d S288CvsRM11Spar --block --id ~/data/alignment/yeast_combine/id2name.csv --dir ~/data/alignment/yeast_combine/S288CvsRM11Spar_mafft --length 5000 --paralle 1

perl d:/wq/Scripts/alignDB/init/init_alignDB.pl -d Acetobacter_pasteurianus
perl d:/wq/Scripts/alignDB/init/gen_alignDB_fas.pl -d Acetobacter_pasteurianus --block --id d:\data\alignment\bac_new\Acetobacter_pasteurianus\round2\id2name.csv --dir d:\data\alignment\bac_new\Acetobacter_pasteurianus\round2\Acetobacter_pasteurianus_mft\ --length 5000 --paralle 1

d:\data\alignment\bac_new\Acinetobacter_baumannii\round2\Acinetobacter_baumannii_mft\
