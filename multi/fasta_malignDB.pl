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
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# program parameters
my $fasta_dir        = '';
my $length_thredhold = 5000;

# run in parallel mode
my $parallel = $Config->{generate}->{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{feature}->{batch};

my $help = 0;
my $man  = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'server=s'    => \$server,
    'port=i'      => \$port,
    'db=s'        => \$db,
    'username=s'  => \$username,
    'password=s'  => \$password,
    'length=i'    => \$length_thredhold,
    'fasta_dir=s' => \$fasta_dir,
    'parallel=i'  => \$parallel,
    'batch=i'     => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# call init_alignDB.pl
#----------------------------------------------------------#
system "perl $FindBin::Bin/../init/init_alignDB.pl"
    . " -d=$db -i=$FindBin::Bin/../minit.sql";

#----------------------------------------------------------#
# Search for all files and push their paths to @fasta_files
#----------------------------------------------------------#
my @fasta_files = File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )
    ->in($fasta_dir);
@fasta_files = sort @fasta_files;

printf "\n----Total .FAS Files: %4s----\n\n", scalar @fasta_files;

my @jobs;
{
    while ( scalar @fasta_files ) {
        my @batching = splice @fasta_files, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job     = shift;
    my @infiles = @$job;

    my $obj = AlignDB::Multi->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    for my $infile (@infiles) {
        print "process " . basename($infile);
        $obj->parse_fasta_file($infile);
    }
};

#----------------------------------------------------------#
# start insert
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
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
            --fasta_dir         .fas files' directory
            --length            threshold of alignment length
            --parallel          run in parallel mode

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

