#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use Path::Tiny;

use AlignDB::Run;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB;
use AlignDB::Outgroup;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

gen_alignDB.pl - Generate alignDB from .fas files

=head1 SYNOPSIS

    perl gen_alignDB.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --dir_align -da STR     .fas files' directory
        --outgroup  -o          alignments have an outgroup
        --length    -l  INT     threshold of alignment length
        --parallel      INT     run in parallel mode

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'     => \( my $server           = $Config->{database}{server} ),
    'port=i'         => \( my $port             = $Config->{database}{port} ),
    'db|d=s'         => \( my $db               = $Config->{database}{db} ),
    'username|u=s'   => \( my $username         = $Config->{database}{username} ),
    'password|p=s'   => \( my $password         = $Config->{database}{password} ),
    'dir_align|da=s' => \( my $dir_align        = $Config->{generate}{dir_align} ),
    'length|lt|l=i'  => \( my $length_threshold = $Config->{generate}{length_threshold} ),
    'parallel=i'     => \( my $parallel         = $Config->{generate}{parallel} ),
    'outgroup|o'     => \( my $outgroup ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Search for all files and push their paths to @files
#----------------------------------------------------------#
$dir_align = path($dir_align)->stringify;
my @files = sort File::Find::Rule->file->name('*.fas')->in($dir_align);
printf "\n----Total .fas Files: %4s----\n\n", scalar @files;
if ( scalar @files == 0 ) {
    @files = sort File::Find::Rule->file->name('*.fas.gz')->in($dir_align);
    printf "\n----Total .fas.gz Files: %4s----\n\n", scalar @files;
}

my @jobs = map { [$_] } @files;

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
        $obj->parse_fas_file( $infile, $opt );
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
    opt      => { threshold => $length_threshold, },
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
