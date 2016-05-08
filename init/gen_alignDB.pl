#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use File::Find::Rule;

use AlignDB::Run;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB;

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

gen_alignDB.pl - Generate alignDB from axt files

=head1 SYNOPSIS

    perl gen_alignDB.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --dir_align -da STR     .axt files' directory
        --target        STR     target_name
        --query         STR     query_name
        --length    -l  INT     threshold of alignment length
        --parallel      INT     run in parallel mode
        --gzip                  open .axt.gz files

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'         => \( my $server           = $Config->{database}{server} ),
    'port|P=i'           => \( my $port             = $Config->{database}{port} ),
    'db|d=s'             => \( my $db               = $Config->{database}{db} ),
    'username|u=s'       => \( my $username         = $Config->{database}{username} ),
    'password|p=s'       => \( my $password         = $Config->{database}{password} ),
    'dir_align|dir|da=s' => \( my $dir_align        = $Config->{taxon}{dir_align} ),
    'target=s'           => \( my $target_name      = $Config->{taxon}{target_name} ),
    'query=s'            => \( my $query_name       = $Config->{taxon}{query_name} ),
    'length|lt|l=i'      => \( my $length_threshold = $Config->{generate}{length_threshold} ),
    'parallel=i'         => \( my $parallel         = $Config->{generate}{parallel} ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Search for all files and push their paths to @files
#----------------------------------------------------------#
my @files;
@files = sort File::Find::Rule->file->name('*.axt')->in($dir_align);
printf "\n----Total .axt Files: %4s----\n\n", scalar @files;
if ( scalar @files == 0 ) {
    @files = sort File::Find::Rule->file->name('*.axt.gz')->in($dir_align);
    printf "\n----Total .axt.gz Files: %4s----\n\n", scalar @files;
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $infile = shift;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("Process $infile...");

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    die "target_name not defined\n" unless length $target_name;
    die "query_name not defined\n"  unless length $query_name;

    $obj->parse_axt_file(
        $infile,
        {   tname     => $target_name,
            qname     => $query_name,
            threshold => $length_threshold,
        }
    );

    $inner_watch->block_message( "$infile has been processed.", "duration" );

    return;
};

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@files,
    code     => $worker,
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
