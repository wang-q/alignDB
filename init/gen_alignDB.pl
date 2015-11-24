#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;

use AlignDB::Run;
use AlignDB::Stopwatch;

use lib "$FindBin::Bin/../lib";
use AlignDB;

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
        --target        STR     "target_taxon_id,target_name"
        --query         STR     "query_taxon_id,query_name"
        --length    -l  INT     threshold of alignment length
        --parallel      INT     run in parallel mode
        --gzip                  open .axt.gz files

=cut

# target, query init values
my $target_taxon_id = $Config->{taxon}{target_taxon_id};
my $target_name     = $Config->{taxon}{target_name};
my $query_taxon_id  = $Config->{taxon}{query_taxon_id};
my $query_name      = $Config->{taxon}{query_name};

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'         => \( my $server           = $Config->{database}{server} ),
    'port|P=i'           => \( my $port             = $Config->{database}{port} ),
    'db|d=s'             => \( my $db               = $Config->{database}{db} ),
    'username|u=s'       => \( my $username         = $Config->{database}{username} ),
    'password|p=s'       => \( my $password         = $Config->{database}{password} ),
    'dir_align|dir|da=s' => \( my $dir_align        = $Config->{taxon}{dir_align} ),
    'target=s'           => \( my $target           = $target_taxon_id . "," . $target_name ),
    'query=s'            => \( my $query            = $query_taxon_id . "," . $query_name ),
    'length|lt|l=i'      => \( my $length_threshold = $Config->{generate}{length_threshold} ),
    'parallel=i'         => \( my $parallel         = $Config->{generate}{parallel} ),
    'gzip'               => \my $gzip,
) or HelpMessage(1);

#----------------------------------------------------------#
# update names
#----------------------------------------------------------#
{
    my ( $target_taxon_id, $target_name ) = split ",", $target;
    my ( $query_taxon_id,  $query_name )  = split ",", $query;
    $target_name = $target_taxon_id unless $target_name;
    $query_name  = $query_taxon_id  unless $query_name;

    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->update_names( { $target_taxon_id => $target_name, $query_taxon_id => $query_name } );
}

#----------------------------------------------------------#
# Search for all files and push their paths to @files
#----------------------------------------------------------#
my @files;
if ( !$gzip ) {
    @files = sort File::Find::Rule->file->name('*.axt')->in($dir_align);
    printf "\n----Total .axt Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files = sort File::Find::Rule->file->name('*.axt.gz')->in($dir_align);
    printf "\n----Total .axt.gz Files: %4s----\n\n", scalar @files;
    $gzip++;
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

    my ( $target_taxon_id, ) = split ",", $target;
    my ( $query_taxon_id, )  = split ",", $query;

    die "target_taxon_id not defined\n" unless $target_taxon_id;
    die "query_taxon_id not defined\n"  unless $query_taxon_id;

    $obj->parse_axt_file(
        $infile,
        {   target_taxon_id => $target_taxon_id,
            query_taxon_id  => $query_taxon_id,
            threshold       => $length_threshold,
            gzip            => $gzip,
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
