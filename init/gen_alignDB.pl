#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;

use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

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
my $dir_align = $Config->{taxon}{dir_align};

# target, query init values
my $target_taxon_id = $Config->{taxon}{target_taxon_id};
my $target_name     = $Config->{taxon}{target_name};
my $query_taxon_id  = $Config->{taxon}{query_taxon_id};
my $query_name      = $Config->{taxon}{query_name};

my $target = $target_taxon_id . "," . $target_name;    # target sequence
my $query  = $query_taxon_id . "," . $query_name;      # query sequence

# program parameter
my $length_threshold
    = $Config->{generate}{length_threshold};    # legnth threshold of align

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

my $gzip;                                       # open .gz

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
    'target=s'           => \$target,
    'query=s'            => \$query,
    'l|lt|length=i'      => \$length_threshold,
    'parallel=i'         => \$parallel,
    'gzip'               => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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
        )
        ->update_names(
        { $target_taxon_id => $target_name, $query_taxon_id => $query_name } );
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

=head1 NAME

    gen_alignDB.pl - Generate alignDB from axt files

=head1 SYNOPSIS

    gen_alignDB.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --username          username
        --password          password
        --db                database name
        --dir_align         .axt files' directory
        --target            "target_taxon_id,target_name"
        --query             "query_taxon_id,query_name"
        --length            threshold of alignment length
        --parallel          run in parallel mode

=cut
