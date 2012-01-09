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

# axt
my $axt_dir = $Config->{taxon}{axt_dir};

# target, query init values
my $target_taxon_id = $Config->{taxon}{target_taxon_id};
my $target_name     = $Config->{taxon}{target_name};
my $query_taxon_id  = $Config->{taxon}{query_taxon_id};
my $query_name      = $Config->{taxon}{query_name};

my $target = $target_taxon_id . "," . $target_name;    # target sequence
my $query  = $query_taxon_id . "," . $query_name;      # query sequence

# program parameter
my $axt_threshold
    = $Config->{generate}{axt_threshold};    # legnth threshold of align
my $insert_dG = $Config->{generate}{insert_dG};    # dG

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'server=s'    => \$server,
    'port=i'      => \$port,
    'db=s'        => \$db,
    'username=s'  => \$username,
    'password=s'  => \$password,
    'axt_dir=s'   => \$axt_dir,
    'target=s'    => \$target,
    'query=s'     => \$query,
    'length=i'    => \$axt_threshold,
    'insert_dG=s' => \$insert_dG,
    'parallel=i'  => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Search for all files and push their paths to @axt_files
#----------------------------------------------------------#
my @axt_files = sort File::Find::Rule->file->name('*.axt')->in($axt_dir);
printf "\n----Total .AXT Files: %4s----\n\n", scalar @axt_files;

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $infile = shift;

    my $inner_watch = AlignDB::Stopwatch->new;
    $inner_watch->block_message("Process $infile...");

    my $obj = AlignDB->new(
        mysql     => "$db:$server",
        user      => $username,
        passwd    => $password,
        insert_dG => $insert_dG,
    );

    my ( $target_taxon_id, $target_name ) = split ",", $target;
    my ( $query_taxon_id,  $query_name )  = split ",", $query;

    die "target_taxon_id not defined\n" unless $target_taxon_id;
    die "query_taxon_id not defined\n"  unless $query_taxon_id;
    $target_name = $target_taxon_id unless $target_name;
    $query_name  = $query_taxon_id  unless $query_name;

    $obj->parse_axt_file(
        {   axt_file        => $infile,
            target_taxon_id => $target_taxon_id,
            target_name     => $target_name,
            query_taxon_id  => $query_taxon_id,
            query_name      => $query_name,
            threshold       => $axt_threshold,
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
    jobs     => \@axt_files,
    code     => $worker,
);
$run->run;

$stopwatch->end_message( "All files have been processed." );

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
        --db                database name
        --username          username
        --password          password
        --axt_dir           .axt files' directory
        --target            "target_taxon_id,target_name"
        --query             "query_taxon_id,query_name"
        --length            threshold of alignment length
        --insert_dG         do deltaG related processes
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
