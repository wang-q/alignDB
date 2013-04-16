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
use AlignDB::Util qw(:all);

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

# dir
my $dir;

my $target;    # target sequence

my $truncated_length = 500_000;

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'd|db=s'       => \$db,
    'dir=s'        => \$dir,
    'target=s'     => \$target,
    'length=i'     => \$truncated_length,
    'parallel=i'   => \$parallel,
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
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    my ( $target_taxon_id, $target_name ) = split ",", $target;
    $target_name = $target_taxon_id unless $target_name;

    $obj->update_names( { $target_taxon_id => $target_name } );
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

    my ( $target_taxon_id, $target_name ) = split ",", $target;

    die "target_taxon_id not defined\n" unless $target_taxon_id;
    $target_name = $target_taxon_id unless $target_name;

    my $chr_name = basename($infile);
    $chr_name =~ s/\..+?$//;

    my ( $seq_of, $seq_names ) = read_fasta($infile);
    my $chr_seq    = $seq_of->{ $seq_names->[0] };
    my $chr_length = length $chr_seq;

    my $id_hash = $obj->get_chr_id_hash($target_taxon_id);
    my $chr_id  = $id_hash->{$chr_name};
    return unless $chr_id;

    for ( my $pos = 0; $pos < $chr_length; $pos += $truncated_length ) {
        my $seq = substr( $chr_seq, $pos, $truncated_length );

        my $info_of = {
            $target_name => {
                taxon_id   => $target_taxon_id,
                name       => $target_name,
                chr_id     => $chr_id,
                chr_name   => $chr_name,
                chr_start  => $pos + 1,
                chr_end    => $pos + length($seq),
                chr_strand => '+',
                seq        => $seq,
            },
        };

        $obj->add_align(
            $info_of,
            [ $target_name, $target_name ],
            [ $seq,         $seq ],
        );
    }

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

>perl init_alignDB.pl -d athvsself
>perl gen_alignDB_genome.pl -d athvsself -t "3702,Ath" --dir e:\data\alignment\arabidopsis\ath_58\  --parallel 4

>perl init_alignDB.pl -d nipvsself
>perl gen_alignDB_genome.pl -d nipvsself -t "39947,Nip" --dir e:\data\alignment\rice\nip_58\  --parallel 4

>perl init_alignDB.pl -d 9311vsself
>perl gen_alignDB_genome.pl -d 9311vsself -t "39946,9311" --dir e:\data\alignment\rice\9311_58\  --parallel 4

>perl init_alignDB.pl -d S288Cvsself
>perl gen_alignDB_genome.pl -d S288Cvsself -t "4932,S288C" --dir d:\data\alignment\yeast65\S288C\  --parallel 4
$perl gen_alignDB_genome.pl -d S288Cvsself -t "4932,S288C" --dir /home/wangq/data/alignment/yeast65/S288C/  --parallel 4
