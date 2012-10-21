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
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'axt_dir=s'  => \$dir,
    'target=s'   => \$target,
    'parallel=i' => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Search for all files and push their paths to @axt_files
#----------------------------------------------------------#
my @files = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )
    ->in($dir);
printf "\n----Total .fa Files: %4s----\n\n", scalar @files;

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
    my $chr_id = $id_hash->{$chr_name};
    return unless $chr_id;

    for ( my $pos = 0; $pos < $chr_length; $pos += $truncated_length ) {
        my $seq = substr( $chr_seq, $pos, $truncated_length );
        
        my $target_info = {
            taxon_id   => $target_taxon_id,
            name       => $target_name,
            chr_id     => $chr_id,
            chr_name   => $chr_name,
            chr_start  => $pos + 1,
            chr_end    => $pos + length($seq),
            chr_strand => '+',
            seq        => $seq,
        };
        my $query_info = $target_info;
    
        $obj->add_align( $target_info, $query_info );

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

$stopwatch->block_message( "All files have been processed.", "duration" );

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

>perl init_alignDB.pl -d athvsself
>perl gen_alignDB_genome.pl -d athvsself -t "3702,Ath" -a e:\data\alignment\arabidopsis\ath_58\  --parallel 4

>perl init_alignDB.pl -d nipvsself
>perl gen_alignDB_genome.pl -d nipvsself -t "39947,Nip" -a e:\data\alignment\rice\nip_58\  --parallel 4

>perl init_alignDB.pl -d 9311vsself
>perl gen_alignDB_genome.pl -d 9311vsself -t "39946,9311" -a e:\data\alignment\rice\9311_58\  --parallel 4

>perl init_alignDB.pl -d S288Cvsself
>perl gen_alignDB_genome.pl -d S288Cvsself -t "4932,S288C" -a d:\data\alignment\yeast65\S288C\  --parallel 4
$perl gen_alignDB_genome.pl -d S288Cvsself -t "4932,S288C" -a /home/wangq/data/alignment/yeast65/S288C/  --parallel 4