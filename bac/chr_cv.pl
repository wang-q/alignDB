#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use File::Find::Rule;

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::GC;

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

# running options
my $base_dir = $Config->{bac}{base_dir};

my $stat_window_size = $Config->{gc}{stat_window_size};
my $stat_window_step = $Config->{gc}{stat_window_step};

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{bac}{db};

# run in parallel mode
my $parallel = $Config->{feature}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{feature}{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'server=s'     => \$server,
    'port=i'       => \$port,
    'db=s'         => \$db,
    'username=s'   => \$username,
    'password=s'   => \$password,
    'b|base_dir=s' => \$base_dir,
    'parallel=i'   => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Compute chromosome CV");

# read all filenames, then grep
my @fna_files = File::Find::Rule->file->name('*.fna')->in($base_dir);

my @jobs;
{

    # $db is not AlignDB , I just want use the segment_gc_stat methods in
    # AlignDB::GC
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    my $dbh = $obj->dbh;

    my $chr_sth = $dbh->prepare(
        qq{
        SELECT t1.accession, t0.organism_name
          FROM strain t0, seq t1
         WHERE t0.taxonomy_id = t1.taxonomy_id
           AND t1.replicon like 'chr%'
           AND t0.seq_ok = 1
        ORDER BY organism_name
        }
    );
    $chr_sth->execute;
    my $ary_ref = $chr_sth->fetchall_arrayref( [0] );
    my @accessions = map { ( @{$_} ) } @{$ary_ref};
    $chr_sth->finish;

    while ( scalar @accessions ) {
        my @batching = splice @accessions, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# start insert
#----------------------------------------------------------#
my $worker = sub {
    my $job        = shift;
    my @accessions = @$job;

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    AlignDB::GC->meta->apply($obj);
    my %opt = (
        stat_window_size => $stat_window_size,
        stat_window_step => $stat_window_step,
    );
    for my $key ( sort keys %opt ) {
        $obj->$key( $opt{$key} );
    }

    # for each chr
    for my $accession (@accessions) {

        print "Finding fasta file for $accession\n";
        my @files = grep {/$accession/} @fna_files;
        if ( @files == 0 ) {
            warn "Can't find fasta file for $accession\n";
            next;
        }
        elsif ( @files > 1 ) {
            warn "Find more than 1 fasta file for $accession\n";
            next;
        }
        else {
            print "OK, $files[0]\n";
        }

        print "Load $accession\n";
        my $in = Bio::SeqIO->new(
            -file   => $files[0],
            -format => 'Fasta'
        );
        my $seq_obj = $in->next_seq;

        print "Start calc...\n";
        insert_segment( $obj, $accession, $seq_obj );
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub insert_segment {
    my $obj       = shift;
    my $accession = shift;
    my $seq_obj   = shift;

    my $dbh     = $obj->dbh;
    my $seq     = $seq_obj->seq;
    my $length  = $seq_obj->length;
    my $seq_set = AlignDB::IntSpan->new("1-$length");

    my @segment_levels = (
        [ 'A', '',    '' ],
        [ 1,   5000,  5000 ],
        [ 2,   1000,  1000 ],
        [ 3,   500,   500 ],
    );

    for (@segment_levels) {

        my $segment_type = $_->[0];
        my $segment_size = $_->[1];
        my $segment_step = $_->[2];

        my @segment_site
            = $obj->segment( $seq_set, $segment_size, $segment_step );

        # prepare segment_insert
        my $segment_insert = $dbh->prepare(
            qq{
            INSERT INTO segment (
                segment_id, accession, segment_type,
                segment_gc_mean, segment_gc_std,
                segment_gc_cv, segment_gc_mdcw
            )
            VALUES (
                NULL, ?, ?,
                ?, ?,
                ?, ?
            )
            }
        );

        for my $segment_set (@segment_site) {
            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( [$seq], $segment_set );

            $segment_insert->execute( $accession, $segment_type, $gc_mean,
                $gc_std, $gc_cv, $gc_mdcw );
        }

        $segment_insert->finish;
    }

    return;
}

__END__

perl bac_batch.pl --base_dir=
