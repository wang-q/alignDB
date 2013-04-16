#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# AlignDB::GC options
my $wave_window_size = $Config->{gc}{wave_window_size};
my $wave_window_step = $Config->{gc}{wave_window_step};
my $vicinal_size     = $Config->{gc}{vicinal_size};
my $fall_range       = $Config->{gc}{fall_range};
my $gsw_size         = $Config->{gc}{gsw_size};
my $stat_window_size = $Config->{gc}{stat_window_size};
my $stat_window_step = $Config->{gc}{stat_window_step};

my $insert_gc      = $Config->{gc}{insert_gc};
my $insert_segment = $Config->{gc}{insert_segment};

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

# use alternative segment levels 200 .. 900, 1000 .. 5000
my $alt_level = 0;

my $one_level = 0;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'           => \$help,
    'man'              => \$man,
    's|server=s'       => \$server,
    'P|port=i'         => \$port,
    'u|username=s'     => \$username,
    'p|password=s'     => \$password,
    'd|db=s'           => \$db,
    'insert_gc=s'      => \$insert_gc,
    'insert_segment=s' => \$insert_segment,
    'alt_level'        => \$alt_level,
    'one_level'        => \$one_level,
    'parallel=i'       => \$parallel,
    'batch=i'          => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update GC tables of $db...");

my @jobs;
{
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    print "Emptying tables...\n";

    # empty tables: segment, gsw, extreme
    $obj->empty_table( 'segment', 'with_window' );
    $obj->empty_table( 'gsw',     'with_window' );
    $obj->empty_table( 'extreme', 'with_window' );

    my @align_ids = @{ $obj->get_align_ids };

    while ( scalar @align_ids ) {
        my @batching = splice @align_ids, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job       = shift;
    my @align_ids = @$job;

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    AlignDB::GC->meta->apply($obj);
    my %opt = (
        wave_window_size => $wave_window_size,
        wave_window_step => $wave_window_step,
        vicinal_size     => $vicinal_size,
        fall_range       => $fall_range,
        gsw_size         => $gsw_size,
        stat_window_size => $stat_window_size,
        stat_window_step => $stat_window_step,
        alt_level        => $alt_level,
        one_level        => $one_level,
    );
    for my $key ( sort keys %opt ) {
        $obj->$key( $opt{$key} );
    }

    # Database handler
    my $dbh = $obj->dbh;

    # alignments' chromosomal location, target_seq and query_seq
    my $align_seq_query = q{
        SELECT  a.align_length,
                a.align_comparable_runlist
        FROM    align a
        WHERE   a.align_id = ?
    };
    my $align_seq_sth = $dbh->prepare($align_seq_query);

    # for each alignment
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);
        $align_seq_sth->execute($align_id);
        my ( $align_length, $comparable_runlist )
            = $align_seq_sth->fetchrow_array;

        # comparable runlist
        my $comparable_set = AlignDB::IntSpan->new($comparable_runlist);

        if ($insert_gc) {
            $obj->insert_extreme( $align_id, $comparable_set, $align_length );
            $obj->insert_gsw( $align_id, $comparable_set );
        }

        if ($insert_segment) {
            $obj->insert_segment( $align_id, $comparable_set );
        }
    }

    return;
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

    insert_gc.pl - Add GC ralated tables to alignDB

=head1 SYNOPSIS

    insert_gc.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --insert_gc         insert gc or not
        --insert_segment    insert segment or not
        --parallel          run in parallel mode
        --batch             number of alignments process in one child process

=cut

perl /home/wangq/Scripts/alignDB/extra/multi_way_batch.pl -d HumanvsCGOR \
    -e human_65 --block --id 9606 -lt 5000 -st 0 --parallel 12 \
    -f /home/wangq/data/alignment/primates/HumanvsCGOR_mft --run 1

perl /home/wangq/Scripts/alignDB/init/insert_gc.pl -d=HumanvsCGOR --parallel 12

perl /home/wangq/Scripts/alignDB/util/dup_db.pl -d HumanvsCGOR -g HumanvsCGOR_alt_level
perl /home/wangq/Scripts/alignDB/init/insert_gc.pl -d HumanvsCGOR_alt_level --parallel 12 --alt_level
perl /home/wangq/Scripts/alignDB/extra/multi_way_batch.pl -d HumanvsCGOR_alt_level \
    -e human_65 --block --id 9606 -lt 5000 -st 0 --parallel 12 \
    -f /home/wangq/data/alignment/primates/HumanvsCGOR_mft \
    --run 21,30,40
perl /home/wangq/Scripts/alignDB/stat/gc_stat_factory.pl -d HumanvsCGOR_alt_level --alt_level -t 0
