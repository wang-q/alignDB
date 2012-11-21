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

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'parallel=i'   => \$parallel,
    'batch=i'      => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update isw-indel relationship of $db...");

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
    $obj->empty_table( 'ssw', 'with_window' );

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

    # for each alignment
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);
        $obj->insert_ssw($align_id);
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

    insert_ssw.pl - Update snp-sliding windows (for three-way alignments)

=head1 SYNOPSIS

    insert_ssw.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --db                database name
        --username          username
        --password          password

=cut
