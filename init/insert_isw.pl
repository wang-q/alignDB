#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use AlignDB::IntSpan;
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

insert_isw.pl - Update indel-sliding windows

=head1 SYNOPSIS

    perl insert_isw.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --outgroup  -o          alignments have an outgroup
        --parallel      INT     run in parallel mode
        --batch         INT     number of alignments in one child process

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server       = $Config->{database}{server} ),
    'port|P=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'       => \( my $db           = $Config->{database}{db} ),
    'username|u=s' => \( my $username     = $Config->{database}{username} ),
    'password|p=s' => \( my $password     = $Config->{database}{password} ),
    'outgroup|o'   => \my $outgroup,
    'parallel=i'   => \( my $parallel     = $Config->{generate}{parallel} ),
    'batch=i'      => \( my $batch_number = $Config->{generate}{batch} ),
) or Getopt::Long::HelpMessage(1);

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

    print "Emptying tables...\n";

    # empty tables
    $obj->empty_table( 'isw', );

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

    # for each alignment
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);
        $obj->insert_isw($align_id);
        $obj->isw_snp_fk($align_id);
        $obj->update_D_values($align_id) if $outgroup;
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
