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
use AlignDB::Multi::GC;

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
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# AlignDB::GC options
my $window_size = $Config->{gc}->{window_size};
my $window_step = $Config->{gc}->{window_step};

my $insert_segment = $Config->{gc}->{insert_segment};

# run in parallel mode
my $parallel = $Config->{generate}->{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{feature}->{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'           => \$help,
    'man'              => \$man,
    'server=s'         => \$server,
    'port=i'           => \$port,
    'db=s'             => \$db,
    'username=s'       => \$username,
    'password=s'       => \$password,
    'insert_segment=s' => \$insert_segment,
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
    my $obj = AlignDB::Multi::GC->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # empty tables: segment
    $obj->empty_table( 'segment', 'with_window' );

    my @align_ids = @{ $obj->get_align_ids };

    while ( scalar @align_ids ) {
        my @batching = splice @align_ids, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# start insert
#----------------------------------------------------------#
my $worker = sub {
    my $job       = shift;
    my @align_ids = @$job;

    my $obj = AlignDB::Multi::GC->new(
        mysql       => "$db:$server",
        user        => $username,
        passwd      => $password,
        window_size => $window_size,
        window_step => $window_step,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # alignments' chromosomal location, target_seq and query_seq
    my $align_seq_query = q{
        SELECT c.chr_name,
               s.chr_start,
               s.chr_end,
               a.align_length,
               a.align_comparable_runlist
        FROM align a, target t, sequence s, chromosome c
        WHERE a.align_id = s.align_id
        AND t.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        AND a.align_id = ?
    };
    my $align_seq_sth = $dbh->prepare($align_seq_query);

    # for each alignment
    for my $align_id (@align_ids) {
        $align_seq_sth->execute($align_id);
        my ( $chr_name, $chr_start, $chr_end, $align_length,
            $comparable_runlist )
            = $align_seq_sth->fetchrow_array;

        print "prosess align $align_id ",
            "in $chr_name $chr_start - $chr_end\n";

        # comparable runlist
        my $comparable_set = AlignDB::IntSpan->new($comparable_runlist);

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
    AlignDB::Multi::GC->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}
exit;

__END__

=head1 NAME

    insert_gc_multi.pl - Add GC ralated tables to alignDB

=head1 SYNOPSIS

    insert_gc.pl [options]
      Options:
        --help            brief help message
        --man             full documentation
        --server          MySQL server IP/Domain name
        --db              database name
        --username        username
        --password        password
        --insert_gc       insert gc or not
        --insert_segment  insert segment or not
       
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
