#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

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

update_segment.pl - update extreme numbers in segments

=head1 SYNOPSIS

    perl update_segment.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server       = $Config->{database}{server} ),
    'port|P=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'       => \( my $db           = $Config->{database}{db} ),
    'username|u=s' => \( my $username     = $Config->{database}{username} ),
    'password|p=s' => \( my $password     = $Config->{database}{password} ),
) or HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update segment table of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# add a column segment_feature4 to segment
{
    $obj->create_column("segment", "segment_feature4", "DOUBLE");
    print "Table segment_feature4 altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    my @align_ids = @{ $obj->get_align_ids };

    # update
    my $segment_query = q{
        # extreme numbers in segments
        UPDATE  segment,
                (SELECT s.id id, COUNT(*) count
                FROM 
                       (SELECT s.segment_id id, 
                               w.window_start start, 
                               w.window_end end
                        FROM   segment s, window w
                        WHERE  s.window_id = w.window_id
                        AND    w.align_id = ?) s,
                       (SELECT e.extreme_id id, 
                               w.window_start start, 
                               w.window_end end
                        FROM   extreme e, window w
                        WHERE  e.window_id = w.window_id
                        AND    w.align_id = ?) e
                WHERE e.start BETWEEN s.start AND s.end
                OR e.end BETWEEN s.start AND s.end
                GROUP BY s.id) seg
        SET segment.segment_feature4 = seg.count
        WHERE segment.segment_id = seg.id
    };
    my $segment_sth = $dbh->prepare($segment_query);

    for my $align_id (@align_ids) {
        print "Processing align_id $align_id\n";
        $segment_sth->execute( $align_id, $align_id );
    }

    $segment_sth->finish;
}

$stopwatch->end_message;

# store program running meta info to database
END {
    $obj->add_meta_stopwatch($stopwatch);
}
exit;

__END__
