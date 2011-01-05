#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

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
my $server     = $Config->{database}->{server};
my $port       = $Config->{database}->{port};
my $username   = $Config->{database}->{username};
my $password   = $Config->{database}->{password};
my $db         = $Config->{database}->{db};

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
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

=head1 NAME

    update_segment.pl - update extreme numbers in segments

=head1 SYNOPSIS

    update_segment.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       

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

