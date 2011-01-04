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
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'username=s' => \$username,
    'password=s' => \$password,
    'db=s'       => \$db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update isw-indel relationship of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    
    # alter table isw if column isw_indel_id does not exist
    my $check_result = $obj->check_column('isw', 'isw_indel_id');
    if (!defined $check_result) {
        print "Alter table structure\n";
        my $alter_query = q{
            alter table `isw` 
                add column `isw_indel_id` int NULL  after `isw_id`
        };
        $dbh->do($alter_query);
    }

    # update query
    my $isw_query_1 = q{
        UPDATE  isw
        SET isw.isw_indel_id = isw.indel_id
        WHERE isw.isw_type IN ('R', 'S')
    };
    my $isw_query_2 = q{
        UPDATE  isw
        SET isw.isw_indel_id = isw.prev_indel_id
        WHERE isw.isw_type IN ('L')
    };
    my $isw_sth_1 = $dbh->prepare($isw_query_1);
    my $isw_sth_2 = $dbh->prepare($isw_query_2);

    print "Updating...\n";
    $isw_sth_1->execute;
    $isw_sth_2->execute;

    $isw_sth_1->finish;
    $isw_sth_2->finish;
}

$stopwatch->end_message;

# store program running meta info to database
END {
    $obj->add_meta_stopwatch($stopwatch);
}
exit;

__END__


=head1 NAME

    update_isw_indel_id.pl - Update isw-indel relationship

=head1 SYNOPSIS

    update_isw_indel_id.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --db                database name
        --username          username
        --password          password

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
