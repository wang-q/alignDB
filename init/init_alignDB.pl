#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Stopwatch;

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

my $init_sql = "$FindBin::Bin/../init.sql";

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
    'init_sql=s' => \$init_sql,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
my $cmd    = "mysql -h$server -P$port -u$username -p$password ";
my $drop   = "-e \"DROP DATABASE IF EXISTS $db;\"";
my $create = "-e \"CREATE DATABASE $db;\"";

print "#drop\n" . "$cmd $drop\n";
system("$cmd $drop");
print "#create\n" . "$cmd $create\n";
system("$cmd $create");
print "#init\n" . "$cmd $db < $init_sql\n";
system("$cmd $db < $init_sql");

#----------------------------------------------------------#
# init object
#----------------------------------------------------------#
$stopwatch->start_message("Init $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
$obj->index_isw_indel_id;
my $dbh = $obj->dbh;

#----------------------------------------------------------#
# taxon
#----------------------------------------------------------#
{
    my $taxon_file = "$FindBin::Bin/taxon.csv";
    print "Use $taxon_file to Init table taxon\n";

    my $sql_string = qq{
        LOAD DATA LOCAL INFILE '$taxon_file'
        INTO TABLE taxon
        FIELDS TERMINATED BY ','
        IGNORE 1 LINES
    };
    my $load_sth = $dbh->prepare($sql_string);
    $load_sth->execute;
    $load_sth->finish;

    # Now retrieve data from the table.
    my $query_sth = $dbh->prepare(
        q{
        SELECT *
        FROM taxon
        WHERE common_name = "Human"
        }
    );
    $query_sth->execute;
    while ( my $ref = $query_sth->fetchrow_hashref() ) {
        my $taxon_id    = $ref->{'taxon_id'};
        my $genus       = $ref->{'genus'};
        my $species     = $ref->{'species'};
        my $common_name = $ref->{'common_name'};
        print "Found a row: taxon_id = $taxon_id, "
            . "genus = $genus, "
            . "species = $species, "
            . "common_name = $common_name\n";
    }
    $query_sth->finish;
}

#----------------------------------------------------------#
# chromosome
#----------------------------------------------------------#
{
    my $chr_file = "$FindBin::Bin/chr_length.csv";
    print "Use $chr_file to Init table chromosome\n";

    my $sql_string = qq{
        LOAD DATA LOCAL INFILE '$chr_file'
        INTO TABLE chromosome
        FIELDS TERMINATED BY ','
        IGNORE 1 LINES
        (taxon_id, chr_name, chr_length)
    };
    my $insert_sth = $dbh->prepare($sql_string);
    $insert_sth->execute;
    $insert_sth->finish;

    # Now retrieve data from the table.
    my $sth = $dbh->prepare("SELECT * FROM chromosome");
    $sth->execute;
    while ( my $ref = $sth->fetchrow_hashref() ) {
        my $taxon_id   = $ref->{'taxon_id'};
        my $chr_name   = $ref->{'chr_name'};
        my $chr_length = $ref->{'chr_length'};
        if ( defined $taxon_id and $taxon_id == 31033 ) {
            print "Found a row: taxon_id = $taxon_id, "
                . "chr_name = $chr_name, "
                . "chr_length = $chr_length\n";
        }
    }
    $sth->finish;
}

$stopwatch->end_message;

# store program running meta info to database
END {
    $obj->add_meta_stopwatch($stopwatch);
}

__END__


=head1 NAME

    init_alignDB.pl - Initiate alignDB

=head1 SYNOPSIS

    init_alignDB.pl [options]
      Options:
        --help          brief help message
        --man           full documentation
        --server        MySQL server IP/Domain name
        --port          MySQL server port
        --db            database name
        --username      username
        --password      password
        --init_sql      init sql filename
      
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
