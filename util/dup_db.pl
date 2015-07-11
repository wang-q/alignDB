#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use IPC::Cmd qw(can_run);

use FindBin;

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
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

my $goal_db;
my $file_dump;

my $gzip;

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
    'f|file=s'     => \$file_dump,
    'g|goal=s'     => \$goal_db,
    'z|gzip'       => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Check external commands
#----------------------------------------------------------#
unless ( can_run('mysql') ) {
    die "[mysql] not found in your system.\n";
}

if ($file_dump) {
    if ( $file_dump =~ /\.gz$/ ) {
        if ( !$gzip ) {
            print
                "The file [$file_dump] seems like a gzipped file, set --gzip ON.";
            $gzip++;
        }
    }
    else {
        if ($gzip) {
            print "The name of file [$file_dump] doesn't end with '.gz'.\n";
            print "Append '.gz' to it.\n";
            $file_dump += '.gz';
            print "It's [$file_dump] now.\n";
        }
    }
}

if ($gzip) {
    unless ( can_run('gzip') ) {
        die "[gzip] not found in your system.\n";
    }
}

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
$stopwatch->start_message("Operation start...");

my $str = " -h$server -P$port -u$username -p$password ";

# dump db
if ($file_dump) {
    if ( !-e $file_dump ) {
        my $dump;
        if ($gzip) {
            $dump = "mysqldump $str $db | gzip > $file_dump";
        }
        else {
            $dump = "mysqldump $str $db  > $file_dump";
        }

        print "===> dump\n$dump\n\n";
        system($dump);
    }
    else {
        print "[$file_dump] exists!\n";
    }
}
else {
    $file_dump = $db . ".dump.sql";
    $file_dump += '.gz' if $gzip;

    my $dump;
    if ($gzip) {
        $dump = "mysqldump $str $db | gzip > $file_dump";
    }
    else {
        $dump = "mysqldump $str $db  > $file_dump";
    }

    print "===> dump\n$dump\n\n";
    system($dump);
}

# load dump
if ($goal_db) {
    my $drop   = "mysql $str -e \"DROP DATABASE IF EXISTS $goal_db;\"";
    my $create = "mysql $str -e \"CREATE DATABASE $goal_db;\"";
    my $duplicate;

    if ($gzip) {
        $duplicate
            = "gzip --stdout --decompress --force $file_dump | mysql $str $goal_db";
    }
    else {
        $duplicate = "mysql $str $goal_db < $file_dump";
    }

    print "#drop\n$drop\n\n";
    system($drop);
    print "#create\n$create\n\n";
    system($create );
    print "#duplicate\n$duplicate\n\n";
    system($duplicate);
}

$stopwatch->end_message;

__END__

=head1 NAME

dup_db.pl - Duplicate a database or dump to a gzipped file

=head1 SYNOPSIS

    perl dup_db.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --db                database name
        --username          username
        --password          password
        --file              dump file name
        --goal              dup db name
        --gzip              dump file is gzipped

=cut
