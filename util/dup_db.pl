#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use IPC::Cmd qw(can_run);

use AlignDB::Stopwatch;

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

dup_db.pl - Duplicate a database or dump to a gzipped file

=head1 SYNOPSIS

    perl dup_db.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --file                  dump file name
        --goal                  dup db name
        --gzip                  dump file is gzipped
        --yes                   overwrite dump file

    perl dup_db.pl -d S288CvsRM11 -z # S288CvsRM11.dump.sql.gz
    perl dup_db.pl -f S288CvsRM11.dump.sql.gz -g S288CvsRM11_dup
    
=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'file|f=s'     => \( my $file_dump ),
    'goal|g=s'     => \( my $goal_db ),
    'gzip|z'       => \( my $gzip ),
    'yes|y'        => \( my $yes ),
) or HelpMessage(1);

#----------------------------------------------------------#
# Check external commands
#----------------------------------------------------------#
if ( !can_run('mysql') ) {
    die "[mysql] not found in your system.\n";
}

if ($file_dump) {
    if ( $file_dump =~ /\.gz$/ ) {
        if ( !$gzip ) {
            print "The file [$file_dump] seems like a gzipped file, set --gzip ON.\n";
            $gzip++;
        }
    }
    else {
        if ($gzip) {
            print "The name of file [$file_dump] doesn't end with '.gz'.\n";
            print "Append '.gz' to it.\n";
            $file_dump .= '.gz';
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
    elsif ($yes) {
        print "Overwrites [$file_dump]\n";
    }
    else {
        print "[$file_dump] exists!\n";
    }
}
else {
    $file_dump = $db . ".dump.sql";
    $file_dump .= '.gz' if $gzip;

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
        $duplicate = "gzip --stdout --decompress --force $file_dump | mysql $str $goal_db";
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
