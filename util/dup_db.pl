#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use IPC::Cmd qw(can_run);

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf_db   = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini")->{database};
my $stopwatch = AlignDB::Stopwatch->new;

my $description = <<'EOF';
Duplicate a database or dump to a gzipped file.

    perl util/dup_db.pl -d S288CvsRM11_1a -z # S288CvsRM11.dump.sql.gz
    perl util/dup_db.pl -f S288CvsRM11_1a.dump.sql.gz -g S288CvsRM11_1a_dup

Usage: perl %c [options]
EOF

(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    $description,
    [ 'help|h', 'display this message' ],
    [],
    ['Database init values'],
    [ 'server|s=s',   'MySQL IP/Domain', { default => $conf_db->{server} }, ],
    [ 'port=i',       'MySQL port',      { default => $conf_db->{port} }, ],
    [ 'username|u=s', 'username',        { default => $conf_db->{username} }, ],
    [ 'password|p=s', 'password',        { default => $conf_db->{password} }, ],
    [ 'db|d=s',       'database name',   { default => $conf_db->{db} }, ],
    [],
    [ 'file|f=s', 'dump file name', ],
    [ 'goal|g=s', 'dup db name', ],
    [ 'gzip|z',   'dump file is gzipped', ],
    [ 'yes|y',    'overwrite dump file', ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

if ( !can_run('mysql') ) {
    die "[mysql] not found in your system.\n";
}

if ( $opt->{file} ) {
    if ( $opt->{file} =~ /\.gz$/ ) {
        if ( !$opt->{gzip} ) {
            print "The file [$opt->{file}] seems like a gzipped file, set --gzip ON.\n";
            $opt->{gzip}++;
        }
    }
    else {
        if ( $opt->{gzip} ) {
            print "The name of file [$opt->{file}] doesn't end with '.gz'.\n";
            print "Append '.gz' to it.\n";
            $opt->{file} .= '.gz';
            print "It's [$opt->{file}] now.\n";
        }
    }
}

if ( $opt->{gzip} ) {
    unless ( can_run('gzip') ) {
        die "[gzip] not found in your system.\n";
    }
}

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
$stopwatch->start_message("Operation start...");

my $str = " -h$opt->{server} -P$opt->{port} -u$opt->{username} -p$opt->{password} ";

# dump db
if ( $opt->{file} ) {
    if ( !-e $opt->{file} or $opt->{yes} ) {
        my $dump;
        if ( $opt->{gzip} ) {
            $dump = "mysqldump $str $opt->{db} | gzip > $opt->{file}";
        }
        else {
            $dump = "mysqldump $str $opt->{db}  > $opt->{file}";
        }

        print "===> dump\n$dump\n\n";
        system($dump);
    }
    else {
        print "[$opt->{file}] exists!\n";
    }
}
else {
    $opt->{file} = $opt->{db} . ".dump.sql";
    $opt->{file} .= '.gz' if $opt->{gzip};

    my $dump;
    if ( $opt->{gzip} ) {
        $dump = "mysqldump $str $opt->{db} | gzip > $opt->{file}";
    }
    else {
        $dump = "mysqldump $str $opt->{db}  > $opt->{file}";
    }

    print "===> dump\n$dump\n\n";
    system($dump);
}

# load dump
if ( $opt->{goal} ) {
    my $drop   = "mysql $str -e \"DROP DATABASE IF EXISTS $opt->{goal};\"";
    my $create = "mysql $str -e \"CREATE DATABASE $opt->{goal};\"";
    my $duplicate;

    if ( $opt->{gzip} ) {
        $duplicate = "gzip --stdout --decompress --force $opt->{file} | mysql $str $opt->{goal}";
    }
    else {
        $duplicate = "mysql $str $opt->{goal} < $opt->{file}";
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
