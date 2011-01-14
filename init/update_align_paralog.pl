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
use AlignDB::Paralog;

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

my $datalib        = $Config->{paralog}{datalib};
my $use_megablast  = $Config->{paralog}{use_megablast};
my $alignment_view = $Config->{paralog}{alignment_view};

my $man  = 0;
my $help = 0;

$|++;

GetOptions(
    'help|?'           => \$help,
    'man|m'            => \$man,
    'server=s'         => \$server,
    'port=i'           => \$port,
    'db|d=s'           => \$db,
    'username=s'       => \$username,
    'password=s'       => \$password,
    'datalib|da=s'     => \$datalib,
    'megablast|mega=s' => \$use_megablast,
    'view|v=s'         => \$alignment_view,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$Config->{paralog}{datalib}        = $datalib;
$Config->{paralog}{use_megablast}  = $use_megablast;
$Config->{paralog}{alignment_view} = $alignment_view;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update align paralog of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# the first hit is the sequence itself, so we want two hits
my $paralog_obj = AlignDB::Paralog->new(
    use_megablast  => $use_megablast,
    alignment_view => $alignment_view,
    datalib        => $datalib,
    Config         => $Config,
);

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{

    # update align table in the new feature column
    my $align_update = q{
        UPDATE align
        SET align_paralog = ?
        WHERE align_id = ?
    };
    my $align_update_sth = $dbh->prepare($align_update);

    # alignments
    my @align_ids = @{$obj->get_align_ids};

    # for each align
    for my $align_id ( @align_ids ) {
        $obj->process_message($align_id);

        my ( $target_seq ) = @{$obj->get_seqs($align_id)};
        $target_seq =~ tr/-//d;
        my $target_length = length $target_seq;

        my $report_ref  = $paralog_obj->web_blast( \$target_seq );
        my $paralog_set = $paralog_obj->paralog_cover($report_ref);

        my $coverage = $paralog_set->cardinality / $target_length;

        $align_update_sth->execute( $coverage, $align_id );

        print Dump {
            target_length => $target_length,
            coverage      => $coverage,
            paralog_set   => $paralog_set->run_list,
        };
    }
}

$stopwatch->end_message;

# store program running meta info to database
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

    update_align_paralog.pl - Add additional paralog info to alignDB
    
=head1 SYNOPSIS

    update_align_paralog.pl [options]
      Options:
        --help               brief help message
        --man                full documentation
        --server             MySQL server IP/Domain name
        --db                 database name
        --username           username
        --password           password
        --datalib|da         blast database
        --megablast|mega     use megablast or not
        --view|v             blast output format

    update_align_paralog.pl -d=Nipvs9311 -da=nip_chro --mega=1 -v=9

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

