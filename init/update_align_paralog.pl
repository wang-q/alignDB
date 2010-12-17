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
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

my $datalib        = $Config->{paralog}->{datalib};
my $use_megablast  = $Config->{paralog}->{use_megablast};
my $alignment_view = $Config->{paralog}->{alignment_view};

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

$Config->{paralog}->{datalib}        = $datalib;
$Config->{paralog}->{use_megablast}  = $use_megablast;
$Config->{paralog}->{alignment_view} = $alignment_view;

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
# Normal db operations
#----------------------------------------------------------#
# add a column align_feature4 to align_extra
{
    $obj->create_column( "align_extra", "align_feature4", "DOUBLE" );
    print "Table align_extra altered\n";
}

# check values in column
{
    my $sql_query = qq{
        SELECT COUNT(align_extra_id)
        FROM align_extra
    };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute;
    my ($count) = $sth->fetchrow_array;

    unless ($count) {
        $sql_query = qq{
            INSERT INTO align_extra (align_id)
            SELECT align.align_id
            FROM align
        };
        $sth = $dbh->prepare($sql_query);
        $sth->execute;
    }
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{

    # update align table in the new feature column
    my $align_extra = q{
        UPDATE align_extra
        SET align_feature4 = ?
        WHERE align_id = ?
    };
    my $align_extra_sth = $dbh->prepare($align_extra);

    # sequence
    my $seq_query = q{
        SELECT s.seq_length, t.target_seq
        FROM align a, target t, sequence s
        WHERE a.align_id = ?
        AND a.align_id = t.align_id
        AND s.seq_id = t.seq_id
    };
    my $seq_sth = $dbh->prepare($seq_query);

    # alignments
    my $align_ids = $obj->get_align_ids;

    # for each align
    for my $align_id ( @{$align_ids} ) {
        print "\nProcessing align_id $align_id\n";

        $seq_sth->execute($align_id);
        my ( $target_length, $target_seq ) = $seq_sth->fetchrow_array;
        $target_seq =~ tr/-//d;

        my $report_ref  = $paralog_obj->web_blast( \$target_seq );
        my $paralog_set = $paralog_obj->paralog_cover($report_ref);

        my $coverage = $paralog_set->cardinality / $target_length;

        $align_extra_sth->execute( $coverage, $align_id );

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
    $obj->add_meta_stopwatch($stopwatch);
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

