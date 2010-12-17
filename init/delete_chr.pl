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

my $chr_id;

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
    'chr_id=i'   => \$chr_id,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

pod2usage(1) unless $chr_id;

#----------------------------------------------------------#
# Delete from alignDB
#----------------------------------------------------------#
$stopwatch->start_message("Delete from $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# delete align
my $align_id_sth = $dbh->prepare(
    'SELECT t.align_id
    FROM sequence s, target t
    WHERE s.chr_id = ?
    AND t.seq_id = s.seq_id'
);

my $delete_align_extra_sth = $dbh->prepare(
    'DELETE FROM align_extra
    WHERE align_id = ?'
);

my $delete_align_sth = $dbh->prepare(
    'DELETE FROM align
    WHERE align_id = ?'
);

# delete target
my $target_seq_id_sth = $dbh->prepare(
    'SELECT seq_id
    FROM target t
    WHERE align_id = ?'
);

my $delete_target_sth = $dbh->prepare(
    'DELETE FROM target
    WHERE align_id = ?'
);

# delete query
my $query_seq_id_sth = $dbh->prepare(
    'SELECT seq_id
    FROM query
    WHERE align_id = ?'
);

my $delete_query_sth = $dbh->prepare(
    'DELETE FROM query
    WHERE align_id = ?'
);

# delete sequence
my $delete_sequence_sth = $dbh->prepare(
    'DELETE FROM sequence
    WHERE seq_id = ?'
);

# delete snp
my $snp_id_sth = $dbh->prepare(
    'SELECT snp_id
    FROM snp
    WHERE align_id = ?'
);

my $delete_ssw_sth = $dbh->prepare(
    'DELETE FROM ssw
    WHERE snp_id = ?'
);

my $delete_snp_extra_sth = $dbh->prepare(
    'DELETE FROM snp_extra
    WHERE snp_id = ?'
);

my $delete_snp_sth = $dbh->prepare(
    'DELETE FROM snp
    WHERE snp_id = ?'
);

# delete indel
my $indel_id_sth = $dbh->prepare(
    'SELECT indel_id
    FROM indel
    WHERE align_id = ?'
);

my $delete_indel_extra_sth = $dbh->prepare(
    'DELETE FROM indel_extra 
    WHERE indel_id = ?'
);

my $delete_indel_sth = $dbh->prepare(
    'DELETE FROM indel
    WHERE indel_id = ?'
);

# delete isw
my $isw_id_sth = $dbh->prepare(
    'SELECT isw_id
    FROM isw
    WHERE indel_id = ?'
);

my $delete_isw_extra_sth = $dbh->prepare(
    'DELETE FROM isw_extra 
    WHERE isw_id = ?'
);

my $delete_isw_sth = $dbh->prepare(
    'DELETE FROM isw 
    WHERE isw_id = ?'
);

# delete gsw
my $delete_gsw_sth = $dbh->prepare(
    'DELETE gsw, window FROM gsw, window
    WHERE gsw.window_id = window.window_id
    AND window.align_id = ?'
);

# delete extreme
my $delete_extreme_sth = $dbh->prepare(
    'DELETE extreme, window FROM extreme, window
    WHERE extreme.window_id = window.window_id
    AND window.align_id = ?'
);

# delete segment
my $delete_segment_sth = $dbh->prepare(
    'DELETE segment, window FROM segment, window
    WHERE segment.window_id = window.window_id
    AND window.align_id = ?'
);

$align_id_sth->execute($chr_id);
while ( my ($align_id) = $align_id_sth->fetchrow_array ) {
    print "Delete align $align_id\n";

    # delete target
    $target_seq_id_sth->execute($align_id);
    my ($target_seq_id) = $target_seq_id_sth->fetchrow_array;
    $delete_sequence_sth->execute($target_seq_id);
    $delete_target_sth->execute($align_id);

    # delete query
    $query_seq_id_sth->execute($align_id);
    my ($query_seq_id) = $query_seq_id_sth->fetchrow_array;
    $delete_sequence_sth->execute($query_seq_id);
    $delete_query_sth->execute($align_id);

    # delete snp
    $snp_id_sth->execute($align_id);
    while ( my ($snp_id) = $snp_id_sth->fetchrow_array ) {
        $delete_ssw_sth->execute($snp_id);
        $delete_snp_extra_sth->execute($snp_id);
        $delete_snp_sth->execute($snp_id);
    }

    # delete indel
    $indel_id_sth->execute($align_id);
    while ( my ($indel_id) = $indel_id_sth->fetchrow_array ) {

        # delete isw
        $isw_id_sth->execute($indel_id);
        while ( my ($isw_id) = $isw_id_sth->fetchrow_array ) {
            $delete_isw_extra_sth->execute($isw_id);
            $delete_isw_sth->execute($isw_id);
        }

        $delete_indel_extra_sth->execute($indel_id);
        $delete_indel_sth->execute($indel_id);
    }

    # delete gsw
    $delete_gsw_sth->execute($align_id);

    # delete extreme
    $delete_extreme_sth->execute($align_id);

    # delete segment
    $delete_segment_sth->execute($align_id);

    # delete align
    $delete_align_extra_sth->execute($align_id);
    $delete_align_sth->execute($align_id);
}

$stopwatch->end_message();

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    $obj->add_meta_stopwatch($stopwatch);
}
exit;

__END__

=head1 NAME

    delete_chr.pl - removed alignments of a certain chr_id,
                    which is usually an unfinished one

=head1 SYNOPSIS

    delete_chr.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --db                database name
        --username          username
        --password          password
        --chr_id            chr_id of alignments to be removed

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
