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
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# chr_id_runlist
my $target_chr_id_runlist = $Config->{write}->{target_chr_id_runlist};
my $query_chr_id_runlist  = $Config->{write}->{query_chr_id_runlist};
my $chr_flag;    # undef => all, 1 => in same chr, 2 => between two chrs

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
    'target=s'   => \$target_chr_id_runlist,
    'query=s'    => \$query_chr_id_runlist,
    'flag=s'    => \$chr_flag,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Write .axt files from $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh();

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
# chr_ids
my $target_chr_id_set = AlignDB::IntSpan->new($target_chr_id_runlist);
my $query_chr_id_set  = AlignDB::IntSpan->new($query_chr_id_runlist);

# select all target_name and query_name in this database
my ( $target_name, $query_name ) = $obj->get_names();

# alignment
my $align_query = q{
    SELECT a.align_id, a.align_length
    FROM align a, target t, sequence s1, query q, sequence s2
    WHERE a.align_id = t.align_id
    AND t.seq_id = s1.seq_id
    AND s1.chr_id = ?
    AND a.align_id = q.align_id
    AND q.seq_id = s2.seq_id
    AND s2.chr_id = ?
    ORDER BY a.align_id
};
my $align_sth = $dbh->prepare($align_query);

# target's chromosomal location
my $target_query = q{
    SELECT  c.chr_name,
            s.chr_start,
            s.chr_end,
            t.target_seq
    FROM target t, sequence s, chromosome c
    WHERE t.seq_id = s.seq_id
    AND s.chr_id = c.chr_id
    AND t.align_id = ?
};
my $target_sth = $dbh->prepare($target_query);

# query's chromosomal location
my $query_query = q{
    SELECT  c.chr_name,
            s.chr_start,
            s.chr_end,
            q.query_seq,
            q.query_strand
    FROM query q, sequence s, chromosome c
    WHERE q.seq_id = s.seq_id
    AND s.chr_id = c.chr_id
    AND q.align_id = ?
};
my $query_sth = $dbh->prepare($query_query);

my @chr_pairs;
for my $tchr_id ( $target_chr_id_set->elements() ) {
    for my $qchr_id ( $query_chr_id_set->elements() ) {
        if (defined $chr_flag and $chr_flag == 1) {# in same chr
            if ($tchr_id == $qchr_id) {
                push @chr_pairs, [ $tchr_id, $qchr_id ];
            }
        }
        elsif (defined $chr_flag and $chr_flag == 2) {# between two chrs
            if ($tchr_id != $qchr_id) {
                push @chr_pairs, [ $tchr_id, $qchr_id ];
            }
        }
        else {
            push @chr_pairs, [ $tchr_id, $qchr_id ];
        }
    }
}

foreach my $chr_pair (@chr_pairs) {
    my ( $tchr_id, $qchr_id ) = @{$chr_pair};

    my ($tchr_name) = $obj->get_chr_info($tchr_id);
    my ($qchr_name) = $obj->get_chr_info($qchr_id);

    print Dump {
        target_chr_id   => $tchr_id,
        target_chr_name => $tchr_name,
        query_chr_id    => $qchr_id,
        query_chr_name  => $qchr_name,
    };

    # for each align sequence
    $align_sth->execute( $tchr_id, $qchr_id );
    my $axt_serial = 0;
    while ( my @row2 = $align_sth->fetchrow_array ) {
        my ( $align_id, $align_length ) = @row2;

        print "Processing align_id $align_id\n";

        # target
        $target_sth->execute($align_id);
        my ($target_chr_name, $target_chr_start,
            $target_chr_end,  $target_seq,
        ) = $target_sth->fetchrow_array;
        print " " x 4, "target info fetched", " " x 10, "\r";

        # query
        $query_sth->execute($align_id);
        my ($query_chr_name, $query_chr_start, $query_chr_end,
            $query_seq,      $query_strand,
        ) = $query_sth->fetchrow_array;
        print " " x 4, "query info fetched", " " x 10, "\r";

        my $score
            = ( $target_chr_end - $target_chr_start + 1 ) * 100;  # sham score

        # append axt file
        {
            unless ( -e $db ) {
                mkdir $db, 0777
                    or die "Cannot create \"$db\" directory: $!";
            }
            my $outfile
                = "$db/$target_chr_name" . "vs$query_chr_name" . ".axt";

            # append axt file
            open my $outfh, '>>', $outfile;
            print {$outfh} "$axt_serial";
            print {$outfh} " $target_chr_name";
            print {$outfh} " $target_chr_start $target_chr_end";
            print {$outfh} " $query_chr_name";
            print {$outfh} " $query_chr_start $query_chr_end";
            print {$outfh} " $query_strand $score\n";
            print {$outfh} $target_seq, "\n";
            print {$outfh} $query_seq, "\n";
            print {$outfh} "\n";
            close $outfh;

            $axt_serial++;
            print " " x 4, "$outfile\n";
        }
    }

    if ( $axt_serial == 0 ) {
        unless ( -e $db ) {
            mkdir $db, 0777
                or die "Cannot create \"$db\" directory: $!";
        }
        my $outfile = "$db/noresult";

        # append axt file
        open my $outfh, '>>', $outfile;
        print {$outfh} "$tchr_name vs. $qchr_name\n";
        close $outfh;
    }
}

$stopwatch->end_message();
exit;

__END__

=head1 NAME

    write_chr_axt.pl - extract sequence of specific chromosomes from alignDB

=head1 SYNOPSIS

    write_chr_axt.pl [options]
      Options:
        --help          brief help message
        --man           full documentation
        --server        MySQL server IP/Domain name
        --db            database name
        --username      username
        --password      password
        --target        target_chr_id_runlist
        --query         query_chr_id_runlist
        --flag          undef => all, 1 => in same chr, 2 => between two chrs

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

