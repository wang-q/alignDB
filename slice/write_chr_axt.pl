#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB;

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

write_chr_axt.pl - extract sequence of specific chromosomes from alignDB

=head1 SYNOPSIS

    write_chr_axt.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --target        STR     target_chr_id_runlist
        --query         STR     query_chr_id_runlist
        --flag          STR     undef => all, 1 => in same chr, 2 => between two chrs

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server             = $Config->{database}{server} ),
    'port|P=i'     => \( my $port               = $Config->{database}{port} ),
    'db|d=s'       => \( my $db                 = $Config->{database}{db} ),
    'username|u=s' => \( my $username           = $Config->{database}{username} ),
    'password|p=s' => \( my $password           = $Config->{database}{password} ),
    'target=s'     => \( my $target_chr_runlist = "1,5-20" ),
    'query=s'      => \( my $query_chr_runlist  = "2,21-37" ),
    'flag=s'       => \my $chr_flag,
) or HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Write .axt files from $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
# chr_ids
my $target_chr_set = AlignDB::IntSpan->new($target_chr_runlist);
my $query_chr_set  = AlignDB::IntSpan->new($query_chr_runlist);

# select all target_name and query_name in this database
my ( $target_name, $query_name ) = $obj->get_names;

# alignment
my $align_query = q{
    SELECT 
        a.align_id,   a.align_length
    FROM
        align a INNER JOIN sequence st
            ON a.align_id = st.align_id
        INNER JOIN target t
            ON st.seq_id = t.seq_id
        INNER JOIN sequence sq
            ON a.align_id = sq.align_id
        INNER JOIN query q
            ON sq.seq_id = q.seq_id
    WHERE
        1 = 1 
        AND st.chr_id = ? 
        AND sq.chr_id = ?
    ORDER BY
        a.align_id
};
my $align_sth = $dbh->prepare($align_query);

my @chr_pairs;
for my $tchr_id ( $target_chr_set->elements ) {
    for my $qchr_id ( $query_chr_set->elements ) {
        if ( defined $chr_flag and $chr_flag == 1 ) {    # in same chr
            if ( $tchr_id == $qchr_id ) {
                push @chr_pairs, [ $tchr_id, $qchr_id ];
            }
        }
        elsif ( defined $chr_flag and $chr_flag == 2 ) {    # between two chrs
            if ( $tchr_id != $qchr_id ) {
                push @chr_pairs, [ $tchr_id, $qchr_id ];
            }
        }
        else {
            push @chr_pairs, [ $tchr_id, $qchr_id ];
        }
    }
}

for my $chr_pair (@chr_pairs) {
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
        my $target_info      = $obj->get_target_info($align_id);
        my $target_chr_name  = $target_info->{chr_name};
        my $target_chr_start = $target_info->{chr_start};
        my $target_chr_end   = $target_info->{chr_end};

        # query
        my ($query_info)    = $obj->get_queries_info($align_id);
        my $query_chr_name  = $query_info->{chr_name};
        my $query_chr_start = $query_info->{chr_start};
        my $query_chr_end   = $query_info->{chr_end};
        my $query_strand    = $query_info->{query_strand};

        my ( $target_seq, $query_seq ) = @{ $obj->get_seqs($align_id) };

        my $score = ( $target_chr_end - $target_chr_start + 1 ) * 100;    # sham score

        # append axt file
        {
            unless ( -e $db ) {
                mkdir $db, 0777
                    or die "Cannot create \"$db\" directory: $!";
            }
            my $outfile = "$db/$target_chr_name" . "vs$query_chr_name" . ".axt";

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

        open my $outfh, '>>', $outfile;
        print {$outfh} "$tchr_name vs. $qchr_name\n";
        close $outfh;
    }
}

$stopwatch->end_message;
exit;

__END__
