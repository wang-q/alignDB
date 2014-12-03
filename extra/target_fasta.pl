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

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

my $outfile;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=s'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'output=s'   => \$outfile,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
$outfile = "$db.fasta" unless $outfile;

open my $outfh, '>', $outfile;

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

my $dbh = $obj->dbh;

my $target_query = q{
    SELECT  tx.common_name,
            c.chr_name,
            s.chr_strand,
            s.chr_start,
            s.chr_end,
            t.target_id,
            s.seq_seq
    FROM sequence s, target t, chromosome c, taxon tx
    WHERE 1 = 1
    AND s.seq_id = t.seq_id
    AND s.chr_id = c.chr_id
    AND c.taxon_id = tx.taxon_id
    ORDER BY common_name, chr_name, chr_start
};
my $target_query_sth = $dbh->prepare($target_query);

# for each sequences
$target_query_sth->execute();
while ( my @row = $target_query_sth->fetchrow_array ) {
    my $target_seq = pop @row;
    $target_seq =~ s/-//g;
    print {$outfh} ">" . $row[0];
    print {$outfh} "." . $row[1];
    print {$outfh} "(" . $row[2] . ")";
    print {$outfh} ":" . $row[3];
    print {$outfh} "-" . $row[4];
    print {$outfh} "|species=" . $row[0];
    print {$outfh} ";target_id=" . $row[5];
    print {$outfh} "\n";
    print {$outfh} $target_seq, "\n";
}

close $outfh;

print "all done!!!\n";

exit;

__END__

=head1 NAME

    target_fasta.pl - Generate multi-fasta file of target sequences

=head1 SYNOPSIS

    perl target_fasta.pl -d db_name [-o db_name.fasta] 

=cut
