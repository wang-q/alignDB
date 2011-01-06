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

my $target_id_query = q{
    SELECT  target_id
    FROM target t
};
my $target_id_query_sth = $dbh->prepare($target_id_query);

my $target_query = q{
    SELECT  CONCAT(">target_id|", t.target_id),
            c.chr_name,
            CONCAT(s.chr_start, "-", s.chr_end),
            s.seq_seq
    FROM sequence s
    INNER JOIN target t ON s.seq_id = t.seq_id
    INNER JOIN chromosome c ON s.chr_id = c.chr_id
    WHERE t.target_id = ?
};
my $target_query_sth = $dbh->prepare($target_query);

# for each chromosome
$target_id_query_sth->execute;
while ( my ($target_id) = $target_id_query_sth->fetchrow_array ) {
    $target_query_sth->execute($target_id);
    while ( my @row = $target_query_sth->fetchrow_array ) {
        my $target_seq = pop @row;
        $target_seq =~ s/-//g;
        my $head = join "|", @row;
        print {$outfh} $head,       "\n";
        print {$outfh} $target_seq, "\n";
    }
}
close $outfh;

print "all done!!!\n";

exit;

__END__

=head1 NAME

    target_fasta.pl - Generate multi-fasta file of target sequences

=head1 SYNOPSIS

    target_fasta.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --output          output filename
       

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
