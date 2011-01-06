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

my $outdir;
my $align_id_runlist = 1;

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
    'output=s'   => \$outdir,
    'align_id=s' => \$align_id_runlist,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
$outdir ||= "./$db";
unless ( -e $outdir ) {
    mkdir $outdir, 0777
        or die "Cannot create \"$outdir\" directory: $!";
}
my $align_id_set = AlignDB::IntSpan->new($align_id_runlist);

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
my $dbh = $obj->dbh;

for my $align_id ( $align_id_set->elements ) {
    my $outfile = $outdir . "/align-$align_id.fasta";
    open my $out_fh, '>', $outfile;
    my ( $target_name, $query_name, $ref_name ) = $obj->get_names($align_id);

    my ( $target_seq, $query_seq ) = @{ $obj->get_seqs($align_id) };

    # reference sequences
    my $ref_seq_query = q{
        SELECT s.seq_seq
        FROM sequence s
        INNER JOIN reference r ON s.seq_id = r.seq_id
        WHERE s.align_id = ?
    };
    my $ref_seq_sth = $dbh->prepare($ref_seq_query);
    $ref_seq_sth->execute($align_id);
    my ($ref_seq) = $ref_seq_sth->fetchrow_array;

    if ( defined $ref_seq ) {
        print {$out_fh} ">$ref_name\n";
        print {$out_fh} $ref_seq, "\n";
    }
    print {$out_fh} ">$target_name\n";
    print {$out_fh} $target_seq, "\n";
    print {$out_fh} ">$query_name\n";
    print {$out_fh} $query_seq, "\n";

    close $out_fh;
}

print "all done!!!\n";

exit;

__END__

=head1 NAME

    align_fasta.pl - Generate multi-fasta alignment file of an align_id

=head1 SYNOPSIS

    align_fasta.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --output          output filename
       --align_id        align_id
       

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
