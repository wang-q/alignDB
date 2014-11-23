#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::KaKs;
use AlignDB::Position;

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

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'username=s' => \$username,
    'password=s' => \$password,
    'db=s'       => \$db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update KaKs $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# add columns to gene
{
    $obj->create_column( "gene", "gene_ka",   "DOUBLE" );
    $obj->create_column( "gene", "gene_ks",   "DOUBLE" );
    $obj->create_column( "gene", "gene_kaks", "DOUBLE" );
    print "Table gene altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    my $gene_query = q{
        SELECT 
                g.gene_stable_id,   g.gene_tl_runlist, g.gene_strand
        FROM
            gene g INNER JOIN window w
                ON g.window_id = w.window_id
        WHERE
            g.gene_stable_id IS NOT NULL 
            AND g.gene_biotype = 'protein_coding' 
            AND g.gene_tl_runlist IS NOT NULL 
            AND g.gene_tl_runlist != '-' 
            AND g.gene_is_full = 1
            AND w.window_pi != 0
            AND w.align_id = ?
    };
    my $gene_sth = $dbh->prepare($gene_query);

    my $update_query = q{
        UPDATE  gene
        SET gene_ka = ?,
            gene_ks = ?,
            gene_kaks = ?
        WHERE gene_stable_id = ?
    };
    my $update_sth = $dbh->prepare($update_query);

    my @align_ids = @{ $obj->get_align_ids };
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);

        $gene_sth->execute($align_id);
    GENE: while ( my @row = $gene_sth->fetchrow_array ) {
            my ( $stable_id, $tl_runlist, $strand ) = @row;
            my $tl_set = AlignDB::IntSpan->new($tl_runlist);

            my @seqs = @{ $obj->get_seqs($align_id) };
            my %seq_of;
            for my $i ( 0 .. $#seqs ) {
                my $sub_seq = $tl_set->substr_span( $seqs[$i] );
                $sub_seq =~ s/\-//g;
                $sub_seq = revcom($sub_seq) if $strand eq '-';

                $sub_seq =~ /^(ATG|GTG|TTG)/ or next GENE;
                $sub_seq =~ /(TAA|TAG|TGA)$/ or next GENE;
                length($sub_seq) % 3 == 0    or next GENE;

                $seq_of{ "seq" . $i } = $sub_seq;
            }
            DumpFile( "$stable_id.yml", \%seq_of );

            print " " x 4, "gene $stable_id\n";

            #open my $out_fh, '>', "$stable_id.fasta";
            #for my $key (sort keys %seq_of) {
            #    print {$out_fh} ">$key\n";
            #    print {$out_fh} $seq_of{$key}, "\n";
            #}
            #close $out_fh;

            my $kaks = AlignDB::KaKs->new;
            $kaks->seq_of( {%seq_of} );
            eval { $kaks->run; };
            if ($@) {
                warn $@;
                next;
            }
            my ( @kas, @kss, @omegas );
            for my $result ( @{ $kaks->results } ) {
                push @kas,    $result->[2];
                push @kss,    $result->[3];
                push @omegas, $result->[4];
            }

            $update_sth->execute( average(@kas), average(@kss),
                average(@omegas), $stable_id );
        }
    }

    $update_sth->finish;
    $gene_sth->finish;
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    update_gene_kaks.pl - Add Ka/gene_ka, Ks/8, Ka/Ks/9

=head1 SYNOPSIS

    update_gene_kaks.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password

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

