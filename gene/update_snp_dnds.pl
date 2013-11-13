#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Codon;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
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
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update dnds info of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

my $pos_obj = AlignDB::Position->new( dbh => $dbh );
my $codon_obj = AlignDB::Codon->new;

# add columns
{
    $obj->create_column( "snp", "exon_id",       "INT" );
    $obj->create_column( "snp", "snp_codon_pos", "INT" );
    $obj->create_column( "snp", "snp_codons",    "CHAR(128)" );
    $obj->create_column( "snp", "snp_syn",       "DOUBLE" );
    $obj->create_column( "snp", "snp_nsy",       "DOUBLE" );
    $obj->create_column( "snp", "snp_stop",      "DOUBLE" );
    print "Table snp altered\n";

    $obj->create_column( "isw", "isw_syn",  "DOUBLE" );
    $obj->create_column( "isw", "isw_nsy",  "DOUBLE" );
    $obj->create_column( "isw", "isw_stop", "DOUBLE" );
    print "Table isw altered\n";

    $obj->create_column( "gene", "gene_syn",  "DOUBLE" );
    $obj->create_column( "gene", "gene_nsy",  "DOUBLE" );
    $obj->create_column( "gene", "gene_stop", "DOUBLE" );
    print "Table gene altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
# select all exons which contain the snp
my $exon_query = q{
    SELECT  e.exon_id, e.exon_strand, 
            e.exon_tl_runlist, e.exon_seq, e.exon_peptide
    FROM exon e, window w
    WHERE e.window_id = w.window_id
    AND w.align_id = ?
    AND e.exon_tl_runlist != '-'
};
my $exon_sth = $dbh->prepare($exon_query);

# select all coding snps in this alignment
my $snp_query = q{
    SELECT  s.snp_id, s.snp_pos
    FROM snp s
    WHERE s.align_id = ?
    AND s.snp_coding = 1
};
my $snp_sth = $dbh->prepare($snp_query);

# update snp table in the new feature column
my $snp_update = q{
    UPDATE snp
    SET exon_id = ?,
        snp_codon_pos = ?,
        snp_codons = ?,
        snp_syn = ?,
        snp_nsy = ?,
        snp_stop = ?
    WHERE snp_id = ?
};
my $snp_update_sth = $dbh->prepare($snp_update);

my @align_ids = @{ $obj->get_align_ids };
ALIGN: for my $align_id (@align_ids) {

    my $target_info    = $obj->get_target_info($align_id);
    my $chr_name       = $target_info->{chr_name};
    my $target_runlist = $target_info->{seq_runlist};
    my $align_length   = $target_info->{align_length};

    next ALIGN if $chr_name =~ /rand|un|contig|hap|scaf/i;

    $obj->process_message($align_id);
    my ( $target_seq, @query_seqs ) = @{ $obj->get_seqs($align_id) };
    $_ = uc $_ for ( $target_seq, @query_seqs ) ;

    # target runlist
    my $target_set = AlignDB::IntSpan->new($target_runlist);

    # get all exons
    my @exons;
    $exon_sth->execute($align_id);
EXON: while ( my @row = $exon_sth->fetchrow_array ) {
        my ( $exon_id, $exon_strand, $exon_tl_runlist, $exon_seq,
            $exon_peptide, )
            = @row;
        
        $exon_seq = uc $exon_seq;

        # extract exon seq in this align
        my $exon_tl_set    = AlignDB::IntSpan->new($exon_tl_runlist);
        my $exon_align_seq = $exon_tl_set->substr_span($target_seq);
        if ( $exon_strand eq '-' ) {
            $exon_align_seq = revcom($exon_align_seq);
        }
        if ( index( $exon_seq, $exon_align_seq ) == -1 ) {
            print " " x 4, "Exon sequence does not match alignment\n";

            print Dump {
                exon_seq       => $exon_seq,
                exon_align_seq => $exon_align_seq
            };
            next EXON;
        }

        # determine exon (may be truncated) frame in this align
        my $exon_frame;
        for ( 0 .. 2 ) {
            my $exon_align_peptide
                = $codon_obj->translate( $exon_align_seq, $_ );
            $exon_align_peptide =~ s/\*//g;
            if ( index( $exon_peptide, $exon_align_peptide ) != -1 ) {
                $exon_frame = $_;
                last;
            }
        }
        if ( !defined $exon_frame ) {
            print " " x 4, "Can't judge exon frame\n";
            print Dump { exon_pep => $exon_peptide };
            next EXON;
        }

        my $exon_info = {
            exon_id     => $exon_id,
            exon_strand => $exon_strand,
            exon_frame  => $exon_frame,
            exon_tl_set => $exon_tl_set,
            exon_length => $exon_tl_set->cardinality,
        };
        push @exons, $exon_info;
    }
    next ALIGN if scalar @exons == 0;

    $snp_sth->execute($align_id);
SNP: while ( my @row = $snp_sth->fetchrow_array ) {
        my ( $snp_id, $snp_pos, ) = @row;

        # only analysis the first exon match this snp
        my $exon;
        for (@exons) {
            if ( $_->{exon_tl_set}->contain($snp_pos) ) {
                $exon = $_;
                last;
            }
        }
        if ( !defined $exon ) {
            print " " x 4, "Can't locate an exon\n";
            next SNP;
        }

        # get all exon infos
        my $exon_id     = $exon->{exon_id};
        my $exon_strand = $exon->{exon_strand};
        my $exon_tl_set = $exon->{exon_tl_set};
        my $exon_frame  = $exon->{exon_frame};
        my $exon_length = $exon->{exon_length};

        my $snp_codon_pos;    # 0, 1, 2
        my $codon_set;
        if ( $exon_strand eq '+' ) {

            # locate the snp in a codon
            # determine snp position in codon
            my $snp_tl_index = $exon_tl_set->index($snp_pos);
            if ( !defined $snp_tl_index ) {
                print " " x 4, "Can't calc snp translate index\n";
                next SNP;
            }
            $snp_codon_pos = ( $snp_tl_index - $exon_frame - 1 ) % 3;

            my $codon_start_index = $snp_tl_index - $snp_codon_pos;
            my $codon_start       = $exon_tl_set->at($codon_start_index);
            my $codon_end         = $exon_tl_set->at( $codon_start_index + 2 );
            if ( !defined $codon_start or !defined $codon_end ) {
                print " " x 4, "Can't calc codon positions\n";
                next SNP;
            }
            elsif ( $codon_start >= $codon_end ) {
                print " " x 4, "Codon start-end error\n";
                next SNP;
            }
            $codon_set = AlignDB::IntSpan->new("$codon_start-$codon_end");
        }
        else {

            # locate the snp in a codon
            # determine snp position in codon
            my $snp_tl_index = $exon_tl_set->index($snp_pos);
            if ( !defined $snp_tl_index ) {
                print " " x 4, "Can't calc snp translate index\n";
                next SNP;
            }
            my $snp_tl_backindex = $exon_length - $snp_tl_index + 1;
            $snp_codon_pos = ( $snp_tl_backindex - $exon_frame - 1 ) % 3;

            my $codon_start_index = $snp_tl_index - ( 2 - $snp_codon_pos );
            my $codon_start       = $exon_tl_set->at($codon_start_index);
            my $codon_end         = $exon_tl_set->at( $codon_start_index + 2 );
            if ( !defined $codon_start or !defined $codon_end ) {
                print " " x 4, "Can't calc codon positions\n";
                next SNP;
            }
            elsif ( $codon_start >= $codon_end ) {
                print " " x 4, "Codon start-end error\n";
                next SNP;
            }
            $codon_set = AlignDB::IntSpan->new("$codon_start-$codon_end");
        }
        if ( $codon_set->cardinality > 3 ) {
            print " " x 4, "Indels in this codon\n";
            next SNP;
        }

        # target codon and query codon
        my $codon_t = $codon_set->substr_span($target_seq);
        if ( $exon_strand eq '-' ) {
            $codon_t = revcom($codon_t);
        }
        my @codons = ($codon_t);
        for my $query_seq (@query_seqs) {
            my $codon_q = $codon_set->substr_span($query_seq);
            if ( $exon_strand eq '-' ) {
                $codon_q = revcom($codon_q);
            }
            push @codons, $codon_q;
        }

        for (@codons) {
            if (/\-/) {
                print " " x 4, "Indels in this codon\n";
                next SNP;
            }
        }

        my ( @syns, @nsys, @stops );
        for ( my $i = 0; $i <= $#codons; $i++ ) {
            for ( my $j = $i + 1; $j <= $#codons; $j++ ) {
                my $codon1 = $codons[$i];
                my $codon2 = $codons[$j];

                my ( $syn, $nsy, $stop ) = (0) x 3;
                if (   $codon_obj->is_ter_codon($codon1)
                    or $codon_obj->is_ter_codon($codon2) )
                {
                    $stop = 1;
                }
                else {
                    ( $syn, $nsy )
                        = $codon_obj->comp_codons( $codon1, $codon2,
                        $snp_codon_pos );
                }
                push @syns,  $syn;
                push @nsys,  $nsy;
                push @stops, $stop;
            }
        }

        $snp_update_sth->execute( $exon_id, $snp_codon_pos,
            join( "|", @codons ),
            average(@syns), average(@nsys), average(@stops), $snp_id );
    }
}
$snp_update_sth->finish;
$snp_sth->finish;

#----------------------------------------------------------#
# isw
#----------------------------------------------------------#
{
    print "Processing isw_syn, nsy, stop\n";

    my $isw_query = q{
        SELECT  i.isw_id id,
                SUM(s.snp_syn) /i.isw_length syn,
                SUM(s.snp_nsy) /i.isw_length nsy,
                SUM(s.snp_stop) /i.isw_length stop
        FROM isw i, snp s
        WHERE i.isw_id = s.isw_id
        AND s.exon_id IS NOT NULL
        GROUP BY i.isw_id
    };
    my $isw_sth = $dbh->prepare($isw_query);

    # update isw table in the new feature column
    my $isw_update = q{
        UPDATE isw
        SET isw_syn = ?,
            isw_nsy = ?,
            isw_stop = ?
        WHERE isw_id = ?
    };
    my $isw_update_sth = $dbh->prepare($isw_update);

    # for isw
    $isw_sth->execute;
    while ( my @row = $isw_sth->fetchrow_array ) {
        my ( $isw_id, $syn, $nsy, $stop ) = @row;
        $isw_update_sth->execute( $syn, $nsy, $stop, $isw_id );
    }
}

#----------------------------------------------------------#
# gene
#----------------------------------------------------------#
{
    print "Processing gene_syn, nsy, stop\n";

    my $gene_query = q{
        SELECT  g.gene_id id,
                g.gene_tl_runlist runlist,
                SUM(s.snp_syn) syn,
                SUM(s.snp_nsy) nsy,
                SUM(s.snp_stop) stop
        FROM gene g 
        inner join exon e on g.gene_id = e.gene_id
        inner join snp s on s.exon_id = e.exon_id
        GROUP BY g.gene_id
    };
    my $gene_sth = $dbh->prepare($gene_query);

    # update gene table in the new feature column
    my $gene_update = q{
        UPDATE gene
        SET gene_syn = ?,
            gene_nsy = ?,
            gene_stop = ?
        WHERE gene_id = ?
    };
    my $gene_update_sth = $dbh->prepare($gene_update);

    # for gene
    $gene_sth->execute;
    while ( my @row = $gene_sth->fetchrow_array ) {
        my ( $gene_id, $runlist, $syn, $nsy, $stop ) = @row;
        my $set    = AlignDB::IntSpan->new($runlist);
        my $length = $set->cardinality;
        $gene_update_sth->execute(
            $syn / $length,
            $nsy / $length,
            $stop / $length,
            $gene_id
        );
    }
}

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
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

    update_snp_dnds.pl - Add additional synonymous/non-synonymous/stop info
                         to alignDB

=head1 SYNOPSIS

    update_snp_dnds.pl [options]
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

B<update_snp_cpg.pl> will Add additional CpG info to alignDB,
1 for CpG and 0 for non.

=cut

