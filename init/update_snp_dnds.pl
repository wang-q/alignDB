#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck;

use AlignDB::Codon;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use App::Fasops::Common;

use FindBin;
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

update_snp_dnds.pl - Add synonymous/non-synonymous/stop info

=head1 SYNOPSIS

    update_snp_dnds.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db_name  = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update dnds info of $db_name...");

my $obj = AlignDB->new(
    mysql  => "$db_name:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my DBI $dbh = $obj->dbh;

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
my DBI $exon_sth = $dbh->prepare(
    q{
    SELECT  e.exon_id, e.exon_strand,
            e.exon_tl_runlist, e.exon_seq, e.exon_peptide
    FROM exon e, window w
    WHERE e.window_id = w.window_id
    AND w.align_id = ?
    AND e.exon_tl_runlist != '-'
    }
);

# select all coding snps in this alignment
my DBI $snp_sth = $dbh->prepare(
    q{
    SELECT  s.snp_id, s.snp_pos
    FROM snp s
    WHERE s.align_id = ?
    AND s.snp_coding = 1
    }
);

# update snp table in the new feature column
my DBI $snp_update_sth = $dbh->prepare(
    q{
    UPDATE snp
    SET exon_id = ?,
        snp_codon_pos = ?,
        snp_codons = ?,
        snp_syn = ?,
        snp_nsy = ?,
        snp_stop = ?
    WHERE snp_id = ?
    }
);

my @align_ids = @{ $obj->get_align_ids };
ALIGN: for my $align_id (@align_ids) {

    my $target_info = $obj->get_target_info($align_id);
    my $chr_name    = $target_info->{chr_name};

    next ALIGN if $chr_name =~ /rand|un|contig|hap|scaf/i;

    $obj->process_message($align_id);
    my ( $target_seq, @query_seqs ) = @{ $obj->get_seqs($align_id) };
    $_ = uc $_ for ( $target_seq, @query_seqs );

    # get all exons
    my @exons;
    $exon_sth->execute($align_id);
EXON: while ( my @row = $exon_sth->fetchrow_array ) {
        my ( $exon_id, $exon_strand, $exon_tl_runlist, $exon_seq, $exon_peptide, ) = @row;

        $exon_seq = uc $exon_seq;

        # extract exon seq in this align
        my $exon_tl_set    = AlignDB::IntSpan->new($exon_tl_runlist);
        my $exon_align_seq = $exon_tl_set->substr_span($target_seq);
        if ( $exon_strand eq '-' ) {
            $exon_align_seq = App::Fasops::Common::revcom($exon_align_seq);
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
            my $exon_align_peptide = $codon_obj->translate( $exon_align_seq, $_ );
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
            if ( $_->{exon_tl_set}->contains($snp_pos) ) {
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
            $codon_t = App::Fasops::Common::revcom($codon_t);
        }
        my @codons = ($codon_t);
        for my $query_seq (@query_seqs) {
            my $codon_q = $codon_set->substr_span($query_seq);
            if ( $exon_strand eq '-' ) {
                $codon_q = App::Fasops::Common::revcom($codon_q);
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
                    ( $syn, $nsy ) = $codon_obj->comp_codons( $codon1, $codon2, $snp_codon_pos );
                }
                push @syns,  $syn;
                push @nsys,  $nsy;
                push @stops, $stop;
            }
        }

        $snp_update_sth->execute(
            $exon_id, $snp_codon_pos,
            join( "|", @codons ),
            App::Fasops::Common::mean(@syns),
            App::Fasops::Common::mean(@nsys),
            App::Fasops::Common::mean(@stops), $snp_id
        );
    }
}
$snp_update_sth->finish;
$snp_sth->finish;

#----------------------------------------------------------#
# isw
#----------------------------------------------------------#
{
    print "Processing isw_syn, nsy, stop\n";

    my DBI $isw_sth = $dbh->prepare(
        q{
        SELECT  i.isw_id id,
                SUM(s.snp_syn) /i.isw_length syn,
                SUM(s.snp_nsy) /i.isw_length nsy,
                SUM(s.snp_stop) /i.isw_length stop
        FROM isw i, snp s
        WHERE i.isw_id = s.isw_id
        AND s.exon_id IS NOT NULL
        GROUP BY i.isw_id
        }
    );

    # update isw table in the new feature column
    my DBI $isw_update_sth = $dbh->prepare(
        q{
        UPDATE isw
        SET isw_syn = ?,
            isw_nsy = ?,
            isw_stop = ?
        WHERE isw_id = ?
        }
    );

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

    my DBI $gene_sth = $dbh->prepare(
        q{
        SELECT  g.gene_id id,
                g.gene_tl_runlist runlist,
                SUM(s.snp_syn) syn,
                SUM(s.snp_nsy) nsy,
                SUM(s.snp_stop) stop
        FROM gene g
        inner join exon e on g.gene_id = e.gene_id
        inner join snp s on s.exon_id = e.exon_id
        GROUP BY g.gene_id
        }
    );

    # update gene table in the new feature column
    my DBI $gene_update_sth = $dbh->prepare(
        q{
        UPDATE gene
        SET gene_syn = ?,
            gene_nsy = ?,
            gene_stop = ?
        WHERE gene_id = ?
        }
    );

    # for gene
    $gene_sth->execute;
    while ( my @row = $gene_sth->fetchrow_array ) {
        my ( $gene_id, $runlist, $syn, $nsy, $stop ) = @row;
        my $set    = AlignDB::IntSpan->new($runlist);
        my $length = $set->cardinality;
        $gene_update_sth->execute( $syn / $length, $nsy / $length, $stop / $length, $gene_id );
    }
}

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB->new(
        mysql  => "$db_name:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}
exit;

__END__
