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
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Update dnds info of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh();

my $pos_obj = AlignDB::Position->new( dbh => $dbh );
my $codon_obj = AlignDB::Codon->new;

# add columns
{
    $obj->create_column( "snp_extra", "snp_feature5", "INT" );     # exon_id
    $obj->create_column( "snp_extra", "snp_feature6", "INT" );     # codon_pos
    $obj->create_column( "snp_extra", "snp_feature7", "CHAR(3)" ); # codon_t
    $obj->create_column( "snp_extra", "snp_feature8", "CHAR(3)" ); # codon_q
    $obj->create_column( "snp_extra", "snp_feature9", "DOUBLE" );  # syn
    $obj->create_column( "snp_extra", "snp_feature10", "DOUBLE" ); # non-syn
    $obj->create_column( "snp_extra", "snp_feature11", "DOUBLE" ); # stop
    print "Table snp_extra altered\n";

    $obj->create_column( "isw_extra", "isw_feature9",  "DOUBLE" );   # syn
    $obj->create_column( "isw_extra", "isw_feature10", "DOUBLE" );   # non-syn
    $obj->create_column( "isw_extra", "isw_feature11", "DOUBLE" );   # stop
    print "Table isw_extra altered\n";
}

# check rows in snp_extra
{
    my $sql_query = qq{
        SELECT COUNT(snp_extra_id)
        FROM snp_extra
    };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my ($count) = $sth->fetchrow_array;

    unless ($count) {
        $sql_query = qq{
            INSERT INTO snp_extra (snp_id)
            SELECT snp.snp_id
            FROM snp
        };
        $sth = $dbh->prepare($sql_query);
        $sth->execute;
    }
}

# check rows in isw_extra
{
    my $sql_query = qq{
        SELECT COUNT(isw_extra_id)
        FROM isw_extra
    };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute;
    my ($count) = $sth->fetchrow_array;

    unless ($count) {
        $sql_query = qq{
            INSERT INTO isw_extra (isw_id)
            SELECT isw.isw_id
            FROM isw
        };
        $sth = $dbh->prepare($sql_query);
        $sth->execute;
    }
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
# alignments
my $align_query = q{
    SELECT align_id
    FROM align 
};
my $align_sth = $dbh->prepare($align_query);

# alignments' chromosomal location, target_seq and query_seq
my $align_seq_query = q{
    SELECT c.chr_name,
           a.align_length,
           s.chr_start,
           s.chr_end,
           t.target_seq,
           t.target_runlist,
           q.query_seq
    FROM align a, target t, query q, sequence s, chromosome c
    WHERE a.align_id = t.align_id
    AND t.seq_id = s.seq_id
    AND a.align_id = q.align_id
    AND s.chr_id = c.chr_id
    AND a.align_id = ?
};
my $align_seq_sth = $dbh->prepare($align_seq_query);

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
    SELECT  s.snp_id, s.snp_pos, s.target_base, s.query_base
    FROM snp s, snp_extra e
    WHERE s.snp_id = e.snp_id
    AND s.align_id = ?
    AND e.snp_feature1 = 1
};
my $snp_sth = $dbh->prepare($snp_query);

# update snp table in the new feature column
my $snp_extra = q{
    UPDATE snp_extra
    SET snp_feature5 = ?,
        snp_feature6 = ?,
        snp_feature7 = ?,
        snp_feature8 = ?,
        snp_feature9 = ?,
        snp_feature10 = ?,
        snp_feature11 = ?
    WHERE snp_id = ?
};
my $snp_extra_sth = $dbh->prepare($snp_extra);

$align_sth->execute;

# for snp
ALIGN: while ( my @row = $align_sth->fetchrow_array ) {
    my ($align_id) = @row;
    print "Processing align_id $align_id\n";

    $align_seq_sth->execute($align_id);
    my ($chr_name,   $align_length,   $chr_start, $chr_end,
        $target_seq, $target_runlist, $query_seq
    ) = $align_seq_sth->fetchrow_array;
    next ALIGN if $chr_name =~ /rand|un|contig|hap|scaf/i;

    print "Prosess align $align_id in $chr_name $chr_start - $chr_end\n";

    $chr_name =~ s/chr0?//i;

    # target runlist
    my $target_set = AlignDB::IntSpan->new($target_runlist);

    # get all exons
    my @exons;
    $exon_sth->execute($align_id);
EXON: while ( my @row = $exon_sth->fetchrow_array ) {
        my ( $exon_id, $exon_strand, $exon_tl_runlist, $exon_seq,
            $exon_peptide )
            = @row;

        # extract exon seq in this align
        my $exon_tl_set    = AlignDB::IntSpan->new($exon_tl_runlist);
        my $exon_align_seq = $exon_tl_set->substr_span($target_seq);
        if ( $exon_strand eq '-' ) {
            $exon_align_seq = revcom($exon_align_seq);
        }

        if ( index( $exon_seq, $exon_align_seq ) == -1 ) {
            print " " x 4, "Exon sequence does not match alignment\n";

            #print Dump {
            #    exon_seq       => $exon_seq,
            #    exon_align_seq => $exon_align_seq
            #};
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
        my ( $snp_id, $snp_pos, $target_base, $query_base ) = @row;

        print "Processing SNP $snp_id, POS $snp_pos\n";

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
            my $codon_end = $exon_tl_set->at( $codon_start_index + 2 );
            if ( !defined $codon_start or !defined $codon_end ) {
                print " " x 4, "Can't calc codon positions\n";
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
            my $codon_start = $exon_tl_set->at($codon_start_index);
            my $codon_end = $exon_tl_set->at( $codon_start_index + 2 );
            if ( !defined $codon_start or !defined $codon_end ) {
                print " " x 4, "Can't calc codon positions\n";
                next SNP;
            }
            $codon_set = AlignDB::IntSpan->new("$codon_start-$codon_end");
        }

        #print Dump {
        #    exon_strand   => $exon_strand,
        #    exon_frame   => $exon_frame,
        #    exon_tl_set   => "$exon_tl_set",
        #    exon_length   => $exon_tl_set->cardinality,
        #    snp_pos       => $snp_pos,
        #    snp_codon_pos => $snp_codon_pos,
        #    codon_set     => "$codon_set",
        #};

        # target codon and query codon
        my $codon_t = $codon_set->substr_span($target_seq);
        my $codon_q = $codon_set->substr_span($query_seq);
        if ( $exon_strand eq '-' ) {
            $codon_t = revcom($codon_t);
            $codon_q = revcom($codon_q);
        }

        my ( $syn, $nsy, $stop ) = (0) x 3;
        if (   $codon_obj->is_ter_codon($codon_t)
            or $codon_obj->is_ter_codon($codon_q) )
        {
            $stop = 1;
        }
        else {
            ( $syn, $nsy )
                = $codon_obj->comp_codons( $codon_t, $codon_q,
                $snp_codon_pos );
        }

        #print Dump {
        #    codon_t => $codon_t,
        #    codon_q => $codon_q,
        #    base_t  => $target_base,
        #    base_q  => $query_base,
        #    syn     => $syn,
        #    nsy     => $nsy,
        #    stop    => $stop,
        #};

        $snp_extra_sth->execute( $exon_id, $snp_codon_pos, $codon_t, $codon_q,
            $syn, $nsy, $stop, $snp_id );
    }

}
$snp_extra_sth->finish;
$snp_sth->finish;

$align_sth->finish;

#----------------------------------------------------------#
# isw_extra
#----------------------------------------------------------#
{
    print "Processing isw_feature9, 10, 11\n";

    my $isw_query = q{
        SELECT  i.isw_id id,
                SUM(e.snp_feature9) /i.isw_length syn,
                SUM(e.snp_feature10) /i.isw_length nsy,
                SUM(e.snp_feature11) /i.isw_length stop
        FROM isw i, snp s, snp_extra e
        WHERE i.isw_id = s.isw_id
        AND s.snp_id = e.snp_id
        AND e.snp_feautre1 = 1
        AND e.snp_feature5 IS NOT NULL
        GROUP BY i.isw_id
    };
    my $isw_sth = $dbh->prepare($isw_query);

    # update isw table in the new feature column
    my $isw_extra = q{
        UPDATE isw_extra
        SET isw_feature9 = ?,
            isw_feature10 = ?,
            isw_feature11 = ?
        WHERE isw_id = ?
    };
    my $isw_extra_sth = $dbh->prepare($isw_extra);

    # for isw
    $isw_sth->execute;
    while ( my @row = $isw_sth->fetchrow_array ) {
        my ( $isw_id, $syn, $nsy, $stop ) = @row;
        $isw_extra_sth->execute( $syn, $nsy, $stop, $isw_id );
    }

    ## update NULL value of isw_feature3 to 0
    #my $isw_null = q{
    #    UPDATE isw_extra
    #    SET isw_feature9 = 0,
    #        isw_feature10 = 0,
    #        isw_feature11 = 0
    #    WHERE isw_feature9 IS NULL
    #};
    #my $isw_null_sth = $dbh->prepare($isw_null);
    #$isw_null_sth->execute();

}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    update_snp_dnds.pl - Add additional synonymous/non-synonymous info
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

