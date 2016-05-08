#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use Bio::EnsEMBL::Registry;

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use AlignDB;
use AlignDB::Position;

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

insert_gene.pl - Add annotation info to alignDB

=head1 SYNOPSIS

    insert_gene.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --ensembl           ensembl database name
        --parallel          run in parallel mode
        --batch             number of alignments process in one child process

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server       = $Config->{database}{server} ),
    'port=i'       => \( my $port         = $Config->{database}{port} ),
    'db|d=s'       => \( my $db           = $Config->{database}{db} ),
    'username|u=s' => \( my $username     = $Config->{database}{username} ),
    'password|p=s' => \( my $password     = $Config->{database}{password} ),
    'ensembl|e=s'  => \( my $ensembl_db   = $Config->{database}{ensembl} ),
    'reg_conf=s'   => \( my $reg_conf     = "$FindBin::RealBin/../ensembl.initrc.pm" ),
    'parallel=i'   => \( my $parallel     = $Config->{generate}{parallel} ),
    'batch=i'      => \( my $batch_number = $Config->{generate}{batch} ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update Gene-related tables of $db...");

# Configure the Bio::EnsEMBL::Registry
Bio::EnsEMBL::Registry->load_all($reg_conf);

my @jobs;
{
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    print "Emptying tables...\n";

    # empty tables: gene, exon, codingsw
    $obj->empty_table( 'gene', 'with_window' );
    $obj->empty_table( 'exon', 'with_window' );

    my @align_ids = @{ $obj->get_align_ids };

    while ( scalar @align_ids ) {
        my @batching = splice @align_ids, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job       = shift;
    my @align_ids = @$job;

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    my DBI $dbh          = $obj->dbh;
    my $pos_obj      = AlignDB::Position->new( dbh => $dbh );
    my $window_maker = $obj->window_maker;

    my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $ensembl_db, 'core', 'Slice' );

    # insert into gene
    my DBI $gene_sth = $dbh->prepare(
        q{
        INSERT INTO gene (
            gene_id, window_id, gene_stable_id,
            gene_external_name, gene_biotype, gene_strand,
            gene_is_full, gene_is_known, gene_multitrans, gene_multiexons,
            gene_tc_runlist, gene_tl_runlist, gene_description
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?, ?
        )
    }
    );

    # insert into exon
    my DBI $exon_sth = $dbh->prepare(
        q{
        INSERT INTO exon (
            exon_id, prev_exon_id, window_id, gene_id, exon_stable_id,
            exon_strand, exon_phase, exon_end_phase,
            exon_frame, exon_is_full,
            exon_tl_runlist,
            exon_seq, exon_peptide)
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?,
            ?,
            ?, ?
        )
    }
    );

    # for each alignment
    for my $align_id (@align_ids) {
        my $target_info    = $obj->get_target_info($align_id);
        my $chr_name       = $target_info->{chr_name};
        my $chr_start      = $target_info->{chr_start};
        my $chr_end        = $target_info->{chr_end};
        my $target_runlist = $target_info->{seq_runlist};

        next if $chr_name =~ /rand|un|contig|hap|scaf/i;

        $obj->process_message($align_id);

        $chr_name =~ s/chr0?//i;

        my $target_set = AlignDB::IntSpan->new($target_runlist);

        # obtain a slice
        my $slice
            = $slice_adaptor->fetch_by_region( 'chromosome', $chr_name, $chr_start, $chr_end );

        my $slice_chr_set = AlignDB::IntSpan->new("$chr_start-$chr_end");

        # insert internal indels, that are, indels in target_set
        # indels in query_set is equal to spans of target_set minus one
        my $internal_indel_flag = 1;

        #----------------------------#
        # INSERT INTO gene
        #----------------------------#
        for my $gene ( @{ $slice->get_all_Genes } ) {
            $gene = $gene->transform('chromosome');

            # gene information
            my %gene_info;
            my @gene_methods = qw{
                start end stable_id external_name biotype strand is_known
                description
            };

            for (@gene_methods) {
                my $result = $gene->$_;
                $result = $result ? $result : undef;
                $gene_info{$_} = $result;
            }

            # convert strand symbol from "1" or "-1" to "+" or "-"
            $gene_info{strand} = $gene_info{strand} > 0 ? "+" : "-";

            # Is this gene a subset of align?
            my $gene_is_full = 1;
            if ( $gene_info{start} < $chr_start ) {
                $gene_info{start} = $chr_start;
                $gene_is_full = 0;
            }
            if ( $gene_info{end} > $chr_end ) {
                $gene_info{end} = $chr_end;
                $gene_is_full = 0;
            }

            # gene position set
            my $gene_start = $pos_obj->at_align( $align_id, $gene_info{start} );
            my $gene_end   = $pos_obj->at_align( $align_id, $gene_info{end} );
            if ( $gene_start >= $gene_end ) {
                print "Gene $gene_info{stable_id} wrong, start >= end\n";
                next;
            }
            my $gene_set = AlignDB::IntSpan->new("$gene_start-$gene_end");
            $gene_set = $gene_set->intersect($target_set);

            # window
            my $cur_gene_window_id
                = $obj->insert_window( $align_id, $gene_set, $internal_indel_flag );

            # Has many transcripts?
            # we use the longest one
            my @transcripts     = @{ $gene->get_all_Transcripts };
            my $gene_multitrans = @transcripts;
            my $transcript;
            if ( $gene_multitrans > 1 ) {
                ($transcript) = sort { $b->length <=> $a->length } @transcripts;
            }
            else {
                ($transcript) = @transcripts;
            }

            my @exons;
            if ($transcript) {
                $transcript = $transcript->transform('chromosome');
                @exons      = @{ $transcript->get_all_Exons };
                @exons      = sort { $a->start <=> $b->start } @exons;
                $_          = $_->transform('chromosome') for @exons;
            }
            my $gene_multiexons = @exons;

            my $gene_tc_set = AlignDB::IntSpan->new;
            my $gene_tl_set = AlignDB::IntSpan->new;

            #----------------------------#
            # PUSH INTO @exon_sites
            #----------------------------#
            my @exon_sites;
            for my $exon (@exons) {

                # exon information
                my %exon_info;

                my @exon_methods = qw{
                    start end stable_id strand phase end_phase frame
                };

                for (@exon_methods) {
                    my $result = $exon->$_;
                    $result = $result ? $result : undef;
                    $exon_info{$_} = $result;
                }

                for (qw{ coding_region_start coding_region_end }) {
                    my $result = $exon->$_($transcript);
                    $result = $result ? $result : undef;
                    $exon_info{$_} = $result;
                }

                # convert strand symbol from "1" or "-1" to "+" or "-"
                $exon_info{strand} = $exon_info{strand} > 0 ? "+" : "-";

                # Is this exon in the align?
                my $i_set = $slice_chr_set->intersect("$exon_info{start}-$exon_info{end}");
                next if $i_set->is_empty;

                # Is this exon a subset of align?
                my $exon_is_full = 1;
                if ( $exon_info{start} < $chr_start ) {
                    $exon_info{start} = $chr_start;
                    $exon_is_full = 0;
                }
                if ( $exon_info{end} > $chr_end ) {
                    $exon_info{end} = $chr_end;
                    $exon_is_full = 0;
                }

                my $exon_seq = $exon->seq->seq;
                my $exon_peptide;
                if ($transcript) {
                    $exon_peptide = $exon->peptide($transcript)->seq;
                }

                # exon position set
                my $exon_start = $pos_obj->at_align( $align_id, $exon_info{start} );
                my $exon_end   = $pos_obj->at_align( $align_id, $exon_info{end} );
                if ( $exon_start >= $exon_end ) {
                    print "Exon $exon_info{stable_id} wrong, start >= end\n";
                    next;
                }
                my $exon_set = AlignDB::IntSpan->new("$exon_start-$exon_end");
                $exon_set = $exon_set->intersect($target_set);

                # coding region set
                my $coding_set = AlignDB::IntSpan->new;
                if (    defined $exon_info{coding_region_start}
                    and defined $exon_info{coding_region_end} )
                {
                    if ( $exon_info{coding_region_start} < $chr_start ) {
                        $exon_info{coding_region_start} = $chr_start;
                        $exon_is_full = 0;
                    }
                    if ( $exon_info{coding_region_end} > $chr_end ) {
                        $exon_info{coding_region_end} = $chr_end;
                        $exon_is_full = 0;
                    }
                    my $coding_region_start
                        = $pos_obj->at_align( $align_id, $exon_info{coding_region_start} );
                    my $coding_region_end
                        = $pos_obj->at_align( $align_id, $exon_info{coding_region_end} );
                    if ( $coding_region_start >= $coding_region_end ) {
                        warn "Exon $exon_info{stable_id} coding_region wrong, start >= end\n";
                        next;
                    }
                    $coding_set->add("$coding_region_start-$coding_region_end");
                    $coding_set = $coding_set->intersect($target_set);
                }

                $exon_info{set}        = $exon_set;
                $exon_info{is_full}    = $exon_is_full;
                $exon_info{tl_runlist} = $coding_set->runlist;
                $exon_info{seq}        = $exon_seq;
                $exon_info{peptide}    = $exon_peptide;

                push @exon_sites, \%exon_info;

                $gene_tc_set->add($exon_set);
                $gene_tl_set->add($coding_set);
            }

            # insert to table
            $gene_sth->execute(
                $cur_gene_window_id,   $gene_info{stable_id}, $gene_info{external_name},
                $gene_info{biotype},   $gene_info{strand},    $gene_is_full,
                $gene_info{is_known},  $gene_multitrans,      $gene_multiexons,
                $gene_tc_set->runlist, $gene_tl_set->runlist, $gene_info{description}
            );
            my $gene_id = $obj->last_insert_id;

            #----------------------------#
            # INSERT INTO exon
            #----------------------------#
            # moving this section and @exon_sites out of gene section
            my $prev_exon_id = 0;
            for my $exon_site (@exon_sites) {

                # window
                my ($cur_exon_window_id)
                    = $obj->insert_window( $align_id, $exon_site->{set}, $internal_indel_flag );

                $exon_sth->execute(
                    $prev_exon_id,            $cur_exon_window_id,  $gene_id,
                    $exon_site->{stable_id},  $exon_site->{strand}, $exon_site->{phase},
                    $exon_site->{end_phase},  $exon_site->{frame},  $exon_site->{is_full},
                    $exon_site->{tl_runlist}, $exon_site->{seq},    $exon_site->{peptide}
                );

                ($prev_exon_id) = $obj->last_insert_id;
            }
        }
    }

    return;
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

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
