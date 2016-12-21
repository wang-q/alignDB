#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use Bio::EnsEMBL::Registry;
use MCE;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use App::Fasops::Common;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use AlignDB::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
Update Gene-related tables to alignDB

    perl init/insert_gene.pl -d S288cvsRM11_1a --parallel 2

Usage: perl %c [options]
EOF

(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    $description,
    [ 'help|h', 'display this message' ],
    [],
    ['Database init values'],
    [ 'server|s=s',   'MySQL IP/Domain', { default => $conf->{database}{server} }, ],
    [ 'port=i',       'MySQL port',      { default => $conf->{database}{port} }, ],
    [ 'username|u=s', 'username',        { default => $conf->{database}{username} }, ],
    [ 'password|p=s', 'password',        { default => $conf->{database}{password} }, ],
    [ 'db|d=s',       'database name',   { default => $conf->{database}{db} }, ],
    [],
    [ 'ensembl|e=s', 'ensembl database name', { default => $conf->{database}{ensembl} }, ],
    [ 'reg_conf=s', 'ensembl.initrc.pm', { default => "$FindBin::RealBin/../ensembl.initrc.pm" }, ],
    [ 'parallel=i', 'run in parallel mode',       { default => $conf->{generate}{parallel} }, ],
    [ 'batch=i',    '#alignments in one process', { default => $conf->{generate}{batch} }, ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# record config
$stopwatch->record_conf($opt);

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update Gene-related tables of [$opt->{db}]...");

# Configure the Bio::EnsEMBL::Registry
Bio::EnsEMBL::Registry->load_all( $opt->{reg_conf} );

my @jobs;
{
    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    # empty tables: gene, exon
    print "Emptying tables...\n";
    $alignDB->empty_table( 'gene', 'with_window' );
    $alignDB->empty_table( 'exon', 'with_window' );

    @jobs = @{ $alignDB->get_align_ids };
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my @align_ids = @{$chunk_ref};
    my $wid       = MCE->wid;

    $stopwatch->block_message("Process task [$chunk_id] by worker #$wid");

    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    my DBI $dbh = $alignDB->dbh;

    my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $opt->{ensembl}, 'core', 'Slice' );

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
        my $target_info    = $alignDB->get_target_info($align_id);
        my $chr_name       = $target_info->{chr_name};
        my $chr_start      = $target_info->{chr_start};
        my $chr_end        = $target_info->{chr_end};
        my $chr_strand     = $target_info->{chr_strand};
        my $target_intspan = $target_info->{seq_intspan};

        next if $chr_name =~ /rand|un|contig|hap|scaf/i;

        $alignDB->process_message($align_id);

        $chr_name =~ s/chr0?//i;

        # obtain a slice
        my $slice
            = $slice_adaptor->fetch_by_region( 'chromosome', $chr_name, $chr_start, $chr_end );

        my $slice_chr_set = AlignDB::IntSpan->new()->add_pair( $chr_start, $chr_end );

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
            my $gene_start = App::Fasops::Common::chr_to_align( $target_intspan, $gene_info{start},
                $chr_start, $chr_strand, );
            my $gene_end = App::Fasops::Common::chr_to_align( $target_intspan, $gene_info{end},
                $chr_start, $chr_strand, );
            if ( $gene_start >= $gene_end ) {
                print "Gene $gene_info{stable_id} wrong, start >= end\n";
                next;
            }
            my $gene_set = AlignDB::IntSpan->new()->add_pair( $gene_start, $gene_end );
            $gene_set = $gene_set->intersect($target_intspan);

            # window
            my $cur_gene_window_id
                = $alignDB->insert_window( $align_id, $gene_set, $internal_indel_flag );

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
                my $exon_start
                    = App::Fasops::Common::chr_to_align( $target_intspan, $exon_info{start},
                    $chr_start, $chr_strand, );
                my $exon_end = App::Fasops::Common::chr_to_align( $target_intspan, $exon_info{end},
                    $chr_start, $chr_strand, );
                if ( $exon_start >= $exon_end ) {
                    print "Exon $exon_info{stable_id} wrong, start >= end\n";
                    next;
                }
                my $exon_set = AlignDB::IntSpan->new("$exon_start-$exon_end");
                $exon_set = $exon_set->intersect($target_intspan);

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
                        = App::Fasops::Common::chr_to_align( $target_intspan,
                        $exon_info{coding_region_start},
                        $chr_start, $chr_strand, );
                    my $coding_region_end
                        = App::Fasops::Common::chr_to_align( $target_intspan,
                        $exon_info{coding_region_end},
                        $chr_start, $chr_strand, );
                    if ( $coding_region_start >= $coding_region_end ) {
                        warn "Exon $exon_info{stable_id} coding_region wrong, start >= end\n";
                        next;
                    }
                    $coding_set->add("$coding_region_start-$coding_region_end");
                    $coding_set = $coding_set->intersect($target_intspan);
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
            my $gene_id = $alignDB->last_insert_id;

            #----------------------------#
            # INSERT INTO exon
            #----------------------------#
            # moving this section and @exon_sites out of gene section
            my $prev_exon_id = 0;
            for my $exon_site (@exon_sites) {

                # window
                my ($cur_exon_window_id)
                    = $alignDB->insert_window( $align_id, $exon_site->{set}, $internal_indel_flag );

                $exon_sth->execute(
                    $prev_exon_id,            $cur_exon_window_id,  $gene_id,
                    $exon_site->{stable_id},  $exon_site->{strand}, $exon_site->{phase},
                    $exon_site->{end_phase},  $exon_site->{frame},  $exon_site->{is_full},
                    $exon_site->{tl_runlist}, $exon_site->{seq},    $exon_site->{peptide}
                );

                ($prev_exon_id) = $alignDB->last_insert_id;
            }
        }
    }

    return;
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $opt->{parallel}, chunk_size => $opt->{batch}, );
$mce->forchunk( \@jobs, $worker, );

$stopwatch->end_message;

# store program's meta info to database
AlignDB::Common->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

exit;

__END__
