#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Text::CSV_XS;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use DateTime::Format::Natural;

use FindBin;

use AlignDB::Stopwatch;

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

# running options
my $base_dir    = $Config->{bac}{base_dir};
my $taxon_dir   = $Config->{bac}{taxon_dir};
my $strain_file = $Config->{bac}{strain_file};
my $seq_file    = $Config->{bac}{seq_file};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    'b|base_dir=s'  => \$base_dir,
    'x|taxon_dir=s' => \$taxon_dir,
    'strain_file=s' => \$strain_file,
    'seq_file=s'    => \$seq_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing bac strains summary...");

my $dbh = DBI->connect("DBI:CSV:");

my $taxon_db = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -directory => "$taxon_dir",
    -nodesfile => "$taxon_dir/nodes.dmp",
    -namesfile => "$taxon_dir/names.dmp",
);

#----------------------------#
# load tab sep. txt files
#----------------------------#
$dbh->{csv_tables}->{t0} = {
    eol       => "\n",
    sep_char  => "\t",
    file      => "$base_dir/lproks_0.txt",
    col_names => [
        "refseq_project_id", "project_id",
        "taxonomy_id",       "organism_name",
        "super_kingdom",     "group",
        "sequence_status",   "genome_size",
        "gc_content",        "gram_stain",
        "shape",             "arrangment",
        "endospores",        "motility",
        "salinity",          "oxygen_req",
        "habitat",           "temp_range",
        "optimal_temp",      "pathogenic_in",
        "disease",           "genbank_accessions",
        "refseq_accessions",
    ],
};
$dbh->{csv_tables}->{t1} = {
    eol       => "\n",
    sep_char  => "\t",
    file      => "$base_dir/lproks_1.txt",
    col_names => [
        "refseq_project_id",     "project_id",
        "taxonomy_id",           "organism_name",
        "super_kingdom",         "group",
        "genome_size",           "gc_content",
        "number_of_chromosomes", "number_of_plasmids",
        "released_date",         "modified_date",
        "genbank_accessions",    "refseq_accessions",
        "publications",          "list_of_center",
    ],
};
$dbh->{csv_tables}->{t2} = {
    eol       => "\n",
    sep_char  => "\t",
    file      => "$base_dir/summary.txt",
    col_names => [
        "accession",  "genbankacc",    "length",   "taxonomy_id",
        "project_id", "organism_name", "replicon", "create_date",
        "update_date",
    ],
};

#----------------------------#
# join t0 and t1
#----------------------------#
{
    my $join_sth = $dbh->prepare(
        qq{
        SELECT t0.taxonomy_id,
               t0.organism_name,
               t0.super_kingdom,
               t0.group,
               t0.genome_size,
               t0.gc_content,
               t1.number_of_chromosomes,
               t1.number_of_plasmids,
               t1.refseq_accessions,
               t0.gram_stain,
               t0.shape,
               t0.arrangment,
               t0.endospores,
               t0.motility,
               t0.salinity,
               t0.oxygen_req,
               t0.habitat,
               t0.temp_range,
               t0.optimal_temp,
               t0.pathogenic_in,
               t0.disease,
               t1.released_date,
               t1.modified_date
        FROM   t0, t1
        WHERE  t0.taxonomy_id = t1.taxonomy_id
        }
    );

    # prepare output csv file
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, ">", $strain_file or die "$strain_file: $!";
    my @headers = qw{
        taxonomy_id
        organism_name
        super_kingdom
        group
        genome_size
        gc_content
        number_of_chromosomes
        number_of_plasmids
        refseq_accessions
        gram_stain
        shape
        arrangment
        endospores
        motility
        salinity
        oxygen_req
        habitat
        temp_range
        optimal_temp
        pathogenic_in
        disease
        released_date
        modified_date
        species
        species_id
        genus
        genus_id
        species_member
        genus_species_member
        genus_strain_member
    };
    $csv->print( $csv_fh, \@headers );

    $join_sth->execute;
    while ( my @row = $join_sth->fetchrow_array ) {

        # parse dates
        my $parser = DateTime::Format::Natural->new( format => 'mm/dd/yy' );
        my $released_date = $parser->parse_datetime( $row[-2] )->ymd;
        my $modified_date = $parser->parse_datetime( $row[-1] )->ymd;
        $row[-2] = $released_date;
        $row[-1] = $modified_date;

        # find each strains' species and genus
        my $taxon_id = $row[0];
        my $name     = $row[1];
        my $bac      = $taxon_db->get_taxon( -taxonid => $taxon_id );
        if ( !$bac ) {
            warn "Can't find taxon for $name\n";
            next;
        }

        my $species = find_ancestor( $bac, 'species' );
        if ($species) {
            push @row, ( $species->scientific_name, $species->id );
        }
        else {
            push @row, ( undef, undef );
            warn "Can't find species for $name\n";
        }

        my $genus = find_ancestor( $bac, 'genus' );
        if ($genus) {
            push @row, ( $genus->scientific_name, $genus->id );
        }
        else {
            push @row, ( undef, undef );
            warn "Can't find genus for $name\n";
        }

        push @row, ( undef, undef, undef );    # member numbers

        # write a line
        $csv->print( $csv_fh, \@row );
    }
    $join_sth->finish;
    close $csv_fh;
}

#----------------------------#
# transform t2
#----------------------------#
{
    my $sth = $dbh->prepare(
        qq{
        SELECT  accession,
                genbankacc,
                taxonomy_id,
                organism_name,
                replicon,
                length
        FROM    t2
        }
    );

    # prepare output csv file
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, ">", $seq_file or die "$seq_file: $!";
    my @headers = qw{
        accession
        accession_version
        genbankacc
        taxonomy_id
        organism_name
        replicon
        length
    };
    $csv->print( $csv_fh, \@headers );

    $sth->execute;
    $sth->fetchrow_array;    # omit first line
    while ( my @row = $sth->fetchrow_array ) {
        my $accession = $row[0];
        $accession =~ s/\.\d+//;
        unshift @row, $accession;
        $row[-1] =~ s/\s+$//;
        $csv->print( $csv_fh, \@row );
    }
    $sth->finish;
    close $csv_fh;
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub find_ancestor {
    my $taxon = shift;
    my $rank = shift || 'species';

    return $taxon if $taxon->rank eq $rank;

RANK: while (1) {
        $taxon = $taxon->ancestor;
        last RANK unless defined $taxon;
        return $taxon if $taxon->rank eq $rank;
    }

    return;
}

__END__

perl bac_strains.pl 
