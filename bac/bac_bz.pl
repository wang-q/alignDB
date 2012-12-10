#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use File::Find::Rule;
use File::Spec;
use File::Copy;
use File::Basename;
use Text::Table;
use List::MoreUtils qw(uniq);
use Archive::Extract;

use Bio::Taxon;
use Bio::DB::Taxonomy;
use Template;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;

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
my $base_dir  = $Config->{bac}{base_dir};
my $taxon_dir = $Config->{bac}{taxon_dir};

my $working_dir = ".";

my $parent_id = "562,585054";    # E.coli and E. fergusonii
my $target_id;
my $outgroup_id;
my $exclude_ids = '0';

# use custom name_str
# working dir and goal db name
# mysql restrict db name length 64
my $name_str;

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{bac}{db};

my $gr;
my $scaffold;
my $td_dir  = $Config->{bac}{td_dir};
my $nb_dir  = $Config->{bac}{nb_dir};
my $nbd_dir = $Config->{bac}{nbd_dir};

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'          => \$help,
    'man'             => \$man,
    'server=s'        => \$server,
    'port=i'          => \$port,
    'db=s'            => \$db,
    'username=s'      => \$username,
    'password=s'      => \$password,
    'b|base_dir=s'    => \$base_dir,
    'x|taxon_dir=s'   => \$taxon_dir,
    'w|working_dir=s' => \$working_dir,
    'p|parent_id=s'   => \$parent_id,
    't|target_id=i'   => \$target_id,
    'o|r=i'           => \$outgroup_id,
    'e|exclude=s'     => \$exclude_ids,
    'n|name_str=s'    => \$name_str,
    'gr'              => \$gr,
    'scaffold'        => \$scaffold,
    'parallel=i'      => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Preparing Blastz whole species");

my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password );

my $id_str;
{    # expand $parent_id
    $taxon_dir = $td_dir if $gr;

    my $taxon_db = Bio::DB::Taxonomy->new(
        -source    => 'flatfile',
        -directory => "$taxon_dir",
        -nodesfile => "$taxon_dir/nodes.dmp",
        -namesfile => "$taxon_dir/names.dmp",
    );
    my @parent_ids = split /,/, $parent_id;

    my $sub_id_set = AlignDB::IntSpan->new;
    for my $p_id (@parent_ids) {
        $sub_id_set->add($p_id);
        my $parent = $taxon_db->get_taxon( -taxonid => $p_id );

        my @taxa = $taxon_db->get_all_Descendents($parent);
        for my $taxon (@taxa) {
            $sub_id_set->add( $taxon->id );
        }
    }

    my $db_id_set = AlignDB::IntSpan->new;
    {
        my $query
            = $gr
            ? $scaffold
                ? q{ SELECT taxonomy_id FROM gr WHERE 1 = 1 }
                : q{ SELECT taxonomy_id FROM gr WHERE status = 'Complete' }
            : q{ SELECT taxonomy_id FROM strain WHERE seq_ok = 1 };
        my $sth = $dbh->prepare($query);
        $sth->execute;
        while ( my ($id) = $sth->fetchrow_array ) {
            $db_id_set->add($id);
        }
    }

    my $id_set = $sub_id_set->intersect($db_id_set);
    $id_set->remove( split /,/, $exclude_ids );
    $id_str = '(' . ( join ",", $id_set->as_array ) . ')';

    die "Wrong id_str $id_str\n" unless $id_str =~ /\d+/;
}

{    # making working dir
    if ( !$name_str ) {
        my $query
            = $gr
            ? qq{ SELECT DISTINCT species FROM gr WHERE taxonomy_id IN $id_str }
            : qq{ SELECT DISTINCT species FROM strain st WHERE st.taxonomy_id IN $id_str };
        my $sth = $dbh->prepare($query);
        $sth->execute;

        while ( my ($name) = $sth->fetchrow_array ) {
            $name_str .= "_$name";
        }
        $name_str =~ s/\W/_/g;
        $name_str =~ s/^_+//g;
        $name_str =~ s/\s+/_/g;
    }

    print "Working on $name_str\n";
    $working_dir = File::Spec->catdir( $working_dir, $name_str );
    $working_dir = File::Spec->rel2abs($working_dir);
    mkdir $working_dir unless -e $working_dir;
    print "Working dir is $working_dir\n";
}

my @query_ids;
{    # find all strains' taxon ids

    # select all strains in this species
    my $query
        = $gr
        ? qq{ SELECT taxonomy_id, organism_name, released_date, status, code FROM gr WHERE taxonomy_id IN $id_str ORDER BY released_date, status, code }
        : qq{ SELECT taxonomy_id, organism_name, released_date FROM strain WHERE taxonomy_id IN $id_str ORDER BY released_date };
    my $sth = $dbh->prepare($query);
    $sth->execute;

    # header line
    my @strains;
    my $table = Text::Table->new( @{ $sth->{NAME} } );
    while ( my @row = $sth->fetchrow_array ) {
        push @strains, [@row];
        $table->load( [@row] );
    }

    my $table_file = File::Spec->catfile( $working_dir, "table.txt" );
    open my $fh, '>', $table_file;
    print {$fh} $table, "\n";
    print $table, "\n";

    {
        my $message = "There are " . scalar @strains . " strains\n";
        print {$fh} $message;
        print $message;
    }

    if ($target_id) {
        my ($exist) = grep { $_->[0] == $target_id } @strains;
        if ( defined $exist ) {
            my $message = "Use [$exist->[1]] as target, as you wish.\n";
            print {$fh} $message;
            print $message;
        }
        else {
            print "Taxon $target_id doesn't exist, please check.\n";
            exit;
        }
    }
    else {
        $target_id = $strains[0]->[0];
        my $message
            = "Use [$strains[0]->[1]] as target, the oldest strain on NCBI.\n";
        print {$fh} $message;
        print $message;
    }

    @query_ids = map { $_->[0] == $target_id ? () : $_->[0] } @strains;

    if ($outgroup_id) {
        my ($exist) = grep { $_ == $outgroup_id } @query_ids;
        if ( defined $exist ) {
            my $message = "Use [$exist] as reference, as you wish.\n";
            print {$fh} $message;
            print $message;

            @query_ids = map { $_ == $outgroup_id ? () : $_ } @query_ids;
            unshift @query_ids, $outgroup_id;
        }
        else {
            print "Taxon $outgroup_id doesn't exist, please check.\n";
        }
    }

    print "\n";
    print {$fh} "perl " . $stopwatch->cmd_line, "\n";

    close $fh;
}

my @new_gff_files;
{    # build fasta files

    $base_dir = $nb_dir if $gr;

    # read all filenames, then grep
    print "Reading file list\n";
    my @fna_files = File::Find::Rule->file->name('*.fna')->in($base_dir);
    my @gff_files = File::Find::Rule->file->name('*.gff')->in($base_dir);
    my @wgs_files;
    if ( $gr and $scaffold ) {
        @wgs_files
            = File::Find::Rule->file->name('*.scaffold.fna.tgz')->in($nbd_dir);
    }

    print "Rewrite seqs for every strains\n";
    for my $taxon_id ( $target_id, @query_ids ) {
        print "taxon_id $taxon_id\n";
        my $id_dir = File::Spec->catdir( $working_dir, $taxon_id );
        mkdir $id_dir unless -e $id_dir;

        my @accs;    # complete accessions

        if ( !$gr ) {
            my $query
                = qq{ SELECT accession FROM seq WHERE taxonomy_id = ? AND replicon like "%chr%" };
            my $sth = $dbh->prepare($query);
            $sth->execute($taxon_id);
            while ( my ($acc) = $sth->fetchrow_array ) {
                push @accs, $acc;
            }
        }
        else {
            my $query = qq{ SELECT chr_refseq FROM gr WHERE taxonomy_id = ? };
            my $sth   = $dbh->prepare($query);
            $sth->execute($taxon_id);
            my ($acc) = $sth->fetchrow_array;
            push @accs,
                ( map { s/\.\d+$//; $_ } grep {defined} ( split /,/, $acc ) );
        }

        for my $acc ( grep {defined} @accs ) {
            my $rc = prep_fa( \@fna_files, $acc, $id_dir );
            if ($rc) {
                warn $rc;
                next;
            }

            if ( $taxon_id eq $target_id ) {
                print "Copy target gff\n";
                my ($gff_file) = grep {/$acc/} @gff_files;
                my $new_gff_file
                    = File::Spec->catfile( $working_dir, "$acc.gff" );
                push @new_gff_files, $new_gff_file;
                copy( $gff_file, $new_gff_file );
            }
        }

        if ($scaffold) {
            my ($wgs) = get_taxon_wgs( $dbh, $taxon_id );

            next unless $wgs;

            $wgs =~ s/\d+$//;
            my $rc = prep_wgs( \@wgs_files, $wgs, $id_dir );
            if ($rc) {
                warn $rc;
            }
        }
    }
}

sub prep_fa {
    my $all_files = shift;
    my $acc       = shift;
    my $dir       = shift;

    my ($fna_file) = grep {/$acc/} @{$all_files};
    if ( !$fna_file ) {
        return "Can't find fasta file for $acc\n";
    }

    my $fa_file = File::Spec->catfile( $dir, "$acc.fa" );
    open my $in_fh,  '<', $fna_file;
    open my $out_fh, '>', $fa_file;
    while (<$in_fh>) {
        if (/>/) {
            print {$out_fh} ">$acc\n";
        }
        else {
            print {$out_fh} $_;
        }
    }
    close $out_fh;
    close $in_fh;

    return;
}

sub get_taxon_wgs {
    my $dbh      = shift;
    my $taxon_id = shift;

    my $query = qq{ SELECT wgs FROM gr WHERE taxonomy_id = ? };
    my $sth   = $dbh->prepare($query);
    $sth->execute($taxon_id);
    my ($wgs) = $sth->fetchrow_array;

    return $wgs;
}

sub prep_wgs {
    my $all_files = shift;
    my $wgs       = shift;
    my $dir       = shift;

    my ($wgs_file) = grep {/NZ_$wgs/} @{$all_files};
    if ( !$wgs_file ) {
        return "Can't find fasta file for $wgs\n";
    }

    my $ae = Archive::Extract->new( archive => $wgs_file );
    my $ok = $ae->extract( to => $dir );

    if ( !$ok ) {
        return $ae->error;
    }

    my (@files) = map { File::Spec->rel2abs( $_, $dir ) } @{ $ae->files };

    for my $file (@files) {
        unless ( -e $file ) {
            return "$file not exists!\n";
        }
        if ( ( stat($file) )[7] < 1024 ) {
            next;
        }
        open my $in_fh, '<', $file;

        my $basename = basename( $file, ".fna" );
        my $fa_file = File::Spec->catfile( $dir, "$basename.fa" );
        open my $out_fh, '>', $fa_file;
        while (<$in_fh>) {
            if (/>/) {
                print {$out_fh} ">$basename\n";
            }
            else {
                print {$out_fh} $_;
            }
        }
        close $out_fh;
        close $in_fh;
    }

    unlink $_ for @files;

    return;
}

sub taxon_info {
    my $dbh      = shift;
    my $taxon_id = shift;
    my $dir      = shift;
    my $gr       = shift;

    my $table = $gr ? 'gr' : 'strain';
    my $query
        = qq{ SELECT taxonomy_id, organism_name, genus, species FROM $table WHERE taxonomy_id = ? };
    my $sth = $dbh->prepare($query);
    $sth->execute($taxon_id);
    my ( $taxonomy_id, $organism_name, $genus, $species )
        = $sth->fetchrow_array;
    $species =~ s/^$genus\s+//;
    my $sub_name = $organism_name;
    $sub_name      =~ s/^$genus\s+//;
    $sub_name      =~ s/^$species\s+//;
    $organism_name =~ s/\W/_/g;
    $organism_name =~ s/_+/_/g;

    return {
        taxon   => $taxonomy_id,
        name    => $organism_name,
        genus   => $genus,
        species => $species,
        subname => $sub_name,
        dir     => File::Spec->catdir( $working_dir, $taxon_id ),
    };
}

{
    my $seq_pair_file = File::Spec->catfile( $working_dir, "seq_pair.csv" );
    {    # write seq_pair.csv and left seq_pair_batch.pl to handle other things
        open my $fh, '>', $seq_pair_file;
        for my $query_id (@query_ids) {
            print {$fh} File::Spec->catdir( $working_dir, $target_id ), ",",
                File::Spec->catdir( $working_dir, $query_id ), "\n";
        }
        close $fh;
    }

    my $tt = Template->new;
    my $text;
    my @data = map { taxon_info( $dbh, $_, $working_dir, $gr ) }
        ( $target_id, @query_ids );

    # taxon.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.genus %],[% item.species %],[% item.subname %],[% item.name %],
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $working_dir, "taxon.csv" )
    ) or die Template->error;

    # chr_length.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $working_dir, "chr_length_chrUn.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% working_dir %]

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
faSize -detailed [% item.dir%]/*.fa > [% item.dir%]/chr.sizes
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]}' [% item.dir %]/chr.sizes >> real_chr.csv
[% END -%]

cat chr_length_chrUn.csv real_chr.csv > chr_length.csv
rm real_chr.csv

echo '# Run the following cmds to merge csv files'
echo
echo perl [% findbin %]/../util/merge_csv.pl -t [% findbin %]/../init/taxon.csv -m [% working_dir %]/taxon.csv
echo
echo perl [% findbin %]/../util/merge_csv.pl -t [% findbin %]/../init/chr_length.csv -m [% working_dir %]/chr_length.csv
echo

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            working_dir => $working_dir,
            findbin     => $FindBin::Bin,
        },
        File::Spec->catfile( $working_dir, "real_chr.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
# bac_bz.pl
# perl [% stopwatch.cmd_line %]

cd [% working_dir %]

# seq_pair_batch.pl blastz
perl [% findbin %]/../extra/seq_pair_batch.pl \
    -d 1 --parallel [% parallel %] \
    -f [% seq_pair_file %] \
    -at 1000 -st 0  -r 0

# seq_pair_batch.pl stat
perl [% findbin %]/../extra/seq_pair_batch.pl \
    -d 1 --parallel [% parallel %] \
    -f [% seq_pair_file %] \
    -at 1000 -st 0  -r 1,2,21,40

## join_dbs.pl
#perl [% findbin %]/../extra/join_dbs.pl \
#    --multi --block --trimmed_fasta --length 1000 \
#    --goal_db [% name_str %] --outgroup 0query --target 0target \
#    --queries [% FOREACH i IN [ 1 .. query_ids.max ] %][% i %]query,[% END %] \
#    --dbs [% FOREACH id IN query_ids %][% target_id %]vs[% id %],[% END %]
#
## multi-way batch
#perl [% findbin %]/../extra/multi_way_batch.pl \
#    -d [% name_str %] \
#    -f [% working_dir %]/[% name_str %] \
#    --gff_file [% gff_files.join(',') %] \
#    -lt 1000 -st 0 --parallel [% parallel %] --batch 5 \
#    --run 10,21,30-32,40,41,43

# join_dbs.pl
perl [% findbin %]/../extra/join_dbs.pl \
    --no_insert --block --trimmed_fasta --length 1000 \
    --goal_db [% name_str %] --outgroup 0query --target 0target \
    --queries [% FOREACH i IN [ 1 .. query_ids.max ] %][% i %]query,[% END %] \
    --dbs [% FOREACH id IN query_ids %][% target_id %]vs[% id %],[% END %]

#----------------------------#
# RAxML
#----------------------------#
# raw phylo guiding tree
if [ ! -d [% working_dir %]/rawphylo ]
then
    mkdir [% working_dir %]/rawphylo
fi

cd [% working_dir %]/rawphylo

perl [% findbin %]/../../blastz/concat_fasta.pl \
    -i [% working_dir %]/[% name_str %]  \
    -o [% working_dir %]/rawphylo/[% name_str %].phy \
    -p

rm [% working_dir %]/rawphylo/RAxML*

[% IF query_ids.size > 2 -%]
raxml -T 2 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
    -o [% query_ids.0 %] -n [% name_str %] \
    -s [% working_dir %]/rawphylo/[% name_str %].phy

cp [% working_dir %]/rawphylo/RAxML_best* [% working_dir %]/rawphylo/[% name_str %].nwk
[% ELSE -%]
echo "(([% target_id %],[% query_ids.1 %]),[% query_ids.0 %]);" > [% working_dir %]/rawphylo/[% name_str %].nwk
[% END -%]

cd [% working_dir %]/..

# drop temp databases
[% FOREACH id IN query_ids -%]
[% sql_cmd %] -e "DROP DATABASE IF EXISTS [% target_id %]vs[% id %];"
[% END -%]

EOF
    $tt->process(
        \$text,
        {   stopwatch     => $stopwatch,
            parallel      => $parallel,
            working_dir   => $working_dir,
            findbin       => $FindBin::Bin,
            seq_pair_file => $seq_pair_file,
            name_str      => $name_str,
            target_id     => $target_id,
            query_ids     => \@query_ids,
            gff_files     => \@new_gff_files,
            sql_cmd       => "mysql -h$server -P$port -u$username -p$password ",
        },
        File::Spec->catfile( $working_dir, "cmd.sh" )
    ) or die Template->error;

    # cmd.bat
    $text = <<'EOF';
REM bac_bz.pl\n";
REM perl [% stopwatch.cmd_line %]

cd /d [% working_dir %]

REM basicstat
perl [% findbin %]/../fig/collect_common_basic.pl -d [% working_dir %]

REM multi chart
perl [% findbin %]/../stat/multi_chart_factory.pl -i [% findbin %]/../stat/[% name_str %].multi.xlsx

REM gc chart\n";
perl [% findbin %]/../stat/gc_chart_factory.pl --add_trend 1 -i [% findbin %]/../stat/[% name_str %].gc.xlsx\n\n";

EOF
    $tt->process(
        \$text,
        {   stopwatch     => $stopwatch,
            parallel      => $parallel,
            working_dir   => $working_dir,
            findbin       => $FindBin::Bin,
            seq_pair_file => $seq_pair_file,
            name_str      => $name_str,
            target_id     => $target_id,
            query_ids     => \@query_ids,
            gff_files     => \@new_gff_files,
        },
        File::Spec->catfile( $working_dir, "cmd.bat" )
    ) or die Template->error;
}

{
    my $round2_dir = File::Spec->catdir( $working_dir, 'round2' );
    mkdir $round2_dir unless -e $round2_dir;
    my $seq_pair_file = File::Spec->catfile( $round2_dir, "seq_pair.csv" );
    {    # write seq_pair.csv and left seq_pair_batch.pl to handle other things
        open my $fh, '>', $seq_pair_file;
        for my $query_id (@query_ids) {
            print {$fh} File::Spec->catdir( $round2_dir, $target_id ), ",",
                File::Spec->catdir( $round2_dir, $query_id ), "\n";
        }
        close $fh;
    }

    {
        open my $fh, '>', File::Spec->catfile( $round2_dir, "id2name.csv" );
        for my $id ( $target_id, @query_ids ) {
            print {$fh} "$id,$id\n";
        }
        close $fh;
    }

    my $tt = Template->new;

    # rm.sh
    my $text = <<'EOF';
#!/bin/bash

cd [% working_dir %]
if [ ! -d [% round2_dir %] ]
then
    mkdir [% round2_dir %]
fi

cp -R [% working_dir %]/[% target_id %] [% round2_dir %]
[% FOREACH id IN query_ids -%]
cp -R [% working_dir %]/[% id %] [% round2_dir %]
[% END -%]

cd [% round2_dir %]
#----------------------------#
# repeatmasker on all fasta
#----------------------------#
for f in `find . -name "*.fa"` ; do
    rename 's/fa$/fasta/' $f ;
done

for f in `find . -name "*.fasta"` ; do
    RepeatMasker $f -xsmall --parallel [% parallel %] ;
done

for f in `find . -name "*.fasta.out"` ; do
    rmOutToGFF3.pl $f > `dirname $f`/`basename $f .fasta.out`.rm.gff;
done

for f in `find . -name "*.fasta"` ; do
    if [ -f $f.masked ];
    then
        rename 's/fasta.masked$/fa/' $f.masked;
        find . -type f -name "`basename $f`*" | xargs rm;
    fi;
done;

echo Please check the following files
find [% round2_dir %] -name "*.fasta"

EOF
    $tt->process(
        \$text,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            round2_dir  => $round2_dir,
            target_id   => $target_id,
            query_ids   => \@query_ids,
        },
        File::Spec->catfile( $round2_dir, "file-rm.sh" )
    ) or die Template->error;

    # pair.sh
    $text = <<'EOF';
#!/bin/bash

cd [% round2_dir %]

#----------------------------#
# seq_pair
#----------------------------#
perl [% findbin %]/../extra/seq_pair_batch.pl -d 1 --parallel [% parallel %] \
    -f [% seq_pair_file %]  -at 1000 -st 0 -r 100-102

#perl [% findbin %]/../extra/seq_pair_batch.pl -d 1 --parallel [% parallel %] \
#    -f [% seq_pair_file %]  -at 1000 -st 0 -r 1,2,21,40

#----------------------------#
# mz
#----------------------------#
perl [% findbin %]/../../blastz/mz.pl \
    [% FOREACH id IN query_ids -%]
    -d [% round2_dir %]/[% target_id %]vs[% id %] \
    [% END -%]
    --tree [% working_dir %]/rawphylo/[% name_str %].nwk \
    --out [% round2_dir %]/[% name_str %] \
    -syn -p [% parallel %]

#----------------------------#
# maf2fasta
#----------------------------#
perl [% findbin %]/../../blastz/maf2fasta.pl \
    --has_outgroup -p [% parallel %] --block \
    -i [% round2_dir %]/[% name_str %] \
    -o [% round2_dir %]/[% name_str %]_fasta

#----------------------------#
# mafft
#----------------------------#
perl [% findbin %]/../../blastz/refine_fasta.pl \
    --msa mafft --block -p [% parallel %] \
    -i [% round2_dir %]/[% name_str %]_fasta \
    -o [% round2_dir %]/[% name_str %]_mft

#----------------------------#
# multi_way_batch
#----------------------------#
perl [% findbin %]/../extra/multi_way_batch.pl \
    -d [% name_str %] \
    -f [% round2_dir %]/[% name_str %]_mft \
    --gff_file [% gff_files.join(',') %] \
    --block --id [% round2_dir %]/id2name.csv \
    -lt 1000 -st 0 -ct 0 --parallel [% parallel %] --batch 5 \
    --run 1,10,21,30-32,40,41,43

#----------------------------#
# RAxML
#----------------------------#
# raw phylo guiding tree
if [ ! -d [% working_dir %]/phylo ]
then
    mkdir [% working_dir %]/phylo
fi

cd [% working_dir %]/phylo

perl [% findbin %]/../../blastz/concat_fasta.pl \
    -i [% round2_dir %]/[% name_str %]_mft  \
    -o [% working_dir %]/phylo/[% name_str %].phy \
    -p

rm [% working_dir %]/phylo/RAxML*

raxml -T 2 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
    -o [% query_ids.0 %] -n [% name_str %] \
    -s [% working_dir %]/phylo/[% name_str %].phy

EOF
    $tt->process(
        \$text,
        {   stopwatch     => $stopwatch,
            parallel      => $parallel,
            working_dir   => $working_dir,
            round2_dir    => $round2_dir,
            findbin       => $FindBin::Bin,
            seq_pair_file => $seq_pair_file,
            name_str      => $name_str,
            target_id     => $target_id,
            query_ids     => \@query_ids,
            gff_files     => \@new_gff_files,
        },
        File::Spec->catfile( $round2_dir, "pair-multi.sh" )
    ) or die Template->error;

}

$stopwatch->end_message;
exit;

__END__

perl bac_bz.pl --base_dir d:\bacteria\bacteria_101015 --parent 562
perl d:/wq/Scripts/tool/replace.pl -d d:/wq/Scripts/alignDB/bac -p "cmd.bat" -f /home/wangq -r d:/wq
