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
use Text::Table;
use List::MoreUtils qw(uniq);

use Bio::Taxon;
use Bio::DB::Taxonomy;
use Template;

use AlignDB::IntSpan;
use AlignDB::Run;
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
my $ref_id;
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
    'r|ref_id=i'      => \$ref_id,
    'e|exclude=s'     => \$exclude_ids,
    'n|name_str=s'    => \$name_str,
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
        my $sth = $dbh->prepare(
            qq{
            SELECT st.taxonomy_id
            FROM strain st
            WHERE st.seq_ok = 1
            }
        );
        $sth->execute;
        while ( my ($id) = $sth->fetchrow_array ) {
            $db_id_set->add($id);
        }
    }

    my $id_set = $sub_id_set->intersect($db_id_set);
    $id_set->remove( split /,/, $exclude_ids );
    $id_str = '(' . ( join ",", $id_set->as_array ) . ')';
}

{    # making working dir
    if ( !$name_str ) {
        my $sth = $dbh->prepare(
            qq{
        SELECT DISTINCT species
        FROM strain st
        WHERE st.taxonomy_id IN $id_str
        }
        );
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
    my $sth = $dbh->prepare(
        qq{
        SELECT  st.taxonomy_id,
                st.organism_name,
                st.released_date
        FROM strain st
        WHERE st.taxonomy_id IN $id_str
        ORDER BY st.released_date
        }
    );
    $sth->execute;

    # header line
    my @strains;
    my $table = Text::Table->new( @{ $sth->{NAME} } );
    while ( my @row = $sth->fetchrow_array ) {
        my $taxon_id = $row[0];
        push @strains, [@row];
        $table->load( \@row );
    }

    my $table_file = File::Spec->catfile( $working_dir, "table.txt" );
    open my $fh, '>', $table_file;
    print {$fh} $table;
    print $table;
    print "\n";

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

    if ($ref_id) {
        my ($exist) = grep { $_ == $ref_id } @query_ids;
        if ( defined $exist ) {
            my $message = "Use [$exist] as reference, as you wish.\n";
            print {$fh} $message;
            print $message;

            @query_ids = map { $_ == $ref_id ? () : $_ } @query_ids;
            unshift @query_ids, $ref_id;
        }
        else {
            print "Taxon $ref_id doesn't exist, please check.\n";
        }
    }

    close $fh;
}

my @new_gff_files;
{    # build fasta files

    # read all filenames, then grep
    print "Reading file list\n";
    my @fna_files = File::Find::Rule->file->name('*.fna')->in($base_dir);
    my @gff_files = File::Find::Rule->file->name('*.gff')->in($base_dir);

    print "Rewrite seqs for every strains\n";
    for my $taxon_id ( $target_id, @query_ids ) {
        my $id_dir = File::Spec->catdir( $working_dir, $taxon_id );
        mkdir $id_dir unless -e $id_dir;

        my $sth = $dbh->prepare(
            qq{
            SELECT s.accession, s.replicon, s.length
            FROM seq s
            WHERE s.taxonomy_id = ?
            AND s.replicon like "%chr%"
            }
        );
        $sth->execute($taxon_id);

        while ( my ( $acc, $rep, $length ) = $sth->fetchrow_array ) {
            my $fa_file = File::Spec->catfile( $id_dir, "$acc.fa" );
            my ($fna_file) = grep {/$acc/} @fna_files;
            if ( !$fna_file ) {
                warn "Can't find fasta file for $acc\n";
                next;
            }

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

            if ( $taxon_id eq $target_id ) {
                print "Copy target gff\n";
                my ($gff_file) = grep {/$acc/} @gff_files;
                my $new_gff_file
                    = File::Spec->catfile( $working_dir, "$acc.gff" );
                push @new_gff_files, $new_gff_file;
                copy( $gff_file, $new_gff_file );
            }
        }
    }
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

    my $text = <<'EOF';
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

[% ids = query_ids; ids.unshift(target_id) -%]
[% FOREACH item IN ids -%]
cp -R [% working_dir %]/[% item %] [% round2_dir %]
[% END -%]

cd [% working_dir %]/round2
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
