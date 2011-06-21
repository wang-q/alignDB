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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{bac}{db};

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
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Blastz whole species");

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

my $name_str;
{    # making working dir
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

    {
        my $table_file = File::Spec->catfile( $working_dir, "table.txt" );
        open my $fh, '>', $table_file;
        print {$fh} $table;
        print $table;
        print "\n";
        close $fh;
    }
    print "There are " . scalar @strains . " strains\n";

    if ($target_id) {
        my ($exist) = grep { $_->[0] == $target_id } @strains;
        if ( defined $exist ) {
            print "Use [$exist->[1]] as target, as you wish.\n";
        }
        else {
            print "Taxon $target_id doesn't exist, please check.\n";
            exit;
        }
    }
    else {
        $target_id = $strains[0]->[0];
        print
            "Use [$strains[0]->[1]] as target, the oldest strain on NCBI.\n";
    }

    @query_ids = map { $_->[0] == $target_id ? () : $_->[0] } @strains;
    
    if ($ref_id) {
        my ($exist) = grep { $_ == $ref_id } @query_ids;
        if (defined $exist) {
            print "Use [$exist] as reference, as you wish.\n";
            @query_ids = map { $_ == $ref_id ? () : $_ } @query_ids;
            unshift @query_ids, $ref_id;
        }
        else {
            print "Taxon $ref_id doesn't exist, please check.\n";
        }
    }
}

my %dir_of;
my $new_gff_file;
{    # build fasta files

    # read all filenames, then grep
    print "Reading file list\n";
    my @fna_files = File::Find::Rule->file->name('*.fna')->in($base_dir);
    my @gff_files = File::Find::Rule->file->name('*.gff')->in($base_dir);

    print "Rewrite seqs for every strains\n";
    for my $taxon_id ( $target_id, @query_ids ) {
        my $id_dir = File::Spec->catdir( $working_dir, $taxon_id );
        $dir_of{$taxon_id} = $id_dir;
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
            
            if ($taxon_id eq $target_id) {
                print "Copy target gff\n";
                my ($gff_file) = grep {/$acc/} @gff_files;
                $new_gff_file = File::Spec->catfile( $working_dir, "$acc.gff" );
                copy($gff_file, $new_gff_file);
            }
        }
    }
    
    
}

my $seq_pair_file = File::Spec->catfile( $working_dir, "seq_pair.csv" );
{    # write seq_pair.csv and left seq_pair_batch.pl to handle other things
    open my $fh, '>', $seq_pair_file;
    for my $query_id (@query_ids) {
        print {$fh} $dir_of{$target_id}, ",", $dir_of{$query_id}, "\n";
    }
    close $fh;
}

{
    my $cmd_file = File::Spec->catfile( $working_dir, "cmd.txt" );
    open my $fh, '>', $cmd_file;

    print {$fh} "# bac_bz.pl\n";
    print {$fh} "perl ", $stopwatch->cmd_line, "\n\n";

    print {$fh} "# seq_pair_batch.pl\n";
    print {$fh} "perl $FindBin::Bin/../extra/seq_pair_batch.pl"
        . " -d 1 -p 4"
        . " -f $seq_pair_file" . "\n\n";

    print {$fh} "# join_dbs.pl\n";
    print {$fh} "perl $FindBin::Bin/../extra/join_dbs.pl"
        . " --no_insert 1 --trimmed_fasta 1"
        . " --length 1000 --reduce_end 10"
        . " --goal_db $name_str"
        . " --outgroup 0query --target 0target"
        . " --queries "
        . ( join ",", map { $_ . "query" } ( 1 .. scalar @query_ids - 1 ) )
        . " --dbs "
        . ( join ",", map { $target_id . "vs" . $_ } @query_ids ) . "\n\n";
    
    print {$fh} "# drop temp databases\n";
    my $sql_cmd = "mysql -h$server -P$port -u$username -p$password ";
    for (@query_ids) {
        my $tempdb = $target_id . "vs" . $_;
        print {$fh} $sql_cmd  . " -e \"DROP DATABASE IF EXISTS $tempdb;\"\n";
    }
    
    print {$fh} "\n";
    print {$fh} "# multi-way batch\n";
    print {$fh} "perl $FindBin::Bin/../extra/multi_way_batch.pl"
        . " -d $name_str"
        . " -f $working_dir/$name_str"
        . " -gff_file $new_gff_file"
        . " --all_freq " . scalar @query_ids
        . " -lt 10000 -st 100000 --parallel=4 --run all\n\n";

    close $fh;
}

$stopwatch->end_message;
exit;

__END__

perl bac_bz.pl --base_dir d:\bacteria\bacteria_101015 --parent 562

