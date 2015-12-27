#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all);
use Math::Combinatorics;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Outgroup;
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

# ref parameter
my $length_threshold = $Config->{generate}{length_threshold} / 2;

# Database info
# Normal order: TvsO, TvsQ1, TvsQ2
my $dbs;
my $target;
my $queries;
my $outgroup;
my $goal_db;

my $trimmed_fasta = 0;
my $no_insert     = 0;

my $block = 0;    # output blocked fasta

my $reduce_end = 10;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    's|server=s'    => \$server,
    'P|port=i'      => \$port,
    'u|username=s'  => \$username,
    'p|password=s'  => \$password,
    'dbs=s'         => \$dbs,
    'goal_db=s'     => \$goal_db,
    'outgroup=s'    => \$outgroup,
    'target=s'      => \$target,
    'queries=s'     => \$queries,
    'length=i'      => \$length_threshold,
    'block'         => \$block,
    'trimmed_fasta' => \$trimmed_fasta,
    'no_insert'     => \$no_insert,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# perl init_alignDB.pl
#----------------------------------------------------------#
$stopwatch->start_message("Joining DBs...");

if ( !$no_insert and $goal_db ) {
    my $cmd
        = "perl $FindBin::Bin/../init/init_alignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -d=$goal_db"
        . " -u=$username"
        . " --password=$password";
    print "\n", "=" x 12, "CMD", "=" x 15, "\n";
    print $cmd , "\n";
    print "=" x 30, "\n";
    system($cmd);
}

#----------------------------------------------------------#
# dbs, names
#----------------------------------------------------------#
my ( @all_dbs, @all_names, @ingroup_names, $target_db );
{
    @all_dbs = grep {defined} split ",", $dbs;
    my @queries = grep {defined} split ",", $queries;

    if ( !$target ) {
        die "Target not defined\n";
    }
    elsif ( !scalar @queries ) {
        die "Queries not defined\n";
    }

    if ($outgroup) {
        print "Outgroup defined.\n";
        print "We will produce alignments with outgroup.\n";

        if ( scalar @all_dbs != scalar @queries + 1 ) {
            printf "DB [%d]\tQueries [%d]\n", scalar @all_dbs, scalar @queries;
            die "Numbers of DB doesn't match with numbers of strains\n";
        }
        @all_names = ( $target, @queries, $outgroup );
    }
    else {
        print "Outgroup not defined.\n";
        print "We will produce alignments without outgroup.\n";

        if ( scalar @all_dbs != scalar @queries ) {
            printf "DB [%d]\tQueries [%d]\n", scalar @all_dbs, scalar @queries;
            die "Numbers of DB doesn't match with numbers of strains\n";
        }
        @all_names = ( $target, @queries );
    }

    @ingroup_names = ( $target, @queries );

    $target =~ /^(\d+)(.+)/;
    $target_db = $all_dbs[$1];
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#

#----------------------------#
# info hash
#----------------------------#
my $db_info_of
    = build_db_info( $server, $username, $password, \@all_dbs, $reduce_end );

#----------------------------#
# build intersect chromosome set
#----------------------------#
my $chr_set_of = build_inter_chr( $db_info_of, \@all_dbs, $target_db );

for my $chr_id ( sort keys %{$chr_set_of} ) {
    my $inter_chr_set = $chr_set_of->{$chr_id};
    my $chr_name      = $db_info_of->{$target_db}{chrs}{$chr_id}{name};

    #----------------------------#
    # process each intersects
    #----------------------------#

    my @segments = $inter_chr_set->spans;
SEG: for my $seg (@segments) {
        my $seg_start  = $seg->[0];
        my $seg_end    = $seg->[1];
        my $seg_length = $seg_end - $seg_start + 1;
        next if $seg_length <= $length_threshold;

        print "$chr_name:$seg_start-$seg_end; length:$seg_length\n";

        for my $db_name (@all_dbs) {
            my $pos_obj = $db_info_of->{$db_name}{pos_obj};
            my ( $align_id, $dummy ) = @{
                $pos_obj->positioning_align_chr_id( $chr_id, $seg_start,
                    $seg_end )
            };

            if ( !defined $align_id ) {
                warn " " x 4, "Find no align in $db_name, jump to next\n";
                next SEG;
            }
            elsif ( defined $dummy ) {
                warn " " x 4, "Overlapped alignment in $db_name!\n";
            }
            $db_info_of->{$db_name}{align_id} = $align_id;
        }

        #----------------------------#
        # get seq, use align coordinates
        #----------------------------#
        print " " x 4, "build seqs\n";
        for my $db_name (@all_dbs) {
            my $align_id = $db_info_of->{$db_name}{align_id};

            my $error
                = build_seq( $db_info_of->{$db_name}, $seg_start, $seg_end );
            if ($error) {
                warn $error . " in $db_name $align_id\n";
                next SEG;
            }
        }

        #----------------------------#
        # build $seq_of, $seq_names, $seq_idxs
        #----------------------------#
        my $name_of  = {};    # taxon_id => name
        my $seq_idxs = [];    # store simplified seq_names
        my ( $seq_of, $seq_names ) = ( {}, [] );
        for my $i ( 0 .. @all_names - 1 ) {
            my $name = $all_names[$i];
            $name =~ /^(\d+)(.+)/;
            my $db_name_idx = $1;
            my $torq        = $2;
            if ( not( $torq =~ /^t/i or $torq =~ /^q/i ) ) {
                die "$torq is not equal to target or query\n";
            }
            my $db_name = $all_dbs[$db_name_idx];

            my $info = $db_info_of->{$db_name}{$torq};

            $name_of->{ $info->{taxon_id} } = $info->{name};

            push @{$seq_idxs}, $i;

            my $header = encode_header($info);
            push @{$seq_names}, $header;
            $seq_of->{$i} = $info->{seq};
        }

        #----------------------------#
        # realign
        #----------------------------#
        print " " x 4, "realign regions\n";
        realign_all( $seq_of, $seq_idxs );

        # trim header and footer indels
        trim_pure_dash( $seq_of, $seq_idxs );

        # trim outgroup only sequence
        if ($outgroup) {
            trim_outgroup( $seq_of, $seq_idxs );
        }

        # record complex indels and ingroup indels
        if ($outgroup) {
            trim_complex_indel( $seq_of, $seq_idxs );
        }

        #----------------------------#
        # output a fasta alignment for further use
        #----------------------------#
        if ($trimmed_fasta) {
            mkdir $goal_db, 0777 if !-e $goal_db;

            my $outfile = "./$goal_db/" . $chr_name . ".fas";
            print " " x 4, "$outfile\n";
            open my $out_fh, '>>', $outfile;
            write_fasta_fh( $out_fh, $seq_of, $seq_idxs, $seq_names );
            print {$out_fh} "\n";
            close $out_fh;
        }

        #----------------------------#
        # insert as a three-way alignDB
        #----------------------------#
        if ( !$no_insert ) {
            my $goal_obj;
            if ($outgroup) {
                $goal_obj = AlignDB::Outgroup->new(
                    mysql  => "$goal_db:$server",
                    user   => $username,
                    passwd => $password,
                );
            }
            else {
                $goal_obj = AlignDB->new(
                    mysql  => "$goal_db:$server",
                    user   => $username,
                    passwd => $password,
                );
            }

            my $info_refs = [];
            my $seq_refs  = [];

            for my $i ( 0 .. @{$seq_names} - 1 ) {
                my $header = $seq_names->[$i];
                my $info   = decode_header($header);
                push @{$info_refs}, $info;

                my $seq = $seq_of->{ $seq_idxs->[$i] };
                push @{$seq_refs}, $seq;
            }
            $goal_obj->update_names($name_of);
            $goal_obj->add_align( $info_refs, $seq_refs );
        }
    }
}

$stopwatch->end_message;

# store program running meta info to database
END {
    if ( !$no_insert ) {
        AlignDB->new(
            mysql  => "$goal_db:$server",
            user   => $username,
            passwd => $password,
        )->add_meta_stopwatch($stopwatch);
    }
}
exit;

#----------------------------------------------------------#
# subs
#----------------------------------------------------------#
sub build_db_info {
    my $server     = shift;
    my $username   = shift;
    my $password   = shift;
    my $all_dbs    = shift;
    my $reduce_end = shift;

    my $db_info_of = {};
    for my $db ( @{$all_dbs} ) {
        my $cur_obj = AlignDB->new(
            mysql  => "$db:$server",
            user   => $username,
            passwd => $password,
        );
        my $cur_dbh = $cur_obj->dbh;
        my $cur_pos_obj = AlignDB::Position->new( dbh => $cur_dbh );
        $db_info_of->{$db} = {
            target => {
                taxon_id => '',
                name     => '',
            },
            query => {
                taxon_id => '',
                name     => '',
            },
            obj     => $cur_obj,
            pos_obj => $cur_pos_obj,
        };

        (   $db_info_of->{$db}{target}{taxon_id},
            $db_info_of->{$db}{query}{taxon_id},
        ) = $cur_obj->get_taxon_ids;

        ( $db_info_of->{$db}{target}{name}, $db_info_of->{$db}{query}{name}, )
            = $cur_obj->get_names;

        my $chr_id_set = AlignDB::IntSpan->new;

        my $chr_ref = $cur_obj->get_chrs('target');

        for my $ref ( @{$chr_ref} ) {
            my ( $chr_id, $chr_name, $chr_length ) = @{$ref};
            my $chr_set = build_chr_set( $cur_dbh, $chr_id, $reduce_end );
            $db_info_of->{$db}{chrs}{$chr_id}{set}  = $chr_set;
            $db_info_of->{$db}{chrs}{$chr_id}{name} = $chr_name;
            $chr_id_set->add($chr_id);
        }
        $db_info_of->{$db}{chr_id_set} = $chr_id_set;
    }

    return $db_info_of;
}

sub build_chr_set {
    my $dbh        = shift;
    my $chr_id     = shift;
    my $reduce_end = shift || 0;

    my $chr_set = AlignDB::IntSpan->new;

    my $chr_query = qq{
        SELECT  s.chr_start + $reduce_end,
                s.chr_end - $reduce_end
        FROM sequence s, chromosome c
        WHERE c.chr_id = ?
        AND s.chr_id = c.chr_id
    };

    # build $chr_set
    my $chr_sth = $dbh->prepare($chr_query);
    $chr_sth->execute($chr_id);
    while ( my @row = $chr_sth->fetchrow_array ) {
        my ( $chr_start, $chr_end ) = @row;
        next if $chr_start > $chr_end;
        $chr_set->add_range( $chr_start, $chr_end );
    }

    return $chr_set;
}

sub build_inter_chr {
    my $db_info_of = shift;
    my $all_dbs    = shift;
    my $target_db  = shift;

    my %chr_set_of;
    for my $chr_id ( $db_info_of->{$target_db}{chr_id_set}->elements ) {
        my ($chr_name) = $db_info_of->{$target_db}{obj}->get_chr_info($chr_id);
        print "\nchr_id: $chr_id\tchr_name: $chr_name\n";

        my $inter_chr_set = AlignDB::IntSpan->new;
        for my $db_name ( @{$all_dbs} ) {
            my $cur_chr_set = $db_info_of->{$db_name}{chrs}{$chr_id}{set};
            $cur_chr_set = AlignDB::IntSpan->new unless $cur_chr_set;
            if ( $inter_chr_set->is_empty ) {
                $inter_chr_set = $cur_chr_set;
            }
            else {
                $inter_chr_set = $inter_chr_set->intersect($cur_chr_set);
            }
        }

        $chr_set_of{$chr_id} = $inter_chr_set;
    }

    return \%chr_set_of;
}

# get seq, use align coordinates
sub build_seq {
    my $db_info   = shift;
    my $seg_start = shift;
    my $seg_end   = shift;

    my $obj      = $db_info->{obj};
    my $pos_obj  = $db_info->{pos_obj};
    my $align_id = $db_info->{align_id};

    my $target_info = $obj->get_target_info($align_id);
    $db_info->{target}{chr_id}     = $target_info->{chr_id};
    $db_info->{target}{chr_name}   = $target_info->{chr_name};
    $db_info->{target}{chr_strand} = $target_info->{chr_strand};

    my ($query_info) = $obj->get_queries_info($align_id);
    $db_info->{query}{chr_id}     = $query_info->{chr_id};
    $db_info->{query}{chr_name}   = $query_info->{chr_name};
    $db_info->{query}{chr_strand} = $query_info->{query_strand};

    ( $db_info->{target}{full_seq}, $db_info->{query}{full_seq} )
        = @{ $obj->get_seqs($align_id) };

    my $align_start = $pos_obj->at_align( $align_id, $seg_start );
    my $align_end   = $pos_obj->at_align( $align_id, $seg_end );

    # align_start and align_end should must be available
    unless ( $align_start and $align_end ) {
        return " " x 8 . "align_start or align_end error";
    }

    my $align_length = $align_end - $align_start + 1;

    # target chr position
    $db_info->{target}{chr_start} = $seg_start;
    $db_info->{target}{chr_end}   = $seg_end;

    # query chr position
    $db_info->{query}{chr_start}
        = $pos_obj->at_query_chr( $align_id, $align_start );
    $db_info->{query}{chr_end}
        = $pos_obj->at_query_chr( $align_id, $align_end );

    $db_info->{target}{seq}
        = substr( $db_info->{target}{full_seq}, $align_start - 1,
        $align_length );
    $db_info->{query}{seq}
        = substr( $db_info->{query}{full_seq}, $align_start - 1,
        $align_length );

    unless (length $db_info->{target}{seq} == length $db_info->{query}{seq}
        and length $db_info->{target}{seq} > 0 )
    {
        return " " x 8 . "seq-length error";
    }

    delete $db_info->{target}{full_seq};
    delete $db_info->{query}{full_seq};

    return;
}

#----------------------------#
# realign all seqs
#----------------------------#
sub realign_all {
    my $seq_of    = shift;
    my $seq_names = shift;

    my @seqs;
    for ( @{$seq_names} ) {
        push @seqs, $seq_of->{$_};
    }

    my $realigned_seqs = multi_align( \@seqs, 'mafft' );

    for my $i ( 0 .. scalar @{$seq_names} - 1 ) {
        $seq_of->{ $seq_names->[$i] } = uc $realigned_seqs->[$i];
    }

    return;
}

sub write_fasta_fh {
    my $fh         = shift;
    my $seq_of     = shift;
    my $seq_names  = shift;
    my $real_names = shift;

    for my $i ( 0 .. @{$seq_names} - 1 ) {
        my $seq = $seq_of->{ $seq_names->[$i] };
        my $header;
        if ($real_names) {
            $header = $real_names->[$i];
        }
        else {
            $header = $seq_names->[$i];
        }

        print {$fh} ">" . $header . "\n";
        print {$fh} $seq . "\n";
    }

    return;
}

__END__

=head1 NAME

    join_dbs.pl - join multiple dbs for three-lineage test or maligndb

=head1 SYNOPSIS

    perl join_dbs.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --port              MySQL server port
        --username          username
        --password          password
        --dbs               DB names list seperated by ','
        --goal_db           goal database name
        --outgroup          outgroup identity (0query)
        --target            target identity (0target)
        --queries           query list (1query,2query)
        --length            threshold of alignment length
        --realign           correct pesudo-alignment error
        --trimmed_fasta     save ref-trimmed fasta files
        --reduce_end        reduce align end to avoid some overlaps in
                              BlastZ results (use 10 instead of 0)
                            For two independent datasets, use 10;
                            for two dependent datasets, use 0
        --no_insert         don't insert into goal_db actually
        --indel_expand
        --indel_join

$ perl join_dbs.pl --dbs S288CvsSpar,S288CvsRM11,S288CvsYJM789 \
    --goal_db S288CvsThree --no_insert --trimmed_fasta \
    --outgroup 0query --target 0target --queries 1query,2query

# windows
perl two_way_batch.pl -d S288CvsRM11 -t "4932,S288C" -q "285006,RM11" -da d:\data\alignment\yeast_combine\S288CvsRM11\ -lt 5000 --parallel 4 --run skeleton
perl two_way_batch.pl -d S288CvsSpar -t "4932,S288C" -q "226125,Spar" -da d:\data\alignment\yeast_combine\S288CvsSpar\ -lt 5000 --parallel 4 --run skeleton

perl join_dbs.pl --dbs S288CvsSpar,S288CvsRM11 --goal_db S288CvsRM11refSpar --target 0target --outgroup 0query --queries 1query --length 5000
perl join_dbs.pl --dbs S288CvsRM11,S288CvsSpar --goal_db S288CvsRM11_Spar --target 0target --queries 0query,1query --length 5000

# linux/mac
cd ~/Scripts/alignDB

perl extra/two_way_batch.pl -d S288CvsRM11 -t "559292,S288C" -q "285006,RM11" \
    -da data/S288CvsRM11 -lt 5000 --parallel 4 --run skeleton

perl extra/two_way_batch.pl -d S288CvsSpar -t "559292,S288C" -q "226125,Spar" \
    -da data/S288CvsSpar -lt 5000 --parallel 4 --run skeleton

perl extra/join_dbs.pl --dbs S288CvsRM11,S288CvsSpar --goal_db S288CvsRM11Spar \
    --target 0target --queries 0query,1query \
    --no_insert --trimmed_fasta --length 5000

perl extra/join_dbs.pl --dbs S288CvsSpar,S288CvsRM11 --goal_db S288CvsRM11refSpar \
    --target 0target --outgroup 0query --queries 1query --length 5000

perl extra/join_dbs.pl --dbs S288CvsRM11,S288CvsSpar --goal_db S288CvsRM11_Spar \
    --target 0target --queries 0query,1query \
    --length 5000
