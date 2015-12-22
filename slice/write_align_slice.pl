#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Basename;
use Path::Tiny;

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB;
use AlignDB::Position;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

write_align_slice.pl - extract alignment slices from alignDB

=head1 SYNOPSIS

    perl write_align_slice.pl -f <yaml file> [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --file      -f  STR     yaml filename
        --type      -t  STR     fas or axt. Default is [fas]
        --outgroup              This alignDB has an outgroup.
                                Work with fas.
        --want_equal            An align should be an island in the slice.
                                Work with axt.
        --parallel      INT     run in parallel mode

=head1 YAML format

    cat bin1.yml
    
    ---
    1: 1-25744,815056-817137
  
    output file: bin1/1.fas

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'file|f=s'     => \( my $yaml_file ),
    'type|t=s'     => \( my $type     = 'fas' ),
    'want_equal'   => \( my $want_equal ),
    'outgroup'     => \( my $outgroup ),
    'parallel=i'   => \( my $parallel = 1 ),
) or HelpMessage(1);

if ( !defined $yaml_file ) {
    die HelpMessage(1);
}
elsif ( !path($yaml_file)->is_file ) {
    die "[$yaml_file] doesn't exist\n";
}

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Write slice files from $db...");

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
print "Loading $yaml_file\n";
my ( $base, $dir ) = fileparse( $yaml_file, ".yaml", ".yml" );

my $slice_set_of = LoadFile($yaml_file);

my $worker = sub {
    my $chr_name  = shift;
    my $slice_set = AlignDB::IntSpan->new( $slice_set_of->{$chr_name} );
    my $outfile   = path( $dir, "$chr_name.axt" )->stringify;
    if ( $type eq 'axt' ) {
        write_slice_axt( $chr_name, $slice_set, $outfile );
    }
    else {
        write_slice( $chr_name, $slice_set, $outfile );
    }
};

my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => [ sort keys %{$slice_set_of} ],
    code     => $worker,
);
$run->run;

sub write_slice {
    my $chr_name  = shift;
    my $slice_set = shift;
    my $outfile   = shift;

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # position finder
    my $pos_obj = AlignDB::Position->new( dbh => $dbh );

    print "Write slice from $chr_name\n";
    print "Output file is $outfile\n";

    # alignment
    my @align_ids = @{ $obj->get_align_ids_of_chr_name($chr_name) };

    for my $align_id (@align_ids) {
        local $| = 1;
        print "Processing align_id $align_id\n";

        # target
        my $target_info      = $obj->get_target_info($align_id);
        my $target_chr_name  = $target_info->{chr_name};
        my $target_chr_start = $target_info->{chr_start};
        my $target_chr_end   = $target_info->{chr_end};
        my $target_runlist   = $target_info->{seq_runlist};

        my $align_chr_set = AlignDB::IntSpan->new("$target_chr_start-$target_chr_end");
        my $iset          = $slice_set->intersect($align_chr_set);
        next if $iset->is_empty;

        my ( $target_seq, @query_seqs ) = @{ $obj->get_seqs($align_id) };
        my $ref_seq;
        if ($outgroup) {
            $ref_seq = $obj->get_seq_ref($align_id);
        }

        my ( $target_name, @query_names ) = $obj->get_names($align_id);
        my $ref_name;
        if ($outgroup) {
            $ref_name = pop @query_names;
        }

        if ( scalar @query_seqs != scalar @query_names ) {
            warn "Number of queries names is not equal to one of seqs.\n";
            next;
        }

        my $target_set = AlignDB::IntSpan->new($target_runlist);

        # there may be two or more subslice intersect this alignment
        for my $ss_set ( $iset->sets ) {

            # position set
            my $ss_start = $pos_obj->at_align( $align_id, $ss_set->min );
            my $ss_end   = $pos_obj->at_align( $align_id, $ss_set->max );
            next if $ss_start >= $ss_end;
            $ss_set = AlignDB::IntSpan->new("$ss_start-$ss_end");
            $ss_set = $ss_set->intersect($target_set);

            next if $ss_set->count <= 1;
            my ( $seg_start, $seg_end ) = ( $ss_set->min, $ss_set->max );
            my $seg_length = $seg_end - $seg_start + 1;

            # align coordinates to target & query chromosome coordinates
            my $target_seg_start = $pos_obj->at_target_chr( $align_id, $seg_start );
            my $target_seg_end   = $pos_obj->at_target_chr( $align_id, $seg_end );

            # append fas file
            # S288C.chrI(+):27520-29557|species=S288C
            {
                print " " x 4
                    . "Append to slice files: "
                    . "$target_chr_name:$target_seg_start-$target_seg_end" . "\n";
                open my $outfh, '>>', $outfile;
                print {$outfh}
                    ">$target_name.$target_chr_name(+):$target_chr_start-$target_chr_end|species=$target_name\n";
                print {$outfh} substr( $target_seq, $seg_start - 1, $seg_length ), "\n";

                for my $i ( 0 .. $#query_seqs ) {
                    print {$outfh} ">" . $query_names[$i] . "\n";
                    print {$outfh}
                        substr( $query_seqs[$i], $seg_start - 1, $seg_length ),
                        "\n";
                }
                if ($outgroup) {
                    print {$outfh} ">" . $ref_name . "\n";
                    print {$outfh}
                        substr( $ref_seq, $seg_start - 1, $seg_length ),
                        "\n";
                }
                print {$outfh} "\n";
                close $outfh;
            }
        }
        print " " x 4 . "Finish write slice files\n";
    }
}

sub write_slice_axt {
    my $chr_name  = shift;
    my $slice_set = shift;
    my $outfile   = shift;

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # position finder
    my $pos_obj = AlignDB::Position->new( dbh => $dbh );

    # get target and query names
    my ( $target_name, $query_name ) = $obj->get_names;

    print "Write slice from $chr_name\n";
    print "Output file is $outfile\n";

    # alignment
    my @align_ids = @{ $obj->get_align_ids_of_chr_name($chr_name) };

    my %align_serial;

    for my $align_id (@align_ids) {
        local $| = 1;
        print "Processing align_id $align_id\n";

        # target
        my $target_info      = $obj->get_target_info($align_id);
        my $target_chr_name  = $target_info->{chr_name};
        my $target_chr_start = $target_info->{chr_start};
        my $target_chr_end   = $target_info->{chr_end};
        my $target_runlist   = $target_info->{seq_runlist};

        # query
        my ($query_info)    = $obj->get_queries_info($align_id);
        my $query_chr_name  = $query_info->{chr_name};
        my $query_chr_start = $query_info->{chr_start};
        my $query_chr_end   = $query_info->{chr_end};
        my $query_runlist   = $query_info->{seq_runlist};
        my $query_strand    = $query_info->{query_strand};

        my ( $target_seq, $query_seq ) = @{ $obj->get_seqs($align_id) };

        my $target_set = AlignDB::IntSpan->new($target_runlist);
        my $query_set  = AlignDB::IntSpan->new($query_runlist);

        my $align_chr_set = AlignDB::IntSpan->new("$target_chr_start-$target_chr_end");
        my $iset          = $slice_set->intersect($align_chr_set);
        next if $iset->is_empty;

        if ($want_equal) {
            my $start_island = $slice_set->find_islands($target_chr_start);
            if ( !$start_island->equal($align_chr_set) ) {
                print "The align is not equal to one island\n";
                next;
            }
        }

        # there may be two or more subslice intersect this alignment
        for my $ss_set ( $iset->sets ) {

            # position set
            my $ss_start = $pos_obj->at_align( $align_id, $ss_set->min );
            my $ss_end   = $pos_obj->at_align( $align_id, $ss_set->max );
            next if $ss_start >= $ss_end;
            $ss_set = AlignDB::IntSpan->new("$ss_start-$ss_end");
            $ss_set = $ss_set->intersect($target_set);

            next if $ss_set->count <= 1;
            my ( $seg_start, $seg_end ) = ( $ss_set->min, $ss_set->max );
            my $seg_length = $seg_end - $seg_start + 1;

            # prepare axt summary line
            $align_serial{$target_chr_name}++;
            my $serial = $align_serial{$target_chr_name} - 1;

            # align coordinates to target & query chromosome coordinates
            my $target_seg_start = $pos_obj->at_target_chr( $align_id, $seg_start );
            my $target_seg_end   = $pos_obj->at_target_chr( $align_id, $seg_end );
            my $query_seg_start = $pos_obj->at_query_chr( $align_id, $seg_start );
            my $query_seg_end   = $pos_obj->at_query_chr( $align_id, $seg_end );
            if ( $query_strand eq '-' ) {
                ( $query_seg_start, $query_seg_end ) = ( $query_seg_end, $query_seg_start );
            }
            my $score = $seg_length * 100;    # sham score

            # append axt file
            {
                print " " x 4
                    . "Append to slice files: "
                    . "$target_chr_name:$target_seg_start-$target_seg_end" . "\n";
                open my $outfh, '>>', $outfile;
                print {$outfh} "$serial";
                print {$outfh} " $target_chr_name";
                print {$outfh} " $target_seg_start $target_seg_end";
                print {$outfh} " $query_chr_name";
                print {$outfh} " $query_seg_start $query_seg_end";
                print {$outfh} " $query_strand $score\n";
                print {$outfh} substr( $target_seq, $seg_start - 1, $seg_length ), "\n";
                print {$outfh} substr( $query_seq, $seg_start - 1, $seg_length ), "\n";
                print {$outfh} "\n";
                close $outfh;
            }
        }
        print " " x 4 . "Finish write slice files\n";
    }
}

$stopwatch->end_message;

__END__
