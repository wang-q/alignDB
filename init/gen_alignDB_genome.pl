#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use MCE;
use Path::Tiny;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use App::Fasops::Common;

use lib "$FindBin::RealBin/../lib";
use AlignDB::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
Generate alignDB from genomic fasta files

    perl init/gen_alignDB_genome.pl -d S288Cvsself --parallel 2 \
        -t S288c --da ~/data/alignment/example/scer/Genomes/S288c

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
    [ 'dir_align|da=s', 'dir contains .fa files',    { required => 1 }, ],
    [ 'target|t=s',     'target name',               { required => 1 }, ],
    [ 'length|lt=i',    'threshold of lengths',      { default  => $conf->{generate}{length} }, ],
    [ 'truncate=i',     'truncated length',          { default  => 100_000 }, ],
    [ 'fill=i',         'fill holes less than this', { default  => 50 }, ],
    [ 'parallel=i',     'run in parallel mode',      { default  => $conf->{generate}{parallel} }, ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# record config
$stopwatch->record_conf($opt);

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

die "target_name not defined\n" unless $opt->{target};

#----------------------------------------------------------#
# Search for all files and push their paths to @axt_files
#----------------------------------------------------------#
$opt->{dir_align} = path( $opt->{dir_align} )->stringify;
my @files
    = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )->in( $opt->{dir_align} );
printf "\n----Total .fa Files: %4s----\n\n", scalar @files;

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $infile = $chunk_ref->[0];
    my $wid    = MCE->wid;

    my $inner    = AlignDB::Stopwatch->new;
    my $basename = path($infile)->basename;
    $inner->block_message("Process task [$chunk_id] by worker #$wid. [$basename]");

    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    my $chr_name = path($infile)->basename( '.fasta', '.fas', '.fa' );

    my $seq_of     = App::Fasops::Common::read_fasta($infile);
    my $chr_seq    = $seq_of->{ ( keys %{$seq_of} )[0] };
    my $chr_length = length $chr_seq;

    my $ambiguous_set = AlignDB::IntSpan->new;
    for ( my $pos = 0; $pos < $chr_length; $pos++ ) {
        my $base = substr $chr_seq, $pos, 1;
        if ( $base =~ /[^ACGT-]/i ) {
            $ambiguous_set->add( $pos + 1 );
        }
    }

    printf "Ambiguous chromosome region for [%s]:\n    %s\n", $chr_name, $ambiguous_set->runlist;

    my $valid_set = AlignDB::IntSpan->new("1-$chr_length");
    $valid_set->subtract($ambiguous_set);
    $valid_set = $valid_set->fill( $opt->{fill} - 1 );    # fill small gaps

    printf "Valid chromosome region for [%s]:\n    %s\n", $chr_name, $valid_set->runlist;

    my @regions;                                          # ([start, end], [start, end], ...)
    for my AlignDB::IntSpan $set ( $valid_set->sets ) {
        my $size = $set->size;
        next if $size < $opt->{length};

        my @set_regions;
        my $pos = $set->min;
        my $max = $set->max;
        while ( $max - $pos + 1 > $opt->{truncate} ) {
            push @set_regions, [ $pos, $pos + $opt->{truncate} - 1 ];
            $pos += $opt->{truncate};
        }
        if ( scalar @set_regions > 0 ) {
            $set_regions[-1]->[1] = $max;
        }
        else {
            @set_regions = ( [ $pos, $max ] );
        }
        push @regions, @set_regions;
    }

    #print Dump \@regions;

    for my $region (@regions) {
        my ( $start, $end ) = @{$region};
        my $seq = substr $chr_seq, $start - 1, $end - $start + 1;

        my $info_refs = [
            {   name   => $opt->{target},
                chr    => $chr_name,
                start  => $start,
                end    => $end,
                strand => '+',
                seq    => $seq,
            },
            {   name   => $opt->{target},
                chr    => $chr_name,
                start  => $start,
                end    => $end,
                strand => '+',
                seq    => $seq,
            },
        ];

        $alignDB->add_align($info_refs);
    }

    $inner->block_message( "[$basename] has been processed.", "duration" );

    return;
};

#----------------------------------------------------------#
# start insert
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $opt->{parallel}, chunk_size => 1, );
$mce->foreach( [ sort @files ], $worker );    # foreach also implies chunk_size => 1

$stopwatch->end_message( "All files have been processed.", "duration" );

# store program's meta info to database
AlignDB::Common->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

exit;

__END__
