#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;

use AlignDB::Stopwatch;

use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
Add additional slippage-like info to alignDB. 1 for slippage-like and 0 for non.

    perl init/update_indel_slippage.pl -d S288cvsRM11_1a

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
    [ 'parallel=i', 'run in parallel mode',       { default => $conf->{generate}{parallel} }, ],
    [ 'batch=i',    '#alignments in one process', { default => $conf->{generate}{batch} }, ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# motif-repeat parameters
my $min_reps = {
    1 => 4,    # mononucl. with >= 4 repeats
};

# record config
$stopwatch->record_conf( { opt => $opt, min_reps => $min_reps, } );

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
$stopwatch->start_message("Update indel-slippage of [$opt->{db}]...");

my @jobs;
{
    my $alignDB = AlignDB->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

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

    my $alignDB = AlignDB->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    # Database handler
    my DBI $dbh = $alignDB->dbh;

    # select all indels in this alignment
    my $indel_query = q{
        SELECT indel_id, indel_start, indel_end, indel_length,
               indel_seq,  left_extand, right_extand
        FROM indel
        WHERE align_id = ?
    };
    my DBI $indel_sth = $dbh->prepare($indel_query);

    # update indel table in the new feature column
    my $indel_update = q{
        UPDATE indel
        SET indel_slippage = ?
        WHERE indel_id = ?
    };
    my DBI $indel_update_sth = $dbh->prepare($indel_update);

    # for every align_id
    for my $align_id (@align_ids) {
        $alignDB->process_message($align_id);

        my ( $target_seq, ) = @{ $alignDB->get_seqs($align_id) };

        $indel_sth->execute($align_id);
        while ( my @row = $indel_sth->fetchrow_array ) {
            my ($indel_id,  $indel_start, $indel_end, $indel_length,
                $indel_seq, $left_extand, $right_extand
            ) = @row;

            my $indel_slippage = 0;
            next unless $indel_seq;

            if ( exists $min_reps->{$indel_length} ) {
                my $reps         = $min_reps->{$indel_length};
                my $fland_length = $indel_length * $reps;

                my $left_flank = " ";    # avoid warning from $flank
                if ( $fland_length <= $left_extand ) {
                    $left_flank
                        = substr( $target_seq, $indel_start - $fland_length - 1, $fland_length );
                }

                my $right_flank = " ";
                if ( $fland_length <= $right_extand ) {
                    $right_flank = substr( $target_seq, $indel_end, $fland_length );
                }

                my $flank = $left_flank . $indel_seq . $right_flank;
                my $regex = $indel_seq . "{$reps,}";

                if ( $flank =~ /$regex/ ) {
                    $indel_slippage = 1;
                }
            }
            else {

                # indel 23-28, length 6: substr 17-22
                # seq start at 1 and string start at 0, so minus 1
                # substr(..., 16, 6)
                my $left_flank;
                if ( $indel_length <= $left_extand ) {
                    $left_flank
                        = substr( $target_seq, $indel_start - $indel_length - 1, $indel_length );
                }

                # indel 23-28, length 6: substr 29-34
                # substr(..., 28, 6)
                my $right_flank;
                if ( $indel_length <= $right_extand ) {
                    $right_flank = substr( $target_seq, $indel_end, $indel_length );
                }

                if ( $left_flank and $indel_seq eq $left_flank ) {
                    $indel_slippage = 1;
                }
                elsif ( $right_flank and $indel_seq eq $right_flank ) {
                    $indel_slippage = 1;
                }
            }

            $indel_update_sth->execute( $indel_slippage, $indel_id );
        }
    }

    $indel_update_sth->finish;
    $indel_sth->finish;

    return;
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $opt->{parallel}, chunk_size => $opt->{batch}, );
$mce->forchunk( \@jobs, $worker, );

$stopwatch->end_message;

# store program's meta info to database
AlignDB->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

exit;

__END__
