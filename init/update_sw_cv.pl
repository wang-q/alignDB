#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;

use AlignDB::GC;
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
Update CV for isw, gsw and ofgsw

    perl init/update_sw_cv.pl -d S288cvsRM11_1a --parallel 2

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

# record config
$stopwatch->record_conf( { opt => $opt, gc => $conf->{gc}, } );

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update [$opt->{db}]...");

#----------------------------#
# Create columnas and find all align_ids
#----------------------------#
my @jobs;
{
    my $alignDB = AlignDB->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    # add column
    $alignDB->create_column( "gsw", "gsw_intra_cv", "DOUBLE" );
    print "Table gsw altered\n";

    @jobs = @{ $alignDB->get_align_ids };
}

#----------------------------------------------------------#
# start update
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
    my $obj_gc = AlignDB::GC->new(
        stat_window_size => $conf->{gc}{stat_window_size},
        stat_window_step => $conf->{gc}{stat_window_step},
    );

    # Database handler
    my DBI $dbh = $alignDB->dbh;

    my DBI $ofgsw_sth = $dbh->prepare(
        q{
        SELECT s.ofgsw_id, w.window_runlist
        FROM ofgsw s, window w
        where s.window_id = w.window_id
        and w.align_id = ?
        }
    );

    my DBI $ofgsw_update_sth = $dbh->prepare(
        q{
        UPDATE ofgsw
        SET ofgsw_cv = ?
        WHERE ofgsw_id = ?
        }
    );

    my DBI $isw_sth = $dbh->prepare(
        q{
        SELECT s.isw_id, s.isw_start, s.isw_end
        FROM isw s, indel i
        where s.indel_id = i.indel_id
        and i.align_id = ?
        }
    );

    my DBI $isw_update_sth = $dbh->prepare(
        q{
        UPDATE isw
        SET isw_cv = ?
        WHERE isw_id = ?
        }
    );

    my DBI $gsw_sth = $dbh->prepare(
        q{
        SELECT s.gsw_id, w.window_runlist
        FROM gsw s, window w
        where s.window_id = w.window_id
        and w.align_id = ?
        }
    );

    my DBI $gsw_update_sth = $dbh->prepare(
        q{
        UPDATE gsw
        SET gsw_cv = ?,
            gsw_intra_cv = ?
        WHERE gsw_id = ?
        }
    );

    for my $align_id (@align_ids) {
        my $target_info    = $alignDB->get_target_info($align_id);
        my $target_runlist = $target_info->{seq_runlist};

        $alignDB->process_message($align_id);

        # sliding in target_set
        my $target_set = AlignDB::IntSpan->new($target_runlist);

        $ofgsw_sth->execute($align_id);
        while ( my @row = $ofgsw_sth->fetchrow_array ) {
            my ( $ofgsw_id, $window_runlist ) = @row;
            my $window_set = AlignDB::IntSpan->new($window_runlist);
            my $resize_set
                = center_resize( $window_set, $target_set, $conf->{gc}{stat_segment_size} );

            next unless $resize_set;

            my $seqs_ref = $alignDB->get_seqs($align_id);
            my ( undef, undef, $gc_cv, undef ) = $obj_gc->segment_gc_stat( $seqs_ref, $resize_set );
            $ofgsw_update_sth->execute( $gc_cv, $ofgsw_id );
        }

        $isw_sth->execute($align_id);
        while ( my @row = $isw_sth->fetchrow_array ) {
            my ( $isw_id, $start, $end ) = @row;
            my $window_set = AlignDB::IntSpan->new("$start-$end");
            my $resize_set
                = center_resize( $window_set, $target_set, $conf->{gc}{stat_segment_size} );

            next unless $resize_set;

            my $seqs_ref = $alignDB->get_seqs($align_id);
            my ( undef, undef, $gc_cv, undef ) = $obj_gc->segment_gc_stat( $seqs_ref, $resize_set );
            $isw_update_sth->execute( $gc_cv, $isw_id );
        }

        $gsw_sth->execute($align_id);
        while ( my @row = $gsw_sth->fetchrow_array ) {
            my ( $gsw_id, $window_runlist ) = @row;
            my $window_set = AlignDB::IntSpan->new($window_runlist);
            my $resize_set
                = center_resize( $window_set, $target_set, $conf->{gc}{stat_segment_size} );

            next unless $resize_set;

            my $seqs_ref = $alignDB->get_seqs($align_id);
            my ( undef, undef, $gc_cv, undef ) = $obj_gc->segment_gc_stat( $seqs_ref, $resize_set );
            my ( undef, undef, $gc_intra_cv, undef )
                = $obj_gc->segment_gc_stat( $seqs_ref, $window_set, 20, 20 );
            $gsw_update_sth->execute( $gc_cv, $gc_intra_cv, $gsw_id );
        }
    }
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

sub center_resize {
    my AlignDB::IntSpan $old_set    = shift;
    my AlignDB::IntSpan $parent_set = shift;
    my $resize                      = shift;

    # find the middles of old_set
    my $half_size           = int( $old_set->size / 2 );
    my $midleft             = $old_set->at($half_size);
    my $midright            = $old_set->at( $half_size + 1 );
    my $midleft_parent_idx  = $parent_set->index($midleft);
    my $midright_parent_idx = $parent_set->index($midright);

    return unless $midleft_parent_idx and $midright_parent_idx;

    # map to parent
    my $parent_size  = $parent_set->size;
    my $half_resize  = int( $resize / 2 );
    my $new_left_idx = $midleft_parent_idx - $half_resize + 1;
    $new_left_idx = 1 if $new_left_idx < 1;
    my $new_right_idx = $midright_parent_idx + $half_resize - 1;
    $new_right_idx = $parent_size if $new_right_idx > $parent_size;

    my $new_set = $parent_set->slice( $new_left_idx, $new_right_idx );

    return $new_set;
}

__END__
