#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;
use Path::Tiny;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use App::RL::Common;

use lib "$FindBin::RealBin/../lib";
use AlignDB::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
Add annotations to alignDB

    perl init/update_annotation.pl -d S288cvsRM11_1a --parallel 2

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
    [ 'annotation|a=s', 'YAML file for cds/repeat', { default => $conf->{generate}{file_anno} }, ],
    [ 'parallel=i',     'run in parallel mode',     { default => $conf->{generate}{parallel} }, ],
    [ 'batch=i', '#alignments in one process', { default => $conf->{generate}{batch} }, ],
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
# other objects are initiated in subroutines
$stopwatch->start_message("Update annotations of $opt->{db}...");

my $file_anno = path( $opt->{annotation} )->stringify;
my $yml       = YAML::Syck::LoadFile($file_anno);

die "Invalid annotation YAML. Need cds.\n"    unless defined $yml->{cds};
die "Invalid annotation YAML. Need repeat.\n" unless defined $yml->{repeat};

my $cds_set_of    = App::RL::Common::runlist2set( $yml->{cds} );
my $repeat_set_of = App::RL::Common::runlist2set( $yml->{repeat} );

my @jobs;
{
    my $alignDB = AlignDB::Common->new(
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

    my $alignDB = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    # Database handler
    my DBI $dbh = $alignDB->dbh;

    #----------------------------#
    # SQL query and DBI sths
    #----------------------------#
    # update align table
    my DBI $align_feature_sth = $dbh->prepare(
        q{
        UPDATE align
        SET align_coding = ?,
            align_repeats = ?,
            align_coding_runlist = ?,
            align_repeats_runlist = ?
        WHERE align_id = ?
        }
    );

    # select all indels in this alignment
    my DBI $indel_query_sth = $dbh->prepare(
        q{
        SELECT indel_id, indel_start, indel_end
        FROM indel
        WHERE align_id = ?
        }
    );

    # update indel table in the new feature column
    my DBI $indel_feature_sth = $dbh->prepare(
        q{
        UPDATE indel
        SET indel_coding  = ?,
            indel_repeats = ?
        WHERE indel_id    = ?
        }
    );

    # select all isws for this indel
    my DBI $isw_query_sth = $dbh->prepare(
        q{
        SELECT isw_id, isw_start, isw_end
        FROM isw
        WHERE indel_id = ?
        }
    );

    # update isw table in the new feature column
    my DBI $isw_feature_sth = $dbh->prepare(
        q{
        UPDATE isw
        SET isw_coding  = ?,
            isw_repeats = ?
        WHERE isw_id    = ?
        }
    );

    # select all snps for this indel
    my DBI $snp_query_sth = $dbh->prepare(
        q{
        SELECT snp_id, snp_pos
        FROM snp
        WHERE align_id = ?
        }
    );

    # update snp table in the new feature column
    my DBI $snp_feature_sth = $dbh->prepare(
        q{
        UPDATE snp
        SET snp_coding  = ?,
            snp_repeats = ?
        WHERE snp_id    = ?
        }
    );

    # select all windows for this alignment
    my DBI $window_query_sth = $dbh->prepare(
        q{
        SELECT window_id, window_runlist
        FROM window
        WHERE align_id = ?
        }
    );

    # update window table in the new feature column
    my DBI $window_update_sth = $dbh->prepare(
        q{
        UPDATE window
        SET window_coding = ?, window_repeats = ?
        WHERE window_id = ?
        }
    );

UPDATE: for my $align_id (@align_ids) {

        #----------------------------#
        # for each alignment
        #----------------------------#
        my $target_info  = $alignDB->get_target_info($align_id);
        my $chr_name     = $target_info->{chr_name};
        my $chr_start    = $target_info->{chr_start};
        my $chr_end      = $target_info->{chr_end};
        my $align_length = $target_info->{align_length};

        next UPDATE if $chr_name =~ /rand|contig|hap|scaf|gi_/i;

        my ($target_seq) = @{ $alignDB->get_seqs($align_id) };
        $alignDB->process_message($align_id);

        $chr_name =~ s/chr0?//i;

        # slicing cds_set and repeat_set
        my $slice_cds_set    = AlignDB::IntSpan->new;
        my $slice_repeat_set = AlignDB::IntSpan->new;
        {
            my $align_chr_set = AlignDB::IntSpan->new;
            $align_chr_set->add_pair( $chr_start, $chr_end );
            if ( defined $cds_set_of->{$chr_name} ) {
                $slice_cds_set = $cds_set_of->{$chr_name}->intersect($align_chr_set);
            }
            if ( defined $repeat_set_of->{$chr_name} ) {
                $slice_repeat_set = $repeat_set_of->{$chr_name}->intersect($align_chr_set);
            }
        }

        #----------------------------#
        # construct transforming array and hash
        #----------------------------#
        # align_position to chr_position transforming array
        # the first element will be ignored
        my @chr_pos;
        $chr_pos[$align_length] = $align_length;    # pre-allocating memory

        # chr_position to align_position transforming hash
        my %align_pos;

        my $indel_count = 0;
        for my $i ( 1 .. $align_length ) {
            my $current_base = substr( $target_seq, $i - 1, 1 );
            if ( $current_base eq '-' ) {
                $indel_count++;
            }
            my $chr_position = $chr_start + $i - 1 - $indel_count;
            $chr_pos[$i] = $chr_position;
            if ( !exists $align_pos{$chr_position} ) {
                $align_pos{$chr_position} = $i;
            }
        }

        #----------------------------#
        # process align
        #----------------------------#
        {
            my $align_chr_start   = $chr_pos[1];
            my $align_chr_end     = $chr_pos[$align_length];
            my $align_chr_runlist = "$align_chr_start-$align_chr_end";

            # feature protions and runlists
            my ( $align_coding, $align_repeats, $align_cds_runlist, $align_repeat_runlist );

            $align_coding = feature_portion( $slice_cds_set, $align_chr_runlist );
            $align_cds_runlist
                = $slice_cds_set->map_set( sub { $align_pos{$_} } )->runlist;

            $align_repeats = feature_portion( $slice_repeat_set, $align_chr_runlist );
            $align_repeat_runlist
                = $slice_repeat_set->map_set( sub { $align_pos{$_} } )->runlist;

            $align_feature_sth->execute( $align_coding, $align_repeats,
                $align_cds_runlist, $align_repeat_runlist, $align_id, );

            $align_feature_sth->finish;
        }

        #----------------------------#
        # process each indels
        #----------------------------#
        $indel_query_sth->execute($align_id);
        while ( my @row3 = $indel_query_sth->fetchrow_array ) {
            my ( $indel_id, $indel_start, $indel_end ) = @row3;
            my $indel_chr_start = $chr_pos[$indel_start];
            my $indel_chr_end   = $chr_pos[$indel_end];

            {
                my $indel_chr_runlist = "$indel_chr_start-$indel_chr_end";

                # feature protions
                my ( $indel_coding, $indel_repeats );
                $indel_coding  = feature_portion( $slice_cds_set,    $indel_chr_runlist );
                $indel_repeats = feature_portion( $slice_repeat_set, $indel_chr_runlist );

                $indel_feature_sth->execute( $indel_coding, $indel_repeats, $indel_id, );
                $indel_feature_sth->finish;
            }

            #----------------------------#
            # process each isws
            #----------------------------#
            {
                $isw_query_sth->execute($indel_id);
                while ( my @row4 = $isw_query_sth->fetchrow_array ) {
                    my ( $isw_id, $isw_start, $isw_end ) = @row4;
                    my $isw_chr_start = $chr_pos[$isw_start];
                    my $isw_chr_end   = $chr_pos[$isw_end];

                    my $isw_chr_runlist = "$isw_chr_start-$isw_chr_end";

                    # feature protions
                    my ( $isw_coding, $isw_repeats );
                    $isw_coding  = feature_portion( $slice_cds_set,    $isw_chr_runlist );
                    $isw_repeats = feature_portion( $slice_repeat_set, $isw_chr_runlist );

                    $isw_feature_sth->execute( $isw_coding, $isw_repeats, $isw_id, );
                }

                $isw_feature_sth->finish;
                $isw_query_sth->finish;
            }
        }
        $indel_query_sth->finish;

        #----------------------------#
        # process each snps
        #----------------------------#
        {
            $snp_query_sth->execute($align_id);
            while ( my @row = $snp_query_sth->fetchrow_array ) {
                my ( $snp_id, $snp_pos ) = @row;
                my $snp_chr_pos = $chr_pos[$snp_pos];

                my $snp_coding  = $slice_cds_set->contains($snp_chr_pos);
                my $snp_repeats = $slice_repeat_set->contains($snp_chr_pos);

                $snp_feature_sth->execute( $snp_coding, $snp_repeats, $snp_id, );
            }

            $snp_feature_sth->finish;
            $snp_query_sth->finish;
        }

        #----------------------------#
        # process each windows
        #----------------------------#
        {
            $window_query_sth->execute($align_id);
            while ( my @row5 = $window_query_sth->fetchrow_array ) {
                my ( $window_id, $window_runlist ) = @row5;
                my $window_set = AlignDB::IntSpan->new($window_runlist);
                my $window_chr_set
                    = $window_set->map_set( sub { $chr_pos[$_] } );

                # feature protions
                my ( $window_coding, $window_repeats );
                $window_coding  = feature_portion( $slice_cds_set,    $window_chr_set );
                $window_repeats = feature_portion( $slice_repeat_set, $window_chr_set );

                $window_update_sth->execute( $window_coding, $window_repeats, $window_id, );
            }

            $window_update_sth->finish;
            $window_query_sth->finish;
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

sub feature_portion {
    my AlignDB::IntSpan $feature_set = shift;
    my $pos_set = shift;

    my $set;
    if ( ref $pos_set eq "AlignDB::IntSpan" ) {
        $set = $pos_set;
    }
    else {
        $set = AlignDB::IntSpan->new($pos_set);
    }

    my $pos_n = $set->size;
    return if $pos_n <= 0;

    my $intersect       = $feature_set->intersect($set);
    my $n               = $intersect->size;
    my $feature_portion = $n / $pos_n;

    return $feature_portion;
}

__END__
