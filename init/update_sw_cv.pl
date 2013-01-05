#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Run;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::GC;

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
my $db       = $Config->{database}{db};

my $stat_segment_size = 500;
my $stat_window_size  = $Config->{gc}{stat_window_size};
my $stat_window_step  = $Config->{gc}{stat_window_step};

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update $db...");

#----------------------------#
# Create columnas and find all align_ids
#----------------------------#
my @align_ids;
{
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # add column
    $obj->create_column( "exonsw",   "exonsw_cv",   "DOUBLE" );
    $obj->create_column( "codingsw", "codingsw_cv", "DOUBLE" );
    $obj->create_column( "ofgsw",    "ofgsw_cv",    "DOUBLE" );
    $obj->create_column( "isw",      "isw_cv",      "DOUBLE" );
    $obj->create_column( "gsw",      "gsw_cv",      "DOUBLE" );
    print "Table exonsw, codingsw, ofgsw, isw and gsw altered\n";

    @align_ids = @{ $obj->get_align_ids };
}

my @jobs;
while ( scalar @align_ids ) {
    my @batching = splice @align_ids, 0, $batch_number;
    push @jobs, [@batching];
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

my $worker = sub {
    my $job       = shift;
    my @align_ids = @$job;

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    Moo::Role->apply_roles_to_object( $obj, qw{ AlignDB::GC } );
    my %opt = (
        stat_window_size => $stat_window_size,
        stat_window_step => $stat_window_step,
    );
    for my $key ( sort keys %opt ) {
        $obj->$key( $opt{$key} );
    }

    # Database handler
    my $dbh = $obj->dbh;

    my $exonsw_sth = $dbh->prepare(
        q{
        SELECT s.exonsw_id, w.window_runlist
        FROM exonsw s, window w
        where s.window_id = w.window_id
        and w.align_id = ?
        }
    );

    my $exonsw_update_sth = $dbh->prepare(
        q{
        UPDATE exonsw
        SET exonsw_cv = ?
        WHERE exonsw_id = ?
        }
    );

    my $codingsw_sth = $dbh->prepare(
        q{
        SELECT s.codingsw_id, w.window_runlist
        FROM codingsw s, window w
        where s.window_id = w.window_id
        and w.align_id = ?
        }
    );

    my $codingsw_update_sth = $dbh->prepare(
        q{
        UPDATE codingsw
        SET codingsw_cv = ?
        WHERE codingsw_id = ?
        }
    );

    my $ofgsw_sth = $dbh->prepare(
        q{
        SELECT s.ofgsw_id, w.window_runlist
        FROM ofgsw s, window w
        where s.window_id = w.window_id
        and w.align_id = ?
        }
    );

    my $ofgsw_update_sth = $dbh->prepare(
        q{
        UPDATE ofgsw
        SET ofgsw_cv = ?
        WHERE ofgsw_id = ?
        }
    );

    my $isw_sth = $dbh->prepare(
        q{
        SELECT s.isw_id, s.isw_start, s.isw_end
        FROM isw s, indel i
        where s.indel_id = i.indel_id
        and i.align_id = ?
        }
    );

    my $isw_update_sth = $dbh->prepare(
        q{
        UPDATE isw
        SET isw_cv = ?
        WHERE isw_id = ?
        }
    );

    my $gsw_sth = $dbh->prepare(
        q{
        SELECT s.gsw_id, w.window_runlist
        FROM gsw s, window w
        where s.window_id = w.window_id
        and w.align_id = ?
        }
    );

    my $gsw_update_sth = $dbh->prepare(
        q{
        UPDATE gsw
        SET gsw_cv = ?
        WHERE gsw_id = ?
        }
    );

    for my $align_id (@align_ids) {
        my $target_info    = $obj->get_target_info($align_id);
        my $target_runlist = $target_info->{seq_runlist};

        $obj->process_message($align_id);

        # sliding in target_set
        my $target_set = AlignDB::IntSpan->new($target_runlist);

        $exonsw_sth->execute($align_id);
        while ( my @row = $exonsw_sth->fetchrow_array ) {
            my ( $exonsw_id, $window_runlist ) = @row;
            my $window_set = AlignDB::IntSpan->new($window_runlist);
            my $resize_set
                = center_resize( $window_set, $target_set, $stat_segment_size );

            my $seqs_ref = $obj->get_seqs($align_id);
            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( $seqs_ref, $resize_set );
            $exonsw_update_sth->execute( $gc_cv, $exonsw_id );
        }

        $codingsw_sth->execute($align_id);
        while ( my @row = $codingsw_sth->fetchrow_array ) {
            my ( $codingsw_id, $window_runlist ) = @row;
            my $window_set = AlignDB::IntSpan->new($window_runlist);
            my $resize_set
                = center_resize( $window_set, $target_set, $stat_segment_size );

            my $seqs_ref = $obj->get_seqs($align_id);
            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( $seqs_ref, $resize_set );
            $codingsw_update_sth->execute( $gc_cv, $codingsw_id );
        }

        $ofgsw_sth->execute($align_id);
        while ( my @row = $ofgsw_sth->fetchrow_array ) {
            my ( $ofgsw_id, $window_runlist ) = @row;
            my $window_set = AlignDB::IntSpan->new($window_runlist);
            my $resize_set
                = center_resize( $window_set, $target_set, $stat_segment_size );

            next unless $resize_set;

            my $seqs_ref = $obj->get_seqs($align_id);
            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( $seqs_ref, $resize_set );
            $ofgsw_update_sth->execute( $gc_cv, $ofgsw_id );
        }

        $isw_sth->execute($align_id);
        while ( my @row = $isw_sth->fetchrow_array ) {
            my ( $isw_id, $start, $end ) = @row;
            my $window_set = AlignDB::IntSpan->new("$start-$end");
            my $resize_set
                = center_resize( $window_set, $target_set, $stat_segment_size );

            next unless $resize_set;

            my $seqs_ref = $obj->get_seqs($align_id);
            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( $seqs_ref, $resize_set );
            $isw_update_sth->execute( $gc_cv, $isw_id );
        }

        $gsw_sth->execute($align_id);
        while ( my @row = $gsw_sth->fetchrow_array ) {
            my ( $gsw_id, $window_runlist ) = @row;
            my $window_set = AlignDB::IntSpan->new($window_runlist);
            my $resize_set
                = center_resize( $window_set, $target_set, $stat_segment_size );

            next unless $resize_set;

            my $seqs_ref = $obj->get_seqs($align_id);
            my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
                = $obj->segment_gc_stat( $seqs_ref, $resize_set );
            $gsw_update_sth->execute( $gc_cv, $gsw_id );
        }
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}
exit;

sub center_resize {
    my $old_set    = shift;
    my $parent_set = shift;
    my $resize     = shift;

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

=head1 NAME

    update_sw_cv.pl - CV for exonsw, codingsw, ofgsw and isw

=head1 SYNOPSIS

    update_sw_cv.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut

