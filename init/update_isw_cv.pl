#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

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

# XXX not done yet
# run in parallel mode
my $parallel = $Config->{feature}{parallel};

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
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update $db...");

my $obj = AlignDB::Multi->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
AlignDB::GC->meta->apply($obj);
my %opt = (
    stat_window_size => $stat_window_size,
    stat_window_step => $stat_window_step,
);
for my $key ( sort keys %opt ) {
    $obj->$key( $opt{$key} );
}

# Database handler
my $dbh = $obj->dbh;

{    # add column isw_cv and codingsw_cv
    $obj->create_column( "isw", "isw_cv", "DOUBLE" );
    print "Table isw altered\n";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

my @align_ids = @{ $obj->get_align_ids };

my $isw_sth = $dbh->prepare(
    q{
    SELECT s.isw_id, s.isw_start, s.isw_end
    FROM indel i, isw s
    WHERE i.indel_id = s.indel_id
    AND i.align_id = ?
    }
);

my $isw_update_sth = $dbh->prepare(
    q{
    UPDATE isw
    SET isw_cv = ?
    WHERE isw_id = ?
    }
);

for my $align_id (@align_ids) {
    my $target_info    = $obj->get_target_info($align_id);
    my $chr_name       = $target_info->{chr_name};
    my $chr_start      = $target_info->{chr_start};
    my $chr_end        = $target_info->{chr_end};
    my $target_runlist = $target_info->{seq_runlist};
    my $align_length   = $target_info->{align_length};

    print $obj->process_message;

    # comparable runlist
    my $target_set = AlignDB::IntSpan->new($target_runlist);

    $isw_sth->execute($align_id);
    while ( my @row = $isw_sth->fetchrow_array ) {
        my ( $isw_id, $isw_start, $isw_end ) = @row;
        my $window_set = AlignDB::IntSpan->new("$isw_start-$isw_end");
        my $resize_set
            = center_resize( $window_set, $target_set, $stat_segment_size );

        my $seqs_ref = $obj->get_seqs($align_id);
        my ( $gc_mean, $gc_std, $gc_cv, $gc_mdcw )
            = $obj->segment_gc_stat( $seqs_ref, $resize_set );
        $isw_update_sth->execute( $gc_cv, $isw_id );
    }
}

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB::Multi->new(
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

    update_isw_cv.pl - Add CV of each isw to alignDB

=head1 SYNOPSIS

    update_isw_cv.pl [options]
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

