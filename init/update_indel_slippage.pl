#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

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
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

#motif-repeat parameters
my $min_reps = {
    1 => 4,    # mononucl. with >= 4 repeats
};

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
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update indel-slippage of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# add a column indel_feature3 to indel_extra
# add a column window_ns_indel to window
{
    $obj->create_column( "indel_extra", "indel_feature3", "DOUBLE" );
    print "Table indel_extra altered\n";

    $obj->create_column( "window", "window_ns_indel", "DOUBLE" );
    print "Table window altered\n";
}

# check rows in indel_extra
# prefill table indel_extra
{
    my $sql_query = qq{
        SELECT COUNT(indel_extra_id)
        FROM indel_extra
    };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my ($count) = $sth->fetchrow_array;

    if ( !$count ) {
        $sql_query = qq{
            INSERT INTO indel_extra (indel_id, prev_indel_id)
            SELECT indel.indel_id, indel.prev_indel_id
            FROM indel
        };
        $sth = $dbh->prepare($sql_query);
        $sth->execute;
    }
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{

    # alignments
    my $align_query = q{
        SELECT align_id
        FROM align 
    };
    my $align_sth = $dbh->prepare($align_query);

    # sequence
    my $seq_query = q{
        SELECT t.target_seq, q.query_seq
        FROM align a, target t, query q
        WHERE a.align_id = ?
        AND a.align_id = q.align_id
        AND a.align_id = t.align_id
    };
    my $seq_sth = $dbh->prepare($seq_query);

    # select all indels in this alignment
    my $indel_query = q{
        SELECT indel_id, indel_start, indel_end, indel_length,
               indel_seq, indel_insert, left_extand, right_extand
        FROM indel
        WHERE align_id = ?
    };
    my $indel_sth = $dbh->prepare($indel_query);

    # update indel table in the new feature column
    my $indel_extra = q{
        UPDATE indel_extra
        SET indel_feature3 = ?
        WHERE indel_id = ?
    };
    my $indel_extra_sth = $dbh->prepare($indel_extra);

    # window
    my $window_query = q{
        SELECT window_id
        FROM window
        WHERE align_id = ?
        AND window_indel > 0
        AND window_start < ?
        AND window_end > ?
    };
    my $window_sth = $dbh->prepare($window_query);

    # update window table with slippage-like indel number
    my $window_update = q{
        UPDATE window
        SET window_ns_indel = IFNULL(window_ns_indel, 0) + 1
        WHERE window_id = ?
    };
    my $window_update_sth = $dbh->prepare($window_update);

    my $window_ns = q{
        UPDATE window
        SET window_ns_indel = window_indel - IFNULL(window_ns_indel, 0)
        WHERE align_id = ?
    };
    my $window_ns_sth = $dbh->prepare($window_ns);

    $align_sth->execute();

    # for indel
    while ( my @row = $align_sth->fetchrow_array ) {
        my ($align_id) = @row;
        print "Processing align_id $align_id\n";

        $seq_sth->execute($align_id);
        my ( $target_seq, $query_seq ) = $seq_sth->fetchrow_array;

        $indel_sth->execute($align_id);
        while ( my @row = $indel_sth->fetchrow_array ) {
            my ($indel_id,  $indel_start,  $indel_end,   $indel_length,
                $indel_seq, $indel_insert, $left_extand, $right_extand
            ) = @row;

            my $insert_host;
            if ( $indel_insert eq 'D' ) {
                $insert_host = \$target_seq;
            }
            else {
                $insert_host = \$query_seq;
            }

            my $indel_feature3 = 0;

            if ( exists $min_reps->{$indel_length} ) {
                my $reps         = $min_reps->{$indel_length};
                my $fland_length = $indel_length * $reps;

                my $left_flank = " ";    # avoid warning from $flank
                if ( $fland_length <= $left_extand ) {
                    $left_flank = substr( $$insert_host,
                        $indel_start - $fland_length - 1,
                        $fland_length );
                }

                my $right_flank = " ";
                if ( $fland_length <= $right_extand ) {
                    $right_flank
                        = substr( $$insert_host, $indel_end, $fland_length );
                }

                my $flank = $left_flank . $indel_seq . $right_flank;
                my $regex = $indel_seq . "{$reps,}";

                if ( $flank =~ /$regex/ ) {
                    $indel_feature3 = 1;
                }
            }
            else {

                # indel 23-28, length 6: substr 17-22
                # seq start at 1 and string start at 0, so minus 1
                # substr(..., 16, 6)
                my $left_flank;
                if ( $indel_length <= $left_extand ) {
                    $left_flank = substr( $$insert_host,
                        $indel_start - $indel_length - 1,
                        $indel_length );
                }

                # indel 23-28, length 6: substr 29-34
                # substr(..., 28, 6)
                my $right_flank;
                if ( $indel_length <= $right_extand ) {
                    $right_flank
                        = substr( $$insert_host, $indel_end, $indel_length );
                }

                if ( $left_flank and $indel_seq eq $left_flank ) {
                    $indel_feature3 = 1;
                }
                elsif ( $right_flank and $indel_seq eq $right_flank ) {
                    $indel_feature3 = 1;
                }
            }

            $indel_extra_sth->execute( $indel_feature3, $indel_id );

            if ($indel_feature3) {
                $window_sth->execute( $align_id, $indel_start, $indel_end );
                while ( my @row = $window_sth->fetchrow_array ) {
                    my ($window_id) = @row;
                    $window_update_sth->execute($window_id);
                }
            }
        }
        $window_ns_sth->execute($align_id);
    }

    $window_ns_sth->finish;
    $window_update_sth->finish;
    $window_sth->finish;
    $indel_extra_sth->finish;
    $indel_sth->finish;
    $align_sth->finish;
}

$stopwatch->end_message;

# store program running meta info to database
END {
    $obj->add_meta_stopwatch($stopwatch);
}
exit;

__END__

=head1 NAME

    update_indel_slippage.pl - Add additional slippage-like info to alignDB
                                 1 for slippage-like and 0 for non

=head1 SYNOPSIS

    update_indel_slippage.pl [options]
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

