#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Stopwatch;

use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

update_indel_slippage.pl - Add additional slippage-like info to alignDB
                            1 for slippage-like and 0 for non

=head1 SYNOPSIS

    perl update_indel_slippage.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password

=cut

# motif-repeat parameters
my $min_reps = {
    1 => 4,    # mononucl. with >= 4 repeats
};

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server       = $Config->{database}{server} ),
    'port|P=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'       => \( my $db           = $Config->{database}{db} ),
    'username|u=s' => \( my $username     = $Config->{database}{username} ),
    'password|p=s' => \( my $password     = $Config->{database}{password} ),
) or Getopt::Long::HelpMessage(1);

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
my DBI $dbh = $obj->dbh;

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    my @align_ids = @{ $obj->get_align_ids };

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

    # for indel
    for my $align_id (@align_ids) {
        print "Processing align_id $align_id\n";

        my ( $target_seq, ) = @{ $obj->get_seqs($align_id) };

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
                        = substr( $target_seq, $indel_start - $fland_length - 1,
                        $fland_length );
                }

                my $right_flank = " ";
                if ( $fland_length <= $right_extand ) {
                    $right_flank
                        = substr( $target_seq, $indel_end, $fland_length );
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
                        = substr( $target_seq, $indel_start - $indel_length - 1,
                        $indel_length );
                }

                # indel 23-28, length 6: substr 29-34
                # substr(..., 28, 6)
                my $right_flank;
                if ( $indel_length <= $right_extand ) {
                    $right_flank
                        = substr( $target_seq, $indel_end, $indel_length );
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
}

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

__END__
