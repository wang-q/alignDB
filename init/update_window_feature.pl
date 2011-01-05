#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Ensembl;

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
my $server     = $Config->{database}->{server};
my $port       = $Config->{database}->{port};
my $username   = $Config->{database}->{username};
my $password   = $Config->{database}->{password};
my $db         = $Config->{database}->{db};
my $ensembl_db = $Config->{database}->{ensembl};

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
    'ensembl=s'  => \$ensembl_db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{
    my @align_ids = @{ $obj->get_align_ids };

    # ensembl handler
    my $ensembl = AlignDB::Ensembl->new(
        server => $server,
        db     => $ensembl_db,
        user   => $username,
        passwd => $password,
    );

    # alignments' chromosomal location
    my $align_seq_query = q{
        SELECT c.chr_name,
               a.align_length,
               s.chr_start,
               s.chr_end,
               t.target_runlist,
               t.target_seq
        FROM align a, target t, sequence s, chromosome c
        WHERE a.align_id = t.align_id
        AND t.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        AND a.align_id = ?
    };
    my $align_seq_sth = $dbh->prepare($align_seq_query);

    # select all windows for this alignment
    my $window_query = q{
        SELECT window_id, window_runlist
        FROM window
        WHERE align_id = ?
        AND window_coding IS NULL
    };
    my $window_query_sth = $dbh->prepare($window_query);

    # update window table in the new feature column
    my $window_update_query = q{
        UPDATE window
        SET window_coding = ?, window_repeats = ?
        WHERE window_id = ?
    };
    my $window_update_sth = $dbh->prepare($window_update_query);

    # for each alignment
    for my $align_id (@align_ids) {
        $align_seq_sth->execute($align_id);
        my @row2 = $align_seq_sth->fetchrow_array;

        my ( $chr_name, $align_length, $chr_start, $chr_end, $target_runlist,
            $target_seq )
            = @row2;
        next if $chr_name =~ /rand|un|contig|hap|scaf/i;

        print "Prosess align $align_id in $chr_name $chr_start - $chr_end\n";

        $chr_name =~ s/chr0?//i;

        # align_position to chr_position transforming array
        # the first element will be ignored
        my @chr_pos     = 0 .. $align_length;
        my $indel_count = 0;
        for ( my $i = 1; $i <= $align_length; $i++ ) {
            my $current_base = substr( $target_seq, $i - 1, 1 );
            if ( $current_base eq '-' ) {
                $indel_count++;
            }
            $chr_pos[$i] = $chr_start + $i - 1 - $indel_count;
        }

        # make a new ensembl slice object
        $ensembl->set_slice( $chr_name, $chr_start, $chr_end );

        # for each window
        {
            $window_query_sth->execute($align_id);
            while ( my @row5 = $window_query_sth->fetchrow_array ) {
                my ( $window_id, $window_runlist ) = @row5;
                my $window_set     = AlignDB::IntSpan->new($window_runlist);
                my $window_chr_set = $window_set->map_set( $chr_pos[$_] );

                my $window_coding = $ensembl->feature_portion( '_cds_set',
                    $window_chr_set );
                my $window_repeats
                    = $ensembl->feature_portion( '_repeat_set',
                    $window_chr_set );
                $window_update_sth->execute( $window_coding,
                    $window_repeats, $window_id, );
            }

            $window_update_sth->finish;
            $window_query_sth->finish;
        }

    }
}

$stopwatch->end_message;

# store program running meta info to database
END {
    $obj->add_meta_stopwatch($stopwatch);
}
exit;

__END__

=head1 NAME

    update_window_feature.pl - Add annotation info to alignDB for null windows

=head1 SYNOPSIS

    update_window_feature.pl [options]
      Options:
        --help            brief help message
        --man             full documentation
        --server          MySQL server IP/Domain name
        --db              database name
        --username        username
        --password        password
        --ensembl         ensembl database name

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
