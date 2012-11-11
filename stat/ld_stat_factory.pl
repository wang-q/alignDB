#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::WriteExcel;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# stat parameter
my $run               = $Config->{stat}{run};
my $combine_threshold = $Config->{stat}{combine_threshold};
my $outfile           = "";

my $max_freq;    # count freq one by one to $max_freq

my $help = 0;
my $man  = 0;

GetOptions(
    'help|?'                 => \$help,
    'man'                    => \$man,
    's|server=s'             => \$server,
    'P|port=s'               => \$port,
    'd|db=s'                 => \$db,
    'u|username=s'           => \$username,
    'p|password=s'           => \$password,
    'o|output=s'             => \$outfile,
    'max|max_freq=s'         => \$max_freq,
    'r|run=s'                => \$run,
    'ct|combine_threshold=i' => \$combine_threshold,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 20 );
    $outfile = "$db.multi.xlsx" unless $outfile;
}
else {
    $run =~ s/\"\'//s;
    my $set = AlignDB::IntSpan->new;
    if ( AlignDB::IntSpan->valid($run) ) {
        $set   = $set->add($run);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $run;
        $set->add(@tasks);
    }

    unless ($outfile) {
        my $runlist = $set->runlist;
        $outfile = "$db.multi.$runlist.xlsx";
    }
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Do stat for $db...");

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

my $lib = "$FindBin::Bin/sql.lib";
my $sql_file = AlignDB::SQL::Library->new( lib => $lib );

# auto detect combine threshold
if ( $combine_threshold == 0 ) {
    my $dbh = $write_obj->dbh;

    my $sql_query = q{
        SELECT SUM(align_length)
        FROM align
    };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute;
    my ($total_length) = $sth->fetchrow_array;

    if ( $total_length <= 1_000_000 ) {
        $combine_threshold = 100;
    }
    elsif ( $total_length <= 10_000_000 ) {
        $combine_threshold = 500;
    }
    else {
        $combine_threshold = 1000;
    }
}

#----------------------------#
# count freq
#----------------------------#
my $all_freq;
{
    my $dbh = $write_obj->dbh;

    my $sql_query = q{
        SELECT DISTINCT COUNT(q.query_id) + 1
        FROM  query q, sequence s
        WHERE q.seq_id = s.seq_id
        GROUP BY s.align_id
    };
    my $sth = $dbh->prepare($sql_query);

    my @counts;
    $sth->execute;
    while ( my ($count) = $sth->fetchrow_array ) {
        push @counts, $count;
    }
    if ( scalar @counts > 1 ) {
        die "Database corrupts, freqs are not consistent\n";
    }

    $all_freq = $counts[0];
}

if ( $all_freq < 3 ) {
    die "all_freq is $all_freq, are you sure this is a AlignDB::Multi DB?\n";
}

my @freqs;
if ($max_freq) {
    for ( 1 .. $max_freq - 1 ) {
        my $name = $_ . "of" . $all_freq;
        push @freqs, [ $name, $_, $_ ];
    }
}
else {

    # for 22 flies, 1, 2, low, mid, high, 20, 21
    {
        my @all_freqs = 1 .. $all_freq - 1;
        if ( scalar @all_freqs <= 7 ) {
            for (@all_freqs) {
                my $name = $_ . "of" . $all_freq;
                push @freqs, [ $name, $_, $_ ];
            }
        }
        else {
            for ( 1, 2 ) {
                my $name = $_ . "of" . $all_freq;
                push @freqs, [ $name, $_, $_ ];
            }

            my @to_be_combs = @all_freqs[ 2 .. $all_freq - 4 ];
            my @chunks      = reverse apportion( scalar @to_be_combs, 3 );
            my @chunks_freq = multi_slice( \@to_be_combs, @chunks );
            for my $chunk (@chunks_freq) {
                my $name = join( '_', @{$chunk} ) . "of" . $all_freq;
                push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
            }

            for ( $all_freq - 2, $all_freq - 1 ) {
                my $name = $_ . "of" . $all_freq;
                push @freqs, [ $name, $_, $_ ];
            }
        }
    }
}

    #----------------------------------------------------------#
    # worksheet -- distance_ld
    #----------------------------------------------------------#
my $ld = sub {

    {
        my $sheet_name = 'ld';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_r AVG_r2 AVG_Dprime AVG_Dprime_abs
                COUNT};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write contents
            my $sql_query = q{
                SELECT
                    w.isw_distance distance,
                    AVG(s.snp_r) AVG_r,
                    AVG(POWER(s.snp_r, 2)) AVG_r2,
                    AVG(s.snp_dprime) AVG_Dprime,
                    AVG(ABS(s.snp_dprime)) AVG_Dprime_abs,
                    COUNT(*) COUNT
                FROM indel i, isw w, snp s
                WHERE 1 = 1
                AND i.indel_id = w.isw_indel_id
                AND w.isw_id = s.isw_id
                AND i.indel_freq != 'unknown'
                AND s.snp_freq != 'unknown'
                GROUP BY w.isw_distance 
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

};

my $ld_insdel = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'ld_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_r AVG_r2 AVG_Dprime AVG_Dprime_abs
                COUNT};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write contents
            my $sql_query = q{
                SELECT
                    w.isw_distance distance,
                    AVG(s.snp_r) AVG_r,
                    AVG(POWER(s.snp_r, 2)) AVG_r2,
                    AVG(s.snp_dprime) AVG_Dprime,
                    AVG(ABS(s.snp_dprime)) AVG_Dprime_abs,
                    COUNT(*) COUNT
                FROM indel i, isw w, snp s
                WHERE 1 = 1
                AND i.indel_id = w.isw_indel_id
                AND w.isw_id = s.isw_id
                AND i.indel_freq != 'unknown'
                AND s.snp_freq != 'unknown'
                AND i.indel_type = ?
                GROUP BY w.isw_distance 
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $level->[1] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@type_levels) {
        &$write_sheet($_);
    }
};

my $ld_freq = sub {
    my @freq_levels = @freqs;

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'ld_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_r AVG_r2 AVG_Dprime AVG_Dprime_abs
                COUNT};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write contents
            my $sql_query = q{
                SELECT
                    w.isw_distance distance,
                    AVG(s.snp_r) AVG_r,
                    AVG(POWER(s.snp_r, 2)) AVG_r2,
                    AVG(s.snp_dprime) AVG_Dprime,
                    AVG(ABS(s.snp_dprime)) AVG_Dprime_abs,
                    COUNT(*) COUNT
                FROM indel i, isw w, snp s
                WHERE 1 = 1
                AND i.indel_id = w.isw_indel_id
                AND w.isw_id = s.isw_id
                AND i.indel_freq != 'unknown'
                AND s.snp_freq != 'unknown'
                AND i.indel_freq >= ?
                AND i.indel_freq <= ?
                GROUP BY w.isw_distance 
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $level->[1], $level->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@freq_levels) {
        &$write_sheet($_);
    }
};

my $ld_insdel_freq = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );
    my @freq_levels = @freqs;

    my $write_sheet = sub {
        my ( $type, $freq ) = @_;
        my $sheet_name = 'ld_' . $type->[0] . '_' . $freq->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_r AVG_r2 AVG_Dprime AVG_Dprime_abs
                COUNT};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write contents
            my $sql_query = q{
                SELECT
                    w.isw_distance distance,
                    AVG(s.snp_r) AVG_r,
                    AVG(POWER(s.snp_r, 2)) AVG_r2,
                    AVG(s.snp_dprime) AVG_Dprime,
                    AVG(ABS(s.snp_dprime)) AVG_Dprime_abs,
                    COUNT(*) COUNT
                FROM indel i, isw w, snp s
                WHERE 1 = 1
                AND i.indel_id = w.isw_indel_id
                AND w.isw_id = s.isw_id
                AND i.indel_freq != 'unknown'
                AND s.snp_freq != 'unknown'
                AND i.indel_type = ?
                AND i.indel_freq >= ?
                AND i.indel_freq <= ?
                GROUP BY w.isw_distance 
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $type->[1], $freq->[1], $freq->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for my $type (@type_levels) {
        for my $freq (@freq_levels) {
            &$write_sheet( $type, $freq );
        }
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$ld;             next; }
    if ( $n == 2 ) { &$ld_insdel;      next; }
    if ( $n == 3 ) { &$ld_freq;        next; }
    if ( $n == 4 ) { &$ld_insdel_freq; next; }
}

$stopwatch->end_message;
exit;

# codes come from http://www.perlmonks.org/?node_id=516493
sub apportion {
    my ( $elements, $pieces ) = @_;
    my $small_chunk     = int $elements / $pieces;
    my $oversized_count = $elements % $pieces;
    (   ( 1 + $small_chunk ) x ($oversized_count),
        ($small_chunk) x ( $pieces - $oversized_count )
    );
}

sub multi_slice {
    my ( $aref, @chunk_sizes ) = @_;
    my $hi_i = -1;
    map {
        my $lo_i = $hi_i + 1;
        $hi_i += $_;
        [ @$aref[ $lo_i .. $hi_i ] ]
    } @chunk_sizes;
}
__END__

=head1 NAME

    multi_stat_factory.pl - Generate statistical Excel files from malignDB

=head1 SYNOPSIS

    multi_stat_factory.pl [options]
        Options:
            --help              brief help message
            --man               full documentation
            --server            MySQL server IP/Domain name
            --db                database name
            --username          username
            --password          password
            --outfile            outfile filename
            --run               run special analysis
            --freq              max freq

=cut

