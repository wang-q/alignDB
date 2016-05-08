#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::WriteExcel;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# stat parameter
my $run     = $Config->{stat}{run};
my $outfile = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=s'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'o|output=s'   => \$outfile,
    'r|run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.three.xlsx" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
}
else {
    $run =~ s/\"\'//s;
    my $set = AlignDB::IntSpan->new();
    if ( AlignDB::IntSpan->valid($run) ) {
        $set   = $set->add($run);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $run;
        $set->add(@tasks);
    }

    my $runlist = $set->runlist();
    $outfile =~ s/(\.xlsx)$/.$runlist$1/;
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Do stat for $db...");

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

my $lib = "$FindBin::Bin/sql.lib";
my $sql_file = AlignDB::SQL::Library->new( lib => $lib );

#----------------------------------------------------------#
# worksheet -- snp_distance
#----------------------------------------------------------#
#
my $snp_distance = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ssw', 'ssw_id' ) ) {
        return;
    }

    my @snp_levels = (
        [ 'all',  0,  999 ],
        [ '-1',   -1, -1 ],
        [ '0-5',  0,  5 ],
        [ '6-10', 6,  10 ],

        #[ '0-0',    0,  0 ],
        #[ '1-1',    1,  1 ],
        #[ '2-2',    2,  2 ],
        #[ '3-3',    3,  3 ],
        #[ '4-4',    4,  4 ],
        #[ '5-5',    5,  5 ],
        #[ '6-6',    6,  6 ],
        #[ '7-7',    7,  7 ],
        #[ '8-8',    8,  8 ],
        #[ '9-9',    9,  9 ],
        #[ '10-10',  10, 10 ],
        [ '11-20',  11, 20 ],
        [ '21-30',  21, 30 ],
        [ '31-40',  31, 40 ],
        [ '41-50',  41, 50 ],
        [ '51-999', 51, 999 ],
        [ '21-999', 21, 999 ],
        [ '31-999', 31, 999 ],
        [ '41-999', 41, 999 ],
    );

    my $write_sheet = sub {
        my ($snp_levels) = @_;
        my $sheet_name = 'snp_distance' . "_" . $snp_levels->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q~
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_snp', 'AVG_d_nosnp', 'AVG_d_complex',
                        'COUNT', 'Ds/Dns'
            ~;
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q~
                # distance effect
                SELECT  ssw_distance `distance`,
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM ssw s, window w, snp, isw i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                AND i.isw_distance BETWEEN ? AND ?
                GROUP BY ssw_distance
            ~;
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $snp_levels->[1], $snp_levels->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@snp_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- snp_LR_distance_
#----------------------------------------------------------#
#
my $snp_LR_distance = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ssw', 'ssw_id' ) ) {
        return;
    }

    my @snp_levels = (
        [ 'R', 0, 0 ],
        [ 'L', 0, 0 ],
        [ 'R', 1, 1 ],
        [ 'L', 1, 1 ],
        [ 'R', 2, 2 ],
        [ 'L', 2, 2 ],
        [ 'R', 3, 3 ],
        [ 'L', 3, 3 ],
        [ 'R', 4, 4 ],
        [ 'L', 4, 4 ],
        [ 'R', 5, 5 ],
        [ 'L', 5, 5 ],
    );

    my $write_sheet = sub {
        my ($snp_levels) = @_;
        my $sheet_name = 'snp_distance' . "_" . join "-", @$snp_levels;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q~
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_snp', 'AVG_d_nosnp', 'AVG_d_complex',
                        'COUNT', 'Ds/Dns'
            ~;
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query_R = q~
                # distance effect
                SELECT  CONCAT(s.ssw_type, s.ssw_distance) `distance`,
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM ssw s, window w, snp, isw i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                AND i.isw_type = ?
                AND i.isw_distance BETWEEN ? AND ?
                AND s.ssw_type = 'L'
                GROUP BY CONCAT(s.ssw_type, s.ssw_distance)
                ORDER BY s.ssw_distance DESC
            ~;
            my $sql_query_L = q~
                # distance effect
                SELECT  CONCAT(s.ssw_type, s.ssw_distance) `distance`,
                        AVG(w.window_pi) `AVG_pi`,
                        AVG(ssw_d_snp) `AVG_d_snp`,
                        AVG(ssw_d_nosnp) `AVG_d_nosnp`,
                        AVG(ssw_d_complex) `AVG_d_complex`,
                        COUNT(*) `COUNT`,
                        AVG(ssw_d_snp) / AVG(ssw_d_nosnp)  `Ds/Dns`
                FROM ssw s, window w, snp, isw i
                WHERE s.window_id = w.window_id
                AND snp.snp_id = s.snp_id
                AND snp.isw_id = i.isw_id
                AND i.isw_type = ?
                AND i.isw_distance BETWEEN ? AND ?
                AND s.ssw_type = 'R'
                GROUP BY CONCAT(s.ssw_type, s.ssw_distance)
                ORDER BY s.ssw_distance 
            ~;

            $sheet_row++;    # add a blank line
            my %option = (
                sql_query  => $sql_query_R,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [@$snp_levels],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );

            $sheet_row++;    # add a blank line
            %option = (
                sql_query  => $sql_query_L,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [@$snp_levels],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@snp_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- snp_indel
#----------------------------------------------------------#
#
my $snp_indel = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ssw', 'ssw_id' ) ) {
        return;
    }

    my $sheet_name = 'snp_indel';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q~
            # header of Table distance
            SELECT 'distance', 'AVG_indel', 'COUNT', 'STD_indel'
        ~;
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q~
            # distance effect
            SELECT AVG(s.ssw_distance) distance,
                   AVG(w.window_indel) AVG_indel,
                   COUNT(w.window_indel) COUNT,
                   STD(w.window_indel) STD_indel
            FROM ssw s, window w
            WHERE s.window_id = w.window_id
            GROUP BY s.ssw_distance
        ~;
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write footer
    {
        $sheet_row += 2;
        my $sql_query = q~
            # distance effect
            SELECT AVG(s.ssw_distance) distance,
                   AVG(w.window_indel) AVG_indel,
                   COUNT(w.window_indel) COUNT,
                   STD(w.window_indel) STD_indel
            FROM ssw s, window w
            WHERE s.window_id = w.window_id
        ~;
        my %option = (
            sql_query      => $sql_query,
            sheet_row      => $sheet_row,
            sheet_col      => $sheet_col,
            content_format => 'TOTAL',
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- ds_dns
#----------------------------------------------------------#
my $ds_dns = sub {
    my $sheet_name = 'ds_dns';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => [qw{distance AVG_d_snp AVG_d_nosnp COUNT Ds/Dns}],
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # init objects
    my $dbh  = $write_obj->dbh;
    my $fmt  = $write_obj->format;
    my @cols = @{ $write_obj->columns };

    my $sql_query1 = qq~
        SELECT  i.isw_id, i.indel_id,
                i.isw_distance, i.isw_density,
                (floor((i.isw_density + 1) / 2) - i.isw_distance) reverse_distance,
                s.snp_occured, COUNT(s.snp_id) * 1.0000 / i.isw_length d
        FROM isw i
        LEFT JOIN snp s ON s.isw_id = i.isw_id
        GROUP BY i.isw_id, s.snp_occured
        HAVING ((s.snp_occured != 'N'
        OR s.snp_occured IS NULL)
        AND i.isw_distance > 30
        AND i.isw_density > 70)
    ~;

    my $sth = $dbh->prepare($sql_query1);
    $sth->execute();

    my $each_interval = {};
    while ( my @row = $sth->fetchrow_array ) {
        unless ( $row[5] ) {
            next;
        }
        $each_interval->{ $row[1] }{ $row[4] }{ $row[5] } += $row[-1];
    }

    my $sql_query2 = qq~
        SELECT  distinct i.isw_id, i.indel_id,
                i.isw_distance, i.isw_density,
                (floor((i.isw_density + 1) / 2) - i.isw_distance) reverse_distance
        FROM isw i
        LEFT JOIN snp s ON s.isw_id = i.isw_id
        GROUP BY i.isw_id, s.snp_occured
        HAVING i.isw_distance > 30
        AND i.isw_density > 70
    ~;

    $sth = $dbh->prepare($sql_query2);
    $sth->execute();

    while ( my @row = $sth->fetchrow_array ) {
        $each_interval->{ $row[1] }{ $row[4] }{count}++;
    }

    foreach my $interval ( keys %$each_interval ) {
        foreach my $rd ( keys %{ $each_interval->{$interval} } ) {
            if ( exists $each_interval->{$interval}{$rd}{T} ) {
                $each_interval->{$interval}{$rd}{T}
                    = $each_interval->{$interval}{$rd}{T} / $each_interval->{$interval}{$rd}{count};
            }
            else {
                $each_interval->{$interval}{$rd}{T} = 0;
            }
            if ( exists $each_interval->{$interval}{$rd}{Q} ) {
                $each_interval->{$interval}{$rd}{Q}
                    = $each_interval->{$interval}{$rd}{Q} / $each_interval->{$interval}{$rd}{count};
            }
            else {
                $each_interval->{$interval}{$rd}{Q} = 0;
            }
            delete $each_interval->{$interval}{$rd}{count};
        }
    }

    # Target SNP and Query SNP probability
    my $interval_tq = {};
    foreach my $interval ( keys %$each_interval ) {
        foreach my $rd ( keys %{ $each_interval->{$interval} } ) {
            $interval_tq->{$interval}{T} += $each_interval->{$interval}{$rd}{T};
            $interval_tq->{$interval}{Q} += $each_interval->{$interval}{$rd}{Q};
        }
    }

    my $rd_count = {};
    foreach ( 1 .. 1000 ) {
        srand( time ^ $$ );
        foreach my $interval ( keys %$each_interval ) {
            foreach my $rd ( keys %{ $each_interval->{$interval} } ) {
                my $t_frq     = $interval_tq->{$interval}{T};
                my $q_frq     = $interval_tq->{$interval}{Q};
                my $total_frq = $t_frq + $q_frq;
                unless ($total_frq) {
                    next;
                }
                my ( $Ds, $Dns );
                if ( $t_frq == 0 ) {
                    ( $Ds, $Dns ) = qw/Q T/;
                }
                elsif ( $q_frq == 0 ) {
                    ( $Ds, $Dns ) = qw/T Q/;
                }
                else {
                    my $random = rand($total_frq);
                    $Ds  = $random < $t_frq ? 'T' : 'Q';
                    $Dns = $random < $t_frq ? 'Q' : 'T';
                }
                $rd_count->{$rd}{Ds}  += $each_interval->{$interval}{$rd}{$Ds};
                $rd_count->{$rd}{Dns} += $each_interval->{$interval}{$rd}{$Dns};
                $rd_count->{$rd}{count}++;
            }
        }
    }

    foreach my $rd ( sort { $a <=> $b } keys %{$rd_count} ) {
        $rd_count->{$rd}{Ds}  = $rd_count->{$rd}{Ds} / $rd_count->{$rd}{count};
        $rd_count->{$rd}{Dns} = $rd_count->{$rd}{Dns} / $rd_count->{$rd}{count};
        $sheet->write( $sheet_row, 0 + $sheet_col,     $rd,                     $fmt->{NORMAL} );
        $sheet->write( $sheet_row, 0 + 1 + $sheet_col, $rd_count->{$rd}{Ds},    $fmt->{NORMAL} );
        $sheet->write( $sheet_row, 0 + 2 + $sheet_col, $rd_count->{$rd}{Dns},   $fmt->{NORMAL} );
        $sheet->write( $sheet_row, 0 + 3 + $sheet_col, $rd_count->{$rd}{count}, $fmt->{NORMAL} );
        $sheet->write(
            $sheet_row,
            0 + 4 + $sheet_col,
            $rd_count->{$rd}{Dns}
            ? $rd_count->{$rd}{Ds} / $rd_count->{$rd}{Dns}
            : '#N/A',
            $fmt->{NORMAL}
        );
        $sheet_row++;
    }

    print "Sheet \"$sheet_name\" has been generated.\n";

    #DumpFile("interval.yaml", $each_interval);
    #DumpFile("interval_tq.yaml", $interval_tq);
    #DumpFile("rd.yaml", $rd_count);

};

for my $n (@tasks) {
    if ( $n == 4 ) { &$snp_distance; &$snp_indel; next; }
    if ( $n == 52 ) { &$ds_dns;          next; }
    if ( $n == 56 ) { &$snp_LR_distance; next; }
}

$stopwatch->end_message();
exit;

__END__
