#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Crypt::Random qw( makerandom );

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::WriteExcel;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# stat parameter
my $run           = $Config->{stat}->{run};
my $sum_threshold = $Config->{stat}->{sum_threshold};
my $outfile       = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'server=s'    => \$server,
    'port=s'      => \$port,
    'db=s'        => \$db,
    'username=s'  => \$username,
    'password=s'  => \$password,
    'output=s'    => \$outfile,
    'run=s'       => \$run,
    'threshold=i' => \$sum_threshold,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.four.xlsx" unless $outfile;

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

#----------------------------------------------------------#
# worksheet -- indel_occured
#----------------------------------------------------------#
my $indel_occured = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my $sheet_name = 'indel_occured';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'Item';
        my $sql_query  = q{
            # header of Table summray
            SELECT  'Indel_occured', 'slippage-like',
                    'AVG_lenght', 'STD_length', 'COUNT'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $query_name = 'Original';
        my $sql_query  = q{
            SELECT  i.indel_occured,
                    e.indel_feature3,
                    AVG(i.indel_length) AVG_length,
                    STD(i.indel_length) STD_length,
                    COUNT(*) COUNT
            FROM indel i, indel_extra e
            WHERE i.indel_id = e.indel_id
            GROUP BY i.indel_occured, e.indel_feature3
            ORDER BY i.indel_occured DESC
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'Another_ref';
        my $sql_query  = q{
            SELECT  e.indel_feature4,
                    e.indel_feature3,
                    AVG(i.indel_length) AVG_length,
                    STD(i.indel_length) STD_length,
                    COUNT(*) COUNT
            FROM indel i, indel_extra e
            WHERE i.indel_id = e.indel_id
            GROUP BY e.indel_feature4, e.indel_feature3
            ORDER BY e.indel_feature4 DESC
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $sql_query  = q{
            SELECT CONCAT(i.indel_occured, "->", e.indel_feature4) `change`,
                   e.indel_feature3 `slippage-like`,
                   AVG (i.indel_length) AVG_length,
                   STD(i.indel_length) STD_length,
                   COUNT(*) COUNT
            FROM indel i,
                 indel_extra e
            WHERE i.indel_id = e.indel_id AND
                  e.indel_feature4 IS NOT NULL AND
                  i.indel_occured = e.indel_feature4
            GROUP BY `change`,
                     e.indel_feature3
            ORDER BY `change` DESC
        };
        my $query_name = 'Occured_not_change';
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        
    }
    
    {
        
        $sheet_row++;
        my $sql_query  = q{
            SELECT CONCAT(i.indel_occured, "->", e.indel_feature4) `change`,
                   e.indel_feature3 `slippage-like`,
                   AVG (i.indel_length) AVG_length,
                   STD(i.indel_length) STD_length,
                   COUNT(*) COUNT
            FROM indel i,
                 indel_extra e
            WHERE i.indel_id = e.indel_id AND
                  e.indel_feature4 IS NOT NULL AND
                  i.indel_occured != e.indel_feature4
            GROUP BY `change`,
                     e.indel_feature3
            ORDER BY `change` DESC
        };
        my $query_name = 'Occured_change';
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'Occured_vanish';
        my $sql_query  = q{
            SELECT CONCAT(i.indel_occured, "->") `vanish`,
                   e.indel_feature3 `slippage-like`,
                   AVG (i.indel_length) AVG_length,
                   STD(i.indel_length) STD_length,
                   COUNT(*) COUNT
            FROM indel i,
                 indel_extra e
            WHERE i.indel_id = e.indel_id AND
                  e.indel_feature4 IS NULL
            GROUP BY `vanish`,
                     e.indel_feature3
            ORDER BY `vanish` DESC
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_occured
#----------------------------------------------------------#
my $snp_occured = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'snp_extra', 'snp_feature4' ) ) {
        return;
    }

    my $sheet_name = 'snp_occured';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $query_name = 'Item';
        my $sql_query  = q{
            # header of Table
            SELECT 'SNP_occured', 'COUNT'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $query_name = 'Original';
        my $sql_query  = q{
            SELECT  s.snp_occured,
                    COUNT(*) COUNT
            FROM snp s
            GROUP BY s.snp_occured
            ORDER BY s.snp_occured DESC
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'Another_ref';
        my $sql_query  = q{
            SELECT  e.snp_feature4,
                    COUNT(*) COUNT
            FROM snp_extra e
            GROUP BY e.snp_feature4
            ORDER BY e.snp_feature4 DESC
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $sql_query  = q{
            SELECT CONCAT(s.snp_occured, "->", e.snp_feature4) `change`,
                   COUNT(*) COUNT
            FROM snp s,
                 snp_extra e
            WHERE s.snp_id = e.snp_id AND
                  e.snp_feature4 IS NOT NULL AND
                  s.snp_occured = e.snp_feature4
            GROUP BY `change`
            ORDER BY `change` DESC
        };
        my $query_name = 'Occured_not_change';
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }
    
    {
        
        $sheet_row++;
        my $sql_query  = q{
            SELECT CONCAT(s.snp_occured, "->", e.snp_feature4) `change`,
                   COUNT(*) COUNT
            FROM snp s,
                 snp_extra e
            WHERE s.snp_id = e.snp_id AND
                  e.snp_feature4 IS NOT NULL AND
                  s.snp_occured != e.snp_feature4
            GROUP BY `change`
            ORDER BY `change` DESC
        };
        my $query_name = 'Occured_change';
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # write contents
    {
        $sheet_row++;
        my $query_name = 'Occured_vanish';
        my $sql_query  = q{
            SELECT CONCAT(s.snp_occured, "->") `vanish`,
                   COUNT(*) COUNT
            FROM snp s,
                 snp_extra e
            WHERE s.snp_id = e.snp_id AND
                  e.snp_feature4 IS NULL
            GROUP BY `vanish`
            ORDER BY `vanish` DESC
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_size_occured
#----------------------------------------------------------#
my $indel_size_occured = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my @size_levels = ( [ 1, 50 ],   );

    my $write_sheet = sub {
        my $level      = shift;
        my $lower      = $level->[0];
        my $upper      = $level->[1];
        my $sheet_name = "indel_occured_$lower-$upper";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $query_name = 'Item';
            my $sql_query  = q{
                # header of Table summray
                SELECT  'Indel_occured', 'slippage-like',
                        'AVG_lenght', 'STD_length', 'COUNT'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                query_name => $query_name,
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $query_name = 'Original';
            my $sql_query  = q{
                SELECT i.indel_occured,
                       e.indel_feature3,
                       AVG (i.indel_length) AVG_length,
                       STD(i.indel_length) STD_length,
                       COUNT(*) COUNT
                FROM indel i,
                     indel_extra e
                WHERE i.indel_id = e.indel_id AND
                      i.indel_length >= ? AND
                      i.indel_length <= ?
                GROUP BY i.indel_occured,
                         e.indel_feature3
                ORDER BY i.indel_occured DESC
            };
            my %option = (
                query_name => $query_name,
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        # write contents
        {
            $sheet_row++;
            my $query_name = 'Another_ref';
            my $sql_query  = q{
                SELECT e.indel_feature4,
                       e.indel_feature3,
                       AVG (i.indel_length) AVG_length,
                       STD(i.indel_length) STD_length,
                       COUNT(*) COUNT
                FROM indel i,
                     indel_extra e
                WHERE i.indel_id = e.indel_id AND
                      i.indel_length >= ? AND
                      i.indel_length <= ?
                GROUP BY e.indel_feature4,
                         e.indel_feature3
                ORDER BY e.indel_feature4 DESC
            };
            my %option = (
                query_name => $query_name,
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        # write contents
        {
            $sheet_row++;
            my $sql_query  = q{
                SELECT CONCAT(i.indel_occured, "->", e.indel_feature4) `change`,
                       e.indel_feature3 `slippage-like`,
                       AVG (i.indel_length) AVG_length,
                       STD(i.indel_length) STD_length,
                       COUNT(*) COUNT
                FROM indel i,
                     indel_extra e
                WHERE i.indel_id = e.indel_id AND
                      e.indel_feature4 IS NOT NULL AND
                      i.indel_occured = e.indel_feature4 AND
                      i.indel_length >= ? AND
                      i.indel_length <= ?
                GROUP BY `change`,
                         e.indel_feature3
                ORDER BY `change` DESC
            };
            my $query_name = 'Occured_not_change';
            my %option = (
                query_name => $query_name,
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
                
        }
        
        {
            $sheet_row++;
            my $sql_query  = q{
                SELECT CONCAT(i.indel_occured, "->", e.indel_feature4) `change`,
                       e.indel_feature3 `slippage-like`,
                       AVG (i.indel_length) AVG_length,
                       STD(i.indel_length) STD_length,
                       COUNT(*) COUNT
                FROM indel i,
                     indel_extra e
                WHERE i.indel_id = e.indel_id AND
                      e.indel_feature4 IS NOT NULL AND
                      i.indel_occured != e.indel_feature4 AND
                      i.indel_length >= ? AND
                      i.indel_length <= ?
                GROUP BY `change`,
                         e.indel_feature3
                ORDER BY `change` DESC
            };
            my $query_name = 'Occured_change';
            my %option = (
                query_name => $query_name,
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        # write contents
        {
            $sheet_row++;
            my $query_name = 'Occured_vanish';
            my $sql_query  = q{
                SELECT CONCAT(i.indel_occured, "->") `vanish`,
                       e.indel_feature3 `slippage-like`,
                       AVG (i.indel_length) AVG_length,
                       STD(i.indel_length) STD_length,
                       COUNT(*) COUNT
                FROM indel i,
                     indel_extra e
                WHERE i.indel_id = e.indel_id AND
                      e.indel_feature4 IS NULL AND
                      i.indel_length >= ? AND
                      i.indel_length <= ?
                GROUP BY `vanish`,
                         e.indel_feature3
                ORDER BY `vanish` DESC
            };
            my %option = (
                query_name => $query_name,
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@size_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_onc_distance
#----------------------------------------------------------#
my $indel_onc_distance = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my $sheet_name = 'indel_onc_distance';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q{
            # header of Table distance
            SELECT  'distance', 'AVG_pi',
                    'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                    'COUNT', 'Di/Dn'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q{
            SELECT isw_distance distance,
                   AVG (isw_pi) AVG_pi,
                   AVG (isw_d_indel) AVG_d_indel,
                   AVG (isw_d_noindel) AVG_d_noindel,
                   AVG (isw_d_complex) AVG_d_complex,
                   COUNT(*) COUNT,
                   AVG (isw_d_indel) / AVG (isw_d_noindel) `Di/Dn`
            FROM isw,
                 (
                  SELECT isw_id
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        e.indel_feature4 IS NOT NULL AND
                        i.indel_occured = e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        isw_type = 'R' AND
                        isw.indel_id = i.indel_id
                  UNION
                  SELECT isw_id
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        e.indel_feature4 IS NOT NULL AND
                        i.indel_occured = e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        isw_type = 'L' AND
                        isw.prev_indel_id = i.indel_id
                 ) i
            WHERE isw_distance >= 0 AND
                  isw_d_indel IS NOT NULL AND
                  isw.isw_id = i.isw_id
            GROUP BY isw_distance
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_oc_distance
#----------------------------------------------------------#
#
my $indel_oc_distance = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my $sheet_name = 'indel_oc_distance';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q{
            # header of Table distance
            SELECT  'distance', 'AVG_pi',
                    'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                    'COUNT', 'Di/Dn'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q{
            SELECT isw_distance distance,
                   AVG (isw_pi) AVG_pi,
                   AVG (isw_d_indel) AVG_d_indel,
                   AVG (isw_d_noindel) AVG_d_noindel,
                   AVG (isw_d_complex) AVG_d_complex,
                   COUNT(*) COUNT,
                   AVG (isw_d_indel) / AVG (isw_d_noindel) `Di/Dn`
            FROM isw,
                 (
                  SELECT isw_id
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        e.indel_feature4 IS NOT NULL AND
                        i.indel_occured != e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        isw_type = 'R' AND
                        isw.indel_id = i.indel_id
                  UNION
                  SELECT isw_id
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        e.indel_feature4 IS NOT NULL AND
                        i.indel_occured != e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        isw_type = 'L' AND
                        isw.prev_indel_id = i.indel_id
                 ) i
            WHERE isw_distance >= 0 AND
                  isw_d_indel IS NOT NULL AND
                  isw.isw_id = i.isw_id
            GROUP BY isw_distance
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }
    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_onc_slip
#----------------------------------------------------------#
#
my $indel_onc_slip = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my @occured_levels = ( [ 'nonslip', 0 ], [ 'slip', 1 ], );

    my $write_sheet = sub {
        my $level      = shift;
        my $sheet_name = 'indel_onc_' . $level->[0];
        my $slip       = $level->[1];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                SELECT isw_distance distance,
                       AVG (isw_pi) AVG_pi,
                       AVG (isw_d_indel) AVG_d_indel,
                       AVG (isw_d_noindel) AVG_d_noindel,
                       AVG (isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG (isw_d_indel) / AVG (isw_d_noindel) `Di/Dn`
                FROM isw,
                     (
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured = e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            e.indel_feature3 = ? AND
                            isw_type = 'R' AND
                            isw.indel_id = i.indel_id
                      UNION
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured = e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            e.indel_feature3 = ? AND
                            isw_type = 'L' AND
                            isw.prev_indel_id = i.indel_id
                     ) i
                WHERE isw_distance >= 0 AND
                      isw_d_indel IS NOT NULL AND
                      isw.isw_id = i.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $slip, $slip ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@occured_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_oc_slip
#----------------------------------------------------------#
#
my $indel_oc_slip = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my @occured_levels = ( [ 'nonslip', 0 ], [ 'slip', 1 ], );

    my $write_sheet = sub {
        my $level      = shift;
        my $sheet_name = 'indel_oc_' . $level->[0];
        my $slip       = $level->[1];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                SELECT isw_distance distance,
                       AVG (isw_pi) AVG_pi,
                       AVG (isw_d_indel) AVG_d_indel,
                       AVG (isw_d_noindel) AVG_d_noindel,
                       AVG (isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG (isw_d_indel) / AVG (isw_d_noindel) `Di/Dn`
                FROM isw,
                     (
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured != e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            e.indel_feature3 = ? AND
                            isw_type = 'R' AND
                            isw.indel_id = i.indel_id
                      UNION
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured != e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            e.indel_feature3 = ? AND
                            isw_type = 'L' AND
                            isw.prev_indel_id = i.indel_id
                     ) i
                WHERE isw_distance >= 0 AND
                      isw_d_indel IS NOT NULL AND
                      isw.isw_id = i.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $slip, $slip ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@occured_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_onc_nonslip_dxr
#----------------------------------------------------------#
#
my $indel_onc_nonslip_dxr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature6' ) ) {
        return;
    }

    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my $sheet_name = 'indel_onc_nonslip_dxr';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q{
            # header of Table distance
            SELECT  'distance', 'AVG_pi',
                    'AVG_d_ir', 'AVG_d_nr',
                    'AVG_d_total',
                    'COUNT'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q{
            SELECT isw_distance distance,
                   AVG (isw_pi) avg_pi,
                   AVG (isw_feature4) avg_d_ir,
                   AVG (isw_feature5) avg_d_nr,
                   AVG (isw_feature8) AVG_d_total,
                   count(*) count
            FROM isw_extra e,
                 (
                  SELECT isw_id,
                         isw_distance,
                         isw_pi
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        i.indel_occured IS NOT NULL AND
                        i.indel_occured = e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        e.indel_feature3 = 0 AND
                        isw_type = 'R' AND
                        isw.indel_id = i.indel_id
                  UNION
                  SELECT isw_id,
                         isw_distance,
                         isw_pi
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        i.indel_occured IS NOT NULL AND
                        i.indel_occured = e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        e.indel_feature3 = 0 AND
                        isw_type = 'L' AND
                        isw.prev_indel_id = i.indel_id
                 ) i
            WHERE i.isw_distance >= 0 AND
                  e.isw_feature6 IS NOT NULL AND
                  i.isw_id = e.isw_id
            GROUP BY isw_distance
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_oc_nonslip_dxr
#----------------------------------------------------------#
#
my $indel_oc_nonslip_dxr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature6' ) ) {
        return;
    }

    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my $sheet_name = 'indel_oc_nonslip_dxr';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # write header
    {
        my $sql_query = q{
            # header of Table distance
            SELECT  'distance', 'AVG_pi',
                    'AVG_d_ir', 'AVG_d_nr',
                    'AVG_d_total',
                    'COUNT'
        };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );
    }

    # write contents
    {
        my $sql_query = q{
            SELECT isw_distance distance,
                   AVG (isw_pi) avg_pi,
                   AVG (isw_feature4) avg_d_ir,
                   AVG (isw_feature5) avg_d_nr,
                   AVG (isw_feature8) AVG_d_total,
                   count(*) count
            FROM isw_extra e,
                 (
                  SELECT isw_id,
                         isw_distance,
                         isw_pi
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        i.indel_occured IS NOT NULL AND
                        i.indel_occured != e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        e.indel_feature3 = 0 AND
                        isw_type = 'R' AND
                        isw.indel_id = i.indel_id
                  UNION
                  SELECT isw_id,
                         isw_distance,
                         isw_pi
                  FROM isw,
                       indel i,
                       indel_extra e
                  WHERE i.indel_id = e.indel_id AND
                        i.indel_occured IS NOT NULL AND
                        i.indel_occured != e.indel_feature4 AND
                        i.indel_occured IN ('T', 'Q') AND
                        e.indel_feature4 IN ('T', 'Q') AND
                        e.indel_feature3 = 0 AND
                        isw_type = 'L' AND
                        isw.prev_indel_id = i.indel_id
                 ) i
            WHERE i.isw_distance >= 0 AND
                  e.isw_feature6 IS NOT NULL AND
                  i.isw_id = e.isw_id
            GROUP BY isw_distance
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_size_onc_nonslip
#----------------------------------------------------------#
#
my $indel_size_onc_nonslip = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my @size_levels = ( [ 1, 50 ],   );

    my $write_sheet = sub {
        my $level      = shift;
        my $lower      = $level->[0];
        my $upper      = $level->[1];
        my $sheet_name = "indel_onc_nonslip_$lower-$upper";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                # distance effect
                SELECT isw_distance distance,
                       AVG (isw_pi) AVG_pi,
                       AVG (isw_d_indel) AVG_d_indel,
                       AVG (isw_d_noindel) AVG_d_noindel,
                       AVG (isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG (isw_d_indel) / AVG (isw_d_noindel) `Di/Dn`
                FROM isw,
                     (
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured = e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'R' AND
                            isw.indel_id = i.indel_id
                      UNION
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured = e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'L' AND
                            isw.prev_indel_id = i.indel_id
                     ) i
                WHERE isw_distance >= 0 AND
                      isw_d_indel IS NOT NULL AND
                      isw.isw_id = i.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper, $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@size_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_size_oc_nonslip
#----------------------------------------------------------#
#
my $indel_size_oc_nonslip = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my @size_levels = ( [ 1, 50 ],   );

    my $write_sheet = sub {
        my $level      = shift;
        my $lower      = $level->[0];
        my $upper      = $level->[1];
        my $sheet_name = "indel_oc_nonslip_$lower-$upper";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                # distance effect
                SELECT isw_distance distance,
                       AVG (isw_pi) AVG_pi,
                       AVG (isw_d_indel) AVG_d_indel,
                       AVG (isw_d_noindel) AVG_d_noindel,
                       AVG (isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG (isw_d_indel) / AVG (isw_d_noindel) `Di/Dn`
                FROM isw,
                     (
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured != e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'R' AND
                            isw.indel_id = i.indel_id
                      UNION
                      SELECT isw_id
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured != e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'L' AND
                            isw.prev_indel_id = i.indel_id
                     ) i
                WHERE isw_distance >= 0 AND
                      isw_d_indel IS NOT NULL AND
                      isw.isw_id = i.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper, $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@size_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_size_onc_nonslip_dxr
#----------------------------------------------------------#
#
my $indel_size_onc_nonslip_dxr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature6' ) ) {
        return;
    }
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my @size_levels = ( [ 1, 50 ],   );

    my $write_sheet = sub {
        my $level      = shift;
        my $lower      = $level->[0];
        my $upper      = $level->[1];
        my $sheet_name = "indel_onc_nonslip_dxr_$lower-$upper";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_ir', 'AVG_d_nr',
                        'AVG_d_total',
                        'COUNT'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                # distance effect
                SELECT isw_distance distance,
                       AVG (isw_pi) avg_pi,
                       AVG (isw_feature4) avg_d_ir,
                       AVG (isw_feature5) avg_d_nr,
                       AVG (isw_feature8) AVG_d_total,
                       count(*) count
                FROM isw_extra e,
                     (
                      SELECT isw_id, isw_distance, isw_pi
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured = e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'R' AND
                            isw.indel_id = i.indel_id
                      UNION
                      SELECT isw_id, isw_distance, isw_pi
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured = e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'L' AND
                            isw.prev_indel_id = i.indel_id
                     ) i
                WHERE i.isw_distance >= 0 AND
                      e.isw_feature6 IS NOT NULL AND
                      i.isw_id = e.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper, $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@size_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_size_oc_nonslip_dxr
#----------------------------------------------------------#
#
my $indel_size_oc_nonslip_dxr = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'isw_extra', 'isw_feature6' ) ) {
        return;
    }
    unless ( $write_obj->check_column( 'indel_extra', 'indel_feature4' ) ) {
        return;
    }

    my @size_levels = ( [ 1, 50 ],   );

    my $write_sheet = sub {
        my $level      = shift;
        my $lower      = $level->[0];
        my $upper      = $level->[1];
        my $sheet_name = "indel_oc_nonslip_dxr_$lower-$upper";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_ir', 'AVG_d_nr',
                        'AVG_d_total',
                        'COUNT'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                # distance effect
                SELECT isw_distance distance,
                       AVG (isw_pi) avg_pi,
                       AVG (isw_feature4) avg_d_ir,
                       AVG (isw_feature5) avg_d_nr,
                       AVG (isw_feature8) AVG_d_total,
                       count(*) count
                FROM isw_extra e,
                     (
                      SELECT isw_id, isw_distance, isw_pi
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured != e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'R' AND
                            isw.indel_id = i.indel_id
                      UNION
                      SELECT isw_id, isw_distance, isw_pi
                      FROM isw,
                           indel i,
                           indel_extra e
                      WHERE i.indel_id = e.indel_id AND
                            i.indel_occured IS NOT NULL AND
                            i.indel_occured != e.indel_feature4 AND
                            i.indel_occured IN ('T', 'Q') AND
                            e.indel_feature4 IN ('T', 'Q') AND
                            i.indel_length BETWEEN ? AND ? AND
                            e.indel_feature3 = 0 AND
                            isw_type = 'L' AND
                            isw.prev_indel_id = i.indel_id
                     ) i
                WHERE i.isw_distance >= 0 AND
                      e.isw_feature6 IS NOT NULL AND
                      i.isw_id = e.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper, $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@size_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_size_didn
#----------------------------------------------------------#
#
my $indel_size_didn = sub {

    my @size_levels = ( [ 1, 50 ],   );

    my $write_sheet = sub {
        my $level      = shift;
        my $lower      = $level->[0];
        my $upper      = $level->[1];
        my $sheet_name = "indel_size_didn_$lower-$upper";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $sql_query = q{
                # header of Table distance
                SELECT  'distance', 'AVG_pi',
                        'AVG_d_indel', 'AVG_d_noindel', 'AVG_d_complex',
                        'COUNT', 'Di/Dn'
            };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_sql( $sheet_name, \%option );
        }

        # write contents
        {
            my $sql_query = q{
                SELECT isw_distance distance,
                       AVG (isw_pi) AVG_pi,
                       AVG (isw_d_indel) AVG_d_indel,
                       AVG (isw_d_noindel) AVG_d_noindel,
                       AVG (isw_d_complex) AVG_d_complex,
                       COUNT(*) COUNT,
                       AVG (isw_d_indel) / AVG (isw_d_noindel) `Di/Dn`
                FROM isw,
                     (
                      SELECT isw_id
                      FROM isw,
                           indel i
                      WHERE i.indel_length >= ? AND
                            indel_length <= ? AND
                            isw_type = 'R' AND
                            isw.indel_id = i.indel_id
                      UNION
                      SELECT isw_id
                      FROM isw,
                           indel i
                      WHERE i.indel_length >= ? AND
                            indel_length <= ? AND
                            isw_type = 'L' AND
                            isw.prev_indel_id = i.indel_id
                     ) i
                WHERE isw_distance >= 0 AND
                      isw_d_indel IS NOT NULL AND
                      isw.isw_id = i.isw_id
                GROUP BY isw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $lower, $upper, $lower, $upper ],
            );
            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@size_levels) {
        &$write_sheet($_);
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$indel_occured;      &$snp_occured;       next; }
    if ( $n == 2 ) { &$indel_size_occured; next; }
    if ( $n == 3 ) { &$indel_size_didn; next; }
    if ( $n == 4 ) { &$indel_onc_distance; &$indel_oc_distance; next; }
    if ( $n == 5 ) { &$indel_onc_slip;     &$indel_oc_slip;     next; }
    if ( $n == 6 ) {
        &$indel_size_onc_nonslip;
        &$indel_size_oc_nonslip;
        next;
    }
    if ( $n == 10 ) {
        &$indel_onc_nonslip_dxr;
        &$indel_oc_nonslip_dxr;
        next;
    }
    if ( $n == 12 ) {
        &$indel_size_onc_nonslip_dxr;
        &$indel_size_oc_nonslip_dxr;
        next;
    }
}

$stopwatch->end_message();
exit;

__END__

=head1 NAME

    four_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    three_stat_factory.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --output          output filename
       --run             run special analysis
       

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
