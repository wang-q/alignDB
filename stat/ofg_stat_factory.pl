#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
use Set::Scalar;
use List::MoreUtils qw(first_index);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::ToXLSX;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

ofg_stat_factory.pl - OFG (other features of genome) stats for alignDB

=head1 SYNOPSIS

    perl ofg_stat_factory.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --outfile   -o  STR     outfile filename
        --by            STR     tag, type or tt
        --run       -r  STR     run special analysis
        --replace       STR=STR replace strings in axis names
        --index                 add an index sheet
        --chart                 add charts

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'output|o=s'   => \( my $outfile ),
    'by=s'         => \( my $by       = "tag" ),
    'run|r=s'      => \( my $run      = $Config->{stat}{run} ),
    'replace=s'    => \my %replace,
    'index'        => \my $add_index_sheet,
    'chart'        => \my $add_chart,
) or Getopt::Long::HelpMessage(1);

$outfile = "$db.ofg.xlsx" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 20 );
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

    my $runlist = $set->runlist;
    $outfile =~ s/(\.xlsx)$/.$runlist$1/;
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Do stat for $db...");

my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password )
    or die "Cannot connect to MySQL database at $db:$server";
my $write_obj = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
    replace => \%replace,
);

#----------------------------------------------------------#
# chart -- ofg
#----------------------------------------------------------#
my $chart_ofg = sub {
    my $sheet = shift;
    my $data  = shift;
    my $is_ld = shift;

    my $x_set = Set::Scalar->new( @{ $data->[0] } );

    my %opt = (
        x_column => 0,
        x_title  => "Distance to ofg",
        top      => 1,
        left     => 13,
    );

    if ( $x_set->has(-10) ) {
        my $idx = first_index { $_ == -10 } @{ $data->[0] };
        $opt{first_row}   = $idx + 1;
        $opt{last_row}    = $idx + 26;
        $opt{x_min_scale} = -10;
        $opt{x_max_scale} = 15;
        $opt{cross}       = 0;
    }
    elsif ( $x_set->has(0) ) {
        my $idx = first_index { $_ == 0 } @{ $data->[0] };
        $opt{first_row}   = $idx + 1;
        $opt{last_row}    = $idx + 16;
        $opt{x_min_scale} = 0;
        $opt{x_max_scale} = 15;
        $opt{cross}       = 0;
    }
    elsif ( $x_set->has(1) ) {
        my $idx = first_index { $_ == 1 } @{ $data->[0] };
        $opt{first_row}   = $idx + 1;
        $opt{last_row}    = $idx + 16;
        $opt{x_min_scale} = 0;
        $opt{x_max_scale} = 15;
    }
    else {
        warn "X column errors\n";
        print YAML::Syck::Dump $data;
        return;
    }

    # chart 1
    $opt{y_column} = 1;
    $opt{y_data}   = $data->[1];
    $opt{y_title}  = "Nucleotide diversity";
    $write_obj->draw_y( $sheet, \%opt );

    # chart 2
    $opt{y_column} = 3;
    $opt{y_data}   = $data->[3];
    $opt{y_title}  = "Indel per 100 bp";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    # chart 3
    $opt{y_column}  = 5;
    $opt{y_data}    = $data->[5];
    $opt{y_title}   = "GC proportion";
    $opt{y2_column} = 7;
    $opt{y2_data}   = $data->[7];
    $opt{y2_title}  = "Window CV";
    $opt{top} += 18;
    $write_obj->draw_2y( $sheet, \%opt );
    delete $opt{y2_column};
    delete $opt{y2_data};
    delete $opt{y2_title};

    # chart 4
    $opt{y_column} = 9;
    $opt{y_data}   = $data->[9];
    $opt{y_title}  = $is_ld ? "r^2" : "Repeats proportion";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    if ( $x_set->has(-90) ) {
        my $idx = first_index { $_ == -90 } @{ $data->[0] };

        $opt{first_row}   = $idx + 1;
        $opt{last_row}    = $idx + 16;
        $opt{x_min_scale} = -90;
        $opt{x_max_scale} = -75;
        $opt{cross}       = -90;

        # chart 5
        $opt{y_column} = 1;
        $opt{y_data}   = $data->[1];
        $opt{y_title}  = "Nucleotide diversity";
        $opt{top}      = 1;
        $opt{left}     = 19;
        $write_obj->draw_y( $sheet, \%opt );

        # chart 6
        $opt{y_column} = 3;
        $opt{y_data}   = $data->[3];
        $opt{y_title}  = "Indel per 100 bp";
        $opt{top} += 18;
        $write_obj->draw_y( $sheet, \%opt );

        # chart 7
        $opt{y_column} = 5;
        $opt{y_data}   = $data->[5];
        $opt{y_title}  = "GC proportion";
        $opt{top} += 18;
        $write_obj->draw_y( $sheet, \%opt );

        # chart 8
        $opt{y_column} = 7;
        $opt{y_data}   = $data->[7];
        $opt{y_title}  = "Window CV";
        $opt{top} += 18;
        $write_obj->draw_y( $sheet, \%opt );

        # chart 9
        $opt{y_column} = 9;
        $opt{y_data}   = $data->[9];
        $opt{y_title}  = $is_ld ? "r^2" : "Repeats proportion";
        $opt{top} += 18;
        $write_obj->draw_y( $sheet, \%opt );
    }
};

#----------------------------------------------------------#
# worksheet -- summary
#----------------------------------------------------------#
my $summary_ofg = sub {
    my $sheet_name = 'summary';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{Type COUNT AVG_length SUM_length};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # contents
        my $query_name = 'ofg count';
        my $sql_query  = q{
            SELECT 
                CONCAT(o.ofg_tag, '_', o.ofg_type) Type, COUNT(*) COUNT
            FROM
                ofg o
            GROUP BY o.ofg_tag , o.ofg_type
            ORDER BY Type
        };

        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }
    $write_obj->increase_row;

    {    # contents
        my $query_name = 'ofg tag count';
        my $sql_query  = q{
            SELECT 
                o.ofg_tag, COUNT(*) COUNT
            FROM
                ofg o
            GROUP BY o.ofg_tag
            ORDER BY o.ofg_tag
        };

        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }
    $write_obj->increase_row;

    {    # contents
        my $query_name = 'ofg type count';
        my $sql_query  = q{
            SELECT  o.ofg_type,
                    COUNT(*) COUNT
            FROM ofg o
            GROUP BY o.ofg_type
            ORDER BY o.ofg_type
        };

        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }
    $write_obj->increase_row;

    {    # contents
        my $query_name = 'ofg';
        my $sql_query  = q{
            SELECT 
                CONCAT(o.ofg_tag, '_', o.ofg_type) Type,
                COUNT(*) COUNT,
                AVG(w.window_length) AVG_length,
                SUM(w.window_length) SUM_length
            FROM
                ofg o,
                window w
            WHERE
                w.window_id = o.window_id
            GROUP BY o.ofg_tag , o.ofg_type
        };

        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }
    $write_obj->increase_row;

    {    # contents
        my $query_name = 'ofgsw_outside';
        my $sql_query  = q{
            SELECT 
                CONCAT(o.ofg_tag, '_', o.ofg_type) Type,
                COUNT(*) COUNT,
                AVG(w.window_length) AVG_length,
                SUM(w.window_length) SUM_length
            FROM
                ofg o,
                window w,
                ofgsw s
            WHERE
                w.window_id = s.window_id
                    AND o.ofg_id = s.ofg_id
            GROUP BY o.ofg_tag , o.ofg_type
        };

        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }
    $write_obj->increase_row;

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- ofg_all
#----------------------------------------------------------#
my $ofg_all = sub {

    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my $sheet_name = "ofg_all";
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my $sql_query = q{
        SELECT 
            s.ofgsw_distance `distance`,
            AVG(w.window_pi) `AVG_pi`,
            STD(w.window_pi) `STD_pi`,
            AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
            STD(w.window_indel / w.window_length * 100) `STD_indel`,
            AVG(w.window_gc) `AVG_gc`,
            STD(w.window_gc) `STD_gc`,
            AVG(s.ofgsw_cv) `AVG_cv`,
            STD(s.ofgsw_cv) `STD_cv`,
            AVG(w.window_repeats) `AVG_repeats`,
            STD(w.window_repeats) `STD_repeats`,
            COUNT(*) COUNT
        FROM
            ofg o,
            ofgsw s,
            window w
        WHERE
            o.ofg_id = s.ofg_id
                AND s.window_id = w.window_id
        GROUP BY s.ofgsw_distance
    };

    my @names = $write_obj->sql2names($sql_query);
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    {    # content
        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query => $sql_query,
                data      => 1,
            }
        );
    }

    if ($add_chart) {    # chart
        $chart_ofg->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- ofg_coding
#----------------------------------------------------------#
my $ofg_coding = sub {

    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }
    unless ( $write_obj->check_column( 'window', 'window_coding' ) ) {
        return;
    }

    my @levels = ( [ "coding", 1 ], [ "noncoding", 0 ], );

    my $write_sheet = sub {
        my ( $order, $value ) = @_;
        my $sheet_name = "ofg_$order";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $sql_query = q{
            SELECT
                s.ofgsw_distance `distance`,
                AVG(w.window_pi) `AVG_pi`,
                STD(w.window_pi) `STD_pi`,
                AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
                STD(w.window_indel / w.window_length * 100) `STD_indel`,
                AVG(w.window_gc) `AVG_gc`,
                STD(w.window_gc) `STD_gc`,
                AVG(s.ofgsw_cv) `AVG_cv`,
                STD(s.ofgsw_cv) `STD_cv`,
                AVG(w.window_repeats) `AVG_repeats`,
                STD(w.window_repeats) `STD_repeats`,
                COUNT(*) COUNT
            FROM ofgsw s,
                 window w,
                 (
                  SELECT s.ofgsw_id
                  FROM ofgsw s,
                        ofg o,
                       window w
                  WHERE o.window_id = w.window_id AND
                        w.window_coding = ? AND
                        s.ofg_id = o.ofg_id
                 ) sw
            WHERE s.window_id = w.window_id AND
                  s.ofgsw_id = sw.ofgsw_id
            GROUP BY s.ofgsw_distance
        };

        my @names = $write_obj->sql2names($sql_query);
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );

        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $sql_query,
                    bind_value => [ $value, ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_ofg->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@levels) {
        $write_sheet->(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- ofg_coding_pure
#----------------------------------------------------------#
my $ofg_coding_pure = sub {

    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }
    unless ( $write_obj->check_column( 'window', 'window_coding' ) ) {
        return;
    }

    my @levels = ( [ "coding", 1 ], [ "noncoding", 0 ], );

    my $write_sheet = sub {
        my ( $order, $value ) = @_;
        my $sheet_name = "ofg_$order" . "_pure";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $sql_query = q{
            SELECT
                s.ofgsw_distance `distance`,
                AVG(w.window_pi) `AVG_pi`,
                STD(w.window_pi) `STD_pi`,
                AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
                STD(w.window_indel / w.window_length * 100) `STD_indel`,
                AVG(w.window_gc) `AVG_gc`,
                STD(w.window_gc) `STD_gc`,
                AVG(s.ofgsw_cv) `AVG_cv`,
                STD(s.ofgsw_cv) `STD_cv`,
                AVG(w.window_repeats) `AVG_repeats`,
                STD(w.window_repeats) `STD_repeats`,
                COUNT(*) COUNT
            FROM ofgsw s,
                 window w,
                 (
                  SELECT s.ofgsw_id
                  FROM ofgsw s,
                        ofg o,
                       window w
                  WHERE o.window_id = w.window_id AND
                        w.window_coding = ? AND
                        s.ofg_id = o.ofg_id
                 ) sw
            WHERE s.window_id = w.window_id AND
                  s.ofgsw_id = sw.ofgsw_id AND
                  w.window_coding = ?
            GROUP BY s.ofgsw_distance
        };

        my @names = $write_obj->sql2names($sql_query);
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $sql_query,
                    bind_value => [ $value, $value ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_ofg->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@levels) {
        $write_sheet->(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- ofg_dG
#----------------------------------------------------------#
my $ofg_dG = sub {

    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_dG' ) ) {
        return;
    }

    my $sheet_name = "ofg_dG";
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my $sql_query = q{
        SELECT
            s.ofgsw_distance `distance`,
            AVG(w.window_pi) `AVG_pi`,
            STD(w.window_pi) `STD_pi`,
            AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
            STD(w.window_indel / w.window_length * 100) `STD_indel`,
            AVG(w.window_gc) `AVG_gc`,
            STD(w.window_gc) `STD_gc`,
            AVG(s.ofgsw_cv) `AVG_cv`,
            STD(s.ofgsw_cv) `STD_cv`,
            AVG(s.ofgsw_dG) `avg_dG`,
            STD(s.ofgsw_dG) `std_dG`,
            COUNT(*) COUNT
          FROM ofg o,
               ofgsw s,
               window w
        WHERE     1 = 1
        AND o.ofg_id = s.ofg_id
        AND s.window_id = w.window_id
        AND s.ofgsw_distance != 0
        GROUP BY s.ofgsw_distance
    };

    my @names = $write_obj->sql2names($sql_query);
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    {    # content
        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query => $sql_query,
                data      => 1,
            }
        );
    }

    if ($add_chart) {    # chart
        $chart_ofg->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

my $ofg_tag_type = sub {

    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my $ary_ref;
    if ( $by eq "tag" ) {
        $ary_ref = get_tags($dbh);
    }
    elsif ( $by eq "type" ) {
        $ary_ref = get_types($dbh);
    }
    elsif ( $by eq "tt" ) {
        $ary_ref = get_tts($dbh);
    }

    my $write_sheet = sub {
        my ( $by_this, $bind ) = @_;

        my $sheet_name;
        if ( $by_this eq "tag" ) {
            $sheet_name = "ofg_tag_$bind";
        }
        elsif ( $by_this eq "type" ) {
            $sheet_name = "ofg_type_$bind";
        }
        elsif ( $by_this eq "tt" ) {
            $sheet_name = "ofg_tt_$bind";
        }
        $sheet_name = substr $sheet_name, 0, 31;    # excel sheet name limit
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $sql_query = q{
            SELECT 
                s.ofgsw_distance `distance`,
                AVG(w.window_pi) `AVG_pi`,
                STD(w.window_pi) `STD_pi`,
                AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
                STD(w.window_indel / w.window_length * 100) `STD_indel`,
                AVG(w.window_gc) `AVG_gc`,
                STD(w.window_gc) `STD_gc`,
                AVG(s.ofgsw_cv) `AVG_cv`,
                STD(s.ofgsw_cv) `STD_cv`,
                AVG(w.window_repeats) `AVG_repeats`,
                STD(w.window_repeats) `STD_repeats`,
                COUNT(*) COUNT
            FROM
                ofg o,
                ofgsw s,
                window w
            WHERE
                o.ofg_id = s.ofg_id
                    AND s.window_id = w.window_id
        };
        if ( $by_this eq "tag" ) {
            $sql_query .= "AND o.ofg_tag = ? GROUP BY s.ofgsw_distance";
        }
        elsif ( $by_this eq "type" ) {
            $sql_query .= "AND o.ofg_type = ? GROUP BY s.ofgsw_distance";
        }
        elsif ( $by_this eq "tt" ) {
            $sql_query .= q{AND CONCAT(o.ofg_tag, "_", o.ofg_type) = ? GROUP BY s.ofgsw_distance};
        }

        my @names = $write_obj->sql2names($sql_query);
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $sql_query,
                    bind_value => [ $bind, ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_ofg->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for ( @{$ary_ref} ) {
        $write_sheet->( $by, $_ );
    }
};

my $ofg_ld_tag_type = sub {

    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_r2_s' ) ) {
        return;
    }

    my $ary_ref;
    if ( $by eq "tag" ) {
        $ary_ref = get_tags($dbh);
    }
    elsif ( $by eq "type" ) {
        $ary_ref = get_types($dbh);
    }
    elsif ( $by eq "tt" ) {
        $ary_ref = get_tts($dbh);
    }

    my $write_sheet = sub {
        my ( $by_this, $bind ) = @_;

        my $sheet_name;
        if ( $by_this eq "tag" ) {
            $sheet_name = "ofg_ld_tag_$bind";
        }
        elsif ( $by_this eq "type" ) {
            $sheet_name = "ofg_ld_type_$bind";
        }
        elsif ( $by_this eq "tt" ) {
            $sheet_name = "ofg_ld_tt_$bind";
        }
        $sheet_name = substr $sheet_name, 0, 31;    # excel sheet name limit
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $sql_query = q{
            SELECT
                s.ofgsw_distance `distance`,
                AVG(w.window_pi) `AVG_pi`,
                STD(w.window_pi) `STD_pi`,
                AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
                STD(w.window_indel / w.window_length * 100) `STD_indel`,
                AVG(w.window_gc) `AVG_gc`,
                STD(w.window_gc) `STD_gc`,
                AVG(s.ofgsw_cv) `AVG_cv`,
                STD(s.ofgsw_cv) `STD_cv`,
                AVG(s.ofgsw_r2_s) `AVG_r^2`,
                STD(s.ofgsw_r2_s) `STD_r^2`,
                COUNT(*) COUNT
            FROM
                ofg o,
                ofgsw s,
                window w
            WHERE
                o.ofg_id = s.ofg_id
                    AND s.window_id = w.window_id
        };
        if ( $by_this eq "tag" ) {
            $sql_query .= "AND o.ofg_tag = ? GROUP BY s.ofgsw_distance";
        }
        elsif ( $by_this eq "type" ) {
            $sql_query .= "AND o.ofg_type = ? GROUP BY s.ofgsw_distance";
        }
        elsif ( $by_this eq "tt" ) {
            $sql_query .= q{AND CONCAT(o.ofg_tag, "_", o.ofg_type) = ? GROUP BY s.ofgsw_distance};
        }

        my @names = $write_obj->sql2names($sql_query);
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $sql_query,
                    bind_value => [ $bind, ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_ofg->( $sheet, $data, 1 );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for ( @{$ary_ref} ) {
        $write_sheet->( $by, $_ );
    }
};

for my $n (@tasks) {
    if ( $n == 1 )  { &$summary_ofg;     next; }
    if ( $n == 2 )  { &$ofg_all;         next; }
    if ( $n == 3 )  { &$ofg_coding;      next; }
    if ( $n == 4 )  { &$ofg_coding_pure; next; }
    if ( $n == 5 )  { &$ofg_dG;          next; }
    if ( $n == 6 )  { &$ofg_tag_type;    next; }
    if ( $n == 11 ) { &$ofg_ld_tag_type; next; }
}

if ($add_index_sheet) {
    $write_obj->add_index_sheet;
    print "Sheet [INDEX] has been generated.\n";
}

$stopwatch->end_message;
exit;

sub get_tags {
    my $dbh = shift;

    my $ary_ref = $dbh->selectcol_arrayref(
        q{
        SELECT DISTINCT o.ofg_tag
        FROM ofg o
        ORDER BY o.ofg_tag
        }
    );

    return $ary_ref;
}

sub get_types {
    my $dbh = shift;

    my $ary_ref = $dbh->selectcol_arrayref(
        q{
        SELECT DISTINCT o.ofg_type
        FROM ofg o
        ORDER BY o.ofg_type
        }
    );

    return $ary_ref;
}

sub get_tts {
    my $dbh = shift;

    my $ary_ref = $dbh->selectcol_arrayref(
        q{
        SELECT DISTINCT CONCAT(o.ofg_tag, "_", o.ofg_type)
        FROM ofg o
        ORDER BY CONCAT(o.ofg_tag, "_", o.ofg_type)
        }
    );

    return $ary_ref;
}

__END__
