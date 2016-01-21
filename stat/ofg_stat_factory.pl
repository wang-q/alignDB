#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::ToXLSX;

use lib "$FindBin::RealBin/../lib";
use AlignDB;
use AlignDB::Ofg;

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
        --combine       INT     
        --piece         INT     
        --index                 add an index sheet

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'output|o=s'   => \( my $outfile ),
    'by=s'         => \( my $by       = "tag" ),
    'run|r=s'      => \( my $run      = $Config->{stat}{run} ),
) or HelpMessage(1);

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
);

#----------------------------------------------------------#
# worksheet -- summary
#----------------------------------------------------------#
my $summary_ofg = sub {
    my $sheet_name = 'summary';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my @names = qw{Type COUNT AVG_length SUM_length};

    {    # write header
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@names,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header( $sheet_name, \%option );
    }

    {    # write contents
        my $query_name = 'ofg count';
        my $sql_query  = q{
            SELECT 
                CONCAT(o.ofg_tag, '_', o.ofg_type) Type, COUNT(*) COUNT
            FROM
                ofg o
            GROUP BY o.ofg_tag , o.ofg_type
            ORDER BY Type
        };

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
    }
    $sheet_row++;

    {    # write contents
        my $data = [];
        push @{$data}, [] for @names;

        my $query_name = 'ofg tag count';
        my $sql_query  = q{
            SELECT 
                o.ofg_tag, COUNT(*) COUNT
            FROM
                ofg o
            GROUP BY o.ofg_tag
            ORDER BY o.ofg_tag
        };

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
    }
    $sheet_row++;

    {    # write contents
        my $data = [];
        push @{$data}, [] for @names;

        my $query_name = 'ofg type count';
        my $sql_query  = q{
            SELECT  o.ofg_type,
                    COUNT(*) COUNT
            FROM ofg o
            GROUP BY o.ofg_type
            ORDER BY o.ofg_type
        };

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
    }
    $sheet_row++;

    {    # write contents
        my $data = [];
        push @{$data}, [] for @names;

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

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
    }
    $sheet_row++;

    {    # write contents
        my $data = [];
        push @{$data}, [] for @names;

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

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
    }
    $sheet_row++;

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- ofg_all
#----------------------------------------------------------#
my $ofg_all = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my $sheet_name = "ofg_all";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $sql_query = q{
        SELECT 
            s.ofgsw_distance `distance`,
            AVG(w.window_pi) `AVG_pi`,
            STD(w.window_pi) `STD_pi`,
            AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
            STD(w.window_indel / w.window_length * 100) `STD_indel`,
            AVG(w.window_target_gc) `AVG_gc`,
            STD(w.window_target_gc) `STD_gc`,
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
    my $data  = [];
    push @{$data}, [] for @names;

    {    # header
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@names,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header( $sheet_name, \%option );
    }

    {    # content
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- ofg_coding
#----------------------------------------------------------#
my $ofg_coding = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
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
        my ( $sheet_row, $sheet_col );

        my $sql_query = q{
            SELECT
                s.ofgsw_distance `distance`,
                AVG(w.window_pi) `AVG_pi`,
                STD(w.window_pi) `STD_pi`,
                AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
                STD(w.window_indel / w.window_length * 100) `STD_indel`,
                AVG(w.window_target_gc) `AVG_gc`,
                STD(w.window_target_gc) `STD_gc`,
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
        my $data  = [];
        push @{$data}, [] for @names;

        {    # header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header( $sheet_name, \%option );
        }

        {    # content
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $value, ],
            );
            ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
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

    # if the target column of the target table does not contain
    #   any values, skip this stat
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
        my ( $sheet_row, $sheet_col );

        my $sql_query = q{
            SELECT
                s.ofgsw_distance `distance`,
                AVG(w.window_pi) `AVG_pi`,
                STD(w.window_pi) `STD_pi`,
                AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
                STD(w.window_indel / w.window_length * 100) `STD_indel`,
                AVG(w.window_target_gc) `AVG_gc`,
                STD(w.window_target_gc) `STD_gc`,
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
        my $data  = [];
        push @{$data}, [] for @names;

        {    # header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header( $sheet_name, \%option );
        }

        {    # content
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $value, $value, ],
            );

            ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
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

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_dG' ) ) {
        return;
    }

    my $sheet_name = "ofg_dG";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    my $sql_query = q{
        SELECT
            s.ofgsw_distance `distance`,
            AVG(w.window_pi) `AVG_pi`,
            STD(w.window_pi) `STD_pi`,
            AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
            STD(w.window_indel / w.window_length * 100) `STD_indel`,
            AVG(w.window_target_gc) `AVG_gc`,
            STD(w.window_target_gc) `STD_gc`,
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
    my $data  = [];
    push @{$data}, [] for @names;

    {    # header
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@names,
        );
        ( $sheet, $sheet_row ) = $write_obj->write_header( $sheet_name, \%option );
    }

    {    # content
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

my $ofg_tag_type = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    unless ( $write_obj->check_column( 'ofgsw', 'ofgsw_id' ) ) {
        return;
    }

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    AlignDB::Ofg->meta->apply($obj);

    my $ary_ref;
    if ( $by eq "tag" ) {
        $ary_ref = $obj->get_tags;
    }
    elsif ( $by eq "type" ) {
        $ary_ref = $obj->get_types;
    }
    elsif ( $by eq "tt" ) {
        $ary_ref = $obj->get_tts;
    }

    my $write_sheet = sub {
        my ( $by, $bind ) = @_;

        my $sheet_name;
        if ( $by eq "tag" ) {
            $sheet_name = "ofg_tag_$bind";
        }
        elsif ( $by eq "type" ) {
            $sheet_name = "ofg_type_$bind";
        }
        elsif ( $by eq "tt" ) {
            $sheet_name = "ofg_tt_$bind";
        }
        $sheet_name = substr $sheet_name, 0, 31;    # excel sheet name limit
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $sql_query = q{
            SELECT 
                s.ofgsw_distance `distance`,
                AVG(w.window_pi) `AVG_pi`,
                STD(w.window_pi) `STD_pi`,
                AVG(w.window_indel / w.window_length * 100) `AVG_indel`,
                STD(w.window_indel / w.window_length * 100) `STD_indel`,
                AVG(w.window_target_gc) `AVG_gc`,
                STD(w.window_target_gc) `STD_gc`,
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
        if ( $by eq "tag" ) {
            $sql_query .= "AND o.ofg_tag = ? GROUP BY s.ofgsw_distance";
        }
        elsif ( $by eq "type" ) {
            $sql_query .= "AND o.ofg_type = ? GROUP BY s.ofgsw_distance";
        }
        elsif ( $by eq "tt" ) {
            $sql_query .= q{AND CONCAT(o.ofg_tag, "_", o.ofg_type) = ? GROUP BY s.ofgsw_distance};
        }

        my @names = $write_obj->sql2names($sql_query);
        my $data  = [];
        push @{$data}, [] for @names;

        {    # header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header( $sheet_name, \%option );
        }

        {    # content
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $bind, ],
            );
            ($sheet_row) = $write_obj->write_sql( $sheet, \%option );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for ( @{$ary_ref} ) {
        $write_sheet->( $by, $_ );
    }
};

for my $n (@tasks) {
    if ( $n == 1 ) { &$summary_ofg;     next; }
    if ( $n == 2 ) { &$ofg_all;         next; }
    if ( $n == 3 ) { &$ofg_coding;      next; }
    if ( $n == 4 ) { &$ofg_coding_pure; next; }
    if ( $n == 5 ) { &$ofg_dG;          next; }
    if ( $n == 6 ) { &$ofg_tag_type;    next; }
}

$stopwatch->end_message;
exit;

__END__
