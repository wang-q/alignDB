#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use Tie::IxHash;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::ToXLSX;

use lib "$FindBin::RealBin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

common_stat_factory.pl - Common stats for alignDB

=head1 SYNOPSIS

    perl common_stat_factory.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --outfile   -o  STR     outfile filename
        --run       -r  STR     run special analysis
        --combine       INT     
        --piece         INT
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
    'run|r=s'      => \( my $run      = $Config->{stat}{run} ),
    'combine=i'    => \( my $combine  = 0 ),
    'piece=i'      => \( my $piece    = 0 ),
    'replace=s'    => \my %replace,
    'index'        => \( my $add_index_sheet, ),
    'chart'        => \( my $add_chart, ),
) or Getopt::Long::HelpMessage(1);

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
    $outfile = "$db.common.xlsx" unless $outfile;
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
        $outfile = "$db.common.$runlist.xlsx";
    }
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Do stat for $db...");

my $aligndb_obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
my DBI $dbh = $aligndb_obj->dbh;

my $write_obj = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
    replace => \%replace,
);

# auto detect combine threshold
if ( $combine == 0 ) {
    ($combine) = $write_obj->calc_threshold;
}

# auto detect combine threshold
if ( $piece == 0 ) {
    ( undef, $piece ) = $write_obj->calc_threshold;
}

# count freq
my $seq_count = $aligndb_obj->get_seq_count;

# sql library
my $sql_file = AlignDB::SQL::Library->new( lib => "$FindBin::RealBin/sql.lib" );

#----------------------------------------------------------#
# chart -- pigccv_*
#----------------------------------------------------------#
my $chart_pigccv = sub {
    my $sheet = shift;
    my $data  = shift;

    my $sheet_name = $sheet->get_name;

    my %opt = (
        x_column    => 0,
        y_column    => 1,
        first_row   => 2,
        last_row    => 17,
        x_max_scale => 15,
        y_data      => $data->[1],
        x_title     => "Distance to indels (d1)",
        y_title     => "Nucleotide diversity",
        top         => 1,
        left        => 10,
    );
    if ( $sheet_name =~ /^d2_/ ) {
        $opt{x_title}     = "Reciprocal of indel density (d2)";
        $opt{x_max_scale} = 30;
        $opt{last_row}    = 33;
        $opt{x_max_scale} = 30;
    }
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 3;
    $opt{y_data}   = $data->[3];
    $opt{y_title}  = "GC proportion";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 5;
    $opt{y_data}   = $data->[5];
    $opt{y_title}  = "Window CV";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column}  = 3;
    $opt{y_data}    = $data->[3];
    $opt{y_title}   = "GC proportion";
    $opt{y2_column} = 5;
    $opt{y2_data}   = $data->[5];
    $opt{y2_title}  = "Window CV";
    $opt{top} += 18;
    $write_obj->draw_2y( $sheet, \%opt );
    delete $opt{y2_column};
    delete $opt{y2_data};
    delete $opt{y2_title};
};

my $chart_dd = sub {
    my $sheet      = shift;
    my $data_of    = shift;
    my $sheet_name = $sheet->get_name;

    # write charting data
    my @keys = keys %{$data_of};
    $write_obj->row(2);
    $write_obj->column(7);

    $write_obj->write_column( $sheet, { column => $data_of->{ $keys[-1] }[0], } );
    for my $key (@keys) {
        $write_obj->write_column(
            $sheet,
            {   query_name => $key,
                column     => $data_of->{$key}[1],
            }
        );
    }

    my %opt = (
        x_column      => 7,
        y_column      => 8,
        y_last_column => 8 + @keys - 1,
        first_row     => 2,
        last_row      => 12,
        x_min_scale   => 0,
        x_max_scale   => 10,
        y_data        => [ map { $data_of->{$_}[1] } @keys ],
        x_title       => "Distance to indels (d1)",
        y_title       => "Nucleotide diversity",
        top           => 20,
        left          => 7,
        height        => 480,
        width         => 480,
    );
    if ( $sheet_name =~ /_gc$/ ) {
        $opt{y_title} = "GC proportion";
    }
    $write_obj->draw_dd( $sheet, \%opt );
};

my $chart_series = sub {
    my $sheet      = shift;
    my $data_of    = shift;
    my $sheet_name = $sheet->get_name;

    # write charting data
    my @keys = keys %{$data_of};
    $write_obj->row(2);
    $write_obj->column(7);

    $write_obj->write_column( $sheet, { column => $data_of->{ $keys[-1] }[0], } );
    for my $key (@keys) {
        $write_obj->write_column(
            $sheet,
            {   query_name => $key,
                column     => $data_of->{$key}[1],
            }
        );
    }

    my %opt = (
        x_column      => 7,
        y_column      => 8,
        y_last_column => 8 + @keys - 1,
        first_row     => 2,
        last_row      => 7,
        x_min_scale   => 0,
        x_max_scale   => 5,
        y_data        => [ map { $data_of->{$_}[1] } @keys ],
        x_title       => "Distance to indels (d1)",
        y_title       => "Nucleotide diversity",
        top           => 11,
        left          => 7,
        height        => 480,
        width         => 480,
    );
    if ( $sheet_name =~ /_gc$/ ) {
        $opt{y_title} = "GC proportion";
    }
    $write_obj->draw_dd( $sheet, \%opt );
};

my $chart_snp_indel_ratio = sub {
    my $sheet = shift;
    my $data  = shift;

    my %opt = (
        x_column  => 1,
        y_column  => 2,
        first_row => 1,
        last_row  => scalar @{ $data->[0] },
        x_data    => $data->[0],
        y_data    => $data->[1],
        x_title   => "Nucleotide diversity",
        y_title   => "SNP/Indel ratio",
        top       => 1,
        left      => 10,
    );
    $write_obj->draw_xy( $sheet, \%opt );
};

my $chart_snp = sub {
    my $sheet      = shift;
    my $data_of    = shift;
    my $sheet_name = $sheet->get_name;

    # write charting data
    my @keys = keys %{$data_of};
    $write_obj->row(2);
    $write_obj->column(7);

    $write_obj->write_column( $sheet, { column => $data_of->{ $keys[-1] }[0], } );
    for my $key (@keys) {
        $write_obj->write_column(
            $sheet,
            {   query_name => $key,
                column     => $data_of->{$key}[1],
            }
        );
    }

    my %opt = (
        x_column      => 7,
        y_column      => 8,
        y_last_column => 8 + @keys - 1,
        first_row     => 2,
        last_row      => 18,
        x_min_scale   => 0,
        x_max_scale   => 15,
        y_data        => [ map { $data_of->{$_}[1] } @keys ],
        x_title       => "Distance to indels (d1)",
        y_title       => "Substitutions",
        top           => 20,
        left          => 7,
        height        => 480,
        width         => 480,
    );
    if ( $sheet_name =~ /density/ ) {
        $opt{x_title}     = "Indel density (d2)";
        $opt{last_row}    = 33;
        $opt{x_max_scale} = 30;
    }
    $write_obj->draw_dd( $sheet, \%opt );
};

#----------------------------------------------------------#
# worksheet -- basic
#----------------------------------------------------------#
my $basic = sub {
    my $sheet_name = 'basic';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{VALUE};
    {    # header
        $sheet = $write_obj->write_header(
            $sheet_name,
            {   query_name => 'Item',
                header     => \@names,
            }
        );
    }

    {    # contents
        $write_obj->write_row(
            $sheet,
            {   query_name => 'No. of genomes',
                row        => [$seq_count],
            }
        );
    }

    {    # contents
        my $query_name = 'Genome size on average (Mb)';
        my $sql_query  = q{
            SELECT
                AVG(l.length) / 1000000.00
            FROM
                (SELECT
                    SUM(c.chr_length) length
                FROM
                    chromosome c
                INNER JOIN (SELECT DISTINCT
                    chr_id
                FROM
                    sequence) s ON c.chr_id = s.chr_id
                GROUP BY common_name) l
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'Aligned length (Mb)';
        my $sql_query  = q{
            SELECT  SUM(a.align_length) / 1000000.00
            FROM    align a
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'Indels per 100 bp';
        my $sql_query  = q{
            SELECT  SUM(a.align_indels) / SUM(a.align_comparables) * 100.0
            FROM    align a
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'Substitutions per 100 bp';
        my $sql_query  = q{
            SELECT  SUM(a.align_differences) * 1.0 / SUM(a.align_comparables) * 100.0
            FROM    align a
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'D on average';
        my $sql_query  = q{
            SELECT -0.75 * log2( 1 - ( 4.0 / 3.0 ) * original.Pi )
            FROM (
                SELECT AVG(a.align_pi) Pi
                FROM align a
                ) original
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'GC-content';
        my $sql_query  = q{
            SELECT sum(a.align_length * a.align_gc )
                   / sum(a.align_length)
            FROM align a
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'Coding portions';
        my $sql_query  = q{
            SELECT  SUM(s.seq_length * a.align_coding) / SUM(s.seq_length)
            FROM    sequence s, align a
            WHERE   s.seq_role = "T"
            AND     s.align_id = a.align_id
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'Repeats portions';
        my $sql_query  = q{
            SELECT  SUM(s.seq_length * a.align_repeats) / SUM(s.seq_length)
            FROM    sequence s, align a
            WHERE   s.seq_role = "T"
            AND     s.align_id = a.align_id
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- process
#----------------------------------------------------------#
my $process = sub {
    my $sheet_name = 'process';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my @names = qw{Order Operation Duration Cmd_line};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names, } );
    }

    {    # contents
        my $sql_query = q{
            SELECT  meta_value
            FROM    meta m
            WHERE   meta_key IN ("a_operation","d_duration", "e_cmd_line")
        };

        my $array_ref = $dbh->selectcol_arrayref($sql_query);

        my $order = 1;
        while ( scalar @{$array_ref} ) {
            my @row = splice @{$array_ref}, 0, 3;
            $write_obj->write_row( $sheet, { row => [ $order, @row ], } );
            $order++;
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- summary
#----------------------------------------------------------#
my $summary = sub {
    my $sheet_name = 'summary';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{AVG MIN MAX STD COUNT SUM};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names, } );
    }

    my $column_stat = sub {
        my $query_name = shift;
        my $table      = shift;
        my $column     = shift;
        my $where      = shift;

        my $sql_query = q{
            # summray stat of _COLUMN_
            SELECT AVG(_COLUMN_) AVG,
                   MIN(_COLUMN_) MIN,
                   MAX(_COLUMN_) MAX,
                   STD(_COLUMN_) STD,
                   COUNT(_COLUMN_) COUNT,
                   SUM(_COLUMN_) SUM
            FROM _TABLE_
        };

        $sql_query =~ s/_TABLE_/$table/g;
        $sql_query =~ s/_COLUMN_/$column/g;
        $sql_query .= $where if $where;

        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    };

    {    # contents
        $column_stat->( 'align_length',       'align', 'align_length' );
        $column_stat->( 'indel_length',       'indel', 'indel_length' );
        $column_stat->( 'indel_left_extand',  'indel', 'left_extand' );
        $column_stat->( 'indel_right_extand', 'indel', 'right_extand' );
        $column_stat->( 'indel_windows',      'isw',   'isw_length', 'WHERE isw_distance <= 0' );
        $column_stat->( 'indel_free_windows', 'isw',   'isw_length', 'WHERE isw_distance > 0' );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- pi_gc_cv
#----------------------------------------------------------#
my $pi_gc_cv = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    {
        my $sheet_name = 'd1_pi_gc_cv';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_pi_gc_cv-0');

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query => $thaw_sql->as_sql,
                    data      => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }

    {
        my $sheet_name = 'd2_pi_gc_cv';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_pi_gc_cv-0');
        $thaw_sql->replace( { distance => 'density' } );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query => $thaw_sql->as_sql,
                    data      => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- comb_pi_gc_cv
#----------------------------------------------------------#
my $comb_pi_gc_cv = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    {
        my $sheet_name = 'd1_comb_pi_gc_cv';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        # make combine
        my $combined = $write_obj->make_combine(
            {   sql_query  => $sql_file->retrieve('common-d1_combine-0')->as_sql,
                threshold  => $combine,
                standalone => [ -1, 0 ],
            }
        );

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        for my $comb ( @{$combined} ) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');
            $thaw_sql->add_where( 'isw.isw_distance' => $comb );

            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => $comb,
                    data       => $data,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }

    {
        my $sheet_name = 'd2_comb_pi_gc_cv';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        # make combine
        my $combined = $write_obj->make_combine(
            {   sql_query  => $sql_file->retrieve('common-d2_combine-0')->as_sql,
                threshold  => $combine,
                standalone => [ -1, 0 ],
            }
        );

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');
            $thaw_sql->replace( { distance => 'density' } );
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        for my $comb ( @{$combined} ) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');
            $thaw_sql->add_where( 'isw.isw_distance' => $comb );
            $thaw_sql->replace( { distance => 'density' } );

            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => $comb,
                    data       => $data,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- group_distance
#----------------------------------------------------------#
my $group_distance = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'group_distance';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{AVG_distance AVG_pi COUNT STD_pi SUM_length length_proportion};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # make last portion
    my ( $all_length, $last_portion );
    {
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_combine-0');
        my $portion               = 0.05;
        my %option                = (
            sql_query => $thaw_sql->as_sql,
            portion   => $portion,
        );
        ( $all_length, $last_portion ) = $write_obj->make_last_portion( \%option );
    }

    my @group_distance = (
        [-1],
        [0],
        [1],
        [2],
        [ 3 .. 5 ],
        [ 6 .. 10 ],
        [ 11 .. 20 ],
        [ 21 .. 50 ],
        [ 51 .. 999 ],
        [ 3 .. 20 ],
        [ 21 .. 999 ],
        [ 101 .. 999 ],
        [0],
        [ 1 .. 5 ],
        [ 6 .. 10 ],
        [ 11 .. 15 ],
        [ 16 .. 20 ],
        [ 21 .. 25 ],
        $last_portion,
    );

    for my $group (@group_distance) {
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_pi_avg-0');
        $thaw_sql->add_select( "SUM(isw_length)",                     'SUM_length' );
        $thaw_sql->add_select( "SUM(isw_length) / $all_length * 100", 'length_proportion' );
        $thaw_sql->add_where( 'isw_distance' => $group );

        my $group_name;
        if ( scalar @{$group} > 1 ) {
            $group_name = $group->[0] . "--" . $group->[-1];
        }
        else {
            $group_name = $group->[0];
        }

        $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $group,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- group_density
#----------------------------------------------------------#
my $group_density = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'group_density';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{AVG_density AVG_pi COUNT STD_pi SUM_length length_proportion};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # make last portion
    my ( $all_length, $last_portion );
    {
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_combine-0');
        $thaw_sql->replace( { distance => 'density' } );
        my $portion = 0.05;
        my %option  = (
            sql_query => $thaw_sql->as_sql,
            portion   => $portion,
        );
        ( $all_length, $last_portion ) = $write_obj->make_last_portion( \%option );
    }

    my @group_density = (
        [-1],
        [0],
        [1],
        [2],
        [ 3 .. 5 ],
        [ 6 .. 10 ],
        [ 11 .. 20 ],
        [ 21 .. 50 ],
        [ 51 .. 999 ],
        [ 3 .. 20 ],
        [ 21 .. 999 ],
        [ 101 .. 999 ],
        $last_portion,
    );

    for my $group (@group_density) {
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_pi_avg-0');
        $thaw_sql->add_select( "SUM(isw_length)",                     'SUM_length' );
        $thaw_sql->add_select( "SUM(isw_length) / $all_length * 100", 'length_proportion' );
        $thaw_sql->add_where( 'isw_distance' => $group );
        $thaw_sql->replace( { distance => 'density' } );

        my $group_name;
        if ( scalar @{$group} > 1 ) {
            $group_name = $group->[0] . "--" . $group->[-1];
        }
        else {
            $group_name = $group->[0];
        }

        $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $group,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- comb_coding
#----------------------------------------------------------#
my $comb_coding = sub {

    unless ( $write_obj->check_column( 'isw', 'isw_coding' ) ) {
        return;
    }

    my @type_levels = ( [ 'coding', 1, 1 ], [ 'non_coding', 0, 0 ], );

    my $write_sheet_d1 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        my $sheet_name = "d1_comb_$name";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        # make combine
        my $combined;
        {
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_make_combine_coding-0');
            $combined = $write_obj->make_combine(
                {   sql_query  => $thaw_sql->as_sql,
                    threshold  => $combine,
                    standalone => [ -1, 0 ],
                    bind_value => [ $feature_1, $feature_2 ],
                    merge_last => 1,
                }
            );
        }

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_coding-0');
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        for my $comb ( @{$combined} ) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_coding-0');
            $thaw_sql->add_where( 'isw.isw_distance' => $comb );

            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $feature_1, $feature_2, @{$comb} ],
                    data       => $data,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        my $sheet_name = "d2_comb_$name";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        # make combine
        my $combined;
        {
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_make_combine_coding-0');
            $thaw_sql->replace( { distance => 'density' } );
            $combined = $write_obj->make_combine(
                {   sql_query  => $thaw_sql->as_sql,
                    threshold  => $combine,
                    standalone => [ -1, 0 ],
                    bind_value => [ $feature_1, $feature_2 ],
                    merge_last => 1,
                }
            );
        }

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_coding-0');
            $thaw_sql->replace( { distance => 'density' } );
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        for my $comb ( @{$combined} ) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_coding-0');
            $thaw_sql->add_where( 'isw.isw_distance' => $comb );
            $thaw_sql->replace( { distance => 'density' } );

            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $feature_1, $feature_2, @{$comb} ],
                    data       => $data,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@type_levels) {
        $write_sheet_d1->($_);
    }

    for (@type_levels) {
        $write_sheet_d2->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- comb_slippage
#----------------------------------------------------------#
my $comb_slippage = sub {

    unless ( $write_obj->check_column( 'indel', 'indel_slippage' ) ) {
        return;
    }

    my @type_levels = ( [ 'slippage', 1, 1 ], [ 'non_slippage', 0, 0 ], );

    my $write_sheet_d1 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        my $sheet_name = "d1_comb_$name";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        # make combine
        my $combined;
        {
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_make_combine_slippage-0');
            $combined = $write_obj->make_combine(
                {   sql_query  => $thaw_sql->as_sql,
                    threshold  => $combine,
                    standalone => [ -1, 0 ],
                    bind_value => [ $feature_1, $feature_2 ],
                    merge_last => 1,
                }
            );
        }

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_slippage-0');
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        for my $comb ( @{$combined} ) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_slippage-0');
            $thaw_sql->add_where( 'isw.isw_distance' => $comb );

            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $feature_1, $feature_2, @{$comb} ],
                    data       => $data,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ( $name, $feature_1, $feature_2 ) = @{ $_[0] };

        my $sheet_name = "d2_comb_$name";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        # make combine
        my $combined;
        {
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_make_combine_slippage-0');
            $thaw_sql->replace( { distance => 'density' } );
            $combined = $write_obj->make_combine(
                {   sql_query  => $thaw_sql->as_sql,
                    threshold  => $combine,
                    standalone => [ -1, 0 ],
                    bind_value => [ $feature_1, $feature_2 ],
                    merge_last => 1,
                }
            );
        }

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_slippage-0');
            $thaw_sql->replace( { distance => 'density' } );
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        for my $comb ( @{$combined} ) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-d1_comb_slippage-0');
            $thaw_sql->add_where( 'isw.isw_distance' => $comb );
            $thaw_sql->replace( { distance => 'density' } );

            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $feature_1, $feature_2, @{$comb} ],
                    data       => $data,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@type_levels) {
        $write_sheet_d1->($_);
    }

    for (@type_levels) {
        $write_sheet_d2->($_);
    }
};

my $dd_group = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    #----------------------------------------------------------#
    # worksheet -- dd_group
    #----------------------------------------------------------#
    {
        my $sheet_name = 'dd_group';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(1);

        my @dd_density_group
            = ( [ 1, 2 ], [ 3, 6 ], [ 7, 10 ], [ 11, 18 ], [ 19, 999 ], );

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-dd_group');
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        tie my %data_of, 'Tie::IxHash';
        for my $item (@dd_density_group) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-dd_group');

            my $group_name = $item->[0] . "--" . $item->[1];
            $write_obj->increase_row;

            my $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    query_name => $group_name,
                    bind_value => [ @{$item}, $item->[0] ],
                    data       => 1,
                }
            );
            $data_of{$group_name} = $data;
        }

        if ($add_chart) {    # chart
            $chart_dd->( $sheet, \%data_of );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }

    #----------------------------------------------------------#
    # worksheet -- dd_group_gc
    #----------------------------------------------------------#
    {
        my $sheet_name = 'dd_group_gc';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(1);

        my @dd_density_group
            = ( [ 1, 2 ], [ 3, 6 ], [ 7, 10 ], [ 11, 18 ], [ 19, 999 ], );

        {    # header
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-dd_group_gc');
            my @names = $thaw_sql->as_header;
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        tie my %data_of, 'Tie::IxHash';
        for my $item (@dd_density_group) {    # content
            my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-dd_group_gc');

            my $group_name = $item->[0] . "--" . $item->[1];
            $write_obj->increase_row;

            my $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    query_name => $group_name,
                    bind_value => [ @{$item}, $item->[0] ],
                    data       => 1,
                }
            );
            $data_of{$group_name} = $data;
        }

        if ($add_chart) {    # chart
            $chart_dd->( $sheet, \%data_of );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }
};

#----------------------------------------------------------#
# worksheet -- indel_length_group
#----------------------------------------------------------#
my $indel_length_group = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'indel_length_group';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @groups = ( [ 1, 5 ], [ 6, 10 ], [ 11, 50 ], [ 51, 300 ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    tie my %data_of, 'Tie::IxHash';
    for my $item (@groups) {    # content
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        $thaw_sql->add_where( 'indel_length' => { op => '>=', value => '1' } );
        $thaw_sql->add_where( 'indel_length' => { op => '<=', value => '1' } );

        my $group_name = $item->[0] . "--" . $item->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $item,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_extand_group
#----------------------------------------------------------#
my $indel_extand_group = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'indel_extand_group';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @groups
        = ( [ 0, 0 ], [ 1, 2 ], [ 3, 4 ], [ 5, 9 ], [ 10, 19 ], [ 20, 999 ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    tie my %data_of, 'Tie::IxHash';
    for my $item (@groups) {    # content
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        $thaw_sql->add_where( 'FLOOR(indel.right_extand / 100)' => { op => '>=', value => '0' } );
        $thaw_sql->add_where( 'FLOOR(indel.right_extand / 100)' => { op => '<=', value => '0' } );

        my $group_name = $item->[0] . "--" . $item->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $item,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_position_group
#----------------------------------------------------------#
my $indel_position_group = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }
    unless ( $write_obj->check_column( 'indel', 'indel_coding' ) ) {
        return;
    }

    my $sheet_name = 'indel_position_group';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @groups
        = ( [ 1, 1, 0, 0 ], [ 1, 1, 1, 1 ], [ 0, 0, 0, 0 ], [ 0, 0, 1, 1 ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    tie my %data_of, 'Tie::IxHash';
    for my $item (@groups) {    # content
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        $thaw_sql->add_where( 'indel_coding'  => { op => '>=', value => '0' } );
        $thaw_sql->add_where( 'indel_coding'  => { op => '<=', value => '0' } );
        $thaw_sql->add_where( 'indel_repeats' => { op => '>=', value => '0' } );
        $thaw_sql->add_where( 'indel_repeats' => { op => '<=', value => '0' } );

        my $group_name = $item->[0] . "--" . $item->[2];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $item,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_coding_group
#----------------------------------------------------------#
my $indel_coding_group = sub {

    unless ( $write_obj->check_column( 'indel', 'indel_coding' ) ) {
        return;
    }

    my $sheet_name = 'indel_coding_group';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @groups = ( [ 0, 0 ], [ 1, 1 ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    tie my %data_of, 'Tie::IxHash';
    for my $item (@groups) {    # content
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        $thaw_sql->add_where( 'indel_coding' => { op => '>=', value => '0' } );
        $thaw_sql->add_where( 'indel_coding' => { op => '<=', value => '0' } );

        my $group_name = $item->[0] . "--" . $item->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $item,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_repeat_group
#----------------------------------------------------------#
my $indel_repeat_group = sub {

    unless ( $write_obj->check_column( 'indel', 'indel_repeats' ) ) {
        return;
    }

    my $sheet_name = 'indel_repeat_group';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @groups
        = ( [ 0, 0 ], [ 1, 1 ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    tie my %data_of, 'Tie::IxHash';
    for my $item (@groups) {    # content
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        $thaw_sql->add_where( 'indel_repeats' => { op => '>=', value => '0' } );
        $thaw_sql->add_where( 'indel_repeats' => { op => '<=', value => '0' } );

        my $group_name = $item->[0] . "--" . $item->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $item,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_slip_group
#----------------------------------------------------------#
my $indel_slip_group = sub {

    unless ( $write_obj->check_column( 'indel', 'indel_slippage' ) ) {
        return;
    }

    my $sheet_name = 'indel_slip_group';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @groups = ( [ 0, 0 ], [ 1, 1 ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    tie my %data_of, 'Tie::IxHash';
    for my $item (@groups) {    # content
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        $thaw_sql->add_where( 'indel_slippage' => { op => '>=', value => '0' } );
        $thaw_sql->add_where( 'indel_slippage' => { op => '<=', value => '0' } );

        my $group_name = $item->[0] . "--" . $item->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $item,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_gc_group
#----------------------------------------------------------#
my $indel_gc_group = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'indel_gc_group';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @groups = ( [ 0, 0.2999 ], [ 0.3, 0.4999 ], [ 0.5, 1 ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    tie my %data_of, 'Tie::IxHash';
    for my $item (@groups) {    # content
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_isw');
        $thaw_sql->add_where( 'indel_gc' => { op => '>=', value => '0' } );
        $thaw_sql->add_where( 'indel_gc' => { op => '<=', value => '0' } );
        $thaw_sql->add_where( 'indel_length' => \' >= 5' );

        my $group_name = $item->[0] . "--" . $item->[1];
        $write_obj->increase_row;

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $item,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;
    }

    if ($add_chart) {    # chart
        $chart_series->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- snp_indel_ratio
#----------------------------------------------------------#
my $snp_indel_ratio = sub {
    my $sheet_name = 'snp_indel_ratio';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # create temporary table
    {
        my $sql_query = q{
            DROP TABLE IF EXISTS pi_group
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    {
        my $sql_query = q{
            # create temporary table
            CREATE TABLE pi_group (p_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (p_id))
                ENGINE=MyISAM
                SELECT a.align_pi `pi`,
                       a.align_differences `snp`,
                       i.indel_number `indel`,
                       a.align_gaps `gaps`,
                       a.align_length `align_length`
                FROM (SELECT align_id, COUNT(indel_id) indel_number
                      FROM indel i
                      GROUP BY align_id) i,
                     align a
                WHERE i.align_id = a.align_id
                ORDER BY pi DESC
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    # make group
    my @group_align;
    {
        my $sql_query = q{
            SELECT p_id, align_length
            FROM pi_group
        };
        my %option = (
            sql_query => $sql_query,
            piece     => $piece,
        );
        @group_align = @{ $write_obj->make_combine_piece( \%option ) };
    }

    my @names
        = qw{AVG_pi AVG_SNP/Indel COUNT AVG_align_length SUM_align_length AVG_SNP/kb AVG_Indel/kb};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data = [];
    push @{$data}, [] for @names;
    {    # align query
        my $sql_query = q{
            SELECT AVG(p.pi) `AVG_pi`,
                   AVG(p.snp / p.indel) `AVG_SNP/Indel`,
                   COUNT(*) COUNT,
                   AVG(p.align_length) `AVG_align_length`,
                   SUM(p.align_length) `SUM_align_length`,
                   AVG(p.snp / p.align_length * 1000) `AVG_SNP/kb`,
                   AVG(p.indel / p.align_length * 1000) `AVG_Indel/kb`
            FROM pi_group p
            WHERE p_id IN
        };

        my @group_names;
        for (@group_align) {
            my @range        = @$_;
            my $in_list      = '(' . join( ',', @range ) . ')';
            my $sql_query_in = $sql_query . $in_list;
            my $group_name;
            if ( scalar @range > 1 ) {
                $group_name = $range[0] . "--" . $range[-1];
            }
            else {
                $group_name = $range[0];
            }
            push @group_names, $group_name;

            my DBI $sth = $dbh->prepare($sql_query_in);
            $sth->execute;
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
        }

        $sheet->write( $write_obj->row, 0, [ [@group_names] ], $write_obj->format->{NAME} );
        $sheet->write( $write_obj->row, 1, $data, $write_obj->format->{NORMAL} );
    }

    {    # drop temporary table
        my $sql_query = q{
            DROP TABLE IF EXISTS pi_group
        };
        my %option = ( sql_query => $sql_query, );
        $write_obj->excute_sql( \%option );
    }

    if ($add_chart) {    # chart
        $chart_snp_indel_ratio->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_length
#----------------------------------------------------------#
my $indel_length = sub {
    my $sheet_name = 'indel_length';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_length-0');

    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # content
        my DBI $sth = $dbh->prepare( $thaw_sql->as_sql );
        $sth->execute;

        my $last_number;
        while ( my @row = $sth->fetchrow_array ) {

            # Highlight 'special' indels
            my $style = 'NORMAL';
            if ( defined $last_number ) {
                if ( $row[1] > $last_number ) {
                    $style = 'HIGHLIGHT';
                }
            }
            $last_number = $row[1];

            for ( my $i = 0; $i < scalar @row; $i++ ) {
                $sheet->write(
                    $write_obj->row, $i + $write_obj->column,
                    $row[$i],        $write_obj->format->{$style}
                );
            }
            $write_obj->increase_row;
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_length_100
#----------------------------------------------------------#
my $indel_length_100 = sub {
    my $sheet_name = 'indel_length_100';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-indel_length-0');
    $thaw_sql->add_where( 'left_extand'  => \'>= 100' );
    $thaw_sql->add_where( 'right_extand' => \'>= 100' );

    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # content
        my DBI $sth = $dbh->prepare( $thaw_sql->as_sql );
        $sth->execute;

        my $last_number;
        while ( my @row = $sth->fetchrow_array ) {

            # Highlight 'special' indels
            my $style = 'NORMAL';
            if ( defined $last_number ) {
                if ( $row[1] > $last_number ) {
                    $style = 'HIGHLIGHT';
                }
            }
            $last_number = $row[1];

            for ( my $i = 0; $i < scalar @row; $i++ ) {
                $sheet->write(
                    $write_obj->row, $i + $write_obj->column,
                    $row[$i],        $write_obj->format->{$style}
                );
            }
            $write_obj->increase_row;
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_snp
#----------------------------------------------------------#
my $distance_snp = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'distance_snp';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # Six base groups
    my @pairs = qw{A<->C A<->G A<->T C<->G C<->T G<->T};

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-distance_snp');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    my %sum_of;
    for (@pairs) {
        my $group_name = $_;
        $write_obj->increase_row;

        my $pair_1 = $_;
        $pair_1 =~ s/\W//g;
        my $pair_2 = reverse $pair_1;
        my $bases = [ $pair_1, $pair_2 ];

        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-distance_snp');

        $thaw_sql->add_where( 'isw_distance'                    => \'>= -1' );
        $thaw_sql->add_where( 'isw_distance'                    => \'<= 15' );
        $thaw_sql->add_where( 'CONCAT(target_base, query_base)' => $bases );

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $bases,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;

        for my $idx ( 0 .. @{ $data->[0] } - 1 ) {
            my $category = $data->[0][$idx];
            my $value    = $data->[1][$idx];
            $sum_of{$category} += $value;
        }
    }

    for my $group (@pairs) {
        my $data = $data_of{$group};
        for my $idx ( 0 .. @{ $data->[0] } - 1 ) {
            my $category = $data->[0][$idx];
            my $value    = $data->[1][$idx];
            $data->[1][$idx] = $value / $sum_of{$category};
        }
    }

    if ($add_chart) {    # chart
        $chart_snp->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- density_snp
#----------------------------------------------------------#
my $density_snp = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'density_snp';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # Six base groups
    my @pairs = qw{A<->C A<->G A<->T C<->G C<->T G<->T};

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-distance_snp');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    my %sum_of;
    for (@pairs) {
        my $group_name = $_;
        $write_obj->increase_row;

        my $pair_1 = $_;
        $pair_1 =~ s/\W//g;
        my $pair_2 = reverse $pair_1;
        my $bases = [ $pair_1, $pair_2 ];

        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-distance_snp');

        $thaw_sql->add_where( 'isw_density'                     => \'>= -1' );
        $thaw_sql->add_where( 'isw_density'                     => \'<= 30' );
        $thaw_sql->add_where( 'CONCAT(target_base, query_base)' => $bases );

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $bases,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;

        for my $idx ( 0 .. @{ $data->[0] } - 1 ) {
            my $category = $data->[0][$idx];
            my $value    = $data->[1][$idx];
            $sum_of{$category} += $value;
        }
    }

    for my $group (@pairs) {
        my $data = $data_of{$group};
        for my $idx ( 0 .. @{ $data->[0] } - 1 ) {
            my $category = $data->[0][$idx];
            my $value    = $data->[1][$idx];
            $data->[1][$idx] = $value / $sum_of{$category};
        }
    }

    if ($add_chart) {    # chart
        $chart_snp->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_tri_trv
#----------------------------------------------------------#
my $distance_tri_trv = sub {
    unless ( $write_obj->check_column( 'isw', 'isw_id' ) ) {
        return;
    }

    my $sheet_name = 'distance_tri_trv';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # base groups
    my @pairs = ( [ 'Tri', qw{AG CT GA TC} ], [ 'Trv', qw{AC AT CA CG GC GT TA TG} ], );

    {    # header
        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-distance_snp');
        my @names = $thaw_sql->as_header;
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    tie my %data_of, 'Tie::IxHash';
    my %sum_of;
    for (@pairs) {
        my $group_name = shift @{$_};
        $write_obj->increase_row;

        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-distance_snp');

        $thaw_sql->add_where( 'isw_distance'                    => \'>= -1' );
        $thaw_sql->add_where( 'isw_distance'                    => \'<= 15' );
        $thaw_sql->add_where( 'CONCAT(target_base, query_base)' => $_ );

        my $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $group_name,
                bind_value => $_,
                data       => 1,
            }
        );
        $data_of{$group_name} = $data;

        for my $idx ( 0 .. @{ $data->[0] } - 1 ) {
            my $category = $data->[0][$idx];
            my $value    = $data->[1][$idx];
            $sum_of{$category} += $value;
        }
    }

    for my $group ( keys %data_of ) {
        my $data = $data_of{$group};
        for my $idx ( 0 .. @{ $data->[0] } - 1 ) {
            my $category = $data->[0][$idx];
            my $value    = $data->[1][$idx];
            $data->[1][$idx] = $value / $sum_of{$category};
        }
    }

    if ($add_chart) {    # chart
        $chart_snp->( $sheet, \%data_of );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- align_coding
#----------------------------------------------------------#
my $align_coding = sub {

    unless ( $write_obj->check_column( 'align', 'align_coding' ) ) {
        return;
    }

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT align_coding
            FROM align
            WHERE align_coding IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @coding_levels = (
        [ 1, $quartiles->[0], $quartiles->[1] ],
        [ 2, $quartiles->[1], $quartiles->[2] ],
        [ 3, $quartiles->[2], $quartiles->[3] ],
        [ 4, $quartiles->[3], $quartiles->[4] ],
        [ 9, 0.4,             0.6 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;

        my $sheet_name = "align_coding_$order";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-align-0');

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        {    # contents
            $thaw_sql->add_where( 'align.align_coding' => { op => '>=', value => '1' } );
            $thaw_sql->add_where( 'align.align_coding' => { op => '<=', value => '1' } );
            $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $low_border, $high_border ],
                }
            );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for (@coding_levels) {
        &$write_sheet(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- align_repeat
#----------------------------------------------------------#
my $align_repeat = sub {

    unless ( $write_obj->check_column( 'align', 'align_repeats' ) ) {
        return;
    }

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT align_repeats
            FROM align
            WHERE align_repeats IS NOT NULL
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @repeat_levels = (
        [ 1, $quartiles->[0], $quartiles->[1] ],
        [ 2, $quartiles->[1], $quartiles->[2] ],
        [ 3, $quartiles->[2], $quartiles->[3] ],
        [ 4, $quartiles->[3], $quartiles->[4] ],
        [ 9, 0.4,             0.6 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;

        my $sheet_name = "align_repeat_$order";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my AlignDB::SQL $thaw_sql = $sql_file->retrieve('common-align-0');

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        {    # contents
            $thaw_sql->add_where( 'align.align_repeats' => { op => '>=', value => '1' } );
            $thaw_sql->add_where( 'align.align_repeats' => { op => '<=', value => '1' } );
            $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $low_border, $high_border ],
                }
            );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for (@repeat_levels) {
        &$write_sheet(@$_);
    }
};

for my $n (@tasks) {
    if ( $n == 1 ) { &$basic; &$process; &$summary; next; }
    if ( $n == 2 )  { &$pi_gc_cv;             next; }
    if ( $n == 3 )  { &$comb_pi_gc_cv;        next; }
    if ( $n == 4 )  { &$group_distance;       &$group_density; next; }
    if ( $n == 5 )  { &$comb_coding;          next; }
    if ( $n == 6 )  { &$comb_slippage;        next; }
    if ( $n == 8 )  { &$dd_group;             next; }
    if ( $n == 9 )  { &$indel_length_group;   next; }
    if ( $n == 10 ) { &$indel_extand_group;   next; }
    if ( $n == 11 ) { &$indel_position_group; next; }
    if ( $n == 12 ) { &$indel_coding_group;   &$indel_repeat_group; next; }
    if ( $n == 13 ) { &$indel_slip_group;     &$indel_gc_group; next; }
    if ( $n == 14 ) { &$snp_indel_ratio;      next; }
    if ( $n == 15 ) { &$indel_length;         &$indel_length_100; next; }
    if ( $n == 17 ) { &$distance_snp;         &$density_snp; next; }
    if ( $n == 18 ) { &$distance_tri_trv;     next; }

    if ( $n == 51 ) { &$align_coding; next; }
    if ( $n == 52 ) { &$align_repeat; next; }
}

if ($add_index_sheet) {
    $write_obj->add_index_sheet;
    print "Sheet [INDEX] has been generated.\n";
}

$stopwatch->end_message;
exit;

__END__
