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
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::ToXLSX;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

gene_stat_factory.pl - Gene stats for alignDB

=head1 SYNOPSIS

    perl gene_stat_factory.pl [options]
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
        --replace       STR=STR replace strings in axis names
        --index                 add an index sheet
        --chart                 add charts

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'output|o=s'   => \( my $outfile ),
    'run|r=s'      => \( my $run      = $Config->{stat}{run} ),
    'combine=i'    => \( my $combine  = 0 ),
    'replace=s'    => \my %replace,
    'index'        => \my $add_index_sheet,
    'chart'        => \my $add_chart,
) or HelpMessage(1);

$outfile = "$db.gene.xlsx" unless $outfile;

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
$stopwatch->start_message("Do gene stat for $db...");

my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password )
    or die "Cannot connect to MySQL database at $db:$server";
my $write_obj = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
    replace => \%replace,
);

my $sql_file = AlignDB::SQL::Library->new( lib => "$FindBin::Bin/sql.lib" );

# auto detect combine threshold
if ( $combine == 0 ) {
    ($combine) = $write_obj->calc_threshold;
}

#----------------------------------------------------------#
# chart -- coding
#----------------------------------------------------------#
my $chart_coding = sub {
    my $sheet = shift;
    my $data  = shift;

    my %opt = (
        x_column    => 0,
        y_column    => 1,
        first_row   => 1,
        last_row    => 16,
        x_min_scale => -5,
        x_max_scale => 10,
        y_data      => $data->[1],
        x_title     => "Distance to coding borders",
        y_title     => "Nucleotide diversity",
        top         => 1,
        left        => 10,
    );
    $write_obj->draw_y( $sheet, \%opt );

    #cross        => -5,

    $opt{y_column} = 2;
    $opt{y_data}   = $data->[2];
    $opt{y_title}  = "Indel per 100 bp";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 3;
    $opt{y_data}   = $data->[3];
    $opt{y_title}  = "GC proportion";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 4;
    $opt{y_data}   = $data->[4];
    $opt{y_title}  = "Window CV";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 5;
    $opt{y_data}   = $data->[5];
    $opt{y_title}  = "Repeats proportion";
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );
};

#----------------------------------------------------------#
# chart -- distance_*
#----------------------------------------------------------#
my $chart_distance = sub {
    my $sheet = shift;
    my $data  = shift;

    my %opt = (
        x_column      => 0,
        y_column      => 2,
        y_last_column => 3,
        first_row     => 2,
        last_row      => 17,
        x_min_scale   => 0,
        x_max_scale   => 15,
        y_data        => [ map { $data->[$_] } 2 .. 3 ],
        x_title       => "Distance to indels (d1)",
        y_title       => "Syn - Nonsyn",
        top           => 1,
        left          => 10,
    );
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column}      = 6;
    $opt{y_last_column} = 6;
    $opt{y_title}       = "dn/ds";
    $opt{y_data}        = $data->[6];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );
};

#----------------------------------------------------------#
# worksheet -- summary_gene
#----------------------------------------------------------#
my $summary_gene = sub {
    my $sheet_name = 'summary';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{ COUNT AVG_length SUM_length AVG_pi indel INDEL/100bp };
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # contents
        my $query_name = 'distinct gene';
        my $sql_query  = q{
            SELECT COUNT(DISTINCT g.gene_stable_id) COUNT
            FROM gene g
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'gene';
        my $sql_query  = q{
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`
            FROM gene g, window w
            WHERE w.window_id = g.window_id
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'full gene';
        my $sql_query  = q{
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`
            FROM gene g, window w
            WHERE w.window_id = g.window_id
            AND g.gene_is_full = 1
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
        my $query_name = 'distinct exon';
        my $sql_query  = q{
            SELECT COUNT(DISTINCT e.exon_stable_id) COUNT
            FROM exon e
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'exon';
        my $sql_query  = q{
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`
            FROM exon e, window w
            WHERE w.window_id = e.window_id
        };
        $write_obj->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # contents
        my $query_name = 'full exon';
        my $sql_query  = q{
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100 `INDEL/100bp`
            FROM exon e, window w
            WHERE w.window_id = e.window_id
            AND e.exon_is_full = 1
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
# worksheet -- coding_all
#----------------------------------------------------------#
my $coding_all = sub {

    # if the target column of the target table does not contain any values,
    # skip this stat
    return unless $write_obj->check_column( 'codingsw', 'codingsw_id' );

    my $sheet_name = "coding_all";
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my @names = qw{ distance AVG_pi AVG_Indel/100bp AVG_gc AVG_CV AVG_repeats COUNT };
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    {    # query
        my $sql_query = q{
            SELECT sw.codingsw_distance `distance`,
                   AVG (w.window_pi) `avg_pi`,
                   AVG (w.window_indel / w.window_length * 100)
                   `avg_indel/100bp`,
                   AVG (w.window_target_gc) `avg_gc`,
                   AVG (sw.codingsw_cv) `avg_cv`,
                   AVG (w.window_repeats) `avg_repeats`,
                   count(*) count
            FROM codingsw sw,
                 window w
            WHERE sw.window_id = w.window_id 
              AND sw.codingsw_distance < 0
              AND w.window_coding = 1
            GROUP BY sw.codingsw_distance
            ORDER BY sw.codingsw_distance ASC
        };
        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query => $sql_query,
                data      => $data,
            }
        );
    }

    {    # query
        my $sql_query = q{
            SELECT sw.codingsw_distance `distance`,
                   AVG (w.window_pi) `avg_pi`,
                   AVG (w.window_indel / w.window_length * 100)
                   `avg_indel/100bp`,
                   AVG (w.window_target_gc) `avg_gc`,
                   AVG (sw.codingsw_cv) `avg_cv`,
                   AVG (w.window_repeats) `avg_repeats`,
                   count(*) count
            FROM codingsw sw,
                 window w
            WHERE sw.window_id = w.window_id
              AND sw.codingsw_distance > 0
              AND w.window_coding = 0
            GROUP BY sw.codingsw_distance
            ORDER BY sw.codingsw_distance ASC
        };
        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query => $sql_query,
                data      => $data,
            }
        );
    }

    if ($add_chart) {    # chart
        $chart_coding->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- combined_dnds
#----------------------------------------------------------#
my $combined_dnds = sub {
    my $sheet_name = 'combined_dnds';
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

    my $thaw_sql = $sql_file->retrieve('dnds-d1_comb_dn_ds-0');

    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('dnds-d1_comb_dn_ds-0');
        $thaw_sql->add_where( 'isw.isw_distance' => $comb );

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $thaw_sql->as_sql,
                query_name => $_->[0],
                bind_value => $comb,
                data       => $data,
            }
        );
    }

    if ($add_chart) {    # chart
        $chart_distance->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance(frequecy)
#----------------------------------------------------------#
my $frequency_dnds = sub {
    unless ( $write_obj->check_column( 'indel', 'indel_freq' ) ) {
        return;
    }

    my @freq_levels = ( [ 1, 1, 1 ] );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'dnds_freq_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('dnds-d1_dn_ds-0');

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # contents
            $thaw_sql->from(  [] );
            $thaw_sql->joins( [] );
            $thaw_sql->add_join(
                isw => {
                    type      => 'inner',
                    table     => 'indel',
                    condition => 'isw.isw_indel_id = indel.indel_id'
                }
            );
            $thaw_sql->add_where( 'indel.indel_freq' => { op => '>=', value => '1' } );
            $thaw_sql->add_where( 'indel.indel_freq' => { op => '<=', value => '1' } );

            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $level->[1], $level->[2] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_distance->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@freq_levels) {
        &$write_sheet($_);
    }
};

for my $n (@tasks) {
    if ( $n == 1 )  { &$summary_gene;   next; }
    if ( $n == 2 )  { &$coding_all;     next; }
    if ( $n == 10 ) { &$combined_dnds;  next; }
    if ( $n == 11 ) { &$frequency_dnds; next; }
}

$stopwatch->end_message;
exit;

__END__
