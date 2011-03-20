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
use AlignDB::WriteExcel;

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
my $run     = $Config->{stat}{run};
my $outfile = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=s'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'output=s'   => \$outfile,
    'run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

#----------------------------------------------------------#
# worksheet -- summary_gene
#----------------------------------------------------------#
my $summary_gene = sub {
    my $sheet_name = 'summary';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{
            COUNT AVG_length SUM_length AVG_pi
            indel INDEL/100bp ns_indel ns_INDEL/100bp
        };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            query_name => $query_name,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $query_name = 'distinct gene';
        my $sql_query  = q{
            SELECT COUNT(DISTINCT g.gene_stable_id) COUNT
            FROM gene g
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'gene';
        my $sql_query  = q{
            # gene
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
            FROM gene g, window w
            WHERE w.window_id = g.window_id
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'full gene';
        my $sql_query  = q{
            # gene
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
            FROM gene g, window w
            WHERE w.window_id = g.window_id
            AND g.gene_is_full = 1
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    {    # write contents
        my $query_name = 'distinct exon';
        my $sql_query  = q{
            # distinct exon
            SELECT COUNT(DISTINCT e.exon_stable_id) COUNT
            FROM exon e
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'exon';
        my $sql_query  = q{
            # exon
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
            FROM exon e, window w
            WHERE w.window_id = e.window_id
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'full exon';
        my $sql_query  = q{
            # full exon
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
            FROM exon e, window w
            WHERE w.window_id = e.window_id
            AND e.exon_is_full = 1
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    # add a blank row
    $sheet_row++;

    {    # write contents
        my $query_name = 'gene_ess';
        my $sql_query  = q{
            # gene_ess
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
            FROM gene g, window w
            WHERE w.window_id = g.window_id
            AND g.gene_feature4 = 1
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'gene_non_ess';
        my $sql_query  = q{
            # gene_non_ess
            SELECT COUNT(*) COUNT,
                   AVG(w.window_length) AVG_length,
                   SUM(w.window_length) SUM_length,
                   AVG(w.window_pi) AVG_pi,
                   SUM(w.window_indel) indel,
                   SUM(w.window_indel) / SUM(w.window_length) * 100
                   `INDEL/100bp`,
                   SUM(w.window_ns_indel) `ns_indel`,
                   SUM(w.window_ns_indel) / SUM(w.window_length) * 100
                   `ns_INDEL/100bp`
            FROM gene g, window w
            WHERE w.window_id = g.window_id
            AND g.gene_feature4 = 0
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
# worksheet -- coding_all
#----------------------------------------------------------#
my $coding_all = sub {

    # if the target column of the target table does not contain any values,
    # skip this stat
    return unless $write_obj->check_column( 'codingsw', 'codingsw_id' );

    my $sheet_name = "coding_all";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ distance AVG_pi AVG_Indel/100bp AVG_gc AVG_CV COUNT };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # query
        my $sql_query = q{
            SELECT s.codingsw_distance `distance`,
                   AVG (w.window_pi) `avg_pi`,
                   AVG (w.window_indel / w.window_length * 100)
                   `avg_indel/100bp`,
                   AVG (w.window_target_gc) `avg_gc`,
                   AVG (s.codingsw_cv) `avg_cv`,
                   count(*) count
            FROM codingsw s,
                 window w,
                 (
                  SELECT s.codingsw_id
                  FROM codingsw s,
                       window w
                  WHERE s.window_id = w.window_id AND
                        ASCII (s.codingsw_type) IN (ASCII ('L'), ASCII ('R'))
                  UNION
                  SELECT s.codingsw_id
                  FROM codingsw s,
                       exon e,
                       window w
                  WHERE s.exon_id = e.exon_id AND
                        e.window_id = w.window_id AND
                        ASCII (s.codingsw_type) IN (ASCII ('l'), ASCII ('r'))
                 ) sw
            WHERE s.window_id = w.window_id AND
                  s.codingsw_id = sw.codingsw_id
            GROUP BY s.codingsw_distance
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
# worksheet -- exon_D
#----------------------------------------------------------#
my $exon_D = sub {

    # if the target column of the target table does not contain any values,
    # skip this stat
    return unless $write_obj->check_column( 'gene',   'gene_feature5' );
    return unless $write_obj->check_column( 'exonsw', 'exonsw_id' );

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT AVG (gene_feature5)
            FROM gene
            WHERE gene_feature5 IS NOT NULL
            GROUP BY gene_stable_id
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @D_levels = (
        [ 1,     $quartiles->[0], $quartiles->[1] ],
        [ 2,     $quartiles->[1], $quartiles->[2] ],
        [ 3,     $quartiles->[2], $quartiles->[3] ],
        [ 4,     $quartiles->[3], $quartiles->[4] ],
        [ "all", 0,               1 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;
        my $sheet_name = "exon_D_$order";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ distance AVG_pi AVG_Indel/100bp AVG_gc AVG_CV COUNT };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT s.exonsw_distance `distance`,
                       AVG (w.window_pi) `avg_pi`,
                       AVG (w.window_indel / w.window_length * 100)
                       `avg_indel/100bp`,
                       AVG (w.window_target_gc) `avg_gc`,
                       AVG (s.exonsw_cv) `avg_cv`,
                       COUNT(*) count
                FROM exonsw s,
                     window w,
                     (
                      SELECT s.exonsw_id
                      FROM gene g,
                           exon e, 
                           exonsw s
                      WHERE g.gene_id = e.gene_id AND
                            g.gene_feature5 BETWEEN ? AND ? AND
                            e.exon_id = s.prev_exon_id AND
                            s.exonsw_type IN ('L', 'l')
                      UNION
                      SELECT s.exonsw_id
                      FROM gene g,
                           exon e, 
                           exonsw s
                      WHERE g.gene_id = e.gene_id AND
                            g.gene_feature5 BETWEEN ? AND ? AND
                            e.exon_id = s.exon_id AND
                            s.exonsw_type IN ('R', 'r')
                     ) sw
                WHERE s.window_id = w.window_id AND
                      s.exonsw_id = sw.exonsw_id
                GROUP BY s.exonsw_distance
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                bind_value =>
                    [ $low_border, $high_border, $low_border, $high_border ],
            );

            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@D_levels) {
        &$write_sheet(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- exon_D_null
#----------------------------------------------------------#
my $exon_D_null = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    return unless $write_obj->check_column( 'gene',   'gene_feature5' );
    return unless $write_obj->check_column( 'exonsw', 'exonsw_id' );

    my $sheet_name = "exon_D_null";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers
            = qw{ distance AVG_pi AVG_Indel/100bp AVG_gc AVG_CV COUNT };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    {    # query
        my $sql_query = q{
            SELECT s.exonsw_distance `distance`,
                   AVG (w.window_pi) `avg_pi`,
                   AVG (w.window_indel / w.window_length * 100)
                   `avg_indel/100bp`,
                   AVG (w.window_target_gc) `avg_gc`,
                   AVG (s.exonsw_cv) `avg_cv`,
                   COUNT(*) count
            FROM exonsw s,
                 window w,
                 (
                  SELECT s.exonsw_id
                  FROM gene g,
                       exon e, 
                       exonsw s
                  WHERE g.gene_id = e.gene_id AND
                        g.gene_feature5 IS NULL AND
                        e.exon_id = s.prev_exon_id AND
                        s.exonsw_type IN ('L', 'l')
                  UNION
                  SELECT s.exonsw_id
                  FROM gene g,
                       exon e, 
                       exonsw s
                  WHERE g.gene_id = e.gene_id AND
                        g.gene_feature5 IS NULL AND
                        e.exon_id = s.exon_id AND
                        s.exonsw_type IN ('R', 'r')
                 ) sw
            WHERE s.window_id = w.window_id AND
                  s.exonsw_id = sw.exonsw_id
            GROUP BY s.exonsw_distance
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
# worksheet -- exon_gc
#----------------------------------------------------------#
my $exon_gc = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    return unless $write_obj->check_column( 'exonsw', 'exonsw_id' );

    # find quartiles
    my $quartiles;
    {
        my $sql_query = q{
            SELECT window_target_gc
            FROM exonsw e,
                 window w
            WHERE e.window_id = w.window_id
        };
        my %option = ( sql_query => $sql_query, );
        $quartiles = $write_obj->quantile_sql( \%option, 4 );
    }

    my @exon_levels = (
        [ 1,     $quartiles->[0], $quartiles->[1] ],
        [ 2,     $quartiles->[1], $quartiles->[2] ],
        [ 3,     $quartiles->[2], $quartiles->[3] ],
        [ 4,     $quartiles->[3], $quartiles->[4] ],
        [ "all", 0,               1 ],
    );

    my $write_sheet = sub {
        my ( $order, $low_border, $high_border ) = @_;
        my $sheet_name = "exon_gc_$order";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ distance AVG_pi AVG_Indel/100bp AVG_gc AVG_CV COUNT };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG (s.exonsw_distance) `avg_distance`,
                       AVG (w.window_pi) `avg_pi`,
                       AVG (w.window_indel / w.window_length * 100)
                       `avg_indel/100bp`,
                       AVG (w.window_target_gc) `avg_gc`,
                       AVG (s.exonsw_cv) `avg_cv`,
                       count(*) count
                FROM exonsw s,
                     window w,
                     (
                      SELECT s.exonsw_id
                      FROM exon e,
                           window w,
                           exonsw s
                      WHERE e.window_id = w.window_id AND
                            w.window_target_gc BETWEEN ? AND ? AND
                            e.exon_id = s.exon_id AND
                            s.exonsw_type IN ('R', 'r')
                      UNION
                      SELECT s.exonsw_id
                      FROM exon e,
                           window w,
                           exonsw s
                      WHERE e.window_id = w.window_id AND
                            w.window_target_gc BETWEEN ? AND ? AND
                            e.exon_id = s.prev_exon_id AND
                            s.exonsw_type IN ('L', 'l')
                     ) sw
                WHERE s.window_id = w.window_id AND
                      s.exonsw_id = sw.exonsw_id
                GROUP BY s.exonsw_distance
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                bind_value =>
                    [ $low_border, $high_border, $low_border, $high_border ],
            );

            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@exon_levels) {
        &$write_sheet(@$_);
    }
};

#----------------------------------------------------------#
# worksheet -- exon_ess
#----------------------------------------------------------#
my $exon_ess = sub {

    # if the target column of the target table does not contain
    #   any values, skip this stat
    return unless $write_obj->check_column( 'exonsw', 'exonsw_id' );
    return unless $write_obj->check_column( 'gene',   'gene_feature4' );

    my @exon_levels = ( [ 'ess', 1 ], [ 'non_ess', 0 ], );

    my $write_sheet = sub {
        my ($exon_type) = @_;
        my $sheet_name = 'exon_ess' . "_$exon_type->[0]";
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ distance AVG_pi AVG_Indel/100bp AVG_gc AVG_CV COUNT };
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # query
            my $sql_query = q{
                SELECT s.exonsw_distance `distance`,
                       AVG (w.window_pi) `avg_pi`,
                       AVG (w.window_indel / w.window_length * 100)
                       `avg_indel/100bp`,
                       AVG (w.window_target_gc) `avg_gc`,
                       AVG (s.exonsw_cv) `avg_cv`,
                       count(*) count
                FROM exonsw s,
                     window w,
                    (
                     SELECT s.exonsw_id
                     FROM gene g,
                          exon e, 
                          exonsw s
                     WHERE g.gene_id = e.gene_id AND
                           g.gene_feature4 = ? AND
                           g.gene_biotype = 'protein_coding' AND
                           e.exon_id = s.prev_exon_id AND
                           s.exonsw_type IN ('L', 'l')
                     UNION
                     SELECT s.exonsw_id
                     FROM gene g,
                          exon e, 
                          exonsw s
                     WHERE g.gene_id = e.gene_id AND
                           g.gene_feature4 = ? AND
                           g.gene_biotype = 'protein_coding' AND
                           e.exon_id = s.exon_id AND
                           s.exonsw_type IN ('R', 'r')
                    ) sw
                WHERE s.window_id = w.window_id AND
                      s.exonsw_id = sw.exonsw_id AND
                      s.exonsw_type != 'S'
                GROUP BY s.exonsw_distance
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $exon_type->[1], $exon_type->[1] ],
            );

            ($sheet_row)
                = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been exonrated.\n";
    };

    foreach (@exon_levels) {
        my $exon_type = $_;
        &$write_sheet($exon_type);
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$summary_gene; next; }
    if ( $n == 2 ) { &$coding_all;   next; }
    if ( $n == 3 ) { &$exon_D;       &$exon_D_null; next; }
    if ( $n == 4 ) { &$exon_gc;      next; }
    if ( $n == 5 ) { &$exon_ess;     next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    gene_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    gene_stat_factory.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --output            output filename
        --run               run special analysis

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
