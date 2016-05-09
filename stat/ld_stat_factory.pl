#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
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

ld_stat_factory.pl - LD stats for alignDB

=head1 SYNOPSIS

    perl ld_stat_factory.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --outfile   -o  STR     outfile filename
        --freq          INT     count freq one by one to $max_freq
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
    'freq=i'       => \( my $max_freq ),
    'run|r=s'      => \( my $run      = $Config->{stat}{run} ),
    'combine=i'    => \( my $combine  = 0 ),
    'piece=i'      => \( my $piece    = 0 ),
    'replace=s'    => \my %replace,
    'index'        => \my $add_index_sheet,
    'chart'        => \my $add_chart,
) or Getopt::Long::HelpMessage(1);

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
    $outfile = "$db.ld.xlsx" unless $outfile;
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
        $outfile = "$db.ld.$runlist.xlsx";
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
my $dbh = $aligndb_obj->dbh;

my $write_obj = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
    replace => \%replace,
);

my $sql_file = AlignDB::SQL::Library->new( lib => "$FindBin::RealBin/sql.lib" );

# auto detect combine threshold
if ( $combine == 0 ) {
    ($combine) = $write_obj->calc_threshold;
}

# auto detect combine threshold
if ( $piece == 0 ) {
    ( undef, $piece ) = $write_obj->calc_threshold;
}

#----------------------------#
# count freq
#----------------------------#
my $seq_count      = $aligndb_obj->get_seq_count;
my $max_indel_freq = $aligndb_obj->get_max_indel_freq;

my @freqs;
if ($max_freq) {
    for ( 1 .. $max_freq - 1 ) {
        my $name = $_ . "of" . $seq_count;
        push @freqs, [ $name, $_, $_ ];
    }
}
else {

    # for 22 flies, low, mid, high
    {
        my @all_freqs = 1 .. $max_indel_freq;
        if ( scalar @all_freqs <= 3 ) {
            for (@all_freqs) {
                my $name = $_ . "of" . $seq_count;
                push @freqs, [ $name, $_, $_ ];
            }
        }
        else {
            my @to_be_combs = @all_freqs[ 0 .. $max_indel_freq - 1 ];
            my @chunks      = reverse apportion( scalar @to_be_combs, 3 );
            my @chunks_freq = multi_slice( \@to_be_combs, @chunks );
            for my $chunk (@chunks_freq) {
                if ( $chunk->[0] == $chunk->[-1] ) {
                    my $name = $chunk->[0] . "of" . $seq_count;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
                else {
                    my $name = join( '_', $chunk->[0], $chunk->[-1] ) . "of" . $seq_count;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
            }
        }
    }
}

print YAML::Syck::Dump [
    {   combine => $combine,
        piece   => $piece,
    },
    {   all_freq => $seq_count,
        freq     => \@freqs,
    }
];

#----------------------------------------------------------#
# chart -- d1_indel_ld
#----------------------------------------------------------#
my $chart_d1_indel_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %opt = (
        x_column    => 0,
        y_column    => 2,
        first_row   => 2,
        last_row    => 17,
        x_max_scale => 15,
        y_data      => $data->[2],
        x_title     => "Distance to indels (d1)",
        y_title     => "r^2",
        top         => 1,
        left        => 6,
    );
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 4;
    $opt{y_title}  = "|Dprime|";
    $opt{y_data}   = $data->[4];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column}  = 1;
    $opt{y_title}   = "r";
    $opt{y_data}    = $data->[1];
    $opt{y2_column} = 3;
    $opt{y2_data}   = $data->[3];
    $opt{y2_title}  = "Dprime";
    $opt{top}       = 1;
    $opt{left}      = 12;
    $write_obj->draw_2y( $sheet, \%opt );
};

#----------------------------------------------------------#
# chart -- d2_indel_ld
#----------------------------------------------------------#
my $chart_d2_indel_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %opt = (
        x_column    => 0,
        y_column    => 2,
        first_row   => 2,
        last_row    => 32,
        x_max_scale => 30,
        y_data      => $data->[2],
        x_title     => "Reciprocal of indel density (d2)",
        y_title     => "r^2",
        top         => 1,
        left        => 6,
    );
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 4;
    $opt{y_title}  = "|Dprime|";
    $opt{y_data}   = $data->[4];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column}  = 1;
    $opt{y_title}   = "r";
    $opt{y_data}    = $data->[1];
    $opt{y2_column} = 3;
    $opt{y2_data}   = $data->[3];
    $opt{y2_title}  = "Dprime";
    $opt{top}       = 1;
    $opt{left}      = 12;
    $write_obj->draw_2y( $sheet, \%opt );
};

#----------------------------------------------------------#
# chart -- d1_snps_ld
#----------------------------------------------------------#
my $chart_d1_snps_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %opt = (
        x_column    => 0,
        y_column    => 1,
        first_row   => 2,
        last_row    => 17,
        x_max_scale => 15,
        y_data      => $data->[1],
        x_title     => "Distance to indels (d1)",
        y_title     => "r^2 to nearest indel ",
        top         => 1,
        left        => 11,
    );
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 3;
    $opt{y_title}  = "|Dprime| to nearest indel ";
    $opt{y_data}   = $data->[3];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 2;
    $opt{y_title}  = "r^2 to near snps";
    $opt{y_data}   = $data->[2];
    $opt{top}      = 1;
    $opt{left}     = 17;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 5;
    $opt{y_title}  = "indel group snps r^2";
    $opt{y_data}   = $data->[5];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 6;
    $opt{y_title}  = "nonindel group snps r^2";
    $opt{y_data}   = $data->[6];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );
};

#----------------------------------------------------------#
# chart -- d2_snps_ld
#----------------------------------------------------------#
my $chart_d2_snps_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %opt = (
        x_column    => 0,
        y_column    => 1,
        first_row   => 2,
        last_row    => 32,
        x_max_scale => 30,
        y_data      => $data->[1],
        x_title     => "Reciprocal of indel density (d2)",
        y_title     => "r^2 to nearest indel ",
        top         => 1,
        left        => 11,
    );
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 3;
    $opt{y_title}  = "|Dprime| to nearest indel ";
    $opt{y_data}   = $data->[3];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 2;
    $opt{y_title}  = "r^2 to near snps";
    $opt{y_data}   = $data->[2];
    $opt{top}      = 1;
    $opt{left}     = 17;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 5;
    $opt{y_title}  = "indel group snps r^2";
    $opt{y_data}   = $data->[5];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column} = 6;
    $opt{y_title}  = "nonindel group snps r^2";
    $opt{y_data}   = $data->[6];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );
};

#----------------------------------------------------------#
# worksheet -- indel_ld
#----------------------------------------------------------#
my $indel_ld = sub {

    {
        my $sheet_name = 'd1_indel_ld';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->limit(20);

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
            $chart_d1_indel_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }

    {
        my $sheet_name = 'd2_indel_ld';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->limit(35);
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
            $chart_d2_indel_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }
};

my $indel_ld_insdel = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet_d1 = sub {
        my ($level) = @_;
        my $sheet_name = 'd1_indel_ld_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->limit(20);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $level->[1] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_d1_indel_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_indel_ld_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $level->[1] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_d2_indel_ld->( $sheet, $data );
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

my $snps_ld = sub {

    {
        my $sheet_name = 'd1_snps_ld';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->limit(20);

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
            $chart_d1_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }

    {
        my $sheet_name = 'd2_snps_ld';
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);

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
            $chart_d2_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }
};

my $snps_ld_insdel = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet_d1 = sub {
        my ($level) = @_;
        my $sheet_name = 'd1_snps_ld_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->limit(20);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $level->[1] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_d1_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_snps_ld_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $level->[1] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_d2_snps_ld->( $sheet, $data );
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

my $snps_ld_freq = sub {
    my @freq_levels = @freqs;

    my $write_sheet_d1 = sub {
        my ($level) = @_;
        my $sheet_name = 'd1_snps_ld_freq_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->limit(20);
        $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
        $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $level->[1], $level->[2] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_d1_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_snps_ld_freq_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);
        $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
        $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $level->[1], $level->[2] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_d2_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@freq_levels) {
        $write_sheet_d1->($_);
    }

    for (@freq_levels) {
        $write_sheet_d2->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- segment_gc_indel
#----------------------------------------------------------#
my $segment_gc_indel = sub {

    my @segment_levels = (3);

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_gc_indel' . "_$segment_type";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(1);

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %opt = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%opt );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_gc `gc`,
                           s.segment_gc_CV `cv`,
                           w.window_coding `coding`,
                           s.segment_r2_s `r2`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY gc DESC, r2 DESC, pi, indel
            };
            my %opt = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],
            );
            $write_obj->excute_sql( \%opt );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                SELECT t_id, length
                FROM tmp_group
            };
            my %opt = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%opt ) };
        }

        my @names
            = qw{AVG_gc AVG_pi AVG_Indel/100bp AVG_CV AVG_coding AVG_r2 AVG_length COUNT SUM_length};
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data = [];
        push @{$data}, [] for @names;
        {    # content
            my $sql_query = q{
                SELECT AVG(t.gc) `AVG_gc`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.cv) `AVG_CV`,
                       AVG(t.coding) `AVG_coding`,
                       AVG(t.r2) `AVG_r2`,
                       AVG(t.length) `AVG_length`,
                       COUNT(*) COUNT,
                       SUM(t.length) `SUM_length`
                FROM tmp_group t
                WHERE t_id IN
            };

            my @group_names;
            for (@combined_segment) {
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

                my $sth = $dbh->prepare($sql_query_in);
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
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %opt = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%opt );
        }

        if ($add_chart) {    # chart
            my %opt = (
                x_column  => 1,
                y_column  => 6,
                first_row => 1,
                last_row  => scalar @combined_segment,
                x_data    => $data->[0],
                y_data    => $data->[5],
                x_title   => "GC proportion",
                y_title   => "near snps r^2",
                top       => 1,
                left      => 11,
            );
            $write_obj->draw_xy( $sheet, \%opt );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@segment_levels) {
        &$write_sheet($_);
    }
};

my $segment_cv_indel = sub {

    my @segment_levels = (3);

    my $write_sheet = sub {
        my ($segment_type) = @_;
        my $sheet_name = 'segment_cv_indel' . "_$segment_type";
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(1);

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %opt = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%opt );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_gc `gc`,
                           s.segment_gc_CV `cv`,
                           w.window_coding `coding`,
                           s.segment_r2_s `r2`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY cv DESC, r2 DESC, pi, indel
            };
            my %opt = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],

            );
            $write_obj->excute_sql( \%opt );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                SELECT t_id, length
                FROM tmp_group
            };
            my %opt = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%opt ) };
        }

        my @names
            = qw{AVG_CV AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_r2 AVG_length COUNT SUM_length Range_gc};
        my $data = [];
        push @{$data}, [] for @names;

        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        {    # query
            my $sql_query = q{
                SELECT AVG(t.CV) `AVG_CV`,
                       AVG(t.pi) `AVG_pi`,
                       AVG(t.indel / t.length * 100) `AVG_Indel/100bp`,
                       AVG(t.gc) `AVG_gc`,
                       AVG(t.coding) `AVG_coding`,
                       AVG(t.r2) `AVG_r2`,
                       AVG(t.length) `AVG_length`,
                       COUNT(*) COUNT,
                       SUM(t.length) `SUM_length`,
                       MAX(t.gc) - MIN(t.gc) `Range_gc`
                FROM tmp_group t
                WHERE t_id IN
            };

            my @group_names;
            for (@combined_segment) {
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

                my $sth = $dbh->prepare($sql_query_in);
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
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %opt = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%opt );
        }

        if ($add_chart) {    # chart
            my %opt = (
                x_column  => 1,
                y_column  => 6,
                first_row => 1,
                last_row  => scalar @combined_segment,
                x_data    => $data->[0],
                y_data    => $data->[5],
                x_title   => "Segment CV",
                y_title   => "near snps r^2",
                top       => 1,
                left      => 12,
            );
            $write_obj->draw_xy( $sheet, \%opt );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@segment_levels) {
        &$write_sheet($_);
    }
};

for my $n (@tasks) {
    if ( $n == 1 ) { &$indel_ld;        next; }
    if ( $n == 2 ) { &$indel_ld_insdel; next; }

    if ( $n == 11 ) { &$snps_ld;        next; }
    if ( $n == 12 ) { &$snps_ld_insdel; next; }
    if ( $n == 13 ) { &$snps_ld_freq;   next; }

    if ( $n == 21 ) { &$segment_gc_indel; next; }
    if ( $n == 22 ) { &$segment_cv_indel; next; }
}

if ($add_index_sheet) {
    $write_obj->add_index_sheet;
    print "Sheet [INDEX] has been generated.\n";
}

$stopwatch->end_message;
exit;

# codes come from http://www.perlmonks.org/?node_id=516493
sub apportion {
    my ( $elements, $pieces ) = @_;
    my $small_chunk     = int $elements / $pieces;
    my $oversized_count = $elements % $pieces;
    ( ( 1 + $small_chunk ) x ($oversized_count), ($small_chunk) x ( $pieces - $oversized_count ) );
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
