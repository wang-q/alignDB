#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Statistics::R;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;
use AlignDB::ToXLSX;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

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
    'freq=i'       => \( my $max_freq ),
    'run|r=s'      => \( my $run      = $Config->{stat}{run} ),
    'combine=i'    => \( my $combine  = 0 ),
    'index'        => \( my $add_index_sheet, ),
    'chart'        => \( my $add_chart, ),
) or HelpMessage(1);

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
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

my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password )
    or die "Cannot connect to MySQL database at $db:$server";
my $write_obj = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
);

my $sql_file = AlignDB::SQL::Library->new( lib => "$FindBin::Bin/sql.lib" );

# auto detect combine threshold
if ( $combine == 0 ) {
    ($combine) = $write_obj->calc_threshold;
}

#----------------------------#
# count freq
#----------------------------#
my $largest    = get_freq($dbh);
my $seq_number = get_seq_number($dbh);

my @freqs;
if ($max_freq) {
    for ( 1 .. $max_freq - 1 ) {
        my $name = $_ . "of" . $largest;
        push @freqs, [ $name, $_, $_ ];
    }
}
else {

    # for 22 flies, low, mid, high
    {
        my @all_freqs = 1 .. $largest;
        if ( scalar @all_freqs <= 3 ) {
            for (@all_freqs) {
                my $name = $_ . "of" . $seq_number;
                push @freqs, [ $name, $_, $_ ];
            }
        }
        else {
            my @to_be_combs = @all_freqs[ 0 .. $largest - 1 ];
            my @chunks      = reverse apportion( scalar @to_be_combs, 3 );
            my @chunks_freq = multi_slice( \@to_be_combs, @chunks );
            for my $chunk (@chunks_freq) {
                if ( $chunk->[0] == $chunk->[-1] ) {
                    my $name = $chunk->[0] . "of" . $seq_number;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
                else {
                    my $name = join( '_', $chunk->[0], $chunk->[-1] ) . "of" . $seq_number;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
            }
        }
    }
}

print Dump [
    { combine => $combine, },
    {   largest    => $largest,
        seq_number => $seq_number,
        freq       => \@freqs,
    }
];

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
        last_row      => 7,
        x_max_scale   => 5,
        y_data        => [ map { $data->[$_] } 2 .. 3 ],
        x_title       => "Distance to indels (d1)",
        y_title       => "Nucleotide diversity",
        top           => 1,
        left          => 10,
    );
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column}      = 4;
    $opt{y_last_column} = 5;
    $opt{y_data}        = [ map { $data->[$_] } 4 .. 5 ], $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );

    $opt{y_column}      = 7;
    $opt{y_last_column} = 7;
    $opt{y_title}       = "Di/Dn";
    $opt{y_data}        = $data->[7];
    $opt{top} += 18;
    $write_obj->draw_y( $sheet, \%opt );
};

#----------------------------------------------------------#
# chart -- pigccv_*
#----------------------------------------------------------#
my $chart_pigccv = sub {
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
        y_title     => "Nucleotide diversity",
        top         => 1,
        left        => 10,
    );
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

#----------------------------------------------------------#
# worksheet -- summary_indel
#----------------------------------------------------------#
my $summary_indel = sub {
    my $sheet_name = 'summary_indel';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @names = qw{AVG MIN MAX STD COUNT SUM};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
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
        my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ] );

        # all indels
        $column_stat->( 'all', 'indel i', 'i.indel_length', );

        for my $level (@freq_levels) {
            my ( $name, $lower, $upper ) = @{$level};

            $column_stat->(
                'all_' . $name, 'indel i',
                'i.indel_length',
                qq{WHERE i.indel_freq >= $lower
                    AND i.indel_freq <= $upper}
            );
        }
        $write_obj->increase_row;

        # insertions
        $column_stat->( 'ins', 'indel i', 'i.indel_length', q{WHERE i.indel_type = 'I'} );

        for my $level (@freq_levels) {
            my ( $name, $lower, $upper ) = @{$level};

            $column_stat->(
                'ins_' . $name, 'indel i',
                'i.indel_length',
                qq{WHERE i.indel_freq >= $lower
                    AND i.indel_freq <= $upper
                    AND i.indel_type = 'I'}
            );
        }
        $write_obj->increase_row;

        # deletions
        $column_stat->( 'del', 'indel i', 'i.indel_length', q{WHERE i.indel_type = 'D'} );

        for my $level (@freq_levels) {
            my ( $name, $lower, $upper ) = @{$level};

            $column_stat->(
                'del_' . $name, 'indel i',
                'i.indel_length',
                qq{WHERE i.indel_freq >= $lower
                    AND i.indel_freq <= $upper
                    AND i.indel_type = 'D'}
            );
        }
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance
#----------------------------------------------------------#
my $distance = sub {
    my $sheet_name = 'distance';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my $thaw_sql = $sql_file->retrieve('multi-distance-0');

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
        $chart_distance->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance2
#----------------------------------------------------------#
my $distance2 = sub {
    my $sheet_name = 'distance2';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my $thaw_sql = $sql_file->retrieve('multi-distance2-0');

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
        $chart_distance->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance3
#----------------------------------------------------------#
my $distance3 = sub {
    my $sheet_name = 'distance3';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my $thaw_sql = $sql_file->retrieve('multi-distance3-0');

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
        $chart_distance->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance_length
#----------------------------------------------------------#
my $distance_length = sub {
    my @length_levels
        = ( [ '1-5', 1, 5 ], [ '6-10', 6, 10 ], [ '11-50', 11, 50 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'distance_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('multi-distance-0');
        $thaw_sql->add_where( 'indel.indel_length' => { op => '>=', value => '1' } );
        $thaw_sql->add_where( 'indel.indel_length' => { op => '<=', value => '1' } );

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
            $chart_distance->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@length_levels) {
        $write_sheet->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- distance(frequecy)
#----------------------------------------------------------#
my $frequency_distance = sub {
    my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'distance_freq_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('multi-distance-0');
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '>=', value => '1' } );
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '<=', value => '1' } );

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
            $chart_distance->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@freq_levels) {
        $write_sheet->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- distance(ins)
#----------------------------------------------------------#
my $distance_insdel = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'distance_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('multi-distance-0');
        $thaw_sql->add_where( 'indel.indel_type' => { op => '=', value => '1' } );

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
            $chart_distance->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@type_levels) {
        $write_sheet->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- distance(low frequecy)
#----------------------------------------------------------#
my $distance_insdel_freq = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );
    my @freq_levels = @freqs;

    my $write_sheet = sub {
        my ( $type, $freq ) = @_;
        my $sheet_name = 'distance_' . $type->[0] . '_' . $freq->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('multi-distance-0');
        $thaw_sql->add_where( 'indel.indel_type' => { op => '=',  value => '1' } );
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '>=', value => '1' } );
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '<=', value => '1' } );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        my $data;
        {    # content
            $data = $write_obj->write_sql(
                $sheet,
                {   sql_query  => $thaw_sql->as_sql,
                    bind_value => [ $type->[1], $freq->[1], $freq->[2] ],
                    data       => 1,
                }
            );
        }

        if ($add_chart) {    # chart
            $chart_distance->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for my $type (@type_levels) {
        for my $freq (@freq_levels) {
            $write_sheet->( $type, $freq );
        }
    }
};

#----------------------------------------------------------#
# worksheet -- indel_length
#----------------------------------------------------------#
my $indel_length = sub {
    my $sheet_name = 'indel_length';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(0);

    my $thaw_sql = $sql_file->retrieve('multi-indel_length-0');

    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # content
        my $sth = $dbh->prepare( $thaw_sql->as_sql );
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
# worksheet -- indel_length_freq
#----------------------------------------------------------#
my $indel_length_freq = sub {

    my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'indel_length_freq_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('multi-indel_length-0');
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '>=', value => '1' } );
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '<=', value => '1' } );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        {    # content
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1], $level->[2], );

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

    for (@freq_levels) {
        $write_sheet->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_length_insdel
#----------------------------------------------------------#
my $indel_length_insdel = sub {

    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'indel_length_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('multi-indel_length-0');
        $thaw_sql->add_where( 'indel.indel_type' => { op => '=', value => 'I' } );

        my @names = $thaw_sql->as_header;
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        {    # content
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1] );

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

    for (@type_levels) {
        $write_sheet->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_type_gc_10
#----------------------------------------------------------#
my $indel_type_gc_10 = sub {
    my $sheet_name = 'indel_type_gc_10';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # indel_type groups
    my @indel_types = ( [ 'Insertion', 'I' ], [ 'Deletion', 'D' ], );

    my $sql_query = q{
        SELECT 
            indel_length indel_length,
            AVG(indel_gc) AVG_indel_gc,
            COUNT(*) `COUNT`,
            STD(indel_gc) STD_indel_gc
        FROM
            indel
        WHERE
            indel_type = ? AND indel_length <= 10
        GROUP BY indel_length
    };

    my @names = $write_obj->sql2names( $sql_query, { bind_value => [ $indel_types[0]->[1] ] } );
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    my $data;
    for (@indel_types) {
        $write_obj->increase_row;

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $sql_query,
                query_name => $_->[0],
                bind_value => [ $_->[1] ],
                data       => $data,
            }
        );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_type_gc_100
#----------------------------------------------------------#
my $indel_type_gc_100 = sub {
    my $sheet_name = 'indel_type_gc_100';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    # indel_type groups
    my @indel_types = ( [ 'Insertion', 'I' ], [ 'Deletion', 'D' ], );

    my $sql_query = q{
        SELECT 
            CEIL(indel_length / 10) * 10,
            AVG(indel_gc) AVG_indel_gc,
            COUNT(*),
            STD(indel_gc) STD_indel_gc
        FROM
            indel
        WHERE
            indel_type = ? AND indel_length <= 100
        GROUP BY CEIL(indel_length / 10)
    };

    my @names = $write_obj->sql2names( $sql_query, { bind_value => [ $indel_types[0]->[1] ] } );
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    # contents
    my $data;
    for (@indel_types) {
        $write_obj->increase_row;

        $data = $write_obj->write_sql(
            $sheet,
            {   sql_query  => $sql_query,
                query_name => $_->[0],
                bind_value => [ $_->[1] ],
                data       => $data,
            }
        );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- pigccv_combined
#----------------------------------------------------------#
my $combined_pigccv = sub {
    my $sheet_name = 'pigccv_combined';
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

    my $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');

    my @names = $thaw_sql->as_header;
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    my $data;
    for my $comb ( @{$combined} ) {    # content
        my $thaw_sql = $sql_file->retrieve('common-d1_comb_pi_gc_cv-0');
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
        $chart_pigccv->( $sheet, $data );
    }

    print "Sheet [$sheet_name] has been generated.\n";

};

#----------------------------------------------------------#
# worksheet -- pigccv_freq
#----------------------------------------------------------#
my $frequency_pigccv = sub {
    my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'pigccv_freq_' . $level->[0];
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my $thaw_sql = $sql_file->retrieve('common-d1_pi_gc_cv-0');
        $thaw_sql->from( [] );
        $thaw_sql->add_join(
            isw => [
                {   type      => 'inner',
                    table     => 'indel',
                    condition => 'isw.isw_indel_id = indel.indel_id'
                },
            ]
        );
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '>=', value => '1' } );
        $thaw_sql->add_where( 'indel.indel_freq' => { op => '<=', value => '1' } );

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
            $chart_pigccv->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    for (@freq_levels) {
        $write_sheet->($_);
    }
};

#----------------------------------------------------------#
# worksheet -- di_dn_ttest
#----------------------------------------------------------#
my $di_dn_ttest = sub {
    my $sheet_name = 'di_dn_ttest';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @group_distance = ( [0], [1], [2], [3], [4], [ 5 .. 10 ], [ 2 .. 5 ] );

    my @names = qw{AVG_distance AVG_D AVG_Di AVG_Dn Di/Dn COUNT P_value};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # content
        my $base_query = qq{
            SELECT  AVG(isw_distance) AVG_distance,
                    AVG(isw_pi) AVG_D,
                    AVG(isw_d_indel) AVG_Di,
                    AVG(isw_d_noindel) AVG_Dn,
                    AVG(isw_d_indel)/AVG(isw_d_noindel) `Di/Dn`,
                    COUNT(*) COUNT
            FROM    isw s, indel i
            WHERE s.isw_indel_id = i.indel_id
            AND s.isw_distance >= 0
            AND s.isw_d_indel IS NOT NULL 
            AND s.isw_distance IN 
        };

        for (@group_distance) {
            my @range     = @$_;
            my $in_list   = '(' . join( ',', @range ) . ')';
            my $sql_query = $base_query . $in_list;
            my $group_name;
            if ( scalar @range > 1 ) {
                $group_name = $range[0] . "--" . $range[-1];
            }
            else {
                $group_name = $range[0];
            }
            $write_obj->write_sql(
                $sheet,
                {   sql_query  => $sql_query,
                    query_name => $group_name,
                }
            );
        }
    }

    {    # t-test
        my $base_query = qq{
            SELECT s.isw_d_indel, s.isw_d_noindel
            FROM    isw s, indel i
            WHERE i.indel_id = s.isw_indel_id
            AND s.isw_distance >= 0
            AND s.isw_d_indel IS NOT NULL 
            AND (s.isw_d_indel + s.isw_d_noindel != 0)
            AND s.isw_distance IN 
        };

        my @p_values;
        for (@group_distance) {
            my @range     = @$_;
            my $in_list   = '(' . join( ',', @range ) . ')';
            my $sql_query = $base_query . $in_list;

            my $sth = $dbh->prepare($sql_query);
            $sth->execute;

            my ( @sample1, @sample2 );
            while ( my @row = $sth->fetchrow_array ) {
                push @sample1, $row[0];
                push @sample2, $row[1];
            }

            my $p_value = r_t_test( \@sample1, \@sample2 );
            if ( !defined $p_value ) {
                $p_value = "NA";
            }
            push @p_values, [$p_value];
        }

        $sheet->write( 1, scalar @names, [ [@p_values] ], $write_obj->format->{NORMAL} );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- di_dn_ttest
#----------------------------------------------------------#
my $di_dn_ttest_ns = sub {
    my $sheet_name = 'di_dn_ttest_ns';
    my $sheet;
    $write_obj->row(0);
    $write_obj->column(1);

    my @group_distance = ( [0], [1], [2], [3], [4], [ 5 .. 10 ], [ 2 .. 5 ] );

    my @names = qw{AVG_distance AVG_D AVG_Di AVG_Dn Di/Dn COUNT P_value};
    {    # header
        $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
    }

    {    # content
        my $base_query = qq{
            SELECT  AVG(isw_distance) AVG_distance,
                    AVG(isw_pi) AVG_D,
                    AVG(isw_d_indel) AVG_Di,
                    AVG(isw_d_noindel) AVG_Dn,
                    AVG(isw_d_indel)/AVG(isw_d_noindel) `Di/Dn`,
                    COUNT(*) COUNT
            FROM    isw s, indel i
            WHERE s.isw_indel_id = i.indel_id
            AND s.isw_distance >= 0
            AND s.isw_d_indel IS NOT NULL
            AND i.indel_slippage = 0
            AND (s.isw_d_indel + s.isw_d_noindel != 0)
            AND s.isw_distance IN 
        };

        for (@group_distance) {
            my @range     = @$_;
            my $in_list   = '(' . join( ',', @range ) . ')';
            my $sql_query = $base_query . $in_list;
            my $group_name;
            if ( scalar @range > 1 ) {
                $group_name = $range[0] . "--" . $range[-1];
            }
            else {
                $group_name = $range[0];
            }
            $write_obj->write_sql(
                $sheet,
                {   sql_query  => $sql_query,
                    query_name => $group_name,
                }
            );
        }
    }

    {    # t-test
        my $base_query = qq{
            SELECT s.isw_d_indel, s.isw_d_noindel
            FROM    isw s, indel i
            WHERE i.indel_id = s.isw_indel_id
            AND s.isw_distance >= 0
            AND s.isw_d_indel IS NOT NULL
            AND i.indel_slippage = 0
            AND (s.isw_d_indel + s.isw_d_noindel != 0)
            AND s.isw_distance IN 
        };

        my @p_values;
        for (@group_distance) {
            my @range     = @$_;
            my $in_list   = '(' . join( ',', @range ) . ')';
            my $sql_query = $base_query . $in_list;

            my $sth = $dbh->prepare($sql_query);
            $sth->execute;

            my ( @sample1, @sample2 );
            while ( my @row = $sth->fetchrow_array ) {
                push @sample1, $row[0];
                push @sample2, $row[1];
            }

            my $p_value = r_t_test( \@sample1, \@sample2 );
            if ( !defined $p_value ) {
                $p_value = "NA";
            }
            push @p_values, [$p_value];
        }

        $sheet->write( 1, scalar @names, [ [@p_values] ], $write_obj->format->{NORMAL} );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

for my $n (@tasks) {
    if ( $n == 1 )  { &$summary_indel;        next; }
    if ( $n == 2 )  { &$distance;             &$distance_length; next; }
    if ( $n == 3 )  { &$frequency_distance;   next; }
    if ( $n == 4 )  { &$distance2;            next; }
    if ( $n == 5 )  { &$distance3;            next; }
    if ( $n == 6 )  { &$distance_insdel;      next; }
    if ( $n == 7 )  { &$distance_insdel_freq; next; }
    if ( $n == 8 )  { &$indel_length;         next; }
    if ( $n == 9 )  { &$indel_length_freq;    next; }
    if ( $n == 10 ) { &$indel_length_insdel;  next; }
    if ( $n == 11 ) { &$indel_type_gc_10;     &$indel_type_gc_100; next; }
    if ( $n == 21 ) { &$combined_pigccv;      next; }
    if ( $n == 22 ) { &$frequency_pigccv;     next; }
    if ( $n == 52 ) { &$di_dn_ttest;          &$di_dn_ttest_ns; next; }
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

sub get_seq_number {
    my $dbh = shift;

    my $sql_query = q{
        SELECT DISTINCT
            COUNT(q.query_id) + 1
        FROM
            query q,
            sequence s
        WHERE
            q.seq_id = s.seq_id
        GROUP BY s.align_id
    };
    my $sth = $dbh->prepare($sql_query);

    my @counts;
    $sth->execute;
    while ( my ($count) = $sth->fetchrow_array ) {
        push @counts, $count;
    }
    if ( scalar @counts > 1 ) {
        warn "Database with non-consistent query numbers!\n";
    }

    return $counts[0];
}

sub get_freq {
    my $dbh = shift;

    my $sql_query = q{
        SELECT 
            MAX(indel_freq)
        FROM
            indel
    };
    my $sth = $dbh->prepare($sql_query);

    my @counts;
    $sth->execute;
    while ( my ($count) = $sth->fetchrow_array ) {
        push @counts, $count;
    }

    return $counts[0];
}

# t-test using R
sub r_t_test {
    my $x = shift;
    my $y = shift;

    die "Give two array-refs to me\n" if ref $x ne 'ARRAY';
    die "Give two array-refs to me\n" if ref $y ne 'ARRAY';
    die "Variable lengths differ\n"   if @$x != @$y;
    return                            if @$x <= 2;

    # Create a communication bridge with R and start R
    my $R = Statistics::R->new;

    $R->set( 'x', $x );
    $R->set( 'y', $y );
    $R->run(q{ result = t.test(x, y, alternative = "greater", paired = TRUE) });
    $R->run(q{ p_value <- result$p.value });

    my $p_value = $R->get('p_value');
    $R->stop;

    return $p_value;
}

__END__
