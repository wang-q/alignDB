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
my $run     = 'all';
my $combine = 0;
my $outfile = "";

my $max_freq;    # count freq one by one to $max_freq

my $help = 0;
my $man  = 0;

GetOptions(
    'help|?'         => \$help,
    'man'            => \$man,
    's|server=s'     => \$server,
    'P|port=s'       => \$port,
    'u|username=s'   => \$username,
    'p|password=s'   => \$password,
    'd|db=s'         => \$db,
    'o|output=s'     => \$outfile,
    'max|max_freq=s' => \$max_freq,
    'r|run=s'        => \$run,
    'cb|combine=i'   => \$combine,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
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
my $all_freq = $write_obj->get_freq;

if ( $all_freq < 2 ) {
    die "all_freq is $all_freq, are you sure this is a AlignDB database?\n";
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
                if ( $chunk->[0] == $chunk->[-1] ) {
                    my $name = $chunk->[0] . "of" . $all_freq;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
                else {
                    my $name = join( '_', $chunk->[0], $chunk->[-1] ) . "of"
                        . $all_freq;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
            }

            for ( $all_freq - 2, $all_freq - 1 ) {
                my $name = $_ . "of" . $all_freq;
                push @freqs, [ $name, $_, $_ ];
            }
        }
    }
}

#----------------------------------------------------------#
# worksheet -- summary_indel
#----------------------------------------------------------#
my $summary_indel = sub {
    my $sheet_name = 'summary_indel';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers    = qw{AVG MIN MAX STD COUNT SUM};
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

        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );

        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    };

    # write contents
    {

        my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ] );

        # all indels
        $column_stat->( 'all', 'indel i', 'i.indel_length', );

        foreach my $level (@freq_levels) {
            my ( $name, $lower, $upper ) = @$level;

            $column_stat->(
                'all_' . $name, 'indel i',
                'i.indel_length',
                qq{WHERE i.indel_freq >= $lower
                    AND i.indel_freq <= $upper}
            );
        }

        $sheet_row++;

        # insertions
        $column_stat->(
            'ins', 'indel i', 'i.indel_length', q{WHERE i.indel_type = 'I'}
        );

        foreach my $level (@freq_levels) {
            my ( $name, $lower, $upper ) = @$level;

            $column_stat->(
                'ins_' . $name, 'indel i',
                'i.indel_length',
                qq{WHERE i.indel_freq >= $lower
                    AND i.indel_freq <= $upper
                    AND i.indel_type = 'I'}
            );
        }

        $sheet_row++;

        # deletions
        $column_stat->(
            'del', 'indel i', 'i.indel_length', q{WHERE i.indel_type = 'D'}
        );

        foreach my $level (@freq_levels) {
            my ( $name, $lower, $upper ) = @$level;

            $column_stat->(
                'del_' . $name, 'indel i',
                'i.indel_length',
                qq{WHERE i.indel_freq >= $lower
                    AND i.indel_freq <= $upper
                    AND i.indel_type = 'D'}
            );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance
#----------------------------------------------------------#
my $distance = sub {
    my $sheet_name = 'distance';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_D AVG_Di AVG_Dn
            AVG_Dbii/2 AVG_Dbnn/2 AVG_Dc Di/Dn COUNT};
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
        my $thaw_sql = $sql_file->retrieve('multi-distance-0');
        my %option   = (
            sql_query => $thaw_sql->as_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance2
#----------------------------------------------------------#
my $distance2 = sub {
    my $sheet_name = 'distance2';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_D AVG_Di2 AVG_Dn2
            AVG_Dbii2/2 AVG_Dbnn2/2 AVG_Dc2 Di2/Dn2 COUNT};
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
        my $thaw_sql = $sql_file->retrieve('multi-distance2-0');
        my %option   = (
            sql_query => $thaw_sql->as_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- distance3
#----------------------------------------------------------#
my $distance3 = sub {
    my $sheet_name = 'distance3';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{distance AVG_D AVG_Di3 AVG_Dn3
            AVG_Dbii3/2 AVG_Dbnn3/2 AVG_Dc3 Di3/Dn3 COUNT};
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
        my $thaw_sql = $sql_file->retrieve('multi-distance3-0');
        my %option   = (
            sql_query => $thaw_sql->as_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
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
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_D AVG_Di AVG_Dn
                AVG_Dbii/2 AVG_Dbnn/2 AVG_Dc Di/Dn COUNT};
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
            my $thaw_sql = $sql_file->retrieve('multi-distance-0');
            $thaw_sql->add_where(
                'indel.indel_length' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'indel.indel_length' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $level->[1], $level->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@length_levels) {
        &$write_sheet($_);
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
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_D AVG_Di AVG_Dn
                AVG_Dbii/2 AVG_Dbnn/2 AVG_Dc Di/Dn COUNT};
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
            my $thaw_sql = $sql_file->retrieve('multi-distance-0');
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
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

#----------------------------------------------------------#
# worksheet -- distance(frequecy)
#----------------------------------------------------------#
my $frequency_distance2 = sub {
    my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'distance2_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_D AVG_Di AVG_Dn
                AVG_Dbii/2 AVG_Dbnn/2 AVG_Dc Di/Dn COUNT};
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
            my $thaw_sql = $sql_file->retrieve('multi-distance2-0');
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
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

#----------------------------------------------------------#
# worksheet -- distance(frequecy)
#----------------------------------------------------------#
my $frequency_distance3 = sub {
    my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'distance3_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_D AVG_Di AVG_Dn
                AVG_Dbii/2 AVG_Dbnn/2 AVG_Dc Di/Dn COUNT};
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
            my $thaw_sql = $sql_file->retrieve('multi-distance3-0');
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
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

#----------------------------------------------------------#
# worksheet -- distance(ins)
#----------------------------------------------------------#
my $distance_insdel = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'distance_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_D AVG_Di AVG_Dn
                AVG_Dbii/2 AVG_Dbnn/2 AVG_Dc Di/Dn COUNT};
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
            my $thaw_sql = $sql_file->retrieve('multi-distance-0');
            $thaw_sql->add_where(
                'indel.indel_type' => { op => '=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
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
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_D AVG_Di AVG_Dn
                AVG_Dbii/2 AVG_Dbnn/2 AVG_Dc Di/Dn COUNT};
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
            my $thaw_sql = $sql_file->retrieve('multi-distance-0');
            $thaw_sql->add_where(
                'indel.indel_type' => { op => '=', value => '1' } );
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '<=', value => '1' } );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $type->[1], $freq->[1], $freq->[2], ]
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach my $type (@type_levels) {
        foreach my $freq (@freq_levels) {
            &$write_sheet( $type, $freq );
        }
    }
};

#----------------------------------------------------------#
# worksheet -- indel_length
#----------------------------------------------------------#
my $indel_length = sub {
    my $sheet_name = 'indel_length';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{indel_length indel_number AVG_gc indel_sum};
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
        my $thaw_sql = $sql_file->retrieve('multi-indel_length-0');

        my %option = (
            sql_query => $thaw_sql->as_sql,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_highlight( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_length
#----------------------------------------------------------#
my $indel_length_freq = sub {

    my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'indel_length_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{indel_length indel_number AVG_gc indel_sum};
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
            my $thaw_sql = $sql_file->retrieve('multi-indel_length-0');
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '>=', value => '1' } );
            $thaw_sql->add_where(
                'indel.indel_freq' => { op => '<=', value => '1' } );

            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $level->[1], $level->[2], ]
            );
            ($sheet_row)
                = $write_obj->write_content_highlight( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@freq_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_length
#----------------------------------------------------------#
my $indel_length_insdel = sub {

    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'indel_length_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{indel_length indel_number AVG_gc indel_sum};
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
            my $thaw_sql = $sql_file->retrieve('multi-indel_length-0');
            $thaw_sql->add_where(
                'indel.indel_type' => { op => '=', value => 'I' } );

            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $level->[1] ]
            );
            ($sheet_row)
                = $write_obj->write_content_highlight( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@type_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- indel_type_gc_10
#----------------------------------------------------------#
my $indel_type_gc_10 = sub {
    my $sheet_name = 'indel_type_gc_10';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # indel_type groups
    my @indel_types = ( [ 'Insertion', ['I'] ], [ 'Deletion', ['D'] ], );

    {    # write header
        my $query_name = 'indel_type_gc_10';
        my @headers    = qw{indel_length AVG_gc COUNT STD_gc };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    for (@indel_types) {
        $sheet_row++;
        my $sql_query = q{
            SELECT  indel_length,
                    AVG(indel_gc) AVG_indel_gc,
                    COUNT(*),
                    STD(indel_gc) STD_indel_gc
            FROM indel
            WHERE indel_type = ?
            AND indel_length <= 10
            GROUP BY indel_length
        };
        my %option = (
            sql_query  => $sql_query,
            query_name => $_->[0],
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            bind_value => $_->[1],
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- indel_type_gc_100
#----------------------------------------------------------#
my $indel_type_gc_100 = sub {
    my $sheet_name = 'indel_type_gc_100';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    # indel_type groups
    my @indel_types
        = ( [ 'Insertion', ['I'] ], [ 'Deletion', ['D'] ], );

    {    # write header
        my $query_name = 'indel_type_gc_100';
        my @headers    = qw{indel_length AVG_gc COUNT STD_gc };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            header     => \@headers,
            query_name => $query_name,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_direct( $sheet_name, \%option );
    }

    # write contents
    for (@indel_types) {
        $sheet_row++;
        my $sql_query = q{
            # indel_length distribution
            SELECT  CEIL(indel_length / 10) * 10,
                    AVG(indel_gc) AVG_indel_gc,
                    COUNT(*),
                    STD(indel_gc) STD_indel_gc
            FROM indel
            WHERE indel_type = ?
            AND indel_length <= 100
            GROUP BY CEIL(indel_length / 10)
        };
        my %option = (
            sql_query  => $sql_query,
            query_name => $_->[0],
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            bind_value => $_->[1],
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

my $combined_pigccv = sub {

    #----------------------------------------------------------#
    # worksheet -- combined_pigccv
    #----------------------------------------------------------#

    # make combine
    my @combined;
    {
        my $thaw_sql   = $sql_file->retrieve('common-d1_combine-0');
        my $standalone = [ -1, 0 ];
        my %option     = (
            sql_query  => $thaw_sql->as_sql,
            threshold  => $combine,
            standalone => $standalone,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }
    {
        my $sheet_name = 'combined_pigccv';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{AVG_distance AVG_pi STD_pi AVG_gc STD_gc
                AVG_cv STD_cv COUNT};
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
                    AVG(isw_distance) AVG_distance,
                    AVG(isw_pi) AVG_pi,
                    STD(isw_pi) STD_pi,
                    AVG(isw_average_gc) AVG_gc,
                    STD(isw_average_gc) STD_gc,
                    AVG(isw_cv) AVG_cv,
                    STD(isw_cv) STD_cv,
                    COUNT(*) COUNT
                FROM isw
                WHERE isw_distance IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                combined  => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    #----------------------------------------------------------#
    # worksheet -- combined_pure_coding
    #----------------------------------------------------------#

    # make combine
    @combined = ();
    {
        my $sql_query = q{
                SELECT  i.isw_distance,
                        COUNT(*) COUNT
                FROM isw i
                inner join indel on i.isw_indel_id = indel.indel_id and indel.indel_coding = 1
                where i.isw_coding = 1
                group by i.isw_distance
            };
        my $standalone = [ -1, 0 ];
        my %option = (
            sql_query  => $sql_query,
            threshold  => $combine,
            standalone => $standalone,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }
    {
        my $sheet_name = 'combined_pure_coding';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{AVG_distance AVG_pi STD_pi AVG_gc STD_gc
                AVG_cv STD_cv COUNT};
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
                    AVG(isw_distance) AVG_distance,
                    AVG(isw_pi) AVG_pi,
                    STD(isw_pi) STD_pi,
                    AVG(isw_average_gc) AVG_gc,
                    STD(isw_average_gc) STD_gc,
                    AVG(isw_cv) AVG_cv,
                    STD(isw_cv) STD_cv,
                    COUNT(*) COUNT
                FROM isw i
                INNER JOIN indel ON i.isw_indel_id = indel.indel_id AND indel.indel_coding = 1
                WHERE i.isw_coding = 1
                AND isw_distance IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                combined  => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    #----------------------------------------------------------#
    # worksheet -- combined_pure_noncoding
    #----------------------------------------------------------#

    # make combine
    @combined = ();
    {
        my $sql_query = q{
                SELECT  i.isw_distance,
                        COUNT(*) COUNT
                FROM isw i
                inner join indel on i.isw_indel_id = indel.indel_id and indel.indel_coding = 0
                where i.isw_coding = 0
                group by i.isw_distance
            };
        my $standalone = [ -1, 0 ];
        my %option = (
            sql_query  => $sql_query,
            threshold  => $combine,
            standalone => $standalone,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }
    {
        my $sheet_name = 'combined_pure_noncoding';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{AVG_distance AVG_pi STD_pi AVG_gc STD_gc
                AVG_cv STD_cv COUNT};
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
                    AVG(isw_distance) AVG_distance,
                    AVG(isw_pi) AVG_pi,
                    STD(isw_pi) STD_pi,
                    AVG(isw_average_gc) AVG_gc,
                    STD(isw_average_gc) STD_gc,
                    AVG(isw_cv) AVG_cv,
                    STD(isw_cv) STD_cv,
                    COUNT(*) COUNT
                FROM isw i
                INNER JOIN indel ON i.isw_indel_id = indel.indel_id AND indel.indel_coding = 0
                WHERE i.isw_coding = 0
                AND isw_distance IN
            };
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                combined  => \@combined,
            );
            ($sheet_row)
                = $write_obj->write_content_combine( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

};

#----------------------------------------------------------#
# worksheet -- di_dn_ttest (frequecy)
#----------------------------------------------------------#
my $di_dn_ttest = sub {
    my @freq_levels = @freqs;

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'di_dn_ttest_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{AVG_distance AVG_D AVG_Di AVG_Dn
                Di-Dn Di/Dn COUNT P_value};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        my @group_distance = ( [0], [1], [2], [3], [4], [ 5 .. 10 ], );

        my @p_values;

        # ttest
        {
            my $sql_query = qq{
                SELECT s.isw_d_indel, s.isw_d_noindel
                FROM    isw s, indel i
                WHERE i.indel_id = s.isw_indel_id
                AND i.indel_freq = ?
                AND s.isw_distance >= 0
                AND s.isw_d_indel IS NOT NULL 
                AND s.isw_distance IN 
            };

            foreach (@group_distance) {
                my @range      = @$_;
                my $in_list    = '(' . join( ',', @range ) . ')';
                my $sql_query2 = $sql_query . $in_list;
                $sql_query2 .= " \nLIMIT 10000"; # prevent to exceed excel limit

                my %option = (
                    sql_query  => $sql_query2,
                    bind_value => [ $level->[1] ],
                );
                my $ttest   = $write_obj->column_ttest( \%option );
                my $p_value = $ttest->{t_prob};

                unless ( defined $p_value ) {
                    $p_value = "NA";

                    my $ttest_sheet_name
                        = "ttest_"
                        . $level->[1] . "of"
                        . $all_freq . "_"
                        . join( '-', @range[ 0, -1 ] );
                    my $stat_sheet;

                    my ( $stat_sheet_row, $stat_sheet_col ) = ( 0, 0 );
                    my %option = (
                        sheet_row => $stat_sheet_row,
                        sheet_col => $stat_sheet_col,
                        header    => [qw{d_indel d_noindel}],
                    );
                    ( $stat_sheet, $stat_sheet_row )
                        = $write_obj->write_header_direct( $ttest_sheet_name,
                        \%option );

                    %option = (
                        sql_query  => $sql_query2,
                        sheet_row  => $stat_sheet_row,
                        sheet_col  => $stat_sheet_col,
                        bind_value => [ $level->[1] ],
                    );
                    $write_obj->write_content_direct( $stat_sheet, \%option );

                }
                push @p_values, [$p_value];
            }
        }

        {
            my $sql_query = qq{
                SELECT  AVG(isw_distance) AVG_distance,
                        AVG(isw_pi) AVG_D,
                        AVG(isw_d_indel) AVG_Di,
                        AVG(isw_d_noindel) AVG_Dn,
                        AVG(isw_d_indel) - AVG(isw_d_noindel)  `Di-Dn`,
                        AVG(isw_d_indel)/AVG(isw_d_noindel) `Di/Dn`,
                        COUNT(*) COUNT
                FROM    isw s, indel i
                WHERE s.isw_indel_id = i.indel_id
                AND i.indel_freq = ?
                AND s.isw_distance >= 0
                AND s.isw_d_indel IS NOT NULL 
                AND s.isw_distance IN 
            };
            my %option = (
                sql_query     => $sql_query,
                sheet_row     => $sheet_row,
                sheet_col     => $sheet_col,
                group         => \@group_distance,
                append_column => \@p_values,
                bind_value    => [ $level->[1] ],
            );
            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@freq_levels) {
        &$write_sheet($_);
    }
};

#----------------------------------------------------------#
# worksheet -- distance(frequecy)
#----------------------------------------------------------#
my $frequency_pigccv = sub {
    my @freq_levels = ( @freqs, [ 'unknown', -1, -1 ], );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'pigccv_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{distance AVG_pi STD_pi AVG_gc STD_gc
                AVG_cv STD_cv COUNT};
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
                    isw_distance,
                    AVG(isw_pi) AVG_D,
                    STD(isw_pi) STD_D,
                    AVG(isw_average_gc) AVG_gc,
                    STD(isw_average_gc) STD_gc,
                    AVG(isw_cv) AVG_cv,
                    STD(isw_cv) STD_cv,
                    COUNT(*) COUNT
                FROM    isw s, indel i
                WHERE s.isw_indel_id = i.indel_id
                AND i.indel_freq >= ?
                AND i.indel_freq <= ?
                GROUP BY isw_distance
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

foreach my $n (@tasks) {
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
    if ( $n == 52 ) { &$di_dn_ttest;          next; }
    if ( $n == 53 ) { &$frequency_distance2;  next; }
    if ( $n == 54 ) { &$frequency_distance3;  next; }
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

