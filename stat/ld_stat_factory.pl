#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::SQL;
use AlignDB::SQL::Library;

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
my $run               = $Config->{stat}{run};
my $sum_threshold     = $Config->{stat}{sum_threshold};
my $combine_threshold = $Config->{stat}{combine_threshold};
my $outfile           = "";

my $max_freq;    # count freq one by one to $max_freq

my $help = 0;
my $man  = 0;

GetOptions(
    'help|?'                 => \$help,
    'man'                    => \$man,
    's|server=s'             => \$server,
    'P|port=s'               => \$port,
    'd|db=s'                 => \$db,
    'u|username=s'           => \$username,
    'p|password=s'           => \$password,
    'o|output=s'             => \$outfile,
    'max|max_freq=s'         => \$max_freq,
    'r|run=s'                => \$run,
    't|st|threshold=i'       => \$sum_threshold,
    'ct|combine_threshold=i' => \$combine_threshold,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

my $sql_file = AlignDB::SQL::Library->new( lib => "$FindBin::Bin/sql.lib" );

# auto detect sum threshold
if ( $sum_threshold == 0 ) {
    ( $sum_threshold, undef ) = $write_obj->calc_threshold;
}

# auto detect combine threshold
if ( $combine_threshold == 0 ) {
    ( undef, $combine_threshold ) = $write_obj->calc_threshold;
}

#----------------------------#
# count freq
#----------------------------#
my $all_freq = $write_obj->get_freq;

if ( $all_freq < 2 ) {
    die "all_freq is $all_freq, are you sure this is an AlignDB DB?\n";
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
# worksheet -- indel_ld
#----------------------------------------------------------#
my $indel_ld = sub {

    {
        my $sheet_name = 'd1_indel_ld';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    {
        my $sheet_name = 'd2_indel_ld';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $indel_ld_insdel = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet_d1 = sub {
        my ($level) = @_;
        my $sheet_name = 'd1_indel_ld_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );
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

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_indel_ld_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );
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
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

    {
        my $sheet_name = 'd2_snps_ld';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            my %option = (
                sql_query => $thaw_sql->as_sql,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    }
};

my $snps_ld_insdel = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );

    my $write_sheet_d1 = sub {
        my ($level) = @_;
        my $sheet_name = 'd1_snps_ld_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );
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

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_snps_ld_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );
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
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
            $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );
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

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_snps_ld_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
            $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );
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

    for (@freq_levels) {
        $write_sheet_d1->($_);
    }

    for (@freq_levels) {
        $write_sheet_d2->($_);
    }
};

my $snps_ld_insdel_freq = sub {
    my @type_levels = ( [ 'ins', 'I' ], [ 'del', 'D' ], );
    my @freq_levels = @freqs;

    my $write_sheet_d1 = sub {
        my ( $type, $freq ) = @_;
        my $sheet_name = 'd1_snps_ld_' . $type->[0] . '_' . $freq->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );
            $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
            $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $type->[1], $freq->[1], $freq->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ( $type, $freq ) = @_;
        my $sheet_name = 'd2_snps_ld_' . $type->[0] . '_' . $freq->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );

        {    # write header
            my @headers = $thaw_sql->as_header;
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
            $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );
            $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
            $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );
            my %option = (
                sql_query  => $thaw_sql->as_sql,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $type->[1], $freq->[1], $freq->[2] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    for my $type (@type_levels) {
        for my $freq (@freq_levels) {
            $write_sheet_d1->( $type, $freq );
        }
    }

    for my $type (@type_levels) {
        for my $freq (@freq_levels) {
            $write_sheet_d2->( $type, $freq );
        }
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
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_average_gc `gc`,
                           s.segment_gc_CV `cv`,
                           w.window_coding `coding`,
                           s.segment_r2_s `r2`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY gc DESC, pi, indel
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],

            );
            $write_obj->excute_sql( \%option );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                SELECT t_id, length
                FROM tmp_group
            };
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $sum_threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
        }

        {    # write header
            my @headers
                = qw{AVG_gc AVG_pi AVG_Indel/100bp AVG_CV AVG_coding AVG_r2 AVG_length COUNT SUM_length};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@combined_segment,
            );

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
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
        my ( $sheet_row, $sheet_col );

        {    # create temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {
            my $sql_query = q{
                # create temporary table
                CREATE TABLE tmp_group (t_id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY (t_id))
                    ENGINE=MyISAM
                    SELECT w.window_pi `pi`,
                           w.window_indel `indel`,
                           w.window_average_gc `gc`,
                           s.segment_gc_CV `cv`,
                           w.window_coding `coding`,
                           s.segment_r2_s `r2`,
                           w.window_length `length`
                    FROM segment s, window w
                    WHERE s.window_id = w.window_id
                    AND s.segment_type = ?
                    ORDER BY cv DESC, pi, indel
            };
            my %option = (
                sql_query  => $sql_query,
                bind_value => [$segment_type],

            );
            $write_obj->excute_sql( \%option );
        }

        # make group
        my @combined_segment;
        {
            my $sql_query = q{
                SELECT t_id, length
                FROM tmp_group
            };
            my $merge_last = 1;
            my %option     = (
                sql_query  => $sql_query,
                threshold  => $sum_threshold,
                merge_last => $merge_last,
            );
            @combined_segment = @{ $write_obj->make_combine( \%option ) };
        }

        {    # write header
            my @headers
                = qw{AVG_CV AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_r2 AVG_length COUNT SUM_length Range_gc};
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
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
            my %option = (
                sql_query => $sql_query,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                group     => \@combined_segment,
            );

            ($sheet_row) = $write_obj->write_content_group( $sheet, \%option );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
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

    #if ( $n == 14 ) { &$snps_ld_insdel_freq; next; }

    if ( $n == 21 ) { &$segment_gc_indel; next; }
    if ( $n == 22 ) { &$segment_cv_indel; next; }
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

    ld_stat_factory.pl - Generate statistical Excel files from malignDB

=head1 SYNOPSIS

    ld_stat_factory.pl [options]
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

