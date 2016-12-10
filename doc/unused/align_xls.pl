#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use Excel::Writer::XLSX;
use List::Util;

use AlignDB::Stopwatch;

use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

align_xls.pl - Generate a colorful excel file for one alignment in a three-way alignDB

=head1 SYNOPSIS

    perl align_xls.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --wrap          INT     wrap length, default is [50] 
        --spacing       INT     wrapped line spacing, default is [0]
        --align_id      INT     align_id
        --outfile       STR     output file name

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'wrap=i'       => \( my $wrap     = 50 ),
    'spacing=i'    => \( my $spacing  = 0 ),
    'align_id=s'   => \( my $align_id = 1 ),
    'output=s'     => \my $outfile,
) or Getopt::Long::HelpMessage(1);

$outfile ||= "$db-align-$align_id.xlsx";

#----------------------------------------------------------#
# Init objects and SQL queries
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Writing xls...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

#@type DBI
my $dbh = $obj->dbh;

# get target, query and reference names via AlignDB methods
my ( $target_name, $query_name, $ref_name ) = $obj->get_names($align_id);
my $max_name_length = List::Util::max( length $target_name, length $query_name, length $ref_name );

if ( !$ref_name or $ref_name eq 'NULL' ) {
    die "$db is not a three-way alignDB\n";
}

# get snp info
my $snp_query = q{
    SELECT s.snp_pos,
           s.snp_occured,
           s.target_base,
           s.query_base,
           s.ref_base
    FROM align a
    INNER JOIN snp s ON a.align_id = s.align_id
    WHERE a.align_id = ?
    ORDER BY s.snp_pos
};
my $snp_query_sth = $dbh->prepare($snp_query);

# select all isws in this alignment
my $indel_query = q{
    SELECT i.indel_start,
           i.indel_end,
           i.indel_length,
           i.indel_occured,
           i.indel_type
    FROM align a
    INNER JOIN indel i ON a.align_id = i.align_id
    WHERE a.align_id = ?
    ORDER BY i.indel_start
};
my $indel_query_sth = $dbh->prepare($indel_query);

# Create workbook and worksheet objects
my $workbook;
unless ( $workbook = Excel::Writer::XLSX->new($outfile) ) {
    die "Cannot create Excel file $outfile\n";
}

#@type Excel::Writer::XLSX::Worksheet
my $sheet = $workbook->add_worksheet;

#----------------------------------------------------------#
# Get data
#----------------------------------------------------------#
print "Get variation data...\n";

# store all variations, including indels and snps
my %variations;

# snp
$snp_query_sth->execute($align_id);
while ( my $hash_ref = $snp_query_sth->fetchrow_hashref ) {
    my $start_pos = $hash_ref->{snp_pos};
    $hash_ref->{var_type} = 'snp';
    $variations{$start_pos} = $hash_ref;
}

# indel
$indel_query_sth->execute($align_id);
while ( my $hash_ref = $indel_query_sth->fetchrow_hashref ) {
    my $start_pos = $hash_ref->{indel_start};
    $hash_ref->{var_type} = 'indel';
    $variations{$start_pos} = $hash_ref;
}

#DumpFile( "var.yaml", \%variations );

#----------------------------------------------------------#
# Excel format
#----------------------------------------------------------#
# species name format
my $name_format = $workbook->add_format(
    font => 'Courier New',
    size => 10,
);

# variation position format
my $pos_format = $workbook->add_format(
    font     => 'Courier New',
    size     => 8,
    align    => 'center',
    valign   => 'vcenter',
    rotation => 90,
);

# snp base format
my $snp_fg_of = {
    A   => { color => 'green', },
    C   => { color => 'blue', },
    G   => { color => 'pink', },
    T   => { color => 'red', },
    N   => { color => 'black' },
    '-' => { color => 'black' },
};

my $snp_bg_of = {
    T => { bg_color => 43, },    # lightyellow
    Q => { bg_color => 42, },    # lightgreen
    N => {},
};

my $snp_format = {};
for my $fg ( keys %{$snp_fg_of} ) {
    for my $bg ( keys %{$snp_bg_of} ) {
        $snp_format->{"$fg$bg"} = $workbook->add_format(
            font   => 'Courier New',
            size   => 10,
            bold   => 1,
            align  => 'center',
            valign => 'vcenter',
            %{ $snp_fg_of->{$fg} },
            %{ $snp_bg_of->{$bg} },
        );
    }
}

# indel format
my $indel_bg_of = {
    T => { bg_color => 'yellow', },
    Q => { bg_color => 11, },         # brightgreen, lime
    N => { bg_color => 'silver', },
    C => { bg_color => 'gray', },
};

my $indel_format = {};
my $merge_format = {};
for my $bg ( keys %{$indel_bg_of} ) {
    $indel_format->{$bg} = $workbook->add_format(
        font   => 'Courier New',
        size   => 10,
        bold   => 1,
        align  => 'center',
        valign => 'vcenter',
        %{ $indel_bg_of->{$bg} },
    );
    $merge_format->{$bg} = $workbook->add_format(
        font   => 'Courier New',
        size   => 10,
        bold   => 1,
        align  => 'center',
        valign => 'vcenter',
        %{ $indel_bg_of->{$bg} },
    );
}

#----------------------------------------------------------#
# write execel
#----------------------------------------------------------#
print "Write excel...\n";

my $col_cursor     = 1;
my $section        = 1;
my $section_height = 4 + $spacing;

for my $pos ( sort { $a <=> $b } keys %variations ) {
    my $var = $variations{$pos};
    my $pos_row = $section_height * ( $section - 1 );

    # write snp
    if ( $var->{var_type} eq 'snp' ) {
        my $snp_pos     = $var->{snp_pos};
        my $snp_occured = $var->{snp_occured};
        my $ref_base    = $var->{ref_base};
        my $target_base = $var->{target_base};
        my $query_base  = $var->{query_base};

        # write position
        $sheet->write( $pos_row, $col_cursor, $snp_pos, $pos_format );

        # write reference
        my $ref_occ = "$ref_base$snp_occured";
        $sheet->write( $pos_row + 1, $col_cursor, $ref_base, $snp_format->{$ref_occ} );

        # write target
        my $target_occ = "$target_base$snp_occured";
        $sheet->write( $pos_row + 2, $col_cursor, $target_base, $snp_format->{$target_occ} );

        # write query
        my $query_occ = "$query_base$snp_occured";
        $sheet->write( $pos_row + 3, $col_cursor, $query_base, $snp_format->{$query_occ} );

        $col_cursor++;
    }

    # write indel
    if ( $var->{var_type} eq 'indel' ) {
        my $indel_start   = $var->{indel_start};
        my $indel_end     = $var->{indel_end};
        my $indel_length  = $var->{indel_length};
        my $indel_occured = $var->{indel_occured};
        my $indel_type    = $var->{indel_type};

        # how many column does this indel take up
        my $col_takeup = List::Util::min( $indel_length, 3 );

        # if exceed the wrap limit, start a new section
        if ( $col_cursor + $col_takeup > $wrap ) {
            $col_cursor = 1;
            $section++;
            $pos_row = $section_height * ( $section - 1 );
        }

        #print Dump($var);

        # offset from the pos_row
        my $indel_offset
            = $indel_occured eq "N" ? 1
            : $indel_occured eq "C" ? 1
            : $indel_occured eq "T" ? 2
            : $indel_occured eq "Q" ? 3
            :                         undef;

        my $indel_string = "$indel_type$indel_length";

        if ( $indel_length == 1 ) {

            # write position
            $sheet->write( $pos_row, $col_cursor, $indel_start, $pos_format );

            # write in indel occured lineage
            $sheet->write( $pos_row + $indel_offset,
                $col_cursor, $indel_string, $indel_format->{$indel_occured} );

            $col_cursor++;
        }
        elsif ( $indel_length == 2 ) {

            # write indel_start position
            my $start_col = $col_cursor;
            $sheet->write( $pos_row, $start_col, $indel_start, $pos_format );
            $col_cursor++;

            # write indel_end position
            my $end_col = $col_cursor;
            $sheet->write( $pos_row, $end_col, $indel_end, $pos_format );
            $col_cursor++;

            # merge two indel position
            $sheet->merge_range(
                $pos_row + $indel_offset,
                $start_col, $pos_row + $indel_offset,
                $end_col, $indel_string, $merge_format->{$indel_occured},
            );
        }
        else {

            # write indel_start position
            my $start_col = $col_cursor;
            $sheet->write( $pos_row, $start_col, $indel_start, $pos_format );
            $col_cursor++;

            # write middle sign
            my $middle_col = $col_cursor;
            $sheet->write( $pos_row, $middle_col, '|', $pos_format );
            $col_cursor++;

            # write indel_end position
            my $end_col = $col_cursor;
            $sheet->write( $pos_row, $end_col, $indel_end, $pos_format );
            $col_cursor++;

            # merge two indel position
            $sheet->merge_range(
                $pos_row + $indel_offset,
                $start_col, $pos_row + $indel_offset,
                $end_col, $indel_string, $merge_format->{$indel_occured},
            );
        }

    }

    if ( $col_cursor > $wrap ) {
        $col_cursor = 1;
        $section++;
    }
}

# write names
for ( 1 .. $section ) {
    my $pos_row = $section_height * ( $_ - 1 );

    $sheet->write( $pos_row + 1, 0, $ref_name,    $name_format );
    $sheet->write( $pos_row + 2, 0, $target_name, $name_format );
    $sheet->write( $pos_row + 3, 0, $query_name,  $name_format );
}

# format column
$sheet->set_column( 0, 0,         $max_name_length + 2 );
$sheet->set_column( 1, $wrap + 3, 1.6 );

$workbook->close;

$stopwatch->end_message;
exit;

__END__
