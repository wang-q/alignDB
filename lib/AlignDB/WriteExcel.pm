package AlignDB::WriteExcel;
use Moose;

use FindBin;
use lib "$FindBin::Bin/../";
extends qw(AlignDB);
use AlignDB::SQL;

use Excel::Writer::XLSX;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw( all any );
use YAML qw(Dump Load DumpFile LoadFile);

use Statistics::PointEstimation;
use Statistics::TTest;
use Statistics::Descriptive;

has 'outfile'  => ( is => 'ro', isa => 'Str' );     # output file, autogenerable
has 'workbook' => ( is => 'rw', isa => 'Object' );  # excel workbook object
has 'format'   => ( is => 'ro', isa => 'HashRef' ); # excel formats
has 'columns'  => ( is => 'ro', isa => 'HashRef' ); # excel column names

sub BUILD {
    my $self = shift;

    # set outfile
    unless ( $self->outfile ) {
        $self->mysql =~ /^(.+):/;
        $self->{outfile} = "$1.auto.xlsx";
    }

    # Create $workbook object
    my $workbook;
    unless ( $workbook = Excel::Writer::XLSX->new( $self->outfile ) ) {
        warn "Cannot create Excel file.\n";
        return;
    }
    $self->{workbook} = $workbook;

    # set $workbook format
    my %font = (
        font => 'Arial',
        size => 10,
    );
    my $format = {
        HEADER => $workbook->add_format(
            align    => 'center',
            bg_color => 42,
            bold     => 1,
            bottom   => 2,
            %font,
        ),
        HIGHLIGHT => $workbook->add_format( color => 'blue',  %font, ),
        NORMAL    => $workbook->add_format( color => 'black', %font, ),
        NAME => $workbook->add_format( bold => 1, color => 57, %font, ),
        TOTAL => $workbook->add_format( bold => 1, top => 2, %font, ),
        DATE  => $workbook->add_format(
            bg_color   => 42,
            bold       => 1,
            num_format => 'yy-m-d hh:mm',
            %font,
        ),
    };
    $self->{format} = $format;

    # set $workbook column names
    my $columns = [
        'A:A',  'B:B',  'C:C',  'D:D',  'E:E',  'F:F',  'G:G',  'H:H',
        'I:I',  'J:J',  'K:K',  'L:L',  'M:M',  'N:N',  'O:O',  'P:P',
        'Q:Q',  'R:R',  'S:S',  'T:T',  'U:U',  'V:V',  'W:W',  'X:X',
        'Y:Y',  'Z:Z',  'AA:A', 'BB:B', 'CC:C', 'DD:D', 'EE:E', 'FF:F',
        'GG:G', 'HH:H', 'II:I', 'JJ:J', 'KK:K', 'LL:L', 'MM:M', 'NN:N',
        'OO:O', 'PP:P', 'QQ:Q', 'RR:R', 'SS:S', 'TT:T', 'UU:U', 'VV:V',
        'WW:W', 'XX:X', 'YY:Y', 'ZZ:Z'
    ];
    $self->{columns} = $columns;

    return;
}

sub write_header_direct {
    my ( $self, $sheet_name, $option ) = @_;

    # init
    my $workbook = $self->workbook;
    my $sheet    = $workbook->add_worksheet($sheet_name);
    my $fmt      = $self->format;
    my @cols     = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};
    my $header    = $option->{header};

    # query name
    my $query_name = $option->{query_name};

    # create table header
    my @cols_name = @$header;
    for ( my $i = 0; $i < $sheet_col; $i++ ) {
        $sheet->set_column( $cols[$i], 12 );
        $sheet->write( $sheet_row, $i, $query_name, $fmt->{HEADER} );
    }
    for ( my $i = 0; $i <= $#cols_name; $i++ ) {
        $sheet->set_column( $cols[ $i + $sheet_col ],
            max( length( $cols_name[$i] ) + 2, 9 ) );
        $sheet->write( $sheet_row, $i + $sheet_col,
            $cols_name[$i], $fmt->{HEADER} );
    }
    $sheet_row++;
    $sheet->freeze_panes( 1, 0 );    # freeze table

    return ( $sheet, $sheet_row );
}

sub write_header_sql {
    my ( $self, $sheet_name, $option ) = @_;

    # init
    my $dbh      = $self->dbh;
    my $workbook = $self->workbook;
    my $sheet    = $workbook->add_worksheet($sheet_name);
    my $fmt      = $self->format;
    my @cols     = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    # query name
    my $query_name = $option->{query_name};

    # init DBI query
    my $sql_query = $option->{sql_query};
    my $sth       = $dbh->prepare($sql_query);
    $sth->execute();

    # create table header
    my @cols_name = @{ $sth->{'NAME'} };
    for ( my $i = 0; $i < $sheet_col; $i++ ) {
        $sheet->set_column( $cols[$i], 16 );
        $sheet->write( $sheet_row, $i, $query_name, $fmt->{HEADER} );
    }
    for ( my $i = 0; $i < scalar @cols_name; $i++ ) {
        $sheet->set_column( $cols[ $i + $sheet_col ],
            max( length( $cols_name[$i] ) + 2, 9 ) );
        $sheet->write( $sheet_row, $i + $sheet_col,
            $cols_name[$i], $fmt->{HEADER} );
    }
    $sheet_row++;
    $sheet->freeze_panes( 1, 0 );    # freeze table

    return ( $sheet, $sheet_row );
}

sub write_row_direct {
    my ( $self, $sheet, $option ) = @_;

    # init
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    # query name
    my $query_name = $option->{query_name};
    if ( defined $query_name ) {
        $sheet->write( $sheet_row, $sheet_col - 1, $query_name, $fmt->{NAME} );
    }

    # array_ref
    my $row = $option->{row};

    # content format
    my $content_format = $option->{content_format};
    unless ( defined $content_format ) {
        $content_format = 'NORMAL';
    }

    # reverse write
    my $write_step = $option->{write_step};
    unless ( defined $write_step ) {
        $write_step = 1;
    }

    # bind value
    my $append_column = $option->{append_column};
    unless ( defined $append_column ) {
        $append_column = [];
    }

    # insert table columns
    for ( my $i = 0; $i < scalar @$row; $i++ ) {
        $sheet->write(
            $sheet_row, $i + $sheet_col,
            $row->[$i], $fmt->{$content_format}
        );
    }
    $sheet_row += $write_step;

    return ($sheet_row);
}

sub write_content_direct {
    my ( $self, $sheet, $option ) = @_;

    # init
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    # query name
    my $query_name = $option->{query_name};
    if ( defined $query_name ) {
        $sheet->write( $sheet_row, $sheet_col - 1, $query_name, $fmt->{NAME} );
    }

    # content format
    my $content_format = $option->{content_format};
    unless ( defined $content_format ) {
        $content_format = 'NORMAL';
    }

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    # reverse write
    my $write_step = $option->{write_step};
    unless ( defined $write_step ) {
        $write_step = 1;
    }

    # append column
    my $append_column = $option->{append_column};
    unless ( defined $append_column ) {
        $append_column = [];
    }

    # init DBI query
    my $sql_query = $option->{sql_query};
    my $sth       = $dbh->prepare($sql_query);
    $sth->execute(@$bind_value);

    # insert table columns
    while ( my @row = $sth->fetchrow_array ) {
        for ( my $i = 0; $i < scalar @row; $i++ ) {
            $sheet->write(
                $sheet_row, $i + $sheet_col,
                $row[$i],   $fmt->{$content_format}
            );
        }
        if ( scalar @$append_column ) {
            my $appand_row = shift @$append_column;
            for ( my $i = 0; $i < scalar @$appand_row; $i++ ) {
                $sheet->write(
                    $sheet_row,        $i + $sheet_col + scalar(@row),
                    $appand_row->[$i], $fmt->{$content_format}
                );
            }
        }
        $sheet_row += $write_step;
    }

    return ($sheet_row);
}

sub write_content_combine {
    my ( $self, $sheet, $option ) = @_;

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};
    my $sql_query = $option->{sql_query};

    my @combined = @{ $option->{combined} };

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    foreach (@combined) {
        my @range      = @$_;
        my $in_list    = '(' . join( ',', @range ) . ')';
        my $sql_query2 = $sql_query . $in_list;
        my %option     = (
            sql_query  => $sql_query2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            bind_value => $bind_value,
        );
        ($sheet_row) = $self->write_content_direct( $sheet, \%option );
    }
    return ($sheet_row);
}

sub write_content_combine_obj {
    my ( $self, $sheet, $option ) = @_;

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};
    my $sql_obj   = $option->{sql_obj};     # an AlignDB::SQL object

    my @combined    = @{ $option->{combined} };
    my $combine_col = $option->{combine_col};

    # this sub does not accept binding values

    foreach (@combined) {
        my @range     = @$_;
        my $sql_clone = $sql_obj->copy;
        $sql_clone->add_where( $combine_col => \@range );
        my %option = (
            sql_query  => $sql_clone->as_sql,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            bind_value => $sql_clone->bind,
        );
        ($sheet_row) = $self->write_content_direct( $sheet, \%option );
    }
    return ($sheet_row);
}

sub write_content_group {
    my ( $self, $sheet, $option ) = @_;

    # init table cursor
    my $sheet_row     = $option->{sheet_row};
    my $sheet_col     = $option->{sheet_col};
    my $sql_query     = $option->{sql_query};
    my $append_column = $option->{append_column};

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    my @group = @{ $option->{group} };

    foreach (@group) {
        my @range      = @$_;
        my $in_list    = '(' . join( ',', @range ) . ')';
        my $sql_query2 = $sql_query . $in_list;
        my $group_name;
        if ( scalar @range > 1 ) {
            $group_name = $range[0] . "--" . $range[-1];
        }
        else {
            $group_name = $range[0];
        }
        my %option = (
            sql_query     => $sql_query2,
            sheet_row     => $sheet_row,
            sheet_col     => $sheet_col,
            query_name    => $group_name,
            append_column => $append_column,
            bind_value    => $bind_value,
        );
        ($sheet_row) = $self->write_content_direct( $sheet, \%option );
    }
    return ($sheet_row);
}

sub write_content_group_obj {
    my ( $self, $sheet, $option ) = @_;

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};
    my $sql_obj   = $option->{sql_obj};     # an AlignDB::SQL object

    my @group     = @{ $option->{group} };
    my $group_col = $option->{group_col};

    # this sub does not accept binding values and appending columns

    foreach (@group) {
        my @range     = @$_;
        my $sql_clone = $sql_obj->copy;
        $sql_clone->add_where( $group_col => \@range );
        my $group_name;
        if ( scalar @range > 1 ) {
            $group_name = $range[0] . "--" . $range[-1];
        }
        else {
            $group_name = $range[0];
        }
        my %option = (
            sql_query  => $sql_clone->as_sql,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $group_name,
            bind_value => $sql_clone->bind,
        );
        ($sheet_row) = $self->write_content_direct( $sheet, \%option );
    }
    return ($sheet_row);
}

sub write_content_indel {
    my ( $self, $sheet, $option ) = @_;

    # init objects
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    my $sql_query_1 = $option->{sql_query_1};
    my $sql_query_2 = $option->{sql_query_2};
    my @group       = @{ $option->{group} };

    foreach (@group) {
        my @range = @$_;
        my $group_name;
        if ( scalar @range > 1 ) {
            $group_name = join "-", @range;
        }
        else {
            $group_name = $range[0];
        }
        $sheet_row++;    # add a blank line
        my %option = (
            sql_query  => $sql_query_1,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            query_name => $group_name,
            bind_value => \@range,
        );
        ($sheet_row) = $self->write_content_direct( $sheet, \%option );

        $sheet_row++;    # add a blank line
        %option = (
            sql_query  => $sql_query_2,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
            bind_value => \@range,
        );

        #print Dump \%option;
        ($sheet_row) = $self->write_content_direct( $sheet, \%option );
    }
}

sub write_content_dd {
    my ( $self, $sheet, $option ) = @_;

    # init objects
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    my $sql_query = $option->{sql_query};
    my @group     = @{ $option->{group} };

    my $group_rows_range = int( $group[-1]->[0] + 1 ) + 2 + 3;

    foreach (@group) {
        my @range = @$_;
        my $group_name;
        if ( scalar @range > 1 ) {
            $group_name = join "-", @range;
        }
        else {
            $group_name = $range[0];
        }
        my $group_row = $sheet_row;
        $group_row++;    # add a blank line
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $group_row,
            sheet_col  => $sheet_col,
            query_name => $group_name,
            bind_value => [ 'L', @range, $range[0] ],
        );
        ($group_row) = $self->write_content_direct( $sheet, \%option );

        $group_row = $sheet_row + $group_rows_range;
        $group_row--;
        %option = (
            sql_query  => $sql_query,
            sheet_row  => $group_row,
            sheet_col  => $sheet_col,
            bind_value => [ 'R', @range, $range[0] ],
            write_step => -1,
        );
        ($group_row) = $self->write_content_direct( $sheet, \%option );

        $sheet_row += $group_rows_range;
    }
}

sub write_content_dd_gc {
    my ( $self, $sheet, $option ) = @_;

    # init objects
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    my $sql_query = $option->{sql_query};
    my @group     = @{ $option->{group} };

    my $group_rows_range = int( $group[-1]->[0] + 1 ) * 2 + 3;

    foreach (@group) {
        my @range = @$_;
        my $group_name;
        if ( scalar @range > 1 ) {
            $group_name = join "-", @range;
        }
        else {
            $group_name = $range[0];
        }
        my $group_row = $sheet_row;
        $group_row++;    # add a blank line
        my %option = (
            sql_query  => $sql_query,
            sheet_row  => $group_row,
            sheet_col  => $sheet_col,
            query_name => $group_name,
            bind_value => [ 'L', @range, $range[0] ],
        );
        ($group_row) = $self->write_content_direct( $sheet, \%option );

        $group_row = $sheet_row + $group_rows_range;
        $group_row--;
        %option = (
            sql_query  => $sql_query,
            sheet_row  => $group_row,
            sheet_col  => $sheet_col,
            bind_value => [ 'R', @range, $range[0] ],
            write_step => -1,
        );
        ($group_row) = $self->write_content_direct( $sheet, \%option );

        $sheet_row += $group_rows_range;
    }
}

sub write_content_highlight {
    my ( $self, $sheet, $option ) = @_;

    # init
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    # init DBI query
    my $sql_query = $option->{sql_query};
    my $sth       = $dbh->prepare($sql_query);
    $sth->execute(@$bind_value);

    # insert table columns
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
            $sheet->write( $sheet_row, $i + $sheet_col,
                $row[$i], $fmt->{$style} );
        }
        $sheet_row++;
    }
}

sub write_content_snp {
    my ( $self, $sheet, $option ) = @_;

    # init objects
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    my $sql_query1 = $option->{sql_query1};
    my $sql_query2 = $option->{sql_query2};
    my @base_pair  = @{ $option->{base_pair} };

    for ( my $i = 0; $i < scalar @base_pair; $i++ ) {
        my $base_pair_1 = $base_pair[$i];
        $base_pair_1 =~ s/\W//g;
        my $base_pair_2 = reverse $base_pair_1;

        my $sub_sheet_row = $sheet_row;

        my $sth = $dbh->prepare($sql_query1);
        $sth->execute;

        my %total_snp = ();
        while ( my @row = $sth->fetchrow_array ) {
            $total_snp{ $row[0] } = $row[1];
        }

        $sth = $dbh->prepare($sql_query2);
        $sth->execute( $base_pair_1, $base_pair_2 );

        # insert table columns
        while ( my @row = $sth->fetchrow_array ) {
            $sheet->write( $sub_sheet_row, 0 + $sheet_col,
                $row[0], $fmt->{NORMAL} );
            $sheet->write(
                $sub_sheet_row,
                $i + 1 + $sheet_col,
                $row[1] / $total_snp{ $row[0] },
                $fmt->{NORMAL}
            );
            $sub_sheet_row++;
        }
        $sth->finish;
    }
}

sub write_content_snp2 {
    my ( $self, $sheet, $option ) = @_;

    # init objects
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    my $sql_query1 = $option->{sql_query1};
    my $sql_query2 = $option->{sql_query2};
    my @base_pair  = @{ $option->{base_pair} };

    for ( my $i = 0; $i < scalar @base_pair; $i++ ) {
        my $base_pair = $base_pair[$i];
        $base_pair =~ s/\W//g;

        my $sub_sheet_row = $sheet_row;

        my $sth = $dbh->prepare($sql_query1);
        $sth->execute;

        my %total_snp = ();
        while ( my @row = $sth->fetchrow_array ) {
            $total_snp{ $row[0] } = $row[1];
        }

        $sth = $dbh->prepare($sql_query2);
        $sth->execute($base_pair);

        # insert table columns
        while ( my @row = $sth->fetchrow_array ) {
            $sheet->write( $sub_sheet_row, 0 + $sheet_col,
                $row[0], $fmt->{NORMAL} );
            $sheet->write(
                $sub_sheet_row,
                $i + 1 + $sheet_col,
                $row[1] / $total_snp{ $row[0] },
                $fmt->{NORMAL}
            );
            $sub_sheet_row++;
        }
        $sth->finish;
    }
}

sub write_content_snp3 {
    my ( $self, $sheet, $option ) = @_;

    # init objects
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    my $sql_query1 = $option->{sql_query1};
    my $sql_query2 = $option->{sql_query2};
    my @base_pair  = @{ $option->{base_pair} };

    for ( my $i = 0; $i < scalar @base_pair; $i++ ) {
        my @bases           = @{ $base_pair[$i] };
        my $base_group_name = shift @bases;
        foreach (@bases) {
            $_ = "\"" . $_ . "\"";
        }
        my $in_list = join ",", @bases;

        my $sub_sheet_row = $sheet_row;

        my $sth = $dbh->prepare($sql_query1);
        $sth->execute;

        my %total_snp = ();
        while ( my @row = $sth->fetchrow_array ) {
            $total_snp{ $row[0] } = $row[1];
        }

        my $base_query = $sql_query2;
        $base_query =~ s/in_list/$in_list/g;
        $sth = $dbh->prepare($base_query);
        $sth->execute();

        # insert table columns
        while ( my @row = $sth->fetchrow_array ) {
            $sheet->write( $sub_sheet_row, 0 + $sheet_col,
                $row[0], $fmt->{NORMAL} );
            $sheet->write(
                $sub_sheet_row,
                $i + 1 + $sheet_col,
                $row[1] / $total_snp{ $row[0] },
                $fmt->{NORMAL}
            );
            $sub_sheet_row++;
        }
        $sth->finish;
    }
}

sub write_content_column {
    my ( $self, $sheet, $option ) = @_;

    # init objects
    my $dbh  = $self->dbh;
    my $fmt  = $self->format;
    my @cols = @{ $self->columns };

    # init table cursor
    my $sheet_row = $option->{sheet_row};
    my $sheet_col = $option->{sheet_col};

    my $sql_query  = $option->{sql_query};
    my @conditions = @{ $option->{conditions} };

    for ( my $i = 0; $i < scalar @conditions; $i++ ) {
        my @bind_values = @{ $conditions[$i] };

        my $sub_sheet_row = $sheet_row;

        my $sth = $dbh->prepare($sql_query);
        $sth->execute(@bind_values);

        # insert table columns
        while ( my @row = $sth->fetchrow_array ) {
            $sheet->write( $sub_sheet_row, 0 + $sheet_col,
                $row[0], $fmt->{NORMAL} );
            $sheet->write(
                $sub_sheet_row, $i + 1 + $sheet_col,
                $row[1],        $fmt->{NORMAL}
            );
            $sub_sheet_row++;
        }
        $sth->finish;
    }
}

sub make_combine {
    my ( $self, $option ) = @_;

    # init objects
    my $dbh = $self->dbh;

    # init parameters
    my $sql_query  = $option->{sql_query};
    my $threshold  = $option->{threshold};
    my $standalone = $option->{standalone};

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    # merge_last
    my $merge_last = $option->{merge_last};
    unless ( defined $merge_last ) {
        $merge_last = 0;
    }

    # init DBI query
    my $sth = $dbh->prepare($sql_query);
    $sth->execute(@$bind_value);

    my @row_count = ();
    while ( my @row = $sth->fetchrow_array ) {
        push @row_count, \@row;
    }

    my @combined;    # return these
    my @temp_combined = ();
    my $temp_count    = 0;
    foreach my $row_ref (@row_count) {
        if ( any { $_ eq $row_ref->[0] } @{$standalone} ) {
            push @combined, [ $row_ref->[0] ];
        }
        elsif ( $temp_count < $threshold ) {
            push @temp_combined, $row_ref->[0];
            $temp_count += $row_ref->[1];

            if ( $temp_count < $threshold ) {
                next;
            }
            else {
                push @combined, [@temp_combined];
                @temp_combined = ();
                $temp_count    = 0;
            }
        }
        else {
            warn "Errors occured in calculating combined distance.\n";
        }
    }

    # Write the last weighted row which COUNT might
    #   be smaller than $threshold
    if ( $temp_count > 0 ) {
        if ($merge_last) {
            if ( @combined == 0 ) {
                @combined = ( [] );
            }
            push @{ $combined[-1] }, @temp_combined;
        }
        else {
            push @combined, [@temp_combined];
        }
    }

    return \@combined;
}

sub make_last_portion {
    my ( $self, $option ) = @_;

    # init objects
    my $dbh = $self->dbh;

    # init parameters
    my $sql_query = $option->{sql_query};
    my $portion   = $option->{portion};

    # init DBI query
    my $sth = $dbh->prepare($sql_query);
    $sth->execute;

    my @row_count = ();
    while ( my @row = $sth->fetchrow_array ) {
        push @row_count, \@row;
    }

    my @last_portion;    # return @last_portion
    my $all_length = 0;  # return $all_length
    foreach (@row_count) {
        $all_length += $_->[2];
    }
    my @rev_row_count = reverse @row_count;
    my $temp_length   = 0;
    foreach (@rev_row_count) {
        push @last_portion, $_->[0];
        $temp_length += $_->[2];
        if ( $temp_length >= $all_length * $portion ) {
            last;
        }
    }

    return ( $all_length, \@last_portion );
}

sub excute_sql {
    my ( $self, $option ) = @_;

    # init
    my $dbh = $self->dbh;

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    # init DBI query
    my $sql_query = $option->{sql_query};
    my $sth       = $dbh->prepare($sql_query);
    $sth->execute(@$bind_value);
}

sub check_column {
    my ( $self, $table, $column ) = @_;

    # init
    my $dbh = $self->dbh;

    {    # check table existing
        my @table_names = $dbh->tables( '', '', '' );

        # table names are quoted by ` (back-quotes) which is the
        #   quote_identifier
        my $table_name = "`$table`";
        unless ( any { $_ =~ /$table_name/i } @table_names ) {
            print " " x 4, "Table $table does not exist\n";
            return 0;
        }
    }

    {    # check column existing
        my $sql_query = qq{
            SHOW FIELDS
            FROM $table
            LIKE "$column"
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute();
        my ($field) = $sth->fetchrow_array;

        if ( not $field ) {
            print " " x 4, "Column $column does not exist\n";
            return 0;
        }
    }

    {    # check values in column
        my $sql_query = qq{
            SELECT COUNT($column)
            FROM $table
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute;
        my ($count) = $sth->fetchrow_array;

        if ( not $count ) {
            print " " x 4, "Column $column has no records\n";
        }

        return $count;
    }

}

sub column_ttest {
    my ( $self, $option ) = @_;

    # init objects
    my $dbh = $self->dbh;

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    # init DBI query
    my $sql_query = $option->{sql_query};
    my $sth       = $dbh->prepare($sql_query);
    $sth->execute(@$bind_value);

    my ( @sample1, @sample2 );

    # insert table columns
    while ( my @row = $sth->fetchrow_array ) {
        push @sample1, $row[0];
        push @sample2, $row[1];
    }

    my $ttest = Statistics::TTest->new;
    $ttest->set_significance(95);
    eval { $ttest->load_data( \@sample1, \@sample2 ); };

    warn $@ if $@;

    return $ttest;
}

sub quantile {
    my ( $self, $data, $part_number ) = @_;

    my $stat = Statistics::Descriptive::Full->new();

    $stat->add_data(@$data);

    my $min = $stat->min;
    my @quantiles;
    my $base = 100 / $part_number;
    for ( 1 .. $part_number - 1 ) {
        my $percentile = $stat->percentile( $_ * $base );
        push @quantiles, $percentile;
    }
    my $max = $stat->max;

    return [ $min, @quantiles, $max, ];
}

sub quantile_sql {
    my ( $self, $option, $part_number ) = @_;

    # init objects
    my $dbh = $self->dbh;

    # bind value
    my $bind_value = $option->{bind_value};
    unless ( defined $bind_value ) {
        $bind_value = [];
    }

    # init DBI query
    my $sql_query = $option->{sql_query};
    my $sth       = $dbh->prepare($sql_query);
    $sth->execute(@$bind_value);

    my @data;

    while ( my @row = $sth->fetchrow_array ) {
        push @data, $row[0];
    }

    return $self->quantile( \@data, $part_number );
}

sub calc_threshold {
    my $self = shift;

    my $dbh = $self->dbh;

    my ( $sum_threshold, $combine_threshold );

    my $sth = $dbh->prepare(q{ SELECT SUM(align_length) FROM align });
    $sth->execute;
    my ($total_length) = $sth->fetchrow_array;

    if ( $total_length <= 1_000_000 ) {
        $sum_threshold = int( $total_length / 10 );
    }
    elsif ( $total_length <= 10_000_000 ) {
        $sum_threshold = int( $total_length / 10 );
    }
    elsif ( $total_length <= 100_000_000 ) {
        $sum_threshold = int( $total_length / 20 );
    }
    elsif ( $total_length <= 1_000_000_000 ) {
        $sum_threshold = int( $total_length / 50 );
    }
    else {
        $sum_threshold = int( $total_length / 100 );
    }

    if ( $total_length <= 1_000_000 ) {
        $combine_threshold = 100;
    }
    elsif ( $total_length <= 10_000_000 ) {
        $combine_threshold = 500;
    }
    else {
        $combine_threshold = 1000;
    }

    return ( $sum_threshold, $combine_threshold );
}

# instance destructor
# invoked only as object method
sub DESTROY {
    my ($self) = shift;

    # close excel objects
    my $workbook = $self->workbook;
    $workbook->close if $workbook;

    # close dbh
    my $dbh = $self->dbh;
    $dbh->disconnect if $dbh;
}

1;

#
# A simple class to use DBI & Excel::Writer::XLSX to make excel files
#
# perl -e "print scalar localtime"
# Tue Mar  7 19:22:28 2006
#
# Version: 0.1
# Author: Wang Qiang
# For Gattaca lab.
