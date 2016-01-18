#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

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
    'piece=i'      => \( my $piece    = 0 ),
    'index'        => \( my $add_index_sheet, ),
) or HelpMessage(1);

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

my $write_obj = AlignDB::ToXLSX->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);
my $dbh = $write_obj->dbh;

my $sql_file = AlignDB::SQL::Library->new( lib => "$FindBin::Bin/sql.lib" );

# auto detect combine threshold
if ( $combine == 0 ) {
    ($combine) = $write_obj->calc_threshold;
}

# auto detect combine threshold
if ( $piece == 0 ) {
    ( undef, $piece ) = $write_obj->calc_threshold;
}

print Dump {
    combine => $combine,
    piece   => $piece,
};

#----------------------------#
# count freq
#----------------------------#
my $all_freq = get_freq($dbh);

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

    # for 22 flies, low, mid, high
    {
        my @all_freqs = 1 .. $all_freq - 1;
        if ( scalar @all_freqs <= 3 ) {
            for (@all_freqs) {
                my $name = $_ . "of" . $all_freq;
                push @freqs, [ $name, $_, $_ ];
            }
        }
        else {
            my @to_be_combs = @all_freqs[ 0 .. $all_freq - 2 ];
            my @chunks      = reverse apportion( scalar @to_be_combs, 3 );
            my @chunks_freq = multi_slice( \@to_be_combs, @chunks );
            for my $chunk (@chunks_freq) {
                if ( $chunk->[0] == $chunk->[-1] ) {
                    my $name = $chunk->[0] . "of" . $all_freq;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
                else {
                    my $name = join( '_', $chunk->[0], $chunk->[-1] ) . "of" . $all_freq;
                    push @freqs, [ $name, $chunk->[0], $chunk->[-1] ];
                }
            }
        }
    }
}

#----------------------------------------------------------#
# chart -- d1_indel_ld
#----------------------------------------------------------#
my $chart_d1_indel_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %option = (
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
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 4;
    $option{y_title}  = "|Dprime|";
    $option{y_data}   = $data->[4];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column}  = 1;
    $option{y_title}   = "r";
    $option{y_data}    = $data->[1];
    $option{y2_column} = 3;
    $option{y2_data}   = $data->[3];
    $option{y2_title}  = "Dprime";
    $option{top}       = 1;
    $option{left}      = 12;
    $write_obj->draw_2y( $sheet, \%option );
};

#----------------------------------------------------------#
# chart -- d2_indel_ld
#----------------------------------------------------------#
my $chart_d2_indel_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %option = (
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
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 4;
    $option{y_title}  = "|Dprime|";
    $option{y_data}   = $data->[4];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column}  = 1;
    $option{y_title}   = "r";
    $option{y_data}    = $data->[1];
    $option{y2_column} = 3;
    $option{y2_data}   = $data->[3];
    $option{y2_title}  = "Dprime";
    $option{top}       = 1;
    $option{left}      = 12;
    $write_obj->draw_2y( $sheet, \%option );
};

#----------------------------------------------------------#
# chart -- d1_snps_ld
#----------------------------------------------------------#
my $chart_d1_snps_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %option = (
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
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 3;
    $option{y_title}  = "|Dprime| to nearest indel ";
    $option{y_data}   = $data->[3];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 2;
    $option{y_title}  = "r^2 to near snps";
    $option{y_data}   = $data->[2];
    $option{top}      = 1;
    $option{left}     = 17;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 5;
    $option{y_title}  = "indel group snps r^2";
    $option{y_data}   = $data->[5];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 6;
    $option{y_title}  = "nonindel group snps r^2";
    $option{y_data}   = $data->[6];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );
};

#----------------------------------------------------------#
# chart -- d2_snps_ld
#----------------------------------------------------------#
my $chart_d2_snps_ld = sub {
    my $sheet = shift;
    my $data  = shift;

    my %option = (
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
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 3;
    $option{y_title}  = "|Dprime| to nearest indel ";
    $option{y_data}   = $data->[3];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 2;
    $option{y_title}  = "r^2 to near snps";
    $option{y_data}   = $data->[2];
    $option{top}      = 1;
    $option{left}     = 17;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 5;
    $option{y_title}  = "indel group snps r^2";
    $option{y_data}   = $data->[5];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );

    $option{y_column} = 6;
    $option{y_title}  = "nonindel group snps r^2";
    $option{y_data}   = $data->[6];
    $option{top} += 18;
    $write_obj->draw_y( $sheet, \%option );
};

#----------------------------------------------------------#
# worksheet -- indel_ld
#----------------------------------------------------------#
my $indel_ld = sub {

    {
        my $sheet_name = 'd1_indel_ld';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->limit(20);

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute;
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
            $chart_d1_indel_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }

    {
        my $sheet_name = 'd2_indel_ld';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->limit(35);
        $thaw_sql->replace( { distance => 'density' } );

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute;
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
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
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->limit(20);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1] );
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
            $chart_d1_indel_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_indel_ld_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-indel_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1] );
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
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
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->limit(20);

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute;
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
            $chart_d1_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    }

    {
        my $sheet_name = 'd2_snps_ld';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute;
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
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
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->limit(20);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1] );
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
            $chart_d1_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_snps_ld_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);
        $thaw_sql->add_where( 'indel.indel_type' => \'= ?' );

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1] );
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
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
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->limit(20);
        $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
        $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1], $level->[2] );
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
            $chart_d1_snps_ld->( $sheet, $data );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    my $write_sheet_d2 = sub {
        my ($level) = @_;
        my $sheet_name = 'd2_snps_ld_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        my $thaw_sql = $sql_file->retrieve('ld-snps_ld-0');
        $thaw_sql->replace( { distance => 'density' } );
        $thaw_sql->limit(35);
        $thaw_sql->add_where( 'indel.indel_freq' => \'>= ?' );
        $thaw_sql->add_where( 'indel.indel_freq' => \'<= ?' );

        my @names = $thaw_sql->as_header;
        my $data  = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
            my $sth = $dbh->prepare( $thaw_sql->as_sql );
            $sth->execute( $level->[1], $level->[2] );
            while ( my @row = $sth->fetchrow_array ) {
                for my $i ( 0 .. $#names ) {
                    push @{ $data->[$i] }, $row[$i];
                }
            }
            $sheet->write( $sheet_row, 0, $data, $write_obj->format->{NORMAL} );
        }

        {    # chart
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
                    ORDER BY gc DESC, r2 DESC, pi, indel
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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        my @names
            = qw{AVG_gc AVG_pi AVG_Indel/100bp AVG_CV AVG_coding AVG_r2 AVG_length COUNT SUM_length};
        my $data = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        {    # write data
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

            $sheet->write( $sheet_row, 0, [ [@group_names] ], $write_obj->format->{NAME} );
            $sheet->write( $sheet_row, 1, $data, $write_obj->format->{NORMAL} );

        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {    # chart
            my %option = (
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
            $write_obj->draw_xy( $sheet, \%option );
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
                    ORDER BY cv DESC, r2 DESC, pi, indel
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
            my %option = (
                sql_query => $sql_query,
                piece     => $piece,
            );
            @combined_segment = @{ $write_obj->make_combine_piece( \%option ) };
        }

        my @names
            = qw{AVG_CV AVG_pi AVG_Indel/100bp AVG_gc AVG_coding AVG_r2 AVG_length COUNT SUM_length Range_gc};
        my $data = [];
        push @{$data}, [] for @names;

        {    # write header
            ( $sheet_row, $sheet_col ) = ( 0, 1 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@names,
            );
            ( $sheet, $sheet_row ) = $write_obj->write_header_direct( $sheet_name, \%option );
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
            $sheet->write( $sheet_row, 0, [ [@group_names] ], $write_obj->format->{NAME} );
            $sheet->write( $sheet_row, 1, $data, $write_obj->format->{NORMAL} );
        }

        {    # drop temporary table
            my $sql_query = q{DROP TABLE IF EXISTS tmp_group};
            my %option = ( sql_query => $sql_query, );
            $write_obj->excute_sql( \%option );
        }

        {    # chart
            my %option = (
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
            $write_obj->draw_xy( $sheet, \%option );
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

sub get_freq {
    my $dbh = shift;

    my $sql_query = q{
        SELECT DISTINCT COUNT(q.query_id) + 1
        FROM  query q, sequence s
        WHERE q.seq_id = s.seq_id
        GROUP BY s.align_id
    };
    my $sth = $dbh->prepare($sql_query);

    my @counts;
    $sth->execute;
    while ( my ($count) = $sth->fetchrow_array ) {
        push @counts, $count;
    }
    if ( scalar @counts > 1 ) {
        warn "Database corrupts, freqs are not consistent\n";
    }

    return $counts[0];
}

__END__
