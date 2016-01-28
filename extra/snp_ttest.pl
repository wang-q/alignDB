#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::MoreUtils qw(firstidx all any uniq );
use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::WriteExcel;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};

my $db      = 'testdata';
my $outfile = '';
my $run     = 'all';

my $help = 0;
my $man  = 0;

GetOptions(
    'help|?'     => \$help,
    'server=s'   => \$server,
    'Port=i'     => \$port,
    'username=s' => \$username,
    'password=s' => \$password,
    'Database=s' => \$db,
    'run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.S.ttest.xls" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 20 );
}
else {
    $run =~ s/\"\'//s;
    my $set = AlignDB::IntSpan->new();
    if ( AlignDB::IntSpan->valid($run) ) {
        $set   = $set->add($run);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $run;
        $set->add(@tasks);
    }

    my $runlist = $set->runlist();
    $outfile =~ s/(\.xls)$/.$runlist$1/;
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Stat $db...");

my $write_excel_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

#----------------------------#
# stat1 -- diff between indel branch and no-indel branch
#----------------------------#
my $stat1 = sub {
    my @isw_levels = (
        [ 'one_third_all', 1, 0, 999 ],
        [ 'one_third_0',   1, 0, 0 ],
        [ 'one_third_1',   1, 1, 1 ],
        [ 'one_third_2',   1, 2, 2 ],
        [ 'one_third_3',   1, 3, 3 ],
        [ 'one_third_4',   1, 4, 4 ],
        [ 'one_third_5',   1, 5, 5 ],
        [ 'one_third_0-1', 1, 0, 1 ],
        [ 'one_third_2-5', 1, 2, 5 ],

        [ 'two_third_all', 2, 0, 999 ],
        [ 'two_third_0',   2, 0, 0 ],
        [ 'two_third_1',   2, 1, 1 ],
        [ 'two_third_2',   2, 2, 2 ],
        [ 'two_third_3',   2, 3, 3 ],
        [ 'two_third_4',   2, 4, 4 ],
        [ 'two_third_5',   2, 5, 5 ],
        [ 'two_third_0-1', 2, 0, 1 ],
        [ 'two_third_2-5', 2, 2, 5 ],
    );

    my $write_sheet = sub {
        my ($isw_level) = @_;
        my ( $isw_name, $indel_freq, $isw_from, $isw_to ) = @$_;

        my $sheet_name = "Stat1_" . $isw_name;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $header
                = [qw{indel_id AVG_S_i AVG_S_n S_ibranch indel_length}];
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                header    => $header,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_excel_obj->write_header_direct( $sheet_name,
                \%option );
        };

        my @sheet_rows;

        my $dbh       = $write_excel_obj->get_dbhandle;
        my $indel_sth = $dbh->prepare(
            qq{
            SELECT  i.indel_id, i.occured, i.indel_start, i.indel_end,
                    i.indel_length
            FROM indel i, indel_anno a
            WHERE i.indel_id = a.indel_id
            AND i.occured != 'unknown'
            AND i.frequency = $indel_freq
            AND a.slippage = 'N'
            AND i.indel_length <= 50
            }
        );

        my $snp_sth = $dbh->prepare(
            q{
            SELECT s.ref_base, s.other_base, s.occured, s.position
            FROM snp s, isw w, indel i
            WHERE s.isw_id = w.isw_id
            AND w.indel_id = i.indel_id
            AND w.distance_to_indel BETWEEN ? AND ?
            AND s.frequency != 3
            AND i.indel_id = ?
            }
        );

        $indel_sth->execute;
    INDEL: while ( my @row = $indel_sth->fetchrow_array ) {
            my ($indel_id,  $indel_occured, $indel_start,
                $indel_end, $indel_length
            ) = @row;

            my $group_i     = AlignDB::IntSpan->new();
            my $group_n     = AlignDB::IntSpan->new();
            my @indel_occur = split '', $indel_occured;
            my $align_cnt   = scalar @indel_occur;

            for my $i ( 0 .. $align_cnt - 1 ) {
                if ( $indel_occur[$i] eq 'o' ) {
                    $group_i->add($i);
                }
                elsif ( $indel_occur[$i] eq 'x' ) {
                    $group_n->add($i);
                }
                else {
                    die "$indel_occur[$i]\n";
                }
            }
            my @group_is = $group_i->elements;
            my @group_ns = $group_n->elements;

            my @snp_cnt          = (0) x $align_cnt;
            my $indel_branch_cnt = 0;

            $snp_sth->execute( $isw_from, $isw_to, $indel_id );
        SNP: while ( my @row = $snp_sth->fetchrow_array ) {
                my ( $ref_base, $other_base, $snp_occured, $snp_position )
                    = @row;

                my @snp_bases = split '', $other_base;
                my $base_cnt = scalar uniq( $ref_base, @snp_bases );
                next SNP if $base_cnt != 2;
                next SNP if any { $_ =~ /n/i } @snp_bases;

                ## drop all snp close to indel (10 bp)
                #if ( $snp_position < $indel_start ) {
                #    if ( $indel_start - $snp_position - 1 <= 10 ) {
                #        next SNP;
                #    }
                #}
                #elsif ( $snp_position > $indel_end ) {
                #    if ( $snp_position - $indel_end - 1 <= 10 ) {
                #        next SNP;
                #    }
                #}
                #else {
                #    die "snp position error\n";
                #}

                for my $i ( 0 .. $align_cnt - 1 ) {
                    if ( $snp_bases[$i] ne $ref_base ) {
                        $snp_cnt[$i]++;
                    }
                }

                if ( $snp_occured eq $indel_occured ) {
                    $indel_branch_cnt++;
                }
            }

            my $S_i       = mean( @snp_cnt[@group_is] );
            my $S_n       = mean( @snp_cnt[@group_ns] );
            my $S_ibranch = $indel_branch_cnt;
            my $cur_row
                = [ $indel_id, $S_i, $S_n, $S_ibranch, $indel_length ];
            push @sheet_rows, $cur_row;
        }

        foreach (@sheet_rows) {
            my %option = (
                row       => $_,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row)
                = $write_excel_obj->write_row_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@isw_levels) {
        &$write_sheet($_);
    }
};

#----------------------------#
# stat1 -- diff between indel branch and no-indel branch
#----------------------------#
my $stat1_all = sub {
    my @isw_levels = (
        [ 'total_all', 0, 999 ],
        [ 'total_0',   0, 0 ],
        [ 'total_1',   1, 1 ],
        [ 'total_2',   2, 2 ],
        [ 'total_3',   3, 3 ],
        [ 'total_4',   4, 4 ],
        [ 'total_5',   5, 5 ],
        [ 'total_0-1', 0, 1 ],
        [ 'total_2-5', 2, 5 ],
    );

    my $write_sheet = sub {
        my ($isw_level) = @_;
        my ( $isw_name, $isw_from, $isw_to ) = @$_;

        my $sheet_name = "Stat1_" . $isw_name;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $header
                = [qw{indel_id AVG_S_i AVG_S_n S_ibranch indel_length}];
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                header    => $header,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_excel_obj->write_header_direct( $sheet_name,
                \%option );
        };

        my @sheet_rows;

        my $dbh       = $write_excel_obj->get_dbhandle;
        my $indel_sth = $dbh->prepare(
            qq{
            SELECT  i.indel_id, i.occured, i.indel_start, i.indel_end,
                    i.indel_length
            FROM indel i, indel_anno a
            WHERE i.indel_id = a.indel_id
            AND i.occured != 'unknown'
            AND i.frequency IN ( 1, 2 )
            AND a.slippage = 'N'
            }
        );

        my $snp_sth = $dbh->prepare(
            q{
            SELECT s.ref_base, s.other_base, s.occured, s.position
            FROM snp s, isw w, indel i
            WHERE s.isw_id = w.isw_id
            AND w.indel_id = i.indel_id
            AND w.distance_to_indel BETWEEN ? AND ?
            AND s.frequency != 3
            AND i.indel_id = ?
            }
        );

        $indel_sth->execute;
    INDEL: while ( my @row = $indel_sth->fetchrow_array ) {
            my ($indel_id,  $indel_occured, $indel_start,
                $indel_end, $indel_length
            ) = @row;

            my $group_i     = AlignDB::IntSpan->new();
            my $group_n     = AlignDB::IntSpan->new();
            my @indel_occur = split '', $indel_occured;
            my $align_cnt   = scalar @indel_occur;

            for my $i ( 0 .. $align_cnt - 1 ) {
                if ( $indel_occur[$i] eq 'o' ) {
                    $group_i->add($i);
                }
                elsif ( $indel_occur[$i] eq 'x' ) {
                    $group_n->add($i);
                }
                else {
                    die "$indel_occur[$i]\n";
                }
            }
            my @group_is = $group_i->elements;
            my @group_ns = $group_n->elements;

            my @snp_cnt          = (0) x $align_cnt;
            my $indel_branch_cnt = 0;

            $snp_sth->execute( $isw_from, $isw_to, $indel_id );
        SNP: while ( my @row = $snp_sth->fetchrow_array ) {
                my ( $ref_base, $other_base, $snp_occured, $snp_position )
                    = @row;

                my @snp_bases = split '', $other_base;
                my $base_cnt = scalar uniq( $ref_base, @snp_bases );
                next SNP if $base_cnt != 2;
                next SNP if any { $_ =~ /n/i } @snp_bases;

                ## drop all snp close to indel (10 bp)
                #if ( $snp_position < $indel_start ) {
                #    if ( $indel_start - $snp_position - 1 <= 10 ) {
                #        next SNP;
                #    }
                #}
                #elsif ( $snp_position > $indel_end ) {
                #    if ( $snp_position - $indel_end - 1 <= 10 ) {
                #        next SNP;
                #    }
                #}
                #else {
                #    die "snp position error\n";
                #}

                for my $i ( 0 .. $align_cnt - 1 ) {
                    if ( $snp_bases[$i] ne $ref_base ) {
                        $snp_cnt[$i]++;
                    }
                }

                if ( $snp_occured eq $indel_occured ) {
                    $indel_branch_cnt++;
                }
            }

            my $S_i       = mean( @snp_cnt[@group_is] );
            my $S_n       = mean( @snp_cnt[@group_ns] );
            my $S_ibranch = $indel_branch_cnt;
            my $cur_row
                = [ $indel_id, $S_i, $S_n, $S_ibranch, $indel_length ];
            push @sheet_rows, $cur_row;
        }

        foreach (@sheet_rows) {
            my %option = (
                row       => $_,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row)
                = $write_excel_obj->write_row_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@isw_levels) {
        &$write_sheet($_);
    }
};

#----------------------------#
# stat2 -- bigbang or something else
#----------------------------#
my $stat2 = sub {
    my @isw_levels = (
        [ 'all', 0, 999 ],
        [ '0',   0, 0 ],
        [ '1',   1, 1 ],
        [ '2',   2, 2 ],
        [ '3',   3, 3 ],
        [ '4',   4, 4 ],
        [ '5',   5, 5 ],
        [ '0-1', 0, 1 ],
        [ '2-5', 2, 5 ],
    );

    my $write_sheet = sub {
        my ($isw_level) = @_;
        my ( $isw_name, $isw_from, $isw_to ) = @$_;

        my $sheet_name = "Stat2_two_third_" . $isw_name;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $header = [qw{indel_id AVG_S_ibranch S_n}];
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                header    => $header,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_excel_obj->write_header_direct( $sheet_name,
                \%option );
        };

        my @sheet_rows;

        my $dbh       = $write_excel_obj->get_dbhandle;
        my $indel_sth = $dbh->prepare(
            q{
            SELECT i.indel_id, i.occured, i.indel_start, i.indel_end
            FROM indel i, indel_anno a
            WHERE i.indel_id = a.indel_id
            AND i.occured != 'unknown'
            AND i.frequency = 2
            AND a.slippage = 'N'
            }
        );

        my $snp_sth = $dbh->prepare(
            q{
            SELECT s.ref_base, s.other_base, s.occured, s.position
            FROM snp s, isw w, indel i
            WHERE s.isw_id = w.isw_id
            AND w.indel_id = i.indel_id
            AND w.distance_to_indel BETWEEN ? AND ?
            AND s.frequency != 3
            AND i.indel_id = ?
            }
        );

        $indel_sth->execute;
    INDEL: while ( my @row = $indel_sth->fetchrow_array ) {
            my ( $indel_id, $indel_occured, $indel_start, $indel_end ) = @row;

            my $group_i     = AlignDB::IntSpan->new();
            my $group_n     = AlignDB::IntSpan->new();
            my @indel_occur = split '', $indel_occured;
            my $align_cnt   = scalar @indel_occur;

            for my $i ( 0 .. $align_cnt - 1 ) {
                if ( $indel_occur[$i] eq 'o' ) {
                    $group_i->add($i);
                }
                elsif ( $indel_occur[$i] eq 'x' ) {
                    $group_n->add($i);
                }
                else {
                    die "$indel_occur[$i]\n";
                }
            }
            my @group_is = $group_i->elements;
            my @group_ns = $group_n->elements;

            my @induced_cnt = (0) x $align_cnt;
            my @noindel_cnt = (0) x $align_cnt;

            $snp_sth->execute( $isw_from, $isw_to, $indel_id );
        SNP: while ( my @row = $snp_sth->fetchrow_array ) {
                my ( $ref_base, $other_base, $snp_occured, $snp_position )
                    = @row;

                my @snp_bases = split '', $other_base;
                my $base_cnt = scalar uniq( $ref_base, @snp_bases );
                next SNP if $base_cnt != 2;
                next SNP if any { $_ =~ /n/i } @snp_bases;

                ## drop all snp close to indel (10 bp)
                #if ( $snp_position < $indel_start ) {
                #    if ( $indel_start - $snp_position - 1 <= 10 ) {
                #        next SNP;
                #    }
                #}
                #elsif ( $snp_position > $indel_end ) {
                #    if ( $snp_position - $indel_end - 1 <= 10 ) {
                #        next SNP;
                #    }
                #}
                #else {
                #    die "snp position error\n";
                #}

                my $ibranch_base_cnt = scalar uniq @snp_bases[@group_is];
                if ( $ibranch_base_cnt > 1 ) {
                    for my $i (@group_is) {
                        if ( $snp_bases[$i] ne $ref_base ) {
                            $induced_cnt[$i]++;
                        }
                    }
                }

                for my $i (@group_ns) {
                    if ( $snp_bases[$i] ne $ref_base ) {
                        $noindel_cnt[$i]++;
                    }
                }
            }

            my $AVG_S_ibranch = mean( @induced_cnt[@group_is] );
            my $S_n           = mean( @noindel_cnt[@group_ns] );
            my $cur_row       = [ $indel_id, $AVG_S_ibranch, $S_n ];
            push @sheet_rows, $cur_row;

        }

        foreach (@sheet_rows) {
            my %option = (
                row       => $_,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row)
                = $write_excel_obj->write_row_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@isw_levels) {
        &$write_sheet($_);
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$stat1; next; }
    if ( $n == 2 ) { &$stat2; next; }
    if ( $n == 3 ) { &$stat1_all; next; }
}

$stopwatch->end_message();
exit;

