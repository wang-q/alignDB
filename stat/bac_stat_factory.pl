#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Template;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::WriteExcel;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

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
my $db       = $Config->{bac}{db};

# stat parameter
my $run     = $Config->{stat}{run};
my $outfile = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=s'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'o|output=s'   => \$outfile,
    'r|run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.bacgc.xlsx" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 50 );
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
$stopwatch->start_message("Do stat for $db...");

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

my $tt = Template->new;

#----------------------------------------------------------#
# worksheet -- taxon_basic
#----------------------------------------------------------#
my $taxon_basic = sub {
    my @columns = qw{
        group genus
    };

    my $write_sheet = sub {
        my $column     = shift;
        my $sheet_name = 'taxon_' . $column;

        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{Type_name AVG_GC STD_GC AVG_CV STD_CV AVG_STD STD_STD AVG_MDCW STD_MDCW COUNT};
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
                header    => \@headers,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        }

        my $template = q{
            my $sql_query = q{
                SELECT  s.`[% column %]`, 
                        avg(t.segment_gc_mean) avg_gc,
                        std(t.segment_gc_mean) std_gc,
                        avg(t.segment_gc_cv) avg_cv,
                        std(t.segment_gc_cv) std_cv,
                        avg(t.segment_gc_std) avg_std,
                        std(t.segment_gc_std) std_std,
                        avg(t.segment_gc_mdcw) avg_mdcw,
                        std(t.segment_gc_mdcw) std_mdcw,
                        count(*) `COUNT`
                FROM    strain s,
                        seq q,
                        segment t
                WHERE   s.project_id = q.project_id
                AND     q.accession = t.accession
                AND     t.segment_type = 1
                GROUP BY s.`[% column %]`
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        };

        {
            my $code;
            $tt->process( \$template, { column => $_, }, \$code )
                or die Template->error;
            eval $code;
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@columns) {
        $write_sheet->($_);
    }

};

#----------------------------------------------------------#
# worksheet -- strain_type
#----------------------------------------------------------#
my $strain_type = sub {
    my $sheet_name = 'strain_type';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $query_name = 'Item';
        my @headers
            = qw{Type_name AVG_GC STD_GC AVG_CV STD_CV AVG_STD STD_STD AVG_MDCW STD_MDCW COUNT};
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

    my $template = q{
        my $query_name = '[% column %]';
        $sheet_row++;
        my $sql_query = q{
            SELECT  s.`[% column %]`, 
                    avg(t.segment_gc_mean) avg_gc,
                    std(t.segment_gc_mean) std_gc,
                    avg(t.segment_gc_cv) avg_cv,
                    std(t.segment_gc_cv) std_cv,
                    avg(t.segment_gc_std) avg_std,
                    std(t.segment_gc_std) std_std,
                    avg(t.segment_gc_mdcw) avg_mdcw,
                    std(t.segment_gc_mdcw) std_mdcw,
                    count(*) `COUNT`
            FROM    strain s,
                    seq q,
                    segment t
            WHERE   s.project_id = q.project_id
            AND     q.accession = t.accession
            AND     t.segment_type = 1
            GROUP BY s.`[% column %]`
            HAVING `COUNT` > 1
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
    };

    my @columns = qw{
        gram_stain	shape	arrangment
        endospores	motility	salinity
        oxygen_req	habitat	temp_range
    };
    for (@columns) {
        my $code;
        $tt->process( \$template, { column => $_, }, \$code )
            or die Template->error;
        eval $code;
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- strain_stat
#----------------------------------------------------------#
my $strain_stat = sub {
    my @levels = ( [ 'full', 1 ], [ '500', 5 ], );

    my $write_sheet = sub {
        my $level      = shift;
        my $sheet_name = 'strain_stat_' . $level->[0];

        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers
                = qw{ organism_name AVG_GC STD_GC AVG_CV STD_CV AVG_STD STD_STD AVG_MDCW STD_MDCW};
            push @headers,
                qw{ super_kingdom group taxonomy_id genome_size chr gram_stain shape arrangment endospores motility salinity oxygen_req habitat temp_range };
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
                SELECT  s.organism_name,
                        avg(t.segment_gc_mean) avg_gc,
                        std(t.segment_gc_mean) std_gc,
                        avg(t.segment_gc_cv) avg_cv,
                        std(t.segment_gc_cv) std_cv,
                        avg(t.segment_gc_std) avg_std,
                        std(t.segment_gc_std) std_std,
                        avg(t.segment_gc_mdcw) avg_mdcw,
                        std(t.segment_gc_mdcw) std_mdcw,
                        s.super_kingdom,
                        s.group,
                        s.taxonomy_id,
                        s.genome_size,
                        s.number_of_chromosomes,
                        s.gram_stain,
                        s.shape,
                        s.arrangment,
                        s.endospores,
                        s.motility,
                        s.salinity,
                        s.oxygen_req,
                        s.habitat,
                        s.temp_range
                FROM    strain s,
                        seq q,
                        segment t
                where s.project_id = q.project_id
                and q.accession = t.accession
                and t.segment_type = ?
                group by s.project_id, t.segment_type
            };
            my %option = (
                sql_query  => $sql_query,
                sheet_row  => $sheet_row,
                sheet_col  => $sheet_col,
                bind_value => [ $level->[1] ],
            );
            ($sheet_row) = $write_obj->write_content_direct( $sheet, \%option );
        };

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    foreach (@levels) {
        $write_sheet->($_);
    }
};

foreach my $n (@tasks) {
    if ( $n == 1 ) { &$taxon_basic; next; }
    if ( $n == 6 ) { &$strain_type; next; }
    if ( $n == 8 ) { &$strain_stat; next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    var_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    three_stat_factory.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --output          output filename
       --run             run special analysis
       

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
