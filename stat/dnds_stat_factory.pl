#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::IntSpan;
use AlignDB::WriteExcel;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# stat parameter
my $run               = $Config->{stat}{run};
my $combine_threshold = $Config->{stat}{combine_threshold};
my $outfile           = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'                 => \$help,
    'man'                    => \$man,
    's|server=s'             => \$server,
    'P|port=s'               => \$port,
    'd|db=s'                 => \$db,
    'u|username=s'           => \$username,
    'p|password=s'           => \$password,
    'o|output=s'             => \$outfile,
    'r|run=s'                => \$run,
    'ct|combine_threshold=i' => \$combine_threshold,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.dnds.xlsx" unless $outfile;

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
    $outfile =~ s/(\.xlsx)$/.$runlist$1/;
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Do stat for $db...");

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

#----------------------------------------------------------#
# worksheet --
#----------------------------------------------------------#
my $combined_dnds = sub {

    # make combine
    my @combined;
    {
        my $sql_query = q{
            SELECT isw_distance distance,
                   COUNT(*) COUNT
            FROM isw i
            GROUP BY isw_distance
        };
        my $standalone = [ -1, 0 ];
        my %option = (
            sql_query  => $sql_query,
            threshold  => $combine_threshold,
            standalone => $standalone,
        );
        @combined = @{ $write_obj->make_combine( \%option ) };
    }

    #----------------------------------------------------------#
    # worksheet -- combined_distance
    #----------------------------------------------------------#
    {
        my $sheet_name = 'combined_dnds';
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{AVG_distance AVG_pi AVG_d_syn AVG_d_nsy AVG_d_stop
                COUNT dn/ds};
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
                SELECT AVG(i.isw_distance) AVG_distance,
                       AVG(i.isw_pi) AVG_pi,
                       AVG(i.isw_syn) AVG_d_syn,
                       AVG(i.isw_nsy) AVG_d_nsy,
                       AVG(i.isw_stop) AVG_d_stop,
                       COUNT(*) COUNT,
                       AVG(i.isw_nsy) / AVG(i.isw_syn)  `dn/ds`
                FROM isw i
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

};

#----------------------------------------------------------#
# worksheet -- distance(frequecy)
#----------------------------------------------------------#
my $frequency_dnds = sub {
    unless ( $write_obj->check_column( 'indel', 'indel_freq' ) ) {
        return;
    }

    my @freq_levels = ( [ 1, 1, 1 ] );

    my $write_sheet = sub {
        my ($level) = @_;
        my $sheet_name = 'dnds_freq_' . $level->[0];
        my $sheet;
        my ( $sheet_row, $sheet_col );

        {    # write header
            my @headers = qw{AVG_distance AVG_pi AVG_d_syn AVG_d_nsy AVG_d_stop
                COUNT dn/ds};
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
                SELECT i.isw_distance distance,
                       AVG(i.isw_pi) AVG_pi,
                       AVG(i.isw_syn) AVG_d_syn,
                       AVG(i.isw_nsy) AVG_d_nsy,
                       AVG(i.isw_stop) AVG_d_stop,
                       COUNT(*) COUNT,
                       AVG(i.isw_nsy) / AVG(i.isw_syn)  `dn/ds`
                FROM    isw i, indel
                WHERE i.isw_indel_id = indel.indel_id
                AND indel.indel_freq >= ?
                AND indel.indel_freq <= ?
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
    if ( $n == 2 ) { &$combined_dnds;  next; }
    if ( $n == 3 ) { &$frequency_dnds; next; }
}

$stopwatch->end_message();
exit;

__END__

=head1 NAME

    gene_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    gene_stat_factory.pl [options]
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
