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
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# stat parameter
my $run     = $Config->{stat}->{run};
my $outfile = "";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=s'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'output=s'   => \$outfile,
    'run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.dnds.xls" unless $outfile;

# prepare to run tasks in @tasks
my @tasks;

if ( $run eq 'all' ) {
    @tasks = ( 1 .. 20 );
}
else {
    $run =~ s/\"\'//s;
    my $set = AlignDB::IntSpan->new();
    if (AlignDB::IntSpan->valid($run)) {
        $set = $set->add($run);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $run;
        $set->add(@tasks) ;
    }
    
    my $runlist = $set->runlist();
    $outfile =~ s/(\.xls)$/.$runlist$1/;
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
my $combined_distance = sub {

    # make combine
    my $sql_query = q~
        # distance_count
        SELECT isw_distance distance,
               COUNT(*) COUNT
        FROM isw i
        GROUP BY isw_distance
    ~;
    my $threshold  = 1000;
    my $standalone = [ -1, 0 ];
    my %option     = (
        sql_query  => $sql_query,
        threshold  => $threshold,
        standalone => $standalone,
    );
    my @combined_distance = @{ $write_obj->make_combine( \%option ) };

    #----------------------------------------------------------#
    # worksheet -- combined_distance
    #----------------------------------------------------------#
    {
        my $sheet_name = 'distance';
        my $sheet;

        # write header
        $sql_query = q~
            # header of Table density
            SELECT 'AVG_distance', 'AVG_pi',
                    'AVG_d_syn', 'AVG_d_nsy', 'AVG_d_stop',
                    'COUNT', 'dn/ds'
        ~;
        my ( $sheet_row, $sheet_col ) = ( 0, 0 );
        %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row )
            = $write_obj->write_header_sql( $sheet_name, \%option );

        # write contents
        $sql_query = q{
            # distance effect
            SELECT AVG(i.isw_distance) AVG_distance,
                   AVG(i.isw_pi) AVG_pi,
                   AVG(e.isw_feature9) AVG_d_syn,
                   AVG(e.isw_feature10) AVG_d_nsy,
                   AVG(e.isw_feature11) AVG_d_stop,
                   COUNT(*) COUNT,
                   AVG(e.isw_feature10) / AVG(e.isw_feature9)  `dn/ds`
            FROM isw i, isw_extra e
            WHERE i.isw_id = e.isw_id
            AND isw_distance IN
        };
        %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            combined  => \@combined_distance,
        );
        ($sheet_row) = $write_obj->write_content_combine( $sheet, \%option );

        print "Sheet \"$sheet_name\" has been generated.\n";
    }

};

foreach my $n (@tasks) {
    #if ( $n == 1 ) { &$summary_gene;  next; }
    if ( $n == 2 ) { &$combined_distance;  next; }
    #if ( $n == 3 ) { &$gene_D;        &$gene_D_null; next; }
    #if ( $n == 4 ) { &$exon_D;        &$exon_D_null; next; }
    #if ( $n == 5 ) { &$exon_gc;       next; }
    #if ( $n == 6 ) { &$gene_ess;      &$gene_pure; next; }
    #if ( $n == 7 ) { &$gene_5;        &$gene_3; next; }
    #if ( $n == 8 ) { &$exon_ess;      next; }
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
