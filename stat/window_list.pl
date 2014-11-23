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
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);
use AlignDB::Position;

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
my $run = $Config->{stat}{run};

my $outfile;

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

$outfile = "$db.window.xlsx" unless $outfile;

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

my $dbh = $write_obj->dbh;
my $pos_obj = AlignDB::Position->new( dbh => $dbh );


#----------------------------------------------------------#
# worksheet -- indel_list
#----------------------------------------------------------#
my $window_list = sub {
    my $sheet_name = 'window_list';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{extreme_id taxon_id chr_name window_start window_end
            extreme_type };
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
                w.align_id, e.extreme_id, a.taxon_id, a.chr_name,
                w.window_start, w.window_end, e.extreme_type
            FROM
                extreme e
                    INNER JOIN
                window w ON e.window_id = w.window_id
                    INNER JOIN
                (SELECT 
                    se.align_id, ta.taxon_id, c.chr_name, se.chr_start
                FROM
                    sequence se, target t, chromosome c, taxon ta
                WHERE
                    se.seq_id = t.seq_id
                        AND se.chr_id = c.chr_id
                        AND c.taxon_id = ta.taxon_id) a ON w.align_id = a.align_id
            ORDER BY a.chr_name , a.chr_start
        };
        my $sth = $dbh->prepare($sql_query);
        $sth->execute();

        while ( my @row = $sth->fetchrow_array ) {
            my $align_id = shift @row;
            for my $i ( 3, 4 ) {
                my $align_pos = $row[$i];
                my $chr_pos = $pos_obj->at_target_chr( $align_id, $align_pos );
                splice @row, $i, 1, $chr_pos;
            }

            ($sheet_row) = $write_obj->write_row_direct(
                $sheet,
                {   row       => \@row,
                    sheet_row => $sheet_row,
                    sheet_col => $sheet_col,
                }
            );
        }
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

foreach my $n (@tasks) {
    if ( $n == 1 )  { &$window_list;    next; }
    #if ( $n == 2 )  { &$snp_basic;      next; }
    #if ( $n == 3 )  { &$indel_list;     next; }
    #if ( $n == 4 )  { &$snp_list;       next; }
    #if ( $n == 5 )  { &$snp_codon_list; next; }
    #if ( $n == 10 ) { &$strain_list;    next; }
}

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    mvar_stat_factory.pl - Generate statistical Excel files from alignDB

=head1 SYNOPSIS

    mvar_stat_factory.pl [options]
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
