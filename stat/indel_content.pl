#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Set::Light;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::WriteExcel;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# indel length and k-nuc length
my $min_length = $Config->{indel}->{min_length};
my $max_length = $Config->{indel}->{max_length};
my $min_k      = $Config->{indel}->{min_k};
my $max_k      = $Config->{indel}->{max_k};

my $outfile;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'server=s'     => \$server,
    'port=i'       => \$port,
    'db=s'         => \$db,
    'username=s'   => \$username,
    'password=s'   => \$password,
    'output=s'     => \$outfile,
    'min_length=i' => \$min_length,
    'max_length=i' => \$max_length,
    'min_k=i'      => \$min_k,
    'max_k=i'      => \$max_k,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$outfile = "$db.indel.length$min_length-$max_length.xlsx" unless $outfile;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis indel content of $db...");

my $write_obj = AlignDB::WriteExcel->new(
    mysql   => "$db:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

# Database handler
my $dbh = $write_obj->dbh;

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

# select all indels in this alignment
my $indel_query = q{
    SELECT indel_seq
    FROM indel
    WHERE align_id = ?
    AND indel_seq NOT LIKE "%N%"
    AND indel_length between ? and ?
};
my $indel_sth = $dbh->prepare($indel_query);

foreach my $k ( $min_k .. $max_k ) {

    # write table
    my $write_sheet = sub {
        my $sheet_name = shift;
        my $rows       = shift;
        my $sum        = shift;
        my $sheet;
        my ( $sheet_row, $sheet_col );

        # write header
        {
            my $header = [qw{k_nuc gc count freq}];
            ( $sheet_row, $sheet_col ) = ( 0, 0 );
            my %option = (
                header    => $header,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ( $sheet, $sheet_row )
                = $write_obj->write_header_direct( $sheet_name, \%option );
        };

        for ( @{$rows} ) {
            my %option = (
                row       => $_,
                sheet_row => $sheet_row,
                sheet_col => $sheet_col,
            );
            ($sheet_row) = $write_obj->write_row_direct( $sheet, \%option );
        }

        # write footer
        {
            $sheet_row += 2;
            my $row = [ 'Total', '', $sum, '1.00' ];
            my %option = (
                row            => $row,
                sheet_row      => $sheet_row,
                sheet_col      => $sheet_col,
                content_format => 'TOTAL',
            );
            ($sheet_row) = $write_obj->write_row_direct( $sheet, \%option );
        }

        print "Sheet \"$sheet_name\" has been generated.\n";
    };

    # total indel number, k_nuc rows, k_nuc_rc rows
    my ( $sum, @rows, @rows2 );
    {
        my %table;

        # for each indel
        my $align_ids = $write_obj->get_align_ids;
        for my $align_id ( @{$align_ids} ) {
            $indel_sth->execute( $align_id, $min_length, $max_length );
            while ( my ($indel_seq) = $indel_sth->fetchrow_array ) {
                k_nuc_incr( \$indel_seq, $k, \%table );
            }
        }

        {    # remove invalid keys
            my @k_nucs    = k_nuc_permu($k);
            my $k_nuc_set = Set::Light->new(@k_nucs);
            for ( keys %table ) {
                unless ( $k_nuc_set->has($_) ) {
                    delete $table{$_};
                }
            }
        }

        # total indel number
        $sum += $table{$_} for keys %table;

        # k_nuc worksheet
        for my $k_nuc ( sort keys %table ) {
            my $row = [
                $k_nuc,         calc_gc_ratio($k_nuc),
                $table{$k_nuc}, $table{$k_nuc} / $sum
            ];
            push @rows, $row;
        }
        @rows = sort { $b->[2] <=> $a->[2] } @rows;

        # k_nuc_rc worksheet, merge revcom
        while (%table) {
            my ($key) = sort keys %table;
            my $value = $table{$key};
            delete $table{$key};

            my $key_rc   = revcom($key);
            my $value_rc = 0;
            if ( exists $table{$key_rc} ) {
                $value_rc = $table{$key_rc};
                delete $table{$key_rc};
            }

            my $k_nuc       = "$key,$key_rc";
            my $k_nuc_count = $value + $value_rc;

            my $row = [
                $k_nuc,       calc_gc_ratio($key),
                $k_nuc_count, $k_nuc_count / $sum
            ];
            push @rows2, $row;
        }
        @rows2 = sort { $b->[2] <=> $a->[2] } @rows2;
    }

    if ( @rows > 0 or @rows2 > 0 ) {

        #$write_sheet->( "k_nuc_$k",    \@rows,  $sum ); # we didn't need this
        $write_sheet->( "k_nuc_rc_$k", \@rows2, $sum );
    }
}

$indel_sth->finish;

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    indel_content.pl - See what's in indels, using the k-nucleotide algorithm

=head1 SYNOPSIS

    update_indel_slippage.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --output          output filename
       --min_length      minimal length of indels
       --max_length      maximal length of indels
       --min_k           minimal number of k-nucleotide
       --max_k           maximal number of k-nucleotide
       

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

