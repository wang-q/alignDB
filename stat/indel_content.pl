#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
use Set::Scalar;

use App::Fasops::Common;

use AlignDB::Stopwatch;
use AlignDB::ToXLSX;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

=head1 NAME

indel_content.pl - See what's in indels, using the k-nucleotide algorithm

=head1 SYNOPSIS

    perl indel_content.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --output        STR     output filename
        --min_length    INT     minimal length of indels, default is [1]
        --max_length    INT     maximal length of indels, default is [50]
        --min_k         INT     minimal number of k-nucleotide, default is [1]
        --max_k         INT     maximal number of k-nucleotide, default is [8]

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server     = $Config->{database}{server} ),
    'port|P=i'     => \( my $port       = $Config->{database}{port} ),
    'db|d=s'       => \( my $db         = $Config->{database}{db} ),
    'username|u=s' => \( my $username   = $Config->{database}{username} ),
    'password|p=s' => \( my $password   = $Config->{database}{password} ),
    'output=s'     => \( my $outfile ),
    'min_length=i' => \( my $min_length = 1 ),
    'max_length=i' => \( my $max_length = 50 ),
    'min_k=i'      => \( my $min_k      = 1 ),
    'max_k=i'      => \( my $max_k      = 6 ),
) or Getopt::Long::HelpMessage(1);

$outfile = "$db.indel.length$min_length-$max_length.xlsx" unless $outfile;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis indel content of $db...");

my DBI $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password )
    or die "Cannot connect to MySQL database at $db:$server";
my $write_obj = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
);

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

# select all indels in this alignment
my $indel_query = q{
    SELECT indel_seq
    FROM indel
    WHERE 1 = 1
    AND indel_seq NOT LIKE "%N%"
    AND indel_length between ? and ?
};
my DBI $indel_sth = $dbh->prepare($indel_query);

for my $k ( $min_k .. $max_k ) {

    # write table
    my $write_sheet = sub {
        my $sheet_name = shift;
        my $rows       = shift;
        my $sum        = shift;
        my $sheet;
        $write_obj->row(0);
        $write_obj->column(0);

        my @names = qw{k_nuc gc count freq};
        {    # header
            $sheet = $write_obj->write_header( $sheet_name, { header => \@names } );
        }

        for ( @{$rows} ) {
            $write_obj->write_row( $sheet, { row => $_ } );
        }

        print "Sheet [$sheet_name] has been generated.\n";
    };

    # total indel number, k_nuc rows, k_nuc_rc rows
    my ( $sum, @rows, @rows2 );
    {
        my %table;

        # for each indel
        $indel_sth->execute( $min_length, $max_length );
        while ( my ($indel_seq) = $indel_sth->fetchrow_array ) {
            k_nuc_incr( \$indel_seq, $k, \%table );
        }

        {    # remove invalid keys
            my @k_nucs    = k_nuc_permu($k);
            my $k_nuc_set = Set::Scalar->new(@k_nucs);
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
                $k_nuc, App::Fasops::Common::calc_gc_ratio( [$k_nuc] ),
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

            my $key_rc   = App::Fasops::Common::revcom($key);
            my $value_rc = 0;
            if ( exists $table{$key_rc} ) {
                $value_rc = $table{$key_rc};
                delete $table{$key_rc};
            }

            my $k_nuc       = "$key,$key_rc";
            my $k_nuc_count = $value + $value_rc;

            my $row = [
                $k_nuc, App::Fasops::Common::calc_gc_ratio( [$key] ),
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

sub k_nuc_permu {
    my $k         = shift;
    my @alphabets = qw{A C G T};

    my %table;
    $table{$_} = '' for @alphabets;
    for ( 2 .. $k ) {
        for my $current_key ( keys %table ) {
            $table{ $current_key . $_ } = '' for @alphabets;
            delete $table{$current_key};
        }
    }

    return sort keys %table;
}

sub k_nuc_count {
    my $seq_ref = shift;
    my $k       = shift;

    my $seq_length = length $$seq_ref;
    my %table;

    for ( 0 .. $seq_length - $k ) {
        $table{ substr( $$seq_ref, $_, $k ) }++;
    }

    return %table;
}

sub k_nuc_incr {
    my $seq_ref  = shift;
    my $k        = shift;
    my $hash_ref = shift;

    my $seq_length = length $$seq_ref;

    for ( 0 .. $seq_length - $k ) {
        $hash_ref->{ substr( $$seq_ref, $_, $k ) }++;
    }
}

__END__
