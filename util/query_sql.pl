#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use YAML::Syck;

use Text::CSV_XS;
use Text::Table;

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf_db = Config::Tiny->read("$FindBin::Bin/../alignDB.ini")->{database};

my ( $opt, $usage ) = Getopt::Long::Descriptive::describe_options(
    "usage: %c %o",
    [ 'help|h', 'display this message' ],
    [],
    ['Database init values'],
    [ 'server|s=s',   'MySQL IP/Domain', { default => $conf_db->{server} } ],
    [ 'port=i',       'MySQL port',      { default => $conf_db->{port} } ],
    [ 'username|u=s', 'username',        { default => $conf_db->{username} } ],
    [ 'password|p=s', 'password',        { default => $conf_db->{password} } ],
    [ 'db|d=s',       'database name',   { default => $conf_db->{db} } ],
    [],
    [ 'query|q=s',  'SQL statement', { default => "SELECT * FROM meta" } ],
    [ 'file|f=s',   'SQL file', ],
    [ 'output|o=s', 'output filename. [stdout] for screen' ],
    [ 'type|t=s', 'output style (csv, neat, table and box)', { default => "csv" } ],
);

$usage->die( { pre_text => "Write sql query results to a file, supporting multiple styles\n" } )
    if $opt->{help};

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#

# Database handler
my $dsn
    = "dbi:mysql:" . "database=" . $opt->{db} . ";host=" . $opt->{server} . ";port=" . $opt->{port};
my DBI $dbh = DBI->connect( $dsn, $opt->{username}, $opt->{password} )
    or die $DBI::errstr;

# Execute sql query
if ( $opt->{file} ) {
    open my $in_fh, '<', $opt->{file};
    my $content = do { local $/; <$in_fh> };
    close $in_fh;
    my @queries = grep {/select/i} split /\;/, $content;

    for my $i ( 0 .. $#queries ) {
        my $outfile = $opt->{output} ? $opt->{output} : "$opt->{db}.$opt->{type}";
        my $index = $i + 1;
        $outfile =~ s/\.(\w+)$/\.$index.$1/ if @queries > 1;
        result( $queries[$i], $opt->{type}, $outfile );
    }
}
elsif ( $opt->{query} ) {
    my $outfile = $opt->{output} ? $opt->{output} : "$opt->{db}.$opt->{type}";
    result( $opt->{query}, $opt->{type}, $outfile );
}
else {
    $usage->die( { pre_text => "Provide a file or a SQL query string\n" } );
}

exit;

sub result {
    my $sql     = shift;
    my $type    = shift;
    my $outfile = shift;

    my DBI $sth = $dbh->prepare($sql)
        or die $dbh->errstr;
    $sth->execute
        or die $sth->errstr;

    my $out_fh;
    if ( lc($outfile) eq "stdout" ) {
        $out_fh = *STDOUT;
    }
    else {
        open $out_fh, ">", $outfile;
    }

    if ( $type eq 'csv' ) {
        my $csv = Text::CSV_XS->new;

        # header line
        my @columns = @{ $sth->{NAME} };
        $csv->combine(@columns);
        print {$out_fh} $csv->string . "\n";

        # all others
        while ( my @row = $sth->fetchrow_array ) {
            $csv->combine(@row);
            print {$out_fh} $csv->string . "\n";
        }
    }
    elsif ( $type eq 'neat' ) {

        # header line
        my @columns = @{ $sth->{NAME} };
        print {$out_fh} DBI::neat_list( \@columns ) . "\n";

        # all others
        while ( my @row = $sth->fetchrow_array ) {
            print {$out_fh} DBI::neat_list( \@row ) . "\n";
        }
    }
    elsif ( $type eq 'box' or $type eq 'table' ) {
        my $is_box = $type eq 'box';

        # header line
        my @columns = map +{ title => $_, align_title => 'center' }, @{ $sth->{NAME} };
        my $c = 0;
        splice @columns, $_ + $c++, 0, \' | ' for 1 .. $#columns;
        my @header_border = ( $is_box ? \' |' : () );
        my $table = Text::Table->new( @header_border, @columns, @header_border );

        # all others
        while ( my @row = $sth->fetchrow_array ) {
            $table->load( \@row );
        }
        my $rule = $table->rule(qw/- +/);
        my @rows_border = ( $is_box ? $rule : () );
        print {$out_fh} join '', @rows_border, $table->title, $rule, $table->body, @rows_border;
    }
    else {
        die "Unknown output style type!\n";
    }

    close $out_fh;

    return;
}

__END__

