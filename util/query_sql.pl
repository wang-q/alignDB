#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
use Text::CSV_XS;
use Text::Table;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

my $conf_db = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini")->{database};

my $description = <<'EOF';
Write results of sql query to file, supporting multiple styles.

Usage: perl %c [options]
EOF

(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    $description,
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
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#

# Database handler
my $dsn
    = "dbi:mysql:" . "database=" . $opt->{db} . ";host=" . $opt->{server} . ";port=" . $opt->{port};

#@type DBI
my $dbh = DBI->connect( $dsn, $opt->{username}, $opt->{password} )
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

    #@type DBI
    my $sth = $dbh->prepare($sql)
        or die $dbh->errstr;
    $sth->execute
        or die $sth->errstr;

    my $out_fh;
    if ( lc($outfile) eq "stdout" ) {
        $out_fh = *STDOUT{IO};
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
