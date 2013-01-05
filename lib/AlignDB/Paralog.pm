package AlignDB::Paralog;
use MooX 'late';

use Config::Tiny;
use Bio::SearchIO;
use WWW::Mechanize;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;

has 'Config' => ( is => 'rw', isa => 'Ref' );

# basic attributes
has 'blast_url'     => ( is => 'rw', isa => 'Str' );
has 'megablast_url' => ( is => 'rw', isa => 'Str' );
has 'datalib'       => ( is => 'rw', isa => 'Str' );    # blast database
has 'program'       => ( is => 'rw', isa => 'Str' );    # e.g. blastn
has 'use_megablast' => ( is => 'rw', isa => 'Bool' );

# blast options
has 'identity' => ( is => 'rw', isa => 'Int' );    # identity percentage
has 'expact'   => ( is => 'rw', isa => 'Num' );    # statistical significance
has 'wordsize' => ( is => 'rw', isa => 'Int' );    # megablast wordsize

has 'alignment_limit'   => ( is => 'rw', isa => 'Int' );  # HSPs limit
has 'description_limit' => ( is => 'rw', isa => 'Int' );  # short descriptions

# result format
has 'alignment_view' => ( is => 'rw', isa => 'Int' );   # bioperl defined this
has 'result_format'  => ( is => 'ro', isa => 'Str' );   # blast result format

# paralog options
has 'para_length'   => ( is => 'rw', isa => 'Int' );  # paralog length limit
has 'para_identity' => ( is => 'rw', isa => 'Num' );  # paralog identity limit

sub BUILD {
    my $self = shift;

    my $arg_ref = $self->Config->{paralog};

    $self->{blast_url}         = $arg_ref->{blast_url};
    $self->{megablast_url}     = $arg_ref->{megablast_url};
    $self->{datalib}           = $arg_ref->{datalib};
    $self->{program}           = $arg_ref->{program};
    $self->{use_megablast}     = $arg_ref->{use_megablast};
    $self->{identity}          = $arg_ref->{identity};
    $self->{expact}            = $arg_ref->{expact};
    $self->{wordsize}          = $arg_ref->{wordsize};
    $self->{alignment_limit}   = $arg_ref->{alignment_limit};
    $self->{description_limit} = $arg_ref->{description_limit};
    $self->{alignment_view}    = $arg_ref->{alignment_view};
    $self->{result_format}     = $arg_ref->{result_format};
    $self->{para_length}       = $arg_ref->{para_length};
    $self->{para_identity}     = $arg_ref->{para_identity};

    return;
}

# Fetch blast results
sub web_blast {
    my ( $self, $seq_ref ) = @_;

    $ENV{'http_proxy'} = '';

    my $megablast = $self->use_megablast;

    my $server_url;
    my $program;
    my $identity;
    if ($megablast) {
        $server_url = $self->megablast_url;
        $identity   = $self->identity;
    }
    else {
        $server_url = $self->blast_url;
        $program    = $self->program;
    }

    my $datalib           = $self->datalib;
    my $alignment_view    = $self->alignment_view;
    my $expact            = $self->expact;
    my $wordsize          = $self->wordsize;
    my $alignment_limit   = $self->alignment_limit;
    my $description_limit = $self->description_limit;

    my %view_name = (
        0 => "blast",         # Pairwise
        7 => "blastxml",      # BLAST XML
        9 => "blasttable",    # Hit Table
    );

    $self->{result_format} = $view_name{$alignment_view};

    # init $mech object
    my $mech = WWW::Mechanize->new;

    print " " x 2, "Go to BLAST form page", " " x 10, "\r";
    $mech->get($server_url);

    print " " x 2, "Submit query", " " x 10, "\r";
    $mech->form_name('MainBlastForm');

    # Program:
    # <select name="PROGRAM">
    $mech->field( "PROGRAM", $program ) if $program;

    # Database:
    # <select name="DATALIB">
    $mech->field( "DATALIB", $datalib );

    # Alignment view:
    # <select name="ALIGNMENT_VIEW">
    $mech->field( "ALIGNMENT_VIEW", $alignment_view );

    # Percent Identity:
    # <select name="PERC_IDENT">
    $mech->field( "PERC_IDENT", $identity ) if $identity;

    # Other advanced options:
    # <input name="OTHER_ADVANCED" maxlength="50">
    #
    # Wordsize
    # -W [integer]
    # Don't use field <select name="WORD_SIZE"> because it is only in
    #   megablast page but not in blast page
    # for blastall:
    #   -W  Word size, default is 11 for blastn, 3 for other programs.
    # for megablast:
    #   -W  Word size, default is 28. The default word size is very high
    #   because sequences aligned by megablast are expected to be nearly
    #   identical.
    #
    # alignment_limit
    # -b [integer]
    # Default: 250 Programs: All
    #   Truncates the report to [integer] number of alignments. There is no
    #   warning when you exceed this limit, so it's generally a good idea
    #   to set [integer] very high unless you're interested only in the top
    #   hits.
    #
    # description_limit
    # -v [integer]
    # Default: 500 Programs: All
    #   Sets the number of database sequences for which to show the one-line
    #   summary descriptions at the top of a BLAST report. You won't be
    #   warned if you exceed [integer]. Also see the -b parameter.
    my $advanced_options;
    $advanced_options .= " -W $wordsize"          if $wordsize >= 3;
    $advanced_options .= " -b $alignment_limit"   if $alignment_limit;
    $advanced_options .= " -v $description_limit" if $description_limit;
    $mech->field( "OTHER_ADVANCED", $advanced_options ) if $advanced_options;

    # Expect:
    # <select name="EXPECT">
    $mech->field( "EXPECT", $expact );

    # Filter:
    # <input name="FILTER">
    # Do not filt sequences,to ensure the first hsp is an entire query one
    $mech->untick( "FILTER", "L" );

    # Graphical Overview :
    # <input name="OVERVIEW" value="on">
    # Do not show graphical overview, to make process faster
    $mech->untick( "OVERVIEW", "on" );

    # Sequence:
    # <textarea name="SEQUENCE">
    $mech->field( "SEQUENCE", ${$seq_ref} );

    print " " x 2, "Waiting for server processing", " " x 10, "\r";
    $mech->click;

    print " " x 2, "Return blast report", " " x 10, "\n";
    my $blast_report = $mech->content( format => "text" );

    # drop the first line of blasttable page
    if ( $alignment_view == 9 ) {
        $blast_report =~ s/^BLAST Search Results.+\n//m;
    }

    return \$blast_report;
}

# Parse blast results
sub paralog_cover {
    my ( $self, $report_ref ) = @_;

    my $result_format = $self->result_format;
    my $para_identity = $self->para_identity;
    my $para_length   = $self->para_length;

    # Now parse blast file
    print " " x 2, "Reread in blast report", " " x 10, "\n";

    # open a scalar reference as a file handle
    open my $blast_fh, '<', $report_ref
        or die "Can't open in-memory file: $!";

    my $searchio = new Bio::SearchIO(
        -format => $result_format,
        -fh     => $blast_fh,
    );

    my $bypass_first = 1;
    my $query_set    = AlignDB::IntSpan->new;

    my $result = $searchio->next_result;
    while ( my $hit = $result->next_hit ) {
        while ( my $hsp = $hit->next_hsp ) {

            # bypass the first hsp, which is the target sequence itself
            if ($bypass_first) {
                $bypass_first = 0;
                next;
            }

            # process the Bio::Search::HSP::HSPI object
            my $hsp_length = $hsp->length( ['query'] );
            my $hsp_identity = $hsp->percent_identity;

            next unless $hsp_length > $para_length;
            next unless $hsp_identity > $para_identity;

            my $qstart = $hsp->query->start;
            my $qend   = $hsp->query->end;
            $query_set->add("$qstart-$qend");
        }
    }

    close $blast_fh;

    return $query_set;
}

# Parse blast results
sub longest_paralog {
    my ( $self, $report_ref ) = @_;

    my $result_format = $self->result_format;
    my $para_identity = $self->para_identity;
    my $para_length   = $self->para_length;

    # Now parse blast file
    print " " x 2, "Reread in blast report", " " x 10, "\n";

    open( my $blast_fh, '<', $report_ref )
        or die "Can't open in-memory file: $!";

    my $searchio = new Bio::SearchIO(
        -format => $result_format,
        -fh     => $blast_fh,
    );

    my $bypass_first   = 1;
    my $longest_length = 1;
    my $longest_align;
    my $longest_hit_name;
    my $longest_hit_strand;

    my $result = $searchio->next_result;
    while ( my $hit = $result->next_hit ) {
        my $hit_name = $hit->name;
        while ( my $hsp = $hit->next_hsp ) {

            # bypass the first hsp, which is the target sequence itself
            if ($bypass_first) {
                $bypass_first = 0;
                next;
            }

            # process the Bio::Search::HSP::HSPI object
            my $hsp_length = $hsp->length( ['query'] );
            my $hsp_identity = $hsp->percent_identity;

            # use "+" for default strand
            # -1 = Minus strand, +1 = Plus strand
            my ( $query_strand, $hit_strand ) = $hsp->strand("list");
            my $hsp_strand = "+";
            if ( $query_strand + $hit_strand == 0 and $query_strand != 0 ) {
                $hsp_strand = "-";
            }

            next unless $hsp_length > $para_length;
            next unless $hsp_identity > $para_identity;

            if ( $hsp_length > $longest_length ) {
                $longest_length   = $hsp_length;
                $longest_align    = $hsp->get_aln; # a Bio::SimpleAlign object
                $longest_hit_name = $hit_name;
                $longest_hit_strand = $hsp_strand;
            }
        }
    }

    close $blast_fh;

    return ( $longest_align, $longest_hit_name, $longest_hit_strand );
}

1;
