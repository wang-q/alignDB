#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Bio::Tools::GFF;

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

# run in parallel mode
my $parallel = $Config->{feature}->{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{feature}->{batch};

my $gff_file;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'gff_file=s' => \$gff_file,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
# other objects are initiated in subroutines
$stopwatch->start_message("Update annotations of $db...");

my $gff_obj = Bio::Tools::GFF->new(
    -file        => $gff_file,
    -gff_version => 3
);

my $cds_set = AlignDB::IntSpan->new;
while ( my $feature = $gff_obj->next_feature ) {
    if ( $feature->primary_tag eq 'CDS' ) {
        $cds_set->add_range( $feature->start, $feature->end );
    }
}

#----------------------------#
# Find all align_ids
#----------------------------#
my @jobs;
{    # create alignDB object for this scope
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    my @align_ids = @{ $obj->get_align_ids };
    while ( scalar @align_ids ) {
        my @batching = splice @align_ids, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job       = shift;
    my @align_ids = @$job;

    #----------------------------#
    # Init objects
    #----------------------------#
    # create alignDB object for this scope
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    #----------------------------#
    # SQL query and DBI sths
    #----------------------------#
    # alignments' chromosomal location
    my $align_seq_query = q{
        SELECT c.chr_name,
               a.align_length,
               s.chr_start,
               s.chr_end,
               s.seq_seq
        FROM align a, target t, sequence s, chromosome c
        WHERE a.align_id = s.align_id
        AND t.seq_id = s.seq_id
        AND s.chr_id = c.chr_id
        AND a.align_id = ?
    };
    my $align_seq_sth = $dbh->prepare($align_seq_query);

    # update align table
    my $align_feature = q{
        UPDATE align
        SET align_coding  = ?,
            align_repeats = ?,
            align_te      = ?
        WHERE align_id    = ?
    };
    my $align_feature_sth = $dbh->prepare($align_feature);

    # select all indels in this alignment
    my $indel_query = q{
        SELECT indel_id, indel_start, indel_end
        FROM indel
        WHERE align_id = ?
    };
    my $indel_query_sth = $dbh->prepare($indel_query);

    # update indel table in the new feature column
    my $indel_feature = q{
        UPDATE indel
        SET indel_coding  = ?,
            indel_repeats = ?
        WHERE indel_id    = ?
    };
    my $indel_feature_sth = $dbh->prepare($indel_feature);

    # select all isws for this indel
    my $isw_query = q{
        SELECT isw_id, isw_start, isw_end
        FROM isw
        WHERE indel_id = ?
    };
    my $isw_query_sth = $dbh->prepare($isw_query);

    # update isw table in the new feature column
    my $isw_feature = q{
        UPDATE isw
        SET isw_coding  = ?,
            isw_repeats = ?
        WHERE isw_id    = ?
    };
    my $isw_feature_sth = $dbh->prepare($isw_feature);

    # select all snps for this indel
    my $snp_query = q{
        SELECT snp_id, snp_pos, mutant_to
        FROM snp
        WHERE align_id = ?
    };
    my $snp_query_sth = $dbh->prepare($snp_query);

    # update snp table in the new feature column
    my $snp_feature = q{
        UPDATE snp
        SET snp_coding  = ?,
            snp_repeats = ?,
            snp_cpg     = ?
        WHERE snp_id    = ?
    };
    my $snp_feature_sth = $dbh->prepare($snp_feature);

    for my $align_id (@align_ids) {

        #----------------------------#
        # for each alignment
        #----------------------------#
        $align_seq_sth->execute($align_id);
        my ( $chr_name, $align_length, $chr_start, $chr_end, $target_seq )
            = $align_seq_sth->fetchrow_array;
        print
            "prosess align $align_id ",
            "in $chr_name $chr_start - $chr_end\n";
        next UPDATE if $chr_name =~ /rand|un|contig|hap|scaf/i;

        $chr_name =~ s/chr0?//i;

        #----------------------------#
        # construct transforming array and hash
        #----------------------------#
        # align_position to chr_position transforming array
        # the first element will be ignored
        my @chr_pos;
        $chr_pos[$align_length] = $align_length;    # pre-allocating memory

        # chr_position to align_position transforming hash
        my %align_pos;

        my $indel_count = 0;
        for my $i ( 1 .. $align_length ) {
            my $current_base = substr( $target_seq, $i - 1, 1 );
            if ( $current_base eq '-' ) {
                $indel_count++;
            }
            my $chr_position = $chr_start + $i - 1 - $indel_count;
            $chr_pos[$i] = $chr_position;
            if ( !exists $align_pos{$chr_position} ) {
                $align_pos{$chr_position} = $i;
            }
        }

        #----------------------------#
        # process align
        #----------------------------#
        {
            my $align_chr_start   = $chr_pos[1];
            my $align_chr_end     = $chr_pos[$align_length];
            my $align_chr_runlist = "$align_chr_start-$align_chr_end";
            my $align_coding
                = feature_portion( $cds_set, $align_chr_runlist );

            $align_feature_sth->execute( $align_coding, undef, undef,
                $align_id, );

            $align_feature_sth->finish;
        }

        #----------------------------#
        # process each indels
        #----------------------------#
        $indel_query_sth->execute($align_id);
        while ( my @row3 = $indel_query_sth->fetchrow_array ) {
            my ( $indel_id, $indel_start, $indel_end ) = @row3;
            my $indel_chr_start = $chr_pos[$indel_start];
            my $indel_chr_end   = $chr_pos[$indel_end];

            {
                my $indel_chr_runlist = "$indel_chr_start-$indel_chr_end";
                my $indel_coding
                    = feature_portion( $cds_set, $indel_chr_runlist );

                $indel_feature_sth->execute( $indel_coding, undef, $indel_id,
                );
                $indel_feature_sth->finish;
            }

            #----------------------------#
            # process each isws
            #----------------------------#
            {
                $isw_query_sth->execute($indel_id);
                while ( my @row4 = $isw_query_sth->fetchrow_array ) {
                    my ( $isw_id, $isw_start, $isw_end ) = @row4;
                    my $isw_chr_start = $chr_pos[$isw_start];
                    my $isw_chr_end   = $chr_pos[$isw_end];

                    my $isw_chr_runlist = "$isw_chr_start-$isw_chr_end";
                    my $isw_coding
                        = feature_portion( $cds_set, $isw_chr_runlist );

                    $isw_feature_sth->execute( $isw_coding, undef, $isw_id, );
                }

                $isw_feature_sth->finish;
                $isw_query_sth->finish;
            }
        }
        $indel_query_sth->finish;

        #----------------------------#
        # process each snps
        #----------------------------#
        {
            $snp_query_sth->execute($align_id);
            while ( my @row = $snp_query_sth->fetchrow_array ) {
                my ( $snp_id, $snp_pos, $snp_mutant_to ) = @row;
                my $snp_chr_pos = $chr_pos[$snp_pos];

                # coding and repeats
                my $snp_coding = $cds_set->member($snp_chr_pos);

                # cpg
                my $snp_cpg = 0;

                my $left_base  = substr( $target_seq, $snp_pos - 2, 1 );
                my $right_base = substr( $target_seq, $snp_pos,     1 );

                # CpG to TpG, C to T transition
                # On the reverse strand, is CpG to CpA
                if ( $snp_mutant_to eq "C->T" ) {    # original base is C
                    if ( $right_base eq "G" ) {
                        $snp_cpg = 1;
                    }
                }
                elsif ( $snp_mutant_to eq "G->A" ) {    # original base is G
                    if ( $left_base eq "C" ) {
                        $snp_cpg = 1;
                    }
                }

                $snp_feature_sth->execute( $snp_coding, undef, $snp_cpg,
                    $snp_id, );
            }

            $snp_feature_sth->finish;
            $snp_query_sth->finish;
        }
    }

    return;
};

my $worker_isw_cpg = sub {
    print "Processing isw_cpg_pi\n";

    # create alignDB object for this scope
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # select all snps in this alignment
    my $isw_query = q{
        SELECT  i.isw_id id,
                COUNT(*) /i.isw_length * 1.0 cpg
        FROM isw i, snp s
        WHERE i.isw_id = s.isw_id
        AND s.snp_cpg = 1
        GROUP BY i.isw_id
    };
    my $isw_sth = $dbh->prepare($isw_query);

    # update isw table in the new feature column
    my $isw_extra = q{
        UPDATE isw
        SET isw_cpg_pi = ?
        WHERE isw_id = ?
    };
    my $isw_extra_sth = $dbh->prepare($isw_extra);

    # for isw
    $isw_sth->execute;
    while ( my @row = $isw_sth->fetchrow_array ) {
        my ( $isw_id, $cpg ) = @row;
        $isw_extra_sth->execute( $cpg, $isw_id );
    }

    # update NULL value of isw_cpg_pi to 0
    my $isw_null = q{
        UPDATE isw
        SET isw_cpg_pi = 0
        WHERE isw_cpg_pi IS NULL
    };
    $obj->execute_sql($isw_null);
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$worker_isw_cpg->();

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}
exit;

sub feature_portion {
    my $feature_set = shift;
    my $pos_set     = shift;

    my $set;
    if ( ref $pos_set eq "AlignDB::IntSpan" ) {
        $set = $pos_set;
    }
    else {
        $set = AlignDB::IntSpan->new($pos_set);
    }

    my $pos_n = $set->cardinality;
    return if $pos_n <= 0;
    my $intersect       = $set->intersect($feature_set);
    my $n               = $intersect->cardinality;
    my $feature_portion = $n / $pos_n;

    return $feature_portion;

}

__END__

=head1 NAME

    update_multi_gff.pl - Add annotation info to alignDB with gff annotations

=head1 SYNOPSIS

    update_multi.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --ensembl           ensembl database name
        --parallel          run in parallel mode
        --batch             number of alignments process in one child process

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
