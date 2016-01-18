#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB;
use AlignDB::Ensembl;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

update_feature.pl - Add annotation info to alignDB

=head1 SYNOPSIS

    perl update_feature.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --ensembl   -e  STR     ensembl database name
        --parallel      INT     run in parallel mode
        --batch         INT     number of alignments in one child process

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server       = $Config->{database}{server} ),
    'port|P=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'       => \( my $db           = $Config->{database}{db} ),
    'username|u=s' => \( my $username     = $Config->{database}{username} ),
    'password|p=s' => \( my $password     = $Config->{database}{password} ),
    'ensembl|e=s'  => \( my $ensembl_db   = $Config->{database}{ensembl} ),
    'parallel=i'   => \( my $parallel     = $Config->{generate}{parallel} ),
    'batch=i'      => \( my $batch_number = $Config->{generate}{batch} ),
) or HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
# other objects are initiated in subroutines
$stopwatch->start_message("Update annotations of $db...");

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

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # ensembl handler
    my $ensembl = AlignDB::Ensembl->new(
        server => $server,
        db     => $ensembl_db,
        user   => $username,
        passwd => $password,
    );

    #----------------------------#
    # SQL query and DBI sths
    #----------------------------#
    # update align table
    my $align_feature = q{
        UPDATE align
        SET align_coding = ?,
            align_repeats = ?,
            align_coding_runlist = ?,
            align_repeats_runlist = ?
        WHERE align_id = ?
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
        SELECT snp_id, snp_pos
        FROM snp
        WHERE align_id = ?
    };
    my $snp_query_sth = $dbh->prepare($snp_query);

    # update snp table in the new feature column
    my $snp_feature = q{
        UPDATE snp
        SET snp_coding  = ?,
            snp_repeats = ?
        WHERE snp_id    = ?
    };
    my $snp_feature_sth = $dbh->prepare($snp_feature);

    # select all windows for this alignment
    my $window_query = q{
        SELECT window_id, window_runlist
        FROM window
        WHERE align_id = ?
    };
    my $window_query_sth = $dbh->prepare($window_query);

    # update window table in the new feature column
    my $window_update_query = q{
        UPDATE window
        SET window_coding = ?, window_repeats = ?
        WHERE window_id = ?
    };
    my $window_update_sth = $dbh->prepare($window_update_query);

UPDATE: for my $align_id (@align_ids) {

        #----------------------------#
        # for each alignment
        #----------------------------#
        my $target_info  = $obj->get_target_info($align_id);
        my $chr_name     = $target_info->{chr_name};
        my $chr_start    = $target_info->{chr_start};
        my $chr_end      = $target_info->{chr_end};
        my $align_length = $target_info->{align_length};
        my ($target_seq) = @{ $obj->get_seqs($align_id) };
        $obj->process_message($align_id);

        next UPDATE if $chr_name =~ /rand|contig|hap|scaf|gi_/i;

        $chr_name =~ s/chr0?//i;

        # make a new ensembl slice object
        $ensembl->set_slice( $chr_name, $chr_start, $chr_end );

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

            # feature portions
            my $align_chr_start   = $chr_pos[1];
            my $align_chr_end     = $chr_pos[$align_length];
            my $align_chr_runlist = "$align_chr_start-$align_chr_end";

            # feature protions
            my $align_coding  = $ensembl->feature_portion( '_cds_set',    $align_chr_runlist );
            my $align_repeats = $ensembl->feature_portion( '_repeat_set', $align_chr_runlist );

            # feature runlists
            my $cds_set    = $ensembl->feature_set_obj('_cds_set');
            my $repeat_set = $ensembl->feature_set_obj('_repeat_set');
            $cds_set    = $cds_set->map_set( sub    { $align_pos{$_} } );
            $repeat_set = $repeat_set->map_set( sub { $align_pos{$_} } );
            $align_feature_sth->execute( $align_coding, $align_repeats,
                $cds_set->runlist, $repeat_set->runlist, $align_id, );

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
                my $indel_coding = $ensembl->feature_portion( '_cds_set', $indel_chr_runlist );
                my $indel_repeats = $ensembl->feature_portion( '_repeat_set', $indel_chr_runlist );
                $indel_feature_sth->execute( $indel_coding, $indel_repeats, $indel_id, );
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
                    my $isw_coding = $ensembl->feature_portion( '_cds_set', $isw_chr_runlist );
                    my $isw_repeats = $ensembl->feature_portion( '_repeat_set', $isw_chr_runlist );
                    $isw_feature_sth->execute( $isw_coding, $isw_repeats, $isw_id, );
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
            my $cds_set    = $ensembl->feature_set_obj('_cds_set');
            my $repeat_set = $ensembl->feature_set_obj('_repeat_set');
            while ( my @row = $snp_query_sth->fetchrow_array ) {
                my ( $snp_id, $snp_pos ) = @row;
                my $snp_chr_pos = $chr_pos[$snp_pos];

                my $snp_coding  = $cds_set->member($snp_chr_pos);
                my $snp_repeats = $repeat_set->member($snp_chr_pos);
                $snp_feature_sth->execute( $snp_coding, $snp_repeats, $snp_id, );
            }
            $snp_feature_sth->finish;
            $snp_query_sth->finish;
        }

        #----------------------------#
        # process each windows
        #----------------------------#
        {
            $window_query_sth->execute($align_id);
            while ( my @row5 = $window_query_sth->fetchrow_array ) {
                my ( $window_id, $window_runlist ) = @row5;
                my $window_set = AlignDB::IntSpan->new($window_runlist);
                my $window_chr_set
                    = $window_set->map_set( sub { $chr_pos[$_] } );

                my $window_coding  = $ensembl->feature_portion( '_cds_set',    $window_chr_set );
                my $window_repeats = $ensembl->feature_portion( '_repeat_set', $window_chr_set );
                $window_update_sth->execute( $window_coding, $window_repeats, $window_id, );
            }

            $window_update_sth->finish;
            $window_query_sth->finish;
        }
    }

    return;
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

__END__
