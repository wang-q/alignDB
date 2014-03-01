#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Bio::Tools::GFF;

use File::Basename;

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
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# support multiply files, seperated by ','
my @gff_files;

# RepeatMasker generated gff files
my @rm_gff_files;

# gff version
my $gff_version = 3;

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'         => \$help,
    'man'            => \$man,
    's|server=s'     => \$server,
    'P|port=i'       => \$port,
    'u|username=s'   => \$username,
    'p|password=s'   => \$password,
    'd|db=s'         => \$db,
    'gff_files=s'    => \@gff_files,
    'rm_gff_files=s' => \@rm_gff_files,
    'gff_version=i'  => \$gff_version,
    'parallel=i'     => \$parallel,
    'batch=i'        => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
# other objects are initiated in subroutines
$stopwatch->start_message("Update annotations of $db...");

my $cds_set_of = {};
for my $file ( grep {defined} split /\,/, join(",", @gff_files) ) {
    next unless -e $file;
    my $basename = basename( $file, '.gff', '.gff3' );
    print "Loading annotations for [$basename]\n";
    my $gff_obj = Bio::Tools::GFF->new(
        -file        => $file,
        -gff_version => $gff_version,
    );

    my $cds_set = AlignDB::IntSpan->new;
    while ( my $feature = $gff_obj->next_feature ) {
        if ( $feature->primary_tag eq 'CDS' ) {
            $cds_set->add_range( $feature->start, $feature->end );
        }
    }
    $cds_set_of->{$basename} = $cds_set;
}

my $repeat_set_of = {};
for my $file ( grep {defined} split /\,/, join(",", @rm_gff_files) ) {
    next unless -e $file;
    my $basename = basename( $file, '.rm.gff', '.rm.gff3', '.gff', '.gff3' );
    print "Loading RepeatMasker annotations for [$basename]\n";
    my $gff_obj = Bio::Tools::GFF->new(
        -file        => $file,
        -gff_version => $gff_version,
    );

    my $repeat_set = AlignDB::IntSpan->new;
    while ( my $feature = $gff_obj->next_feature ) {
        $repeat_set->add_range( $feature->start, $feature->end );
    }
    $repeat_set_of->{$basename} = $repeat_set;
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
    # update align table
    my $align_feature = q{
        UPDATE align
        SET align_coding = ?,
            align_repeats = ?,
            align_te = ?,
            align_coding_runlist = ?,
            align_repeats_runlist = ?,
            align_te_runlist = ?
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

            # feature protions and runlists
            my ($align_coding,      $align_repeats,
                $align_cds_runlist, $align_repeat_runlist
            );

            if ( $cds_set_of->{$chr_name} ) {
                $align_coding = feature_portion( $cds_set_of->{$chr_name},
                    $align_chr_runlist );
                $align_cds_runlist
                    = $cds_set_of->{$chr_name}
                    ->map_set( sub { $align_pos{$_} } )->runlist;
            }

            if ( $repeat_set_of->{$chr_name} ) {
                $align_repeats = feature_portion( $repeat_set_of->{$chr_name},
                    $align_chr_runlist );
                $align_repeat_runlist
                    = $repeat_set_of->{$chr_name}
                    ->map_set( sub { $align_pos{$_} } )->runlist;
            }

            $align_feature_sth->execute( $align_coding, $align_repeats, undef,
                $align_cds_runlist, $align_repeat_runlist, undef, $align_id, );

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

                # feature protions
                my ( $indel_coding, $indel_repeats );
                if ( $cds_set_of->{$chr_name} ) {
                    $indel_coding = feature_portion( $cds_set_of->{$chr_name},
                        $indel_chr_runlist );
                }
                if ( $repeat_set_of->{$chr_name} ) {
                    $indel_repeats
                        = feature_portion( $repeat_set_of->{$chr_name},
                        $indel_chr_runlist );
                }

                $indel_feature_sth->execute( $indel_coding, $indel_repeats,
                    $indel_id, );
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

                    # feature protions
                    my ( $isw_coding, $isw_repeats );
                    if ( $cds_set_of->{$chr_name} ) {
                        $isw_coding = feature_portion( $cds_set_of->{$chr_name},
                            $isw_chr_runlist );
                    }
                    if ( $repeat_set_of->{$chr_name} ) {
                        $isw_repeats
                            = feature_portion( $repeat_set_of->{$chr_name},
                            $isw_chr_runlist );
                    }

                    $isw_feature_sth->execute( $isw_coding, $isw_repeats,
                        $isw_id, );
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
                my ( $snp_id, $snp_pos ) = @row;
                my $snp_chr_pos = $chr_pos[$snp_pos];

                # feature protions
                my ( $snp_coding, $snp_repeats );
                if ( $cds_set_of->{$chr_name} ) {
                    $snp_coding
                        = $cds_set_of->{$chr_name}->member($snp_chr_pos);
                }
                if ( $repeat_set_of->{$chr_name} ) {
                    $snp_repeats
                        = feature_portion( $repeat_set_of->{$chr_name},
                        $snp_chr_pos );
                }

                $snp_feature_sth->execute( $snp_coding, $snp_repeats, $snp_id,
                );
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

                # feature protions
                my ( $window_coding, $window_repeats );
                if ( $cds_set_of->{$chr_name} ) {
                    $window_coding = feature_portion( $cds_set_of->{$chr_name},
                        $window_chr_set );
                }
                if ( $repeat_set_of->{$chr_name} ) {
                    $window_repeats
                        = feature_portion( $repeat_set_of->{$chr_name},
                        $window_chr_set );
                }

                $window_update_sth->execute( $window_coding, $window_repeats,
                    $window_id, );
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

    update_feature_gff.pl - Add annotation info to alignDB with gff annotations

=head1 SYNOPSIS

    update_feature_gff.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --gff               gff files, support multiply file, seperated by ,
        --parallel          run in parallel mode
        --batch             number of alignments process in one child process

=cut
