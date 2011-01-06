#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Ensembl;

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
my $server     = $Config->{database}{server};
my $port       = $Config->{database}{port};
my $username   = $Config->{database}{username};
my $password   = $Config->{database}{password};
my $db         = $Config->{database}{db};
my $ensembl_db = $Config->{database}{ensembl};

my $process_align  = $Config->{feature}{align};
my $process_indel  = $Config->{feature}{indel};
my $process_isw    = $Config->{feature}{isw};
my $process_snp    = $Config->{feature}{snp};
my $process_window = $Config->{feature}{window};

my $align_id_runlist;    # only update align_ids in this set

# run in parallel mode
my $parallel = $Config->{feature}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{feature}{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'           => \$help,
    'man'              => \$man,
    'server=s'         => \$server,
    'port=i'           => \$port,
    'db=s'             => \$db,
    'username=s'       => \$username,
    'password=s'       => \$password,
    'ensembl=s'        => \$ensembl_db,
    'process_align=s'  => \$process_align,
    'process_indel=s'  => \$process_indel,
    'process_isw=s'    => \$process_isw,
    'process_snp=s'    => \$process_snp,
    'process_window=s' => \$process_window,
    'runlist=s'        => \$align_id_runlist,
    'parallel=i'       => \$parallel,
    'batch=i'          => \$batch_number,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
# other objects are initiated in subroutines
$stopwatch->start_message("Update annotations of $db...");

#----------------------------#
# Clean previous data and find all align_ids
#----------------------------#
my @align_ids;
if ($align_id_runlist) {
    @align_ids = AlignDB::IntSpan->new($align_id_runlist)->elements;
}
else {
    print "Clean previous data\n";

    # create alignDB object for this scope
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    if ($process_align) {
        $obj->empty_table('align_extra');
    }

    if ($process_indel) {
        $obj->empty_table('indel_extra');
    }

    if ($process_isw) {
        $obj->empty_table('isw_extra');
    }

    if ($process_snp) {
        $obj->empty_table('snp_extra');
    }

    @align_ids = @{ $obj->get_align_ids };
}

my @jobs;
while ( scalar @align_ids ) {
    my @batching = splice @align_ids, 0, $batch_number;
    push @jobs, [@batching];
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
    # update align_extra table in the new feature column
    my $align_extra = q{
        INSERT INTO align_extra (
            align_extra_id, align_id, 
            align_feature1, align_feature2, align_feature3,
            align_feature5, align_feature6, align_feature7
        )
        VALUES (
            NULL, ?,
            ?, ?, ?,
            ?, ?, ?
        )
    };
    my $align_extra_sth = $dbh->prepare($align_extra);

    # select all indels in this alignment
    my $indel_query = q{
        SELECT indel_id, prev_indel_id, indel_start, indel_end
        FROM indel
        WHERE align_id = ?
    };
    my $indel_query_sth = $dbh->prepare($indel_query);

    # update indel table in the new feature column
    my $indel_extra = q{
        INSERT INTO indel_extra (
            indel_extra_id, indel_id, prev_indel_id,
            indel_feature1, indel_feature2
        )
        VALUES (NULL, ?, ?, ?, ?)
    };
    my $indel_extra_sth = $dbh->prepare($indel_extra);

    # select all isws for this indel
    my $isw_query = q{
        SELECT isw_id, isw_start, isw_end
        FROM isw
        WHERE indel_id = ?
    };
    my $isw_query_sth = $dbh->prepare($isw_query);

    # update isw table in the new feature column
    my $isw_extra = q{
        INSERT INTO isw_extra (
            isw_extra_id, isw_id, isw_feature1, isw_feature2
        )
        VALUES (NULL, ?, ?, ?)
    };
    my $isw_extra_sth = $dbh->prepare($isw_extra);

    # select all snps for this indel
    my $snp_query = q{
        SELECT snp_id, snp_pos
        FROM snp
        WHERE align_id = ?
    };
    my $snp_query_sth = $dbh->prepare($snp_query);

    # update snp table in the new feature column
    my $snp_extra = q{
        INSERT INTO snp_extra   (
            snp_extra_id, snp_id, snp_feature1, snp_feature2
        )
        VALUES (NULL, ?, ?, ?)
    };
    my $snp_extra_sth = $dbh->prepare($snp_extra);

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
        
        next UPDATE if $chr_name =~ /rand|un|contig|hap|scaf|gi_/i;

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
        if ($process_align) {

            # feature portions
            my $align_chr_start   = $chr_pos[1];
            my $align_chr_end     = $chr_pos[$align_length];
            my $align_chr_runlist = "$align_chr_start-$align_chr_end";
            my $align_feature1
                = $ensembl->feature_portion( '_cds_set', $align_chr_runlist );
            my $align_feature2 = $ensembl->feature_portion( '_repeat_set',
                $align_chr_runlist );
            my $align_feature3
                = $ensembl->feature_portion( '_te_set', $align_chr_runlist );

            # feature runlists
            my $cds_set    = $ensembl->feature_set_obj('_cds_set');
            my $repeat_set = $ensembl->feature_set_obj('_repeat_set');
            my $te_set     = $ensembl->feature_set_obj('_te_set');
            $cds_set    = $cds_set->map_set( sub    { $align_pos{$_} } );
            $repeat_set = $repeat_set->map_set( sub { $align_pos{$_} } );
            $te_set     = $te_set->map_set( sub     { $align_pos{$_} } );
            $align_extra_sth->execute(
                $align_id,         $align_feature1,
                $align_feature2,   $align_feature3,
                $cds_set->runlist, $repeat_set->runlist,
                $te_set->runlist,
            );

            $align_extra_sth->finish;
        }

        #----------------------------#
        # process each indels
        #----------------------------#
        $indel_query_sth->execute($align_id);
        while ( my @row3 = $indel_query_sth->fetchrow_array ) {
            my ( $indel_id, $prev_indel_id, $indel_start, $indel_end )
                = @row3;
            my $indel_chr_start = $chr_pos[$indel_start];
            my $indel_chr_end   = $chr_pos[$indel_end];

            if ($process_indel) {

               #my ($indel_feature1, $indel_feature2) =
               #  $ensembl->locate_position($indel_chr_start, $indel_chr_end);
                my $indel_chr_runlist = "$indel_chr_start-$indel_chr_end";
                my $indel_feature1    = $ensembl->feature_portion( '_cds_set',
                    $indel_chr_runlist );
                my $indel_feature2 = $ensembl->feature_portion( '_repeat_set',
                    $indel_chr_runlist );
                $indel_extra_sth->execute(
                    $indel_id,       $prev_indel_id,
                    $indel_feature1, $indel_feature2
                );
                $indel_extra_sth->finish;
            }

            #----------------------------#
            # process each isws
            #----------------------------#
            if ($process_isw) {
                $isw_query_sth->execute($indel_id);
                while ( my @row4 = $isw_query_sth->fetchrow_array ) {
                    my ( $isw_id, $isw_start, $isw_end ) = @row4;
                    my $isw_chr_start = $chr_pos[$isw_start];
                    my $isw_chr_end   = $chr_pos[$isw_end];

                    my $isw_chr_runlist = "$isw_chr_start-$isw_chr_end";
                    my $isw_feature1 = $ensembl->feature_portion( '_cds_set',
                        $isw_chr_runlist );
                    my $isw_feature2
                        = $ensembl->feature_portion( '_repeat_set',
                        $isw_chr_runlist );
                    $isw_extra_sth->execute( $isw_id, $isw_feature1,
                        $isw_feature2 );
                }

                $isw_extra_sth->finish;
                $isw_query_sth->finish;
            }
        }
        $indel_query_sth->finish;

        #----------------------------#
        # process each snps
        #----------------------------#
        if ($process_snp) {
            $snp_query_sth->execute($align_id);
            my $cds_set    = $ensembl->feature_set_obj('_cds_set');
            my $repeat_set = $ensembl->feature_set_obj('_repeat_set');
            while ( my @row = $snp_query_sth->fetchrow_array ) {
                my ( $snp_id, $snp_pos ) = @row;
                my $snp_chr_pos = $chr_pos[$snp_pos];

                my $snp_feature1 = $cds_set->member($snp_chr_pos);
                my $snp_feature2 = $repeat_set->member($snp_chr_pos);
                $snp_extra_sth->execute( $snp_id, $snp_feature1,
                    $snp_feature2 );
            }

            $snp_extra_sth->finish;
            $snp_query_sth->finish;
        }

        #----------------------------#
        # process each windows
        #----------------------------#
        if ($process_window) {
            $window_query_sth->execute($align_id);
            while ( my @row5 = $window_query_sth->fetchrow_array ) {
                my ( $window_id, $window_runlist ) = @row5;
                my $window_set = AlignDB::IntSpan->new($window_runlist);
                my $window_chr_set
                    = $window_set->map_set( sub { $chr_pos[$_] } );

                #my ($window_coding, $window_repeats) =
                #  $ensembl->locate_set_position($window_chr_set);
                my $window_coding = $ensembl->feature_portion( '_cds_set',
                    $window_chr_set );
                my $window_repeats = $ensembl->feature_portion( '_repeat_set',
                    $window_chr_set );
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

__END__

=head1 NAME

    update_feature.pl - Add annotation info to alignDB

=head1 SYNOPSIS

    update_feature.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --ensembl           ensembl database name
        --process_align   
        --process_indel  
        --process_isw
        --process_snp
        --process_window
        --runlist           only update align_ids in this set
                            use this opt to continue previous work
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
