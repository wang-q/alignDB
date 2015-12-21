#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Spec;

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

write_runlist_feature.pl - extract runlists of a certain feature from alignDB

=head1 SYNOPSIS

    perl write_runlist_feature.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --ensembl   -e  STR     ensembl database name
        --length    -l  INT     threshold of alignment length
        --feature       STR     feature name, default is [non_repeat]
        --inset         INT     removing $inset bases from each end of each span of $set.
                                If $inset is negative, then -$inset integers are added to each end of each span.
        --invert                write inverted sets
        --parallel      INT     run in parallel mode

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'    => \( my $server           = $Config->{database}{server} ),
    'port|P=i'      => \( my $port             = $Config->{database}{port} ),
    'db|d=s'        => \( my $db               = $Config->{database}{db} ),
    'username|u=s'  => \( my $username         = $Config->{database}{username} ),
    'password|p=s'  => \( my $password         = $Config->{database}{password} ),
    'ensembl|e'     => \( my $ensembl_db       = $Config->{database}{ensembl} ),
    'length|lt|l=i' => \( my $length_threshold = $Config->{generate}{length_threshold} ),
    'feature=s'  => \( my $feature  = 'non_repeat' ),
    'inset=i'    => \my $inset,
    'invert'     => \my $invert,
    'parallel=i' => \( my $parallel = $Config->{generate}{parallel} ),
) or HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Write slice files from $db...");

# output dir
my $dir = "${db}_${feature}";
$dir .= "_inset$inset" if $inset;
$dir .= "_invert"      if $invert;
mkdir $dir, 0777 unless -e $dir;

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

    # select all target chromosomes in this database
    my @chrs = @{ $obj->get_chrs('target') };
    @jobs = @chrs;
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job = shift;
    my $chr = $job;

    local $| = 1;

    #----------------------------#
    # Init objects
    #----------------------------#
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # ensembl handler
    my $ensembl = AlignDB::Ensembl->new(
        server => $server,
        db     => $ensembl_db,
        user   => $username,
        passwd => $password,
    );

    my ( $chr_id, $chr_name, $chr_length ) = @{$chr};
    print "id => $chr_id, name => $chr_name, length => $chr_length\n";

    # for each align
    my @align_ids     = @{ $obj->get_align_ids_of_chr($chr_id) };
    my $chr_ftr_set   = AlignDB::IntSpan->new;
    my $chr_align_set = AlignDB::IntSpan->new;
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);

        # target
        my $target_info      = $obj->get_target_info($align_id);
        my $target_chr_name  = $target_info->{chr_name};
        my $target_chr_start = $target_info->{chr_start};
        my $target_chr_end   = $target_info->{chr_end};

        $chr_align_set->add("$target_chr_start-$target_chr_end");

        # make a new ensembl slice object
        my $ensembl_chr_name = $target_chr_name;
        $ensembl_chr_name =~ s/chr0?//i;

        #print "ensembl_chr_name $ensembl_chr_name\n";
        eval { $ensembl->set_slice( $ensembl_chr_name, $target_chr_start, $target_chr_end ); };
        if ($@) {
            warn "Can't get annotation\n";
            next;
        }

        my $slice       = $ensembl->slice;
        my $ftr_chr_set = $slice->{"_$feature\_set"};

        next unless $ftr_chr_set;
        next if $ftr_chr_set->is_empty;
        if ($inset) {
            $ftr_chr_set->inset($inset);
        }

        for my $set ( $ftr_chr_set->sets ) {
            next if $set->size < $length_threshold;
            $chr_ftr_set->add($set);
        }
    }

    # $chr_ftr_set should be subset of $chr_align_set
    $chr_ftr_set = $chr_ftr_set->intersect($chr_align_set);

    # the inverted set
    if ($invert) {
        $chr_ftr_set = $chr_align_set->diff($chr_ftr_set);
    }

    my $filename = File::Spec->catfile( $dir, "$chr_name.$feature.yml" );
    if ( $chr_ftr_set->is_not_empty ) {
        DumpFile( $filename, $chr_ftr_set->runlist );
        print "Finish merge $feature of $chr_name\n";
    }
    else {
        print "Finish $chr_name has no $feature\n";
    }
};

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$stopwatch->end_message;

__END__
