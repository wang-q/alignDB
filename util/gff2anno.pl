#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use IO::Zlib;
use Path::Tiny;
use AlignDB::IntSpan;
use App::RL::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

gff2anno.pl - Convert gff3 file to chromosome runlists

=head1 SYNOPSIS

    perl gff2anno.pl <gff files> [options]
      Options:
        --help          -?          brief help message
        --type          -t  STR     primary tag (the third field)
        --remove                    remove chr0 and .\d+ from region name

=cut

GetOptions(
    'help|?'   => sub { Getopt::Long::HelpMessage(0) },
    'type|t=s' => \my $type,
    'remove|r' => \my $remove,
) or Getopt::Long::HelpMessage(1);

my $set_of = {};
for my $infile (@ARGV) {
    my $in_fh = IO::Zlib->new( $infile, "rb" );

    while (1) {
        my $line = <$in_fh>;
        last unless $line;
        next if $line =~ /^#/;

        chomp $line;
        my @array = split( "\t", $line );
        my $feature_type = $array[2];

        if ( defined $type ) {
            next if $type ne $feature_type;
        }

        my $chr_name  = $array[0];
        my $chr_start = $array[3];
        my $chr_end   = $array[4];

        if ($remove) {
            $chr_name =~ s/chr0?//i;
            $chr_name =~ s/\.\d+$//;
        }

        if ( !exists $set_of->{$chr_name} ) {
            $set_of->{$chr_name} = AlignDB::IntSpan->new;
        }

        $set_of->{$chr_name}->add_pair( $chr_start, $chr_end );
    }
}

for my $chr ( keys %{$set_of} ) {
    $set_of->{$chr} = $set_of->{$chr}->runlist;
}

print YAML::Syck::Dump($set_of);

exit;

__END__
