#!/usr/bin/perl
use strict;
use warnings;

my $crebas_file = "crebas.sql";
my $init_file   = "init.sql";

open my $infh, '<', $crebas_file;
my $content = do { local $/; <$infh> };
close $infh;

$content =~ s/type = InnoDB/ENGINE = MyISAM/g;
$content =~ s/^drop.*$//mg;
$content =~ s/^alter.*?\;$//smg;
$content =~ s/\n{3,}/\n\n/g;

open my $out, '>', $init_file;
print $out $content;
close $out;
