#!/usr/bin/perl
use strict;
use warnings;

my $raw_file = "gui3.raw.ui";
my $clean_file   = "gui3.ui";

open my $in_fh, '<', $raw_file;
open my $out_fh, '>', $clean_file;
while (my $line = <$in_fh>) {
    $line =~ /primary_icon_activatable/ and next;
    $line =~ /secondary_icon_activatable/ and next;
    $line =~ /primary_icon_sensitive/ and next;
    $line =~ /secondary_icon_sensitive/ and next;
    $line =~ /invisible_char/ and next;
    $line =~ s/ swapped="no"//;
    print {$out_fh} $line;
}
close $in_fh;
close $out_fh;

