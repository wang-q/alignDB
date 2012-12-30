#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Template;
use YAML qw(Dump Load DumpFile LoadFile);

my $orders   = [qw{ 50 40 20 30 10 9 8 7 6 5 4 3 2 }];
my $level_of = {
    50 => '5k ',
    40 => '4k ',
    20 => '3k ',
    30 => '2k ',
    10 => '1k ',
    9  => '900',
    8  => '800',
    7  => '700',
    6  => '600',
    5  => '500',
    4  => '400',
    3  => '300',
    2  => '200',
};

my $tt = Template->new;
my $text;

$text = <<'EOF';
---
charts:
  HumanvsCGOR_alt_level.gc.chart.xlsx:
[% FOREACH i IN orders -%]
    segment_gc_indel_[% i %]:
      1:
        - 0
        - [% loop.index %]
      2:
        - 1
        - [% loop.index %]
    segment_cv_indel_[% i %]:
      1:
        - 2
        - [% loop.index %]
      2:
        - 3
        - [% loop.index %]
[% END -%]
texts:
[% FOREACH i IN orders -%]
  - text: "[% lables.shift %]"
    size: 12
    bold: 1
    pos:
      - 0
      - [% loop.index - 0.05 %]
[% END -%]
[% FOREACH i IN orders -%]
  - text: "[% level_of.$i %] Human, chimp, gorilla and orangutan"
    size: 10
    bold: 1
    pos:
      - 0.1
      - [% loop.index - 0.05 %]
[% END -%]

EOF
$tt->process( \$text,
    { level_of => $level_of, lables => [ 'a' .. 'z' ], orders => $orders, },
    "Fig_gc_alt_level.yml" )
    or die Template->error;

__END__
