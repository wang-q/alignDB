#!/usr/bin/perl
use strict;
use warnings;

use File::Find::Rule;
use File::Basename;
use File::Spec;

use FindBin;

use Template;
use YAML qw(Dump Load DumpFile LoadFile);

my $dir = shift || '.';

$dir = File::Spec->rel2abs($dir);
my $output = File::Spec->catfile($dir, "basicstat.yml");

my @xls_files = sort grep {/common/i}
    File::Find::Rule->file->maxdepth(1)->name( '*.xls', '*.xlsx' )->in($dir);

my @data;
for my $i ( 0 .. $#xls_files ) {
    my $file = $xls_files[$i];
    my $basename = basename( $file, '.common.xls', '.common.xlsx',
        '.common.chart.xls', '.common.chart.xlsx' );
    my $item = {
        name       => $basename,
        idx        => $i + 2,
        file => basename($file),
    };
    push @data, $item;
}

my $tt = Template->new;
$tt->process( \*DATA, { data => \@data, }, $output ) or die Template->error;
print "$output\n";

system("perl $FindBin::Bin/./excel_table.pl -i $output");

__DATA__
---
autofit: A:J
texts:
  - text: "Name"
    pos: A1
  - text: "Target length (Mb)"
    pos: B1
  - text: "Query length (Mb)"
    pos: C1
  - text: "Aligned length (Mb)"
    pos: D1
  - text: "Indels per 100 bp"
    pos: E1
  - text: "D on average"
    pos: F1
  - text: "GC-content"
    pos: G1
  - text: "File"
    pos: H1
[% FOREACH item IN data -%]
  - text: [% item.name %]
    pos: A[% item.idx %]
  - text: [% item.file %]
    pos: H[% item.idx %]
[% END -%]
borders:
  - range: A1:H1
    top: 1
  - range: A1:H1
    bottom: 1
ranges:
[% FOREACH item IN data -%]
  [% item.file %]:
    basic:
      - copy: B2
        paste: B[% item.idx %]
      - copy: B3
        paste: C[% item.idx %]
      - copy: B4
        paste: D[% item.idx %]
      - copy: B5
        paste: E[% item.idx %]
      - copy: B6
        paste: F[% item.idx %]
      - copy: B7
        paste: G[% item.idx %]
[% END -%]
