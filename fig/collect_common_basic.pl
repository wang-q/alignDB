#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;
use File::Spec;

use FindBin;

use Template;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# running options
my $dir = '.';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'  => \$help,
    'man'     => \$man,
    'd|dir=s' => \$dir,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$dir = File::Spec->rel2abs($dir);
my $output = File::Spec->catfile( $dir, "basicstat.yml" );

my @xls_files = sort grep {/common|multi/i}
    File::Find::Rule->file->maxdepth(1)->name('*.xlsx')->in($dir);

my @data;
for my $i ( 0 .. $#xls_files ) {
    my $file = $xls_files[$i];
    my $basename
        = basename( $file, '.common.xlsx', '.common.chart.xlsx', '.multi.xlsx',
        '.multi.chart.xlsx' );
    my $item = {
        name => $basename,
        file => basename($file),
    };
    push @data, $item;
}

my $tt   = Template->new;
my $text = <<'EOF';
---
autofit: A:J
texts:
  - text: "Name"
    pos: A1
  - text: "No. of strains"
    pos: B1
  - text: "Target length (Mb)"
    pos: C1
  - text: "Aligned length (Mb)"
    pos: D1
  - text: "Indels per 100 bp"
    pos: E1
  - text: "SNVs per 100 bp"
    pos: F1
  - text: "D on average"
    pos: G1
  - text: "GC-content"
    pos: H1
  - text: "avg-align"
    pos: I1
  - text: "File"
    pos: L1
[% FOREACH item IN data -%]
  - text: [% item.name %]
    pos: A[% loop.index + 2%]
  - text: [% item.file %]
    pos: L[% loop.index + 2%]
[% END -%]
borders:
  - range: A1:L1
    top: 1
  - range: A1:L1
    bottom: 1
ranges:
[% FOREACH item IN data -%]
  [% item.file %]:
    basic:
      - copy: B2
        paste: B[% loop.index + 2%]
      - copy: B3
        paste: C[% loop.index + 2%]
      - copy: B4
        paste: D[% loop.index + 2%]
      - copy: B5
        paste: E[% loop.index + 2%]
      - copy: B6
        paste: F[% loop.index + 2%]
      - copy: B7
        paste: G[% loop.index + 2%]
      - copy: B8
        paste: H[% loop.index + 2%]
    summary:
      - copy: B2
        paste: I[% loop.index + 2%]
[% END -%]
EOF

$tt->process( \$text, { data => \@data, }, $output ) or die Template->error;
print "$output\n";

system("perl $FindBin::Bin/./excel_table.pl -i $output");

__DATA__
