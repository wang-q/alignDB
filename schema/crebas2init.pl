#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Path::Tiny;
use FindBin;

my $init_file = "init.sql";

my $content = path( $FindBin::RealBin, "crebas.sql" )->slurp;

$content =~ s/\r//g;
$content =~ s/type = InnoDB/ENGINE = MyISAM/g;
$content =~ s/^drop.*$//mg;
$content =~ s/alter.+?\;//sg;
$content =~ s/\n{3,}/\n\n/g;

$content .= q{
create index common_name_idx on chromosome
(
   common_name
);

create index taxon_id_idx on chromosome
(
   taxon_id
);

};

path( $FindBin::RealBin, "init.sql" )->spew($content);
