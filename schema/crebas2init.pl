#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Path::Tiny;
use FindBin;

my $content = path( $FindBin::RealBin, "crebas.sql" )->slurp;

$content =~ s/\r//g;
$content =~ s/type = InnoDB/ENGINE = MyISAM/g;
$content =~ s/^drop.*$//mg;
$content =~ s/alter.+?\;//sg;
$content =~ s/\n{3,}/\n\n/g;

# indexes
$content .= q{
#----------------------------#
# isw
#----------------------------#
CREATE INDEX indel_isw_id_FK ON isw
(
    isw_indel_id
);

};

path( $FindBin::RealBin, "init.sql" )->spew($content);
