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

# indexes
$content .= q{
CREATE INDEX chr_common_name ON chromosome
(
   common_name, chr_name
);

CREATE INDEX taxon_id_idx ON chromosome
(
   taxon_id
);

CREATE INDEX chr_chr_name ON chromosome
(
    chr_name
);

CREATE INDEX seq_chr_id ON sequence
(
    chr_id, chr_start, chr_end
);

CREATE INDEX seq_chr_name ON sequence
(
    chr_name, chr_start, chr_end
);

CREATE INDEX seq_chr_start ON sequence
(
    chr_start, chr_end
);

CREATE INDEX indel_type ON indel
(
    indel_type
);

CREATE INDEX indel_isw_id_FK ON isw
(
    isw_indel_id
);

CREATE INDEX isw_distance ON isw
(
    isw_distance
);

CREATE INDEX isw_density ON isw
(
    isw_density
);

CREATE INDEX isw_type ON isw
(
    isw_type
);

};

path( $FindBin::RealBin, "init.sql" )->spew($content);
