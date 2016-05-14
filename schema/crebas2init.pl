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
# chromosome
#----------------------------#
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

#----------------------------#
# sequence
#----------------------------#
CREATE INDEX seq_chr_common_name ON sequence
(
    common_name, chr_name, chr_start, chr_end
);

CREATE INDEX seq_chr_name ON sequence
(
    chr_name, chr_start, chr_end
);

CREATE INDEX seq_chr_start ON sequence
(
    chr_start, chr_end
);

CREATE INDEX seq_role ON sequence
(
    seq_role
);

CREATE INDEX seq_position ON sequence
(
    seq_position
);

#----------------------------#
# indel
#----------------------------#
#CREATE INDEX indel_type ON indel
#(
#    indel_type
#);
#
#CREATE INDEX indel_align_id ON indel
#(
#    align_id, indel_start, indel_end
#);
#
#CREATE INDEX indel_start ON indel
#(
#    indel_start, indel_end
#);
#
#CREATE INDEX indel_freq ON indel
#(
#    indel_freq
#);

#----------------------------#
# snp
#----------------------------#
#CREATE INDEX snp_align_id ON snp
#(
#    align_id, snp_pos
#);
#
#CREATE INDEX snp_pos ON snp
#(
#    snp_pos
#);
#
#CREATE INDEX snp_freq ON snp
#(
#    snp_freq
#);

#----------------------------#
# isw
#----------------------------#
CREATE INDEX indel_isw_id_FK ON isw
(
    isw_indel_id
);

#CREATE INDEX isw_align_id ON isw
#(
#    align_id, isw_start, isw_end
#);
#
#CREATE INDEX isw_start ON isw
#(
#    isw_start, isw_end
#);
#
#CREATE INDEX isw_distance ON isw
#(
#    isw_distance
#);
#
#CREATE INDEX isw_density ON isw
#(
#    isw_density
#);
#
#CREATE INDEX isw_type ON isw
#(
#    isw_type
#);

#----------------------------#
# window
#----------------------------#
#CREATE INDEX window_align_id ON window
#(
#    align_id, window_start, window_end
#);
#
#CREATE INDEX window_start ON window
#(
#    window_start, window_end
#);

};

path( $FindBin::RealBin, "init.sql" )->spew($content);
