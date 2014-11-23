#!/usr/bin/perl
use strict;
use warnings;

use YAML::Syck qw(Dump Load DumpFile LoadFile);

use AlignDB::DeltaG;

my $dG = AlignDB::DeltaG->new;

DumpFile( "dG.yml", $dG);
