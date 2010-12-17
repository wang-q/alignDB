#!/usr/bin/perl
use strict;
use warnings;

package MyApp;
use Wx;
use base 'Wx::App';

use FindBin;
use lib "$FindBin::Bin/../lib";
require AlignDB::GUI::Frame;

sub OnInit {
    my $self = shift;

    my $frame = AlignDB::GUI::Frame->new;
    $frame->Show(1);
    $self->SetTopWindow($frame);

    return 1;
}

package main;

my $app = MyApp->new;
$app->MainLoop;

