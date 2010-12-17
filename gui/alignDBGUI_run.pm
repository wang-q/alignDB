package alignDBGUI;
use strict;
use warnings;

require Exporter;
our @ISA       = qw(Exporter);
our @EXPORT_OK = qw( run about);

use Tk::MListbox;
use Parse::CSV;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin qw($RealBin);
use lib $RealBin;

use vars qw(
  $about_dialog
  $load_dialog
  $version
  $top
  $pixmap_main
  $pixmap_open
  $pixmap_info
  $pixmap_auto
  $pixmap_load
  $_button_auto_db_name
  $_button_auto_stat_name
  $_button_auto_graph_name
  $_button_open_axt_dir
  $_button_about
  $_button_load_target
  $_button_load_query
);

sub run {
    &load_icon_data;

    # set main window options
    $top->Icon(-image => $pixmap_main);
    $top->resizable(0, 0);
    $top->deiconify();
    $top->raise();
    $top->title("alignDB GUI");

    ## set toplevel window position
    #$top->geometry("+100+100");

    # init button images
    $_button_auto_db_name->configure(-image    => $pixmap_auto);
    $_button_auto_stat_name->configure(-image  => $pixmap_auto);
    $_button_auto_graph_name->configure(-image => $pixmap_auto);
    $_button_open_axt_dir->configure(-image    => $pixmap_open);
    $_button_about->configure(-image           => $pixmap_info);
    $_button_load_target->configure(-image     => $pixmap_load);
    $_button_load_query->configure(-image      => $pixmap_load);
}

#----------------------------------------------------------#
# Exec cmd and output a standard report
#----------------------------------------------------------#
sub exec_cmd {
    my $cmd = shift;
    print "\n", "=" x 12, "CMD", "=" x 15, "\n";
    print $cmd , "\n";
    print "=" x 30, "\n";

    # run $cmd in background
    if ($^O eq 'MSWin32') {
        $cmd = "start /B " . $cmd;
    }
    else {
        $cmd = $cmd . " &";
    }

    system($cmd);
}

sub about {
    if (Exists($about_dialog)) {
        $about_dialog->deiconify();
        $about_dialog->raise();
        return;
    }

    my $text = <<EOT;
alignDB GUI v$version
Copyright (C) 2004-2006 Wang Qiang\n
EOT
    my $text2 = <<EOT;
A GUI shell for alignDB scripts.\n
This program is free software;
you can redistribute it and/or modify it under the same terms as Perl itself.\n
EOT
    my $address = <<EOT;
Mail me: wang__qiang\@hotmail.com
EOT

    $about_dialog = $top->Toplevel(-title => 'About alignDB GUI ...');
    my $about_dialog_f =
      $about_dialog->Frame()
      ->pack(-padx => 4, -pady => 4, -expand => 1, -fill => 'both');
    my $about_icon =
      $about_dialog_f->Label(-image => $pixmap_main)
      ->pack(-side => 'left', -anchor => 'n');

    # Want the text widget to have the same background colour as the window
    # (so the user can't tell it's a text widget!)  We're using a text widget
    # rather than a frame to get some formatting options (i.e. bold text)
    my $colour     = $about_dialog->cget(-bg);
    my $about_text = $about_dialog_f->Text(
        -relief => 'flat',
        -wrap   => 'word',
        -bg     => $colour,
        -width  => 40,
        -height => 10,
        -font   => '{Courier New} 10',
    )->pack(-expand => 0, -fill => 'both', -anchor => 'nw');
    $about_text->tagConfigure('bold',   -font    => '{Courier New} 10 bold');
    $about_text->tagConfigure('center', -justify => 'center');
    $about_text->insert('end', $text, 'bold');
    $about_text->insert('end', $text2);
    $about_text->insert('end', $address);
    my $about_OK = $about_dialog->Button(
        -text    => 'OK',
        -padx    => 4,
        -padx    => 4,
        -width   => 6,
        -font    => 'Tahoma 9',
        -command => sub { $about_dialog->destroy },
    )->pack(-pady => 4, -side => 'bottom');
    $about_dialog->Icon(-image => $pixmap_main);

    # set $about_dialog window position
    my ($top_width, $top_height, $top_x, $top_y) = &get_position($top);
    my ($width, $height) = &get_position($about_dialog);
    my $x = int(($top_width - $width) / 2 + $top_x);
    my $y = int(($top_height - $height) / 2 + $top_y);
    $about_dialog->geometry("+$x+$y");
}

sub load {
    my ($taxon_id_ref, $name_ref) = @_;

    if (Exists($load_dialog)) {
        $load_dialog->deiconify();
        $load_dialog->raise();
        return;
    }

    $load_dialog = $top->Toplevel(-title => 'Load taxon ...');

    # use MListbox to display taxon info
    #
    my $taxon_Mlistbox = $load_dialog->Scrolled(
        "MListbox",
        -scrollbars       => 'oe',
        -borderwidth      => '3',
        -font             => '{Courier New} 9',
        -height           => '9',
        -relief           => 'ridge',
        -selectforeground => '#ffffff',
        -selectmode       => 'extended',
        -separatorcolor   => '#EEEEEE',
        -width            => '0',
    )->pack(-padx => 4, -pady => 4, -expand => 1, -fill => 'both');

    $taxon_Mlistbox->configure(
        -columns => [
            [ -text => 'taxon_id' ],
            [ -text => 'name' ],
            [ -text => 'genus' ],
            [ -text => 'species' ],
        ],
    );

    $taxon_Mlistbox->columnGet(0)->Subwidget('heading')
      ->configure(-font => '{Tahoma} 9',);
    $taxon_Mlistbox->columnGet(1)->Subwidget('heading')
      ->configure(-font => '{Tahoma} 9',);
    $taxon_Mlistbox->columnGet(2)->Subwidget('heading')
      ->configure(-font => '{Tahoma} 9',);
    $taxon_Mlistbox->columnGet(3)->Subwidget('heading')
      ->configure(-font => '{Tahoma} 9',);
    $taxon_Mlistbox->columnGet(0)->configure(-width => 8,);
    $taxon_Mlistbox->columnGet(1)->configure(-width => 8,);
    $taxon_Mlistbox->columnGet(2)->configure(-width => 15,);
    $taxon_Mlistbox->columnGet(3)->configure(-width => 15,);

    # OK button, set taxon_id and name,
    # close this window to affect nothing
    my $selecte_taxon_cmd = sub {

        # only get the first selectd elemenet
        my @selected = $taxon_Mlistbox->curselection;
        if (@selected) {
            my @taxon = $taxon_Mlistbox->getRow($selected[0]);
            if (@taxon) {
                $$taxon_id_ref = $taxon[0];
                $$name_ref     = $taxon[1];
            }
        }
        $load_dialog->destroy;
    };
    my $load_OK = $load_dialog->Button(
        -text    => 'OK',
        -padx    => 4,
        -padx    => 4,
        -width   => 6,
        -font    => 'Tahoma 9',
        -command => $selecte_taxon_cmd,
    )->pack(-pady => 4, -side => 'bottom');

    # set icon of this window
    $load_dialog->Icon(-image => $pixmap_main);
    
    # read out normal taxons from $taxon_file
    # and then put them into $taxon_Mlistbox
    my $taxon_file = "$FindBin::Bin/normal_taxon.csv";
    my $csv = Parse::CSV->new( file => $taxon_file );
    
    # First record contains list of columns
    my $fields = $csv->fetch;

    my @taxons;
    while (my $record = $csv->fetch) {
        last unless @$record;
        push @taxons, $record;
    }
    print STDERR Dump @taxons;

    $taxon_Mlistbox->delete(0, 'end');
    $taxon_Mlistbox->insert('end', @taxons);

    # set $load_dialog window position
    my ($top_width, $top_height, $top_x, $top_y) = &get_position($top);
    my ($width, $height) = &get_position($load_dialog);
    my $x = int(($top_width - $width) / 2 + $top_x);
    my $y = int(($top_height - $height) / 2 + $top_y);
    $load_dialog->geometry("+$x+$y");
}

#----------------------------------------------------------#
# Get window position
#----------------------------------------------------------#
sub get_position {
    my $window   = shift;
    my $geometry = $window->geometry;
    my ($width, $height, $x, $y) =
      ($geometry =~ /(\d.*)x(\d.*?)([+-]\d+)([+-]\d+)/);
    return ($width, $height, $x, $y);
}

#----------------------------------------------------------#
# All icons were borrowed from perlprimer. Great work!
#----------------------------------------------------------#
sub load_icon_data {
    my $icon_main = <<'end_of_pixmap';
/* XPM */ 
static char * ppl_xpm[] = {
"32 32 180 2",
"  	c white",
"! 	c black",
"# 	c #707070",
"$ 	c #888888",
"% 	c #848484",
"& 	c #808080",
"' 	c #8C8C8C",
"( 	c #949494",
") 	c #989898",
"* 	c #A8A8A8",
"+ 	c #B0B0B0",
", 	c #D0D0D0",
"- 	c #CCCCCC",
". 	c #C8C8C8",
"0 	c #D4D4D4",
"1 	c #DCDCDC",
"2 	c #E0E0E0",
"3 	c #4C5864",
"4 	c #B0B8BC",
"5 	c #E0DCDC",
"6 	c #E0E0DC",
"7 	c #E4E0E0",
"8 	c #E8E4E4",
"9 	c #E8E8E8",
": 	c #ECE8E8",
"; 	c #F0ECEC",
"< 	c #0C202C",
"= 	c #30404C",
"> 	c #98A4A8",
"? 	c #A0A8AC",
"@ 	c #ACB8C4",
"A 	c #B0BCC4",
"B 	c #B0C0C8",
"C 	c #B4C0C8",
"D 	c #B4C0CC",
"E 	c #B4C4CC",
"F 	c #ACBCC4",
"G 	c #748088",
"H 	c #14242C",
"I 	c #141818",
"J 	c #282C30",
"K 	c #44484C",
"L 	c #34383C",
"M 	c #8CA0B0",
"N 	c #B4C4D0",
"O 	c #B0C4D0",
"P 	c #B0C0D0",
"Q 	c #ACC0CC",
"R 	c #A8BCCC",
"S 	c #A4B8C8",
"T 	c #A0B8C8",
"U 	c #809CB0",
"V 	c #040408",
"W 	c #181818",
"X 	c #303030",
"Y 	c #545050",
"Z 	c #4C4C4C",
"[ 	c #E0E8EC",
"] 	c #E8E8F0",
"^ 	c #E4E8F0",
"_ 	c #DCE4EC",
"` 	c #D8E0E8",
"a 	c #D0DCE4",
"b 	c #CCD8E0",
"c 	c #C8D4DC",
"d 	c #C4D0DC",
"e 	c #C0CCD8",
"f 	c #345870",
"g 	c #444444",
"h 	c #D0DCE0",
"i 	c #D0D8E0",
"j 	c #D4DCE4",
"k 	c #C4D0D8",
"l 	c #B8C8D4",
"m 	c #A0B4C4",
"n 	c #A8B8C4",
"o 	c #C0D0D8",
"p 	c #BCCCD4",
"q 	c #A8B8C8",
"r 	c #98B0C0",
"s 	c #404040",
"t 	c #A4B4C0",
"u 	c #90A8BC",
"v 	c #A0B0B8",
"w 	c #A0B0C0",
"x 	c #88A0B0",
"y 	c #90A8B8",
"z 	c #305470",
"{ 	c #98A8B8",
"| 	c #8098B0",
"} 	c #B8C8D0",
"~ 	c #C0C8D0",
" !	c #B0BCC8",
"!!	c #7090A0",
"#!	c #9CB4C4",
"$!	c #F0F0F0",
"%!	c #F8F8F8",
"&!	c #D0D8D8",
"'!	c #7890A0",
"(!	c #88A8B8",
")!	c #88A0B8",
"*!	c #305070",
"+!	c #8098A8",
",!	c #7898A8",
"-!	c #80A0B0",
".!	c #E0E8E8",
"0!	c #285068",
"1!	c #383838",
"2!	c #6888A0",
"3!	c #CCD4DC",
"4!	c #C8D0D8",
"5!	c #608098",
"6!	c #7090A8",
"7!	c #2C2C2C",
"8!	c #484848",
"9!	c #7088A0",
":!	c #A8BCC8",
";!	c white",
"<!	c #6088A0",
"=!	c #000404",
">!	c #181814",
"?!	c #282828",
"@!	c #F0F8F8",
"A!	c #DCE4E8",
"B!	c #E8F0F0",
"C!	c #7890A8",
"D!	c #507898",
"E!	c #588098",
"F!	c #204868",
"G!	c #6890A0",
"H!	c #D8D8D8",
"I!	c #A8B8C0",
"J!	c #6080A0",
"K!	c #507890",
"L!	c #101010",
"M!	c #202020",
"N!	c #587898",
"O!	c #C8D0D0",
"P!	c #487090",
"Q!	c #407090",
"R!	c #204860",
"S!	c #000008",
"T!	c #587890",
"U!	c #D8E0E4",
"V!	c #F8F0F0",
"W!	c #386888",
"X!	c #184060",
"Y!	c #182030",
"Z!	c #507090",
"[!	c #406888",
"]!	c #306080",
"^!	c #104060",
"_!	c #487088",
"`!	c #285880",
"a!	c #103858",
"b!	c #080808",
"c!	c #406880",
"d!	c #C0C8C8",
"e!	c #205878",
"f!	c #386080",
"g!	c #285878",
"h!	c #205078",
"i!	c #083858",
"j!	c #306078",
"k!	c #185070",
"l!	c #104870",
"m!	c #205070",
"n!	c #185078",
"o!	c #184870",
"p!	c #104868",
"q!	c #083058",
"r!	c #001020",
"s!	c #084068",
"t!	c #003050",
"u!	c #081828",
"v!	c #104068",
"w!	c #081018",
"x!	c #001038",
"y!	c #082028",
"z!	c #101818",
"                                                                ",
"          # $ % % % % % % % & % $ $ $ $ $ ' ( ) * *             ",
"          + , , , - - - . - . . , 0 0 1 2 2 2 2 2 1 (           ",
"        3 4 5 1 1 2 5 1 1 2 1 6 5 1 2 7 7 8 9 : ; 9 )           ",
"        < = > ? ? ? @ A B B C C D C D E D D D B D F G           ",
"        H ! I J K L M N O N N N N P P Q Q R R S T R U           ",
"        H V W X Y Z D [ [ [ ] ^ [ _ ` ` a b c d e d Q f         ",
"        H V W X Z g @ h i j ` ` j i c k e l N P R Q m f         ",
"        H V W X Z g n b b h j j i b c o p l N Q q R r f         ",
"        H V W X Z s t c c b b i h b c d p O Q q m T u f         ",
"        H V W X Z s v o e o c p w x x x r Q R m r r y z         ",
"        H V W X Z s { l l p r | w } ~  !x !!y #!u u x z         ",
"        H V W X Z s x T R x r $!%!$!$!$!%!&!'!| (!)!U *!        ",
"        H V W X Z s +!u | { %!b y ,!,!-!x i .!'!,!U ,!0!        ",
"        H V W X Z 1!'!-!U ; p 2!r i c x p { 3!4!5!,!6!0!        ",
"        H V W 7!8!1!9!,!F 9 !!:!;!%!;!%!;!m x 7 !!<!<!0!        ",
"        H =!>!?!g 1!2!6!&!E ,!@!A!,!U B!;!U C!1 U D!E!F!        ",
"        H V >!?!s X 5!G!H!I!x ;!#!2!5!j $!J!C!2 C!K!K!F!        ",
"        H =!L!M!1!7!N!G!O!C x %!} 5!y %!3!<!~ O!D!P!Q!R!        ",
"        H S!L!M!1!?!T!J!w 2 | 4!%!U!9 V!9 H!~ 2!P!Q!W!X!        ",
"        Y!=!L!M!X ?!Z!E!6!, ~ U B k y n E y +!5!W![!]!^!        ",
"        Y!=!L!W 7!M!_!D!N!| 0 O!x 2!5!J!6!t 0 2!]!W!`!a!        ",
"        Y!=!b!W ?!M!c!P!P!K!C!A H!0 , 0 0 d!,!Q!]!]!e!a!        ",
"        Y!=!b!I M!W f!Q![![!P!E!6!U x +!!!N![!]!`!g!h!i!        ",
"        H =!b!L!W W j!W!W!W!W!W![!Q!Q![!W!]!`!g!e!h!k!i!        ",
"        Y!! b!L!W L!g!]!]!]!]!]!]!]!`!`!g!g!e!h!h!k!l!i!        ",
"        < =!V b!L!L!m!`!g!g!`!`!g!g!g!g!h!h!n!k!o!l!p!q!        ",
"        r!! =!b!L!b!o!h!h!h!h!h!h!h!h!h!k!k!o!l!p!s!s!t!        ",
"        u!=!! V b!V v!o!o!o!o!o!o!o!l!l!l!l!s!s!s!s!s!a!        ",
"          w!=!V b!=!v!n!k!k!n!n!k!k!o!l!l!l!l!l!l!l!l!x!        ",
"          y!z!z!z!w!a!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!^!          ",
"                                                                ",
};
end_of_pixmap

    my $icon_open = <<'end_of_pixmap';
/* XPM */
static char * open_20_xpm[] = {
"20 20 27 1",
" 	c None",
".	c #020501",
"+	c #37372E",
"@	c #5B5D5A",
"#	c #8E8F7E",
"$	c #BABCAC",
"%	c #B8BB98",
"&	c #ACB0A3",
"*	c #1A1B18",
"=	c #43433B",
"-	c #CDD0AD",
";	c #737561",
">	c #A1A68A",
",	c #898F73",
"'	c #D7D9C6",
")	c #525548",
"!	c #AEB293",
"~	c #C4C7A4",
"{	c #242521",
"]	c #61624E",
"^	c #808673",
"/	c #C2C5B4",
"(	c #959B7E",
"_	c #D5D8B6",
":	c #B1B39D",
"<	c #999B8E",
"[	c #6C6D59",
"                    ",
"                    ",
"  .....             ",
" +#<<<^=            ",
" @____~[{           ",
"=$_~~~-(+           ",
"=/%!!%%>;======+    ",
"+->(>>>>>>:::::#.   ",
"{/,)+++++++++++{..  ",
"{:]<$$$///$////$$$@ ",
"{$=&__---_-----~~!) ",
"{<<'_---------~%!,+ ",
"{<<'-~--~~~~%~%:>]* ",
"*&/~%%%%%%%!!>>>^+  ",
"*$/!!>>>>>>>((,,].  ",
"{&,^;;;;[[[]]]))+   ",
" ...............    ",
"                    ",
"                    ",
"                    "};
end_of_pixmap

    my $icon_info = <<'end_of_pixmap';
/* XPM */
static char * info_20_xpm[] = {
"20 20 147 2",
"  	c None",
". 	c #ABABB1",
"+ 	c #D6D6DE",
"@ 	c #E5E5EF",
"# 	c #E5E5F2",
"$ 	c #D8D8E6",
"% 	c #BBBBC9",
"& 	c #6F7078",
"* 	c #CECED4",
"= 	c #F4F4F9",
"- 	c #F8F8FB",
"; 	c #F0F0F7",
"> 	c #E9E9F4",
", 	c #E4E4F2",
"' 	c #DEDFF0",
") 	c #D7D8EA",
"! 	c #8B8B97",
"~ 	c #C4C4CA",
"{ 	c #F6F6FB",
"] 	c #FCFCFD",
"^ 	c #F9F9FC",
"/ 	c #F1F1F8",
"( 	c #ECECF6",
"_ 	c #E6E6F3",
": 	c #E0E1F1",
"< 	c #DADDEF",
"[ 	c #D6D9EB",
"} 	c #6C6D77",
"| 	c #EEEEF3",
"1 	c #FBFBFD",
"2 	c #F3F3F9",
"3 	c #EEEEF7",
"4 	c #E8E8F5",
"5 	c #E2E2F2",
"6 	c #D7DBED",
"7 	c #C5C8DA",
"8 	c #BFBFC4",
"9 	c #FAFAFC",
"0 	c #F6F6FA",
"a 	c #F2F3F9",
"b 	c #E7E8F5",
"c 	c #DCE0EF",
"d 	c #D6DCEC",
"e 	c #D5DBED",
"f 	c #55575E",
"g 	c #CBCBD2",
"h 	c #F2F2F9",
"i 	c #F7F7FB",
"j 	c #F9F9FB",
"k 	c #F7F8FB",
"l 	c #F4F4FA",
"m 	c #E0E4F1",
"n 	c #D6DFEC",
"o 	c #D4DCED",
"p 	c #6A6C77",
"q 	c #C6C6CD",
"r 	c #EFEFF7",
"s 	c #F6F6F8",
"t 	c #F9F9FA",
"u 	c #F6F7FA",
"v 	c #ECECF2",
"w 	c #E1E6F1",
"x 	c #D6E3EB",
"y 	c #D4DEEC",
"z 	c #676973",
"A 	c #A8A8AF",
"B 	c #EAEAF4",
"C 	c #EDEEF5",
"D 	c #F1F1F7",
"E 	c #E7E7E9",
"F 	c #F2F2F8",
"G 	c #EBECF3",
"H 	c #E3E6ED",
"I 	c #DCE5EC",
"J 	c #D5E8EA",
"K 	c #D3E0EC",
"L 	c #404248",
"M 	c #DADBE6",
"N 	c #E3EAEF",
"O 	c #E5EBF1",
"P 	c #E0E3E6",
"Q 	c #E5EAF2",
"R 	c #E3E8F2",
"S 	c #D8DDE6",
"T 	c #CED6DC",
"U 	c #D6E9E9",
"V 	c #D3E9E9",
"W 	c #B2BBC8",
"X 	c #8D8F96",
"Y 	c #DEE3ED",
"Z 	c #DCE7EB",
"` 	c #C6D0D0",
" .	c #D8E2E9",
"..	c #D7E0E9",
"+.	c #CED5D4",
"@.	c #D4E4E6",
"#.	c #D3EBE9",
"$.	c #C8D6DF",
"%.	c #43454B",
"&.	c #8D8E98",
"*.	c #D5D8E6",
"=.	c #CBCFD2",
"-.	c #C5CBC5",
";.	c #CFD5CE",
">.	c #C1C9C2",
",.	c #D1DDDE",
"'.	c #C4D0DA",
").	c #55585F",
"!.	c #847F75",
"~.	c #C9C3B4",
"{.	c #C7C5C0",
"].	c #CDC9BE",
"^.	c #D0CAB5",
"/.	c #C1B99B",
"(.	c #2E2E2C",
"_.	c #F0E4BA",
":.	c #EBE0B7",
"<.	c #E1D29A",
"[.	c #D5C58C",
"}.	c #AE9E66",
"|.	c #13110A",
"1.	c #F1E5BC",
"2.	c #F1E7C5",
"3.	c #E6D69E",
"4.	c #D8C890",
"5.	c #AFA06C",
"6.	c #15130C",
"7.	c #EFE4BF",
"8.	c #E7DEBE",
"9.	c #DACD9C",
"0.	c #CDBE8A",
"a.	c #B0A16D",
"b.	c #14120B",
"c.	c #E7DDB6",
"d.	c #E2D8B4",
"e.	c #D1C28C",
"f.	c #C3B47D",
"g.	c #9A8D5F",
"h.	c #0F0D08",
"i.	c #898060",
"j.	c #BCB18A",
"k.	c #AD9F6C",
"l.	c #82764C",
"m.	c #312C1A",
"n.	c #2D2D2D",
"o.	c #404040",
"p.	c #060606",
"                                        ",
"            . + @ # $ % &               ",
"          * = - ; > , ' ) !             ",
"        ~ { ] ^ / ( _ : < [ }           ",
"        | ^ 1 - 2 3 4 5 < 6 7           ",
"      8 = ^ ] 9 0 a 3 b c d e f         ",
"      g h i j j ^ k l 3 m n o p         ",
"      q r 2 s t 9 9 u v w x y z         ",
"      A B C D E F ; G H I J K L         ",
"        M N O P Q R S T U V W           ",
"        X Y Z `  ...+.@.#.$.%.          ",
"          &.*.=.-.;.>.,.'.).            ",
"            !.~.{.].^./.(.              ",
"              _.:.<.[.}.|.              ",
"              1.2.3.4.5.6.              ",
"              7.8.9.0.a.b.              ",
"              c.d.e.f.g.h.              ",
"              i.j.k.l.m.                ",
"                n.o.p.                  ",
"                                        "};
end_of_pixmap

    my $icon_auto = <<'end_of_pixmap';
/* XPM */
static char * preferences_20_xpm[] = {
"20 20 138 2",
"  	c None",
". 	c #201F1E",
"+ 	c #242422",
"@ 	c #1B1B19",
"# 	c #97948D",
"$ 	c #B6B3AC",
"% 	c #5E5C57",
"& 	c #2E2D2C",
"* 	c #9F9C98",
"= 	c #C9C6BE",
"- 	c #20201E",
"; 	c #404040",
"> 	c #676767",
", 	c #888786",
"' 	c #D6D4CE",
") 	c #222120",
"! 	c #272727",
"~ 	c #BBBBBB",
"{ 	c #BEBEBE",
"] 	c #373737",
"^ 	c #7C7B77",
"/ 	c #434341",
"( 	c #000000",
"_ 	c #3E3E3E",
": 	c #CDCCCB",
"< 	c #CECBC4",
"[ 	c #72706B",
"} 	c #262626",
"| 	c #D2D2D2",
"1 	c #696969",
"2 	c #111111",
"3 	c #6C6966",
"4 	c #BCB8B1",
"5 	c #CFCBC4",
"6 	c #E7E6E2",
"7 	c #D6D2CD",
"8 	c #BFBBB1",
"9 	c #BBB6AC",
"0 	c #161514",
"a 	c #A8A8A8",
"b 	c #6B6863",
"c 	c #827F79",
"d 	c #878580",
"e 	c #8E8C87",
"f 	c #C7C5C0",
"g 	c #D0CDC6",
"h 	c #8B8781",
"i 	c #272624",
"j 	c #3D3D3D",
"k 	c #959595",
"l 	c #5B5B5A",
"m 	c #DEDDDC",
"n 	c #D4D2CD",
"o 	c #716F6A",
"p 	c #2C2C2C",
"q 	c #616161",
"r 	c #1A1A1A",
"s 	c #494949",
"t 	c #D2D2D1",
"u 	c #545352",
"v 	c #585858",
"w 	c #131313",
"x 	c #393939",
"y 	c #6B6B6B",
"z 	c #3A3937",
"A 	c #4B4A45",
"B 	c #0D0D0C",
"C 	c #171B20",
"D 	c #191D21",
"E 	c #040505",
"F 	c #7C7C7C",
"G 	c #232323",
"H 	c #717170",
"I 	c #898885",
"J 	c #85827C",
"K 	c #151413",
"L 	c #4B5967",
"M 	c #9CB0C6",
"N 	c #99AFC5",
"O 	c #343F4A",
"P 	c #3C3C3C",
"Q 	c #AAAAA9",
"R 	c #CAC8C3",
"S 	c #242321",
"T 	c #53606E",
"U 	c #A9BACC",
"V 	c #99A8B8",
"W 	c #94ACC3",
"X 	c #6B839D",
"Y 	c #575757",
"Z 	c #D9D9D8",
"` 	c #D6D3CE",
" .	c #87847E",
"..	c #292725",
"+.	c #5A6673",
"@.	c #BBC9D8",
"#.	c #8C9DAF",
"$.	c #7A90A6",
"%.	c #859FBB",
"&.	c #637A93",
"*.	c #484848",
"=.	c #EFEEED",
"-.	c #CBC9C5",
";.	c #74726F",
">.	c #1D1C1B",
",.	c #5F6A76",
"'.	c #D1DAE2",
").	c #899BAF",
"!.	c #748BA4",
"~.	c #87A1BC",
"{.	c #69819B",
"].	c #1D232B",
"^.	c #474747",
"/.	c #C0BFBE",
"(.	c #6A6968",
"_.	c #7D7D7C",
":.	c #4B5661",
"<.	c #D2DCE5",
"[.	c #A3AEBC",
"}.	c #7289A2",
"|.	c #8AA5C0",
"1.	c #5E758D",
"2.	c #1B2128",
"3.	c #3A3A3A",
"4.	c #969695",
"5.	c #8F8F8E",
"6.	c #4D5864",
"7.	c #BFCFDF",
"8.	c #A0ACB9",
"9.	c #92ABC4",
"0.	c #5C728A",
"a.	c #202831",
"b.	c #212932",
"c.	c #6A8097",
"d.	c #93A9C0",
"e.	c #6D8299",
"f.	c #242C36",
"g.	c #1D242C",
"            . +                         ",
"          @ # $ %                       ",
"            & * = -             ; >     ",
"              , ' )           ! ~ { ]   ",
"      ^ / ( _ : < [ (         } | 1 2   ",
"      3 4 5 6 7 8 9 0       } a ( (     ",
"      ( b c d e f g h i   j k (         ",
"          ( ( ( l m n o p q r           ",
"                  s t u v w             ",
"                    x y z A B           ",
"              C D E F G H I J K         ",
"            L M N O (   P Q R h S       ",
"          T U V W X (     Y Z `  ...    ",
"        +.@.#.$.%.&.(       *.=.-.;.>.  ",
"      ,.'.).!.~.{.].          ^./.(._.r ",
"    :.<.[.}.|.1.2.              3.4.5.r ",
"    6.7.8.9.0.a.(                 ( (   ",
"    b.c.d.e.a.(                         ",
"      ].f.g.(                           ",
"                                        "};
end_of_pixmap

    my $icon_load = <<'end_of_pixmap';
/* XPM */
static char * report_20_xpm[] = {
"20 20 95 2",
"  	c None",
". 	c #323232",
"+ 	c #282828",
"@ 	c #292929",
"# 	c #2C2C2C",
"$ 	c #E8E8E8",
"% 	c #ECECEC",
"& 	c #EDEDED",
"* 	c #EEEEEE",
"= 	c #EFEFEF",
"- 	c #F0F0F0",
"; 	c #F1F1F1",
"> 	c #F2F2F2",
", 	c #F3F3F3",
"' 	c #F4F4F4",
") 	c #000000",
"! 	c #E1E1E1",
"~ 	c #E2E2E2",
"{ 	c #E3E3E3",
"] 	c #E5E4E5",
"^ 	c #E6E6E6",
"/ 	c #E7E6E7",
"( 	c #E8E7E8",
"_ 	c #E9E9E9",
": 	c #EAEAEA",
"< 	c #EBECEC",
"[ 	c #ECEDED",
"} 	c #EFEEEF",
"| 	c #F0F0F1",
"1 	c #EDEDEC",
"2 	c #E4E5E5",
"3 	c #777777",
"4 	c #EDEEED",
"5 	c #E4E4E4",
"6 	c #E6E5E5",
"7 	c #E7E7E7",
"8 	c #EEEDEE",
"9 	c #E6E7E6",
"0 	c #EEEEEF",
"a 	c #E5E5E6",
"b 	c #E6E6E7",
"c 	c #E7E7E8",
"d 	c #E9E9E8",
"e 	c #EBEBEB",
"f 	c #EEEDED",
"g 	c #F5F5F5",
"h 	c #E8E8E9",
"i 	c #F7F7F7",
"j 	c #EBECEB",
"k 	c #F6F6F6",
"l 	c #F6F7F7",
"m 	c #F8F8F7",
"n 	c #F8F9F9",
"o 	c #F4F4F5",
"p 	c #F9F9F8",
"q 	c #FAFAFA",
"r 	c #FBFBFB",
"s 	c #EFEFF0",
"t 	c #F4F3F4",
"u 	c #F8F8F8",
"v 	c #F9F8F9",
"w 	c #FAF9F9",
"x 	c #FBFBFC",
"y 	c #FCFCFC",
"z 	c #FDFDFE",
"A 	c #F5F4F4",
"B 	c #F1F0F1",
"C 	c #F1F2F1",
"D 	c #F2F3F3",
"E 	c #F3F3F4",
"F 	c #F7F7F8",
"G 	c #F9F9F9",
"H 	c #FAF9FA",
"I 	c #FEFDFD",
"J 	c #FEFEFE",
"K 	c #F0EFF0",
"L 	c #F1F1F2",
"M 	c #F2F2F3",
"N 	c #F3F4F3",
"O 	c #F4F5F4",
"P 	c #F6F6F5",
"Q 	c #F7F7F6",
"R 	c #FAFAFB",
"S 	c #FCFBFC",
"T 	c #FDFDFD",
"U 	c #9A9A9A",
"V 	c #A2A2A2",
"W 	c #A4A4A4",
"X 	c #A5A5A5",
"Y 	c #A6A6A6",
"Z 	c #A7A7A7",
"` 	c #A8A8A8",
" .	c #A9A9A9",
"..	c #AAAAAA",
"+.	c #979797",
"  . + + + + + @ @ @ @ @ @ @ @ @ #       ",
". $ % & & * * = - - ; > > , ' ' & )     ",
"+ % ! ~ { ] ^ / ( _ : < [ * } - | )     ",
"+ 1 ) ) 2 3 3 3 3 3 3 3 3 3 3 3 > )     ",
"+ 4 { 5 6 7 $ _ : % [ 8 = - ; > , )     ",
"+ * ) ) 9 3 3 3 3 3 3 3 3 3 3 3 ' )     ",
"@ 0 a b c d : e % f 0 - ; > , ' g )     ",
"@ = ) ) h 3 3 3 3 3 3 3 3 3 3 3 i )     ",
"@ - c h : j % & * - ; > , ' k l m )     ",
"@ ; ) ) e 3 3 3 3 3 3 3 3 3 3 3 n )     ",
"@ ; : e % & * - ; > , o g k m p q )     ",
"@ > ) ) & 3 3 3 3 3 3 3 3 3 3 3 r )     ",
"@ , % & * s | > , t g k u v w x y )     ",
"@ ' ) ) = 3 3 3 3 3 3 3 3 3 3 3 z )     ",
"@ A * = B C D E g k F G H r y I J )     ",
"@ % K L M N O P Q u G R S y T J ' )     ",
"# U V W W X X Y Y Z `  .......` +.)     ",
"  ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) )       ",
"                                        ",
"                                        "};
end_of_pixmap

    $pixmap_main = $top->Pixmap(-data => $icon_main);
    $pixmap_open = $top->Pixmap(-data => $icon_open);
    $pixmap_info = $top->Pixmap(-data => $icon_info);
    $pixmap_auto = $top->Pixmap(-data => $icon_auto);
    $pixmap_load = $top->Pixmap(-data => $icon_load);
}

1;
