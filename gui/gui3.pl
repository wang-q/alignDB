#!/usr/bin/perl
use strict;
use warnings;

package AlignDB::GUI;
use Moose;
use MooseX::AttributeHelpers;
use Carp;

use Gtk2 '-init';
use Glib qw(TRUE FALSE);
use Gtk2::GladeXML;
use Gtk2::Helper;

use Config::Tiny;
use Proc::Background;
use File::Spec;
use Text::CSV_XS;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

has 'app'  => ( is => 'ro', isa => 'Object', );
has 'win'  => ( is => 'ro', isa => 'Object', );
has 'text' => ( is => 'ro', isa => 'Object', );

has 'processes' => (
    metaclass => 'Collection::Array',
    is        => 'ro',
    isa       => 'ArrayRef',
    default   => sub { [] },
    provides  => {
        push     => 'add_processes',
        elements => 'all_processes',
        count    => 'count_processes',
    }
);

sub BUILD {
    my $self = shift;

    # Load the UI from the Glade-2 file
    my $app = Gtk2::GladeXML->new("$FindBin::Bin/gui3.glade");
    $self->{app} = $app;

    # Connect signals magically
    $app->signal_autoconnect_from_package($self);

    $self->read_config;
    $self->fill_combobox;

    {    # Init the console textview
        my $textview   = $app->get_widget('textview_console');
        my $textbuffer = $textview->get_buffer;
        $self->{text} = $textbuffer;

        $textbuffer->create_tag( "bold",   font => "Courier Bold 9", );
        $textbuffer->create_tag( "normal", font => "Courier 8", );
        $textbuffer->create_tag( "italic", font => "Courier Italic 8", );

        # create a mark at the end of the buffer, with right gravity,
        # so that when you insert text, the mark always stays on
        # the right (at the end).
        my $end_mark
            = $textbuffer->create_mark( 'end', $textbuffer->get_end_iter,
            FALSE );

        # every time we insert text, scroll to that mark.
        $textbuffer->signal_connect(
            insert_text => sub {
                $textview->scroll_to_mark( $end_mark, 0.0, TRUE, 0.0, 1.0 );
            }
        );
    }

    my $win = $app->get_widget('window_main');
    $self->{win} = $win;
    $win->signal_connect( 'delete-event' => sub { Gtk2->main_quit } );
    $win->show;

    $self->on_togglebutton_growl_send_toggled;

    Gtk2->main;
    return;
}

#----------------------------#
# setter and getter of all entries or checkbuttons
#----------------------------#
sub set_value {
    my $self  = shift;
    my $name  = shift;
    my $value = shift;

    my $widget = $self->{app}->get_widget($name);
    my $class  = ref $widget;

    if ( $class eq 'Gtk2::Entry' ) {
        $widget->set_text($value);
    }
    elsif ( $class eq 'Gtk2::CheckButton' ) {
        $widget->set_active($value);
    }
    elsif ( $class eq 'Gtk2::ToggleButton' ) {
        $widget->set_active($value);
    }
    else {
        carp "Widget type is [$class]\n";
        carp "$name doesn't exist or bad types.\n";
    }

    return;
}

sub get_value {
    my $self = shift;
    my $name = shift;

    my $widget = $self->{app}->get_widget($name);
    my $class  = ref $widget;

    my $value;

    if ( $class eq 'Gtk2::Entry' ) {
        $value = $widget->get_text;
    }
    elsif ( $class eq 'Gtk2::CheckButton' ) {
        $value = $widget->get_active ? 1 : 0;
    }
    elsif ( $class eq 'Gtk2::ToggleButton' ) {
        $value = $widget->get_active ? 1 : 0;
    }
    elsif ( $class eq 'Gtk2::ComboBox' ) {
        $value = $widget->get_active_text;
    }
    else {
        carp "Widget type is [$class]\n";
        carp "$name doesn't exist or bad types.\n";
    }

    return $value;
}

sub append_text {
    my $self   = shift;
    my $string = shift;
    my $tags   = shift || 'normal';

    my $text = $self->text;
    $text->insert_with_tags_by_name( $text->get_end_iter, $string, $tags );
    return;
}

#----------------------------#
# execute system commands
#----------------------------#
sub exec_cmd {
    my $self = shift;
    my $cmd  = shift;

    $self->append_text( "\n" . "=" x 12 . "CMD" . "=" x 15 . "\n", "bold" );
    $self->append_text( $cmd . "\n",                               "bold" );
    $self->append_text( "=" x 30 . "\n",                           "bold" );

    $cmd .= ' 1>&2';    # redirect STDOUT to STDERR so we can see the outputs

    my $proc = Proc::Background->new($cmd);
    $proc->{cmd} = $cmd;
    $self->add_processes($proc);

    return;
}

#----------------------------#
# read-out configs
#----------------------------#
sub read_config {
    my $self = shift;

    my $Config = Config::Tiny->new;
    $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

    # Database init values
    $self->set_value( "entry_server",   $Config->{database}{server} );
    $self->set_value( "entry_port",     $Config->{database}{port} );
    $self->set_value( "entry_username", $Config->{database}{username} );
    $self->set_value( "entry_password", $Config->{database}{password} );

    # target, query init values
    $self->set_value( "entry_target_id",   $Config->{taxon}{target_taxon_id} );
    $self->set_value( "entry_target_name", $Config->{taxon}{target_name} );
    $self->set_value( "entry_query_id",    $Config->{taxon}{query_taxon_id} );
    $self->set_value( "entry_query_name",  $Config->{taxon}{query_name} );
    $self->set_value( "entry_db_name",     $Config->{database}{db} );
    $self->set_value( "entry_ensembl",     $Config->{database}{ensembl} );

    # generate
    $self->set_value( "entry_axt_dir", $Config->{taxon}{axt_dir} );
    $self->set_value( "entry_axt_threshold",
        $Config->{generate}{axt_threshold} );
    $self->set_value( "entry_parallel", $Config->{generate}{parallel} );

    # insert GC
    $self->set_value( "checkbutton_insert_gc", $Config->{gc}{insert_gc} );
    $self->set_value( "checkbutton_insert_segment",
        $Config->{gc}{insert_segment} );

    # insert gene
    $self->set_value( "checkbutton_insert_exonsw",
        $Config->{gene}{insert_exonsw} );
    $self->set_value( "checkbutton_insert_codingsw",
        $Config->{gene}{insert_codingsw} );

    # update feature
    $self->set_value( "checkbutton_process_align", $Config->{feature}{align} );
    $self->set_value( "checkbutton_process_indel", $Config->{feature}{indel} );
    $self->set_value( "checkbutton_process_isw",   $Config->{feature}{isw} );
    $self->set_value( "checkbutton_process_snp",   $Config->{feature}{snp} );
    $self->set_value( "checkbutton_process_window",
        $Config->{feature}{window} );

    # three-way
    $self->set_value( "entry_first_db",  $Config->{ref}{first_db} );
    $self->set_value( "entry_second_db", $Config->{ref}{second_db} );
    $self->set_value( "entry_goal_db",   $Config->{ref}{goal_db} );

    # common stat parameter
    $self->set_value( "entry_common_run", $Config->{stat}{run} );
    $self->set_value( "entry_common_threshold",
        $Config->{stat}{sum_threshold} );
    $self->set_value( "checkbutton_common_jc", $Config->{stat}{jc_correction} );
    $self->set_value( "checkbutton_common_stamp", $Config->{stat}{time_stamp} );
    $self->set_value( "checkbutton_common_add_index",
        $Config->{stat}{add_index_sheet} );

    # gc stat parameter
    $self->set_value( "entry_gc_run",         $Config->{stat}{run} );
    $self->set_value( "entry_gc_threshold",   $Config->{stat}{sum_threshold} );
    $self->set_value( "checkbutton_gc_jc",    $Config->{stat}{jc_correction} );
    $self->set_value( "checkbutton_gc_stamp", $Config->{stat}{time_stamp} );
    $self->set_value( "checkbutton_gc_add_index",
        $Config->{stat}{add_index_sheet} );

    # three stat parameter
    $self->set_value( "entry_three_run",       $Config->{stat}{run} );
    $self->set_value( "entry_three_threshold", $Config->{stat}{sum_threshold} );
    $self->set_value( "checkbutton_three_jc",  $Config->{stat}{jc_correction} );
    $self->set_value( "checkbutton_three_stamp", $Config->{stat}{time_stamp} );
    $self->set_value( "checkbutton_three_add_index",
        $Config->{stat}{add_index_sheet} );

    # growl parameter
    $self->set_value( "entry_growl_appname",  $Config->{growl}{appname} );
    $self->set_value( "entry_growl_host",     $Config->{growl}{host} );
    $self->set_value( "entry_growl_password", $Config->{growl}{password} );
    $self->set_value( "checkbutton_growl_starting",
        $Config->{growl}{starting} );
    $self->set_value( "checkbutton_growl_ending", $Config->{growl}{ending} );
    $self->set_value( "checkbutton_growl_other",  $Config->{growl}{other} );
    $self->set_value( "togglebutton_growl_send",  $Config->{growl}{send} );

    return;
}

sub fill_combobox {
    my $self = shift;

    my $model = Gtk2::ListStore->new('Glib::String');
    for (qw{ 0target 0query 1target 1query }) {
        $model->set( $model->append, 0, $_ );
    }
    for (qw{ combobox_first combobox_second combobox_outgroup }) {
        my $cb = $self->{app}->get_widget($_);
        $cb->set_model($model);
        my $cr = Gtk2::CellRendererText->new;
        $cb->pack_start( $cr, TRUE );
        $cb->add_attribute( $cr, 'text', 0 );
    }

    $self->{app}->get_widget("combobox_first")->set_active(0);
    $self->{app}->get_widget("combobox_second")->set_active(1);
    $self->{app}->get_widget("combobox_outgroup")->set_active(3);

    return;
}

#----------------------------#
# menubar and toolbar events
#----------------------------#
sub on_imagemenuitem_about_activate {
    my $self   = shift;
    my $widget = shift;

    Gtk2->show_about_dialog(
        Gtk2::Window->new,
        program_name => 'AlignDB GUI3',
        version      => '0.7',
        copyright    => "(C) 2004-2012 WANG, Qiang",
        authors      => ['WANG, Qiang <wangq@nju.edu.cn>'],
        documenters  => ['WANG, Qiang <wangq@nju.edu.cn>'],
        artists      => ['WANG, Qiang <wangq@nju.edu.cn>'],
        comments     => "The third generation of GUI interface for AlignDB",
        title        => "About AlignDB GUI3",
        website      => "http://chenlab.nju.edu.cn",
        wrap_license => TRUE,
        license =>
            "This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.\n",
    );

    return;
}

sub on_imagemenuitem_quit_activate {
    my $self   = shift;
    my $widget = shift;

    Gtk2->main_quit;
    return;
}

sub on_toolbutton_default_clicked {
    my $self   = shift;
    my $widget = shift;

    $self->read_config;
    return;
}

sub on_toolbutton_process_clicked {
    my $self   = shift;
    my $widget = shift;

    my $count = $self->count_processes;
    $self->append_text( "There are $count process(es) totally.\n", "italic" );
    return unless $count > 0;

    for my $proc ( $self->all_processes ) {
        $proc->alive;
        $self->append_text(
            Dump {
                start => scalar localtime $proc->start_time,
                end   => $proc->end_time
                ? scalar localtime $proc->end_time
                : undef,
                alive => $proc->alive ? 'yes' : 'no',
                pid   => $proc->pid,
                cmd   => $proc->{cmd},
            }
        );
    }

    return;
}

sub on_togglebutton_growl_send_toggled {
    my $self   = shift;
    my $widget = shift;

    $ENV{growl_send}     = $self->get_value("togglebutton_growl_send");
    $ENV{growl_appname}  = $self->get_value("entry_growl_appname");
    $ENV{growl_host}     = $self->get_value("entry_growl_host");
    $ENV{growl_password} = $self->get_value("entry_growl_password");
    $ENV{growl_starting} = $self->get_value("checkbutton_growl_starting");
    $ENV{growl_ending}   = $self->get_value("checkbutton_growl_ending");
    $ENV{growl_other}    = $self->get_value("checkbutton_growl_other");

    return;
}

#----------------------------#
# dialogs
#----------------------------#
sub dialog_taxon {
    my $self = shift;

    # Create the dialog and init two buttons
    my $dialog = Gtk2::Dialog->new(
        "Load taxon...",
        $self->win, [qw{modal destroy-with-parent}],
        'gtk-ok'    => 'ok',
        'gtk-close' => 'close',
    );

    # This vbox comes with the dialog
    my $vbox = $dialog->vbox;

    my $label = Gtk2::Label->new("Choose a taxon:");
    $label->set_alignment( 0, 0.5 );
    $label->set_justify('left');
    $vbox->pack_start( $label, FALSE, FALSE, 5 );

    # read out normal taxons and put them into listmodel
    my $file = "$FindBin::Bin/normal_taxon.csv";
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $file or die "$file: $!";
    $csv->getline($csv_fh);    # bypass title line

    my $model
        = Gtk2::ListStore->new( 'Glib::Int', 'Glib::String', 'Glib::String' );
    while ( my $row = $csv->getline($csv_fh) ) {
        my $id      = $row->[0];
        my $name    = $row->[1];
        my $species = $row->[2] . ' ' . $row->[3];

        # The iter is a pointer in the treestore. We use to add data.
        my $iter = $model->append;
        $model->set(
            $iter,
            0 => $id,
            1 => $name,
            2 => $species,
        );
    }
    close $csv_fh;

    my $treeview = Gtk2::TreeView->new_with_model($model);

    # Add columns
    $treeview->insert_column_with_attributes( 0, "id",
        Gtk2::CellRendererText->new, text => 0 );
    $treeview->insert_column_with_attributes( 1, "name",
        Gtk2::CellRendererText->new, text => 1 );
    $treeview->insert_column_with_attributes( 2, "species",
        Gtk2::CellRendererText->new, text => 2 );
    $treeview->get_column(0)->set_sort_column_id(0);
    $treeview->get_column(1)->set_sort_column_id(1);
    $treeview->get_column(2)->set_sort_column_id(2);

    # get the Gtk2::TreeSelection of $treeview
    my $treeselection = $treeview->get_selection;
    $treeselection->set_mode('single');

    # add TreeView to a scrolledwindow
    my $sw = Gtk2::ScrolledWindow->new;
    $sw->set_size_request( 360, 240 );
    $sw->set_policy( 'automatic', 'automatic' );
    $sw->add($treeview);
    $vbox->pack_start( $sw, TRUE, TRUE, 0 );
    $vbox->show_all;

    my $response = $dialog->run;
    my ( $id, $name );
    if ( $response eq 'ok' ) {
        my $iter = $treeselection->get_selected;
        if ( defined $iter ) {

            # we want data at the model's columns where the iter is pointing
            $id   = $model->get( $iter, 0 );
            $name = $model->get( $iter, 1 );
        }
    }
    $dialog->destroy;

    return ( $id, $name );
}

sub dialog_choose_db {
    my $self = shift;

    # Create the dialog and init two buttons
    my $dialog = Gtk2::Dialog->new(
        "Choose DB...",
        $self->win, [qw{modal destroy-with-parent}],
        'gtk-ok'    => 'ok',
        'gtk-close' => 'close',
    );

    # This vbox comes with the dialog
    my $vbox = $dialog->vbox;

    my $label = Gtk2::Label->new("Choose an AlignDB database:");
    $label->set_alignment( 0, 0.5 );
    $label->set_justify('left');
    $vbox->pack_start( $label, FALSE, FALSE, 5 );

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");

    my $cmd
        = "mysql"
        . " -h$server"
        . " -P$port"
        . " -u$username"
        . " -p$password"
        . ' -e "SHOW DATABASES"';
    my @dbs = grep {/vs/} split "\n", `$cmd`;

    my $model = Gtk2::ListStore->new('Glib::String');
    for (@dbs) {
        my $iter = $model->append;
        $model->set( $iter, 0 => $_, );
    }

    my $treeview = Gtk2::TreeView->new_with_model($model);

    # Add columns
    $treeview->insert_column_with_attributes( 0, "DB name",
        Gtk2::CellRendererText->new, text => 0 );
    $treeview->get_column(0)->set_sort_column_id(0);

    # get the Gtk2::TreeSelection of $treeview
    my $treeselection = $treeview->get_selection;
    $treeselection->set_mode('single');

    # add TreeView to a scrolledwindow
    my $sw = Gtk2::ScrolledWindow->new;
    $sw->set_size_request( 360, 240 );
    $sw->set_policy( 'automatic', 'automatic' );
    $sw->add($treeview);
    $vbox->pack_start( $sw, TRUE, TRUE, 0 );
    $vbox->show_all;

    my $response = $dialog->run;
    my $db_name;
    if ( $response eq 'ok' ) {
        my $iter = $treeselection->get_selected;
        if ( defined $iter ) {
            $db_name = $model->get( $iter, 0 );
        }
    }

    $dialog->destroy;

    if ( defined $db_name ) {

        my $obj = AlignDB->new(
            mysql  => "$db_name:$server",
            user   => $username,
            passwd => $password,
        );

        my ( $target_id,   $query_id )   = $obj->get_taxon_ids;
        my ( $target_name, $query_name ) = $obj->get_names;

        return {
            target_id   => $target_id,
            query_id    => $query_id,
            target_name => $target_name,
            query_name  => $query_name,
            db_name     => $db_name,
        };
    }

    return;
}

#----------------------------#
# button callback events
#----------------------------#
sub on_button_load_target_clicked {
    my $self   = shift;
    my $widget = shift;

    my ( $id, $name ) = $self->dialog_taxon;

    if ( defined $id and defined $name ) {
        $self->set_value( "entry_target_id",   $id );
        $self->set_value( "entry_target_name", $name );
    }

    return;
}

sub on_button_load_query_clicked {
    my $self   = shift;
    my $widget = shift;

    my ( $id, $name ) = $self->dialog_taxon;

    if ( defined $id and defined $name ) {
        $self->set_value( "entry_query_id",   $id );
        $self->set_value( "entry_query_name", $name );
    }

    return;
}

sub on_button_auto_db_name_clicked {
    my $self   = shift;
    my $widget = shift;

    my $target_name = $self->get_value("entry_target_name");
    my $query_name  = $self->get_value("entry_query_name");
    $self->set_value( "entry_db_name", "$target_name" . "vs" . "$query_name" );

    return;
}

sub on_button_choose_db_clicked {
    my $self   = shift;
    my $widget = shift;

    my $result = $self->dialog_choose_db;
    return unless $result;

    $self->set_value( "entry_target_id",   $result->{target_id} );
    $self->set_value( "entry_query_id",    $result->{query_id} );
    $self->set_value( "entry_target_name", $result->{target_name} );
    $self->set_value( "entry_query_name",  $result->{query_name} );
    $self->set_value( "entry_db_name",     $result->{db_name} );

    return;
}

sub on_button_open_axt_dir_clicked {
    my $self   = shift;
    my $widget = shift;

    my $dia = Gtk2::FileChooserDialog->new(
        'Choose a dir', $self->win, 'select-folder',
        'gtk-cancel' => 'cancel',
        'gtk-ok'     => 'ok'
    );

    $dia->show;
    $dia->signal_connect(
        'response' => sub {
            my ( $dia, $response_id ) = @_;
            if ( $response_id eq 'ok' ) {
                my $dir = $dia->get_filename;
                if ( defined $dir ) {
                    $self->set_value( "entry_axt_dir", $dir );
                }
            }
            $dia->destroy;
            return FALSE;
        }
    );

    return;
}

sub on_button_auto_common_stat_file_clicked {
    my $self   = shift;
    my $widget = shift;

    my $db_name = $self->get_value("entry_db_name");
    my $outfile = "$FindBin::Bin/../stat/$db_name.common.xlsx";
    $outfile = File::Spec->rel2abs($outfile);
    $self->set_value( "entry_common_stat_file", $outfile );

    return;
}

sub on_button_auto_gc_stat_file_clicked {
    my $self   = shift;
    my $widget = shift;

    my $db_name = $self->get_value("entry_db_name");
    my $outfile = "$FindBin::Bin/../stat/$db_name.gc.xlsx";
    $outfile = File::Spec->rel2abs($outfile);
    $self->set_value( "entry_gc_stat_file", $outfile );

    return;
}

sub on_button_auto_three_stat_file_clicked {
    my $self   = shift;
    my $widget = shift;

    my $db_name = $self->get_value("entry_goal_db");
    my $outfile = "$FindBin::Bin/../stat/$db_name.three.xlsx";
    $outfile = File::Spec->rel2abs($outfile);
    $self->set_value( "entry_three_stat_file", $outfile );

    return;
}

sub on_button_choose_first_db_clicked {
    my $self   = shift;
    my $widget = shift;

    my $result = $self->dialog_choose_db;
    return unless $result;

    $self->set_value( "entry_first_db", $result->{db_name} );
    return;
}

sub on_button_choose_second_db_clicked {
    my $self   = shift;
    my $widget = shift;

    my $result = $self->dialog_choose_db;
    return unless $result;

    $self->set_value( "entry_second_db", $result->{db_name} );
    return;
}

sub on_button_choose_goal_db_clicked {
    my $self   = shift;
    my $widget = shift;

    my $result = $self->dialog_choose_db;
    return unless $result;

    $self->set_value( "entry_goal_db", $result->{db_name} );
    return;
}

sub on_button_auto_goal_db_name_clicked {
    my $self   = shift;
    my $widget = shift;

    my $first_db  = $self->get_value("entry_first_db");
    my $second_db = $self->get_value("entry_second_db");
    my $goal_db;

    my $first    = $self->get_value("combobox_first");
    my $second   = $self->get_value("combobox_second");
    my $outgroup = $self->get_value("combobox_outgroup");

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");

    my $first_obj = AlignDB->new(
        mysql  => "$first_db:$server",
        user   => $username,
        passwd => $password,
    );
    my ( $first_target_name, $first_query_name ) = $first_obj->get_names;

    my $second_obj = AlignDB->new(
        mysql  => "$second_db:$server",
        user   => $username,
        passwd => $password,
    );
    my ( $second_target_name, $second_query_name ) = $second_obj->get_names;

    my %name_of = (
        '0target' => $first_target_name,
        '0query'  => $first_query_name,
        '1target' => $second_target_name,
        '1query'  => $second_query_name,
    );
    $goal_db
        = $name_of{$first} . 'vs'
        . $name_of{$second} . 'ref'
        . $name_of{$outgroup};

    $self->set_value( "entry_goal_db", $goal_db );
    return;
}

sub on_button_clear_clicked {
    my $self   = shift;
    my $widget = shift;

    my $text = $self->text;
    $text->delete( $text->get_start_iter, $text->get_end_iter );

    return;
}

sub on_button_copy_clicked {
    my $self   = shift;
    my $widget = shift;

    my $text = $self->text;
    $text->select_range( $text->get_start_iter, $text->get_end_iter );
    my $clip = $self->app->get_widget('textview_console')->get_clipboard;
    $text->copy_clipboard($clip);

    return;
}

#----------------------------------------------------------#
# run *.pl events
#----------------------------------------------------------#

sub on_button_test_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");

    my $cmd
        = "mysql"
        . " -h$server"
        . " -P$port"
        . " -u$username"
        . " -p$password"
        . q{ -e "SELECT 123456789"};
    my $success = grep {/123456789/} split "\n", `$cmd`;

    if ($success) {
        $self->append_text("Test connection successful\n");
    }
    else {
        $self->append_text("Test connection failed\n");
    }

    return;
}

sub on_button_init_aligndb_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $cmd
        = "perl $FindBin::Bin/../init/init_alignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_gen_aligndb_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $target_id     = $self->get_value("entry_target_id");
    my $target_name   = $self->get_value("entry_target_name");
    my $query_id      = $self->get_value("entry_query_id");
    my $query_name    = $self->get_value("entry_query_name");
    my $axt_dir       = $self->get_value("entry_axt_dir");
    my $axt_threshold = $self->get_value("entry_axt_threshold");

    my $cmd
        = "perl $FindBin::Bin/../init/gen_alignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -t=\"$target_id,$target_name\""
        . " -q=\"$query_id,$query_name\""
        . " -a=$axt_dir"
        . " -l=$axt_threshold"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);

    return;
}

sub on_button_insert_gc_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $insert_gc      = $self->get_value("checkbutton_insert_gc");
    my $insert_segment = $self->get_value("checkbutton_insert_segment");

    my $cmd
        = "perl $FindBin::Bin/../init/insert_gc.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --insert_gc=$insert_gc"
        . " --insert_segment=$insert_segment"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_insert_gene_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");
    my $ensembl  = $self->get_value("entry_ensembl");

    my $parallel = $self->get_value("entry_parallel");

    my $insert_exonsw   = $self->get_value("checkbutton_insert_exonsw");
    my $insert_codingsw = $self->get_value("checkbutton_insert_codingsw");

    my $cmd
        = "perl $FindBin::Bin/../gene/insert_gene.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -e=$ensembl"
        . " --insert_exonsw=$insert_exonsw"
        . " --insert_codingsw=$insert_codingsw"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_swcv_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $cmd
        = "perl $FindBin::Bin/../init/update_sw_cv.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_feature_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $ensembl = $self->get_value("entry_ensembl");

    my $parallel = $self->get_value("entry_parallel");

    my $process_align  = $self->get_value("checkbutton_process_align");
    my $process_indel  = $self->get_value("checkbutton_process_indel");
    my $process_isw    = $self->get_value("checkbutton_process_isw");
    my $process_snp    = $self->get_value("checkbutton_process_snp");
    my $process_window = $self->get_value("checkbutton_process_window");

    my $cmd
        = "perl $FindBin::Bin/../init/update_feature.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -e=$ensembl"
        . " --process_align=$process_align"
        . " --process_indel=$process_indel"
        . " --process_isw=$process_isw"
        . " --process_snp=$process_snp"
        . " --process_window=$process_window"
        . " --parallel=$parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_slippage_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $cmd
        = "perl $FindBin::Bin/../init/update_indel_slippage.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_segment_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $cmd
        = "perl $FindBin::Bin/../init/update_segment.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_join_dbs_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");

    my $first_db  = $self->get_value("entry_first_db");
    my $second_db = $self->get_value("entry_second_db");
    my $goal_db   = $self->get_value("entry_goal_db");

    my $first    = $self->get_value("combobox_first");
    my $second   = $self->get_value("combobox_second");
    my $outgroup = $self->get_value("combobox_outgroup");

    my $cmd
        = "perl $FindBin::Bin/../extra/join_dbs.pl"
        . " --server=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " --dbs=$first_db,$second_db"
        . " --goal_db=$goal_db"
        . " --outgroup=$outgroup"
        . " --target=$first"
        . " --queries=$second";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_cpg_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_goal_db");

    my $cmd
        = "perl $FindBin::Bin/../extra/update_snp_cpg.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_common_stat_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $output    = $self->get_value("entry_common_stat_file");
    my $run       = $self->get_value("entry_common_run");
    my $threshold = $self->get_value("entry_common_threshold");

    my $cmd
        = "perl $FindBin::Bin/../stat/common_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$output"
        . " -r=$run"
        . " -t=$threshold";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_common_chart_clicked {
    my $self   = shift;
    my $widget = shift;

    my $stat_file = $self->get_value("entry_common_stat_file");
    return if $^O ne "MSWin32";
    return if !$stat_file;

    my $jc_correction   = $self->get_value("checkbutton_common_jc");
    my $time_stamp      = $self->get_value("checkbutton_common_stamp");
    my $add_index_sheet = $self->get_value("checkbutton_common_add_index");

    my $cmd
        = "perl $FindBin::Bin/../stat/common_chart_factory.pl"
        . " -i=$stat_file"
        . " -j=$jc_correction"
        . " -t=$time_stamp"
        . " -a=$add_index_sheet";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_gc_stat_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $output    = $self->get_value("entry_gc_stat_file");
    my $run       = $self->get_value("entry_gc_run");
    my $threshold = $self->get_value("entry_gc_threshold");

    my $cmd
        = "perl $FindBin::Bin/../stat/gc_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$output"
        . " -r=$run"
        . " -t=$threshold";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_gc_chart_clicked {
    my $self   = shift;
    my $widget = shift;

    my $stat_file = $self->get_value("entry_gc_stat_file");
    return if $^O ne "MSWin32";
    return if !$stat_file;

    my $jc_correction   = $self->get_value("checkbutton_gc_jc");
    my $time_stamp      = $self->get_value("checkbutton_gc_stamp");
    my $add_index_sheet = $self->get_value("checkbutton_gc_add_index");

    my $cmd
        = "perl $FindBin::Bin/../stat/gc_chart_factory.pl"
        . " -i=$stat_file"
        . " -j=$jc_correction"
        . " -t=$time_stamp"
        . " -a=$add_index_sheet";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_three_stat_clicked {
    my $self   = shift;
    my $widget = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $output    = $self->get_value("entry_three_stat_file");
    my $run       = $self->get_value("entry_three_run");
    my $threshold = $self->get_value("entry_three_threshold");

    my $cmd
        = "perl $FindBin::Bin/../stat/three_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$output"
        . " -r=$run"
        . " -t=$threshold";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_three_chart_clicked {
    my $self   = shift;
    my $widget = shift;

    my $stat_file = $self->get_value("entry_three_stat_file");
    return if $^O ne "MSWin32";
    return if !$stat_file;

    my $jc_correction   = $self->get_value("checkbutton_three_jc");
    my $time_stamp      = $self->get_value("checkbutton_three_stamp");
    my $add_index_sheet = $self->get_value("checkbutton_three_add_index");

    my $cmd
        = "perl $FindBin::Bin/../stat/three_chart_factory.pl"
        . " -i=$stat_file"
        . " -j=$jc_correction"
        . " -t=$time_stamp"
        . " -a=$add_index_sheet";

    $self->exec_cmd($cmd);
    return;
}

package main;

AlignDB::GUI->new;

exit;
