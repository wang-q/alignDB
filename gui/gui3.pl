#!/usr/bin/perl
use strict;
use warnings;
use autodie;

package AlignDB::GUI;
use Moose;
use MooseX::AttributeHelpers;

use Config::Tiny;
use FindBin;
use YAML::Syck;

use Gtk3 '-init';
use Glib qw(TRUE FALSE);

use Path::Tiny;
use Proc::Background;
use List::MoreUtils qw(uniq);

use lib "$FindBin::RealBin/../lib";
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

    # Load the UI
    $self->{app} = Gtk3::Builder->new;
    $self->{app}->add_from_file("$FindBin::RealBin/ui.xml");

    # Connect signals magically
    $self->{app}->connect_signals( undef, $self );

    $self->read_config;

    {    # Init the console textview
        my $textview   = $self->{app}->get_object('textview_console');
        my $textbuffer = $textview->get_buffer;
        $self->{text} = $textbuffer;

        $textbuffer->create_tag( "bold",   font       => "Courier Bold 10", );
        $textbuffer->create_tag( "normal", font       => "Courier 10", );
        $textbuffer->create_tag( "italic", font       => "Courier Italic 10", );
        $textbuffer->create_tag( "blue",   foreground => "blue" );

        # create a mark at the end of the buffer, with right gravity,
        # so that when you insert text, the mark always stays on
        # the right (at the end).
        my $end_mark = $textbuffer->create_mark( 'end', $textbuffer->get_end_iter, FALSE );

        # every time we insert text, scroll to that mark.
        $textbuffer->signal_connect(
            insert_text => sub {
                $textview->scroll_to_mark( $end_mark, 0.0, TRUE, 0.0, 1.0 );
            }
        );
    }

    $self->{win} = $self->{app}->get_object('window_main');
    $self->{win}->signal_connect( 'delete-event' => sub { Gtk3->main_quit } );
    $self->{win}->show;

    # set label_db_name color
    $self->{app}->get_object('label_db_name')
        ->set_markup("<span foreground='blue'>db name:</span>");
    Gtk3->main;
    return;
}

#----------------------------#
# setter and getter of all entries or checkbuttons
#----------------------------#
sub set_value {
    my $self  = shift;
    my $name  = shift;
    my $value = shift;

    my $object = $self->{app}->get_object($name);
    my $class  = ref $object;

    if ( $class eq 'Gtk3::Entry' ) {
        $object->set_text($value);
    }
    elsif ( $class eq 'Gtk3::CheckButton' ) {
        $object->set_active($value);
    }
    elsif ( $class eq 'Gtk3::ToggleButton' ) {
        $object->set_active($value);
    }
    else {
        warn "Widget type is [$class]\n";
        warn "$name doesn't exist or bad types.\n";
    }

    return;
}

sub get_value {
    my $self = shift;
    my $name = shift;

    my $object = $self->{app}->get_object($name);
    my $class  = ref $object;

    my $value;

    if ( $class eq 'Gtk3::Entry' ) {
        $value = $object->get_text;
    }
    elsif ( $class eq 'Gtk3::CheckButton' ) {
        $value = $object->get_active ? 1 : 0;
    }
    elsif ( $class eq 'Gtk3::ToggleButton' ) {
        $value = $object->get_active ? 1 : 0;
    }
    elsif ( $class eq 'Gtk3::ComboBox' ) {
        my $iter  = $object->get_active_iter;
        my $model = $object->get_model;
        $value = $model->get( $iter, 0 );
    }
    elsif ( $class eq 'Gtk3::Notebook' ) {
        $value = $object->get_current_page;
    }
    else {
        warn "Widget type is [$class]\n";
        warn "$name doesn't exist or bad types.\n";
    }

    return $value;
}

sub append_text {
    my $self   = shift;
    my $string = shift;
    my (@tags) = @_;

    (@tags) = ('normal') unless @tags;

    my $text = $self->text;
    $text->insert_with_tags_by_name( $text->get_end_iter, $string, @tags );
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
# read configs
#----------------------------#
sub read_config {
    my $self = shift;

    my $Config = Config::Tiny->new;
    $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

    # parallel init values
    $self->set_value( "entry_parallel", $Config->{generate}{parallel} );

    # server init values
    $self->set_value( "entry_server",   $Config->{database}{server} );
    $self->set_value( "entry_port",     $Config->{database}{port} );
    $self->set_value( "entry_username", $Config->{database}{username} );
    $self->set_value( "entry_password", $Config->{database}{password} );

    # database init values
    $self->set_value( "entry_db_name",          $Config->{database}{db} );
    $self->set_value( "entry_ensembl",          $Config->{database}{ensembl} );
    $self->set_value( "entry_length_threshold", $Config->{generate}{length_threshold} );

    # axt
    $self->set_value( "entry_target_name",   $Config->{taxon}{target_name} );
    $self->set_value( "entry_query_name",    $Config->{taxon}{query_name} );
    $self->set_value( "entry_dir_align_axt", path( $Config->{taxon}{dir_align} )->stringify );

    # fas
    $self->set_value( "entry_dir_align_fas", path( $Config->{taxon}{dir_align_fas} )->stringify );

    # insert GC
    $self->set_value( "checkbutton_insert_gc",      $Config->{gc}{insert_gc} );
    $self->set_value( "checkbutton_insert_segment", $Config->{gc}{insert_segment} );

    return;
}

#----------------------------#
# menubar and toolbar events
#----------------------------#
sub on_toolbutton_about_clicked {
    my $self = shift;

    Gtk3->show_about_dialog(
        Gtk3::Window->new,
        program_name => 'AlignDB GUI3',
        version      => '0,9',
        copyright    => "(C) 2005-2015 Qiang Wang",
        authors      => ['Qiang Wang <wangq@nju.edu.cn>'],
        comments     => "The GUI interface for AlignDB",
        title        => "About AlignDB GUI3",
        wrap_license => TRUE,
        license =>
            "This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.\n",
    );

    return;
}

sub on_toolbutton_quit_clicked {
    my $self = shift;

    Gtk3->main_quit;
    return;
}

sub on_toolbutton_default_clicked {
    my $self = shift;

    $self->read_config;
    return;
}

sub on_toolbutton_process_clicked {
    my $self = shift;

    my $count = $self->count_processes;

    for my $proc ( $self->all_processes ) {
        $proc->alive;
        if ( $proc->alive ) {
            $self->append_text(
                YAML::Syck::Dump {
                    start => scalar localtime $proc->start_time,
                    pid   => $proc->pid,
                    cmd   => $proc->{cmd},
                }
            );
        }
        else {
            $self->append_text(
                YAML::Syck::Dump {
                    end => scalar localtime $proc->end_time,
                    cmd => $proc->{cmd},
                }
            );
        }
    }

    $self->append_text( "-" x 50 . "\n" );
    $self->append_text( "There are ", "italic" );
    $self->append_text( $count, "bold", "blue" );
    $self->append_text( " process(es) totally.\n", "italic" );
    $self->append_text( "-" x 50 . "\n" );

    return;
}

#----------------------------#
# dialogs
#----------------------------#
sub dialog_taxon {
    my $self = shift;

    # Create the dialog and init two buttons
    my $dialog = Gtk3::Dialog->new(
        "Load taxon...",
        $self->win, [qw{modal destroy-with-parent}],
        'gtk-ok'    => 'ok',
        'gtk-close' => 'close',
    );

    my $label = Gtk3::Label->new("Choose a taxon:");
    $label->set_alignment( 0, 0.5 );
    $label->set_justify('left');

    # put common names to listmodel
    my $model = Gtk3::ListStore->new('Glib::String');
    my @common_names;
    for my $line ( path( $FindBin::RealBin, "..", "data", "chr_length.csv" )->lines ) {
        my ($name) = split( /,/, $line );
        push @common_names, $name;

    }
    @common_names = uniq(@common_names);
    for (@common_names) {

        # The iter is a pointer in the treestore. We use to add data.
        my $iter = $model->append;
        $model->set( $iter, 0 => $_, );
    }

    my $treeview = Gtk3::TreeView->new_with_model($model);

    # Add columns
    $treeview->insert_column_with_attributes( 0, "name", Gtk3::CellRendererText->new, text => 0 );
    $treeview->get_column(0)->set_sort_column_id(0);

    # get the Gtk3::TreeSelection of $treeview
    my $treeselection = $treeview->get_selection;
    $treeselection->set_mode('single');

    # add TreeView to a scrolledwindow
    my $sw = Gtk3::ScrolledWindow->new;
    $sw->set_size_request( 360, 240 );
    $sw->set_policy( 'automatic', 'automatic' );
    $sw->add($treeview);

    # Create a vbox
    my $vbox = Gtk3::Box->new( 'vertical', 5 );
    $vbox->pack_start( $label, FALSE, FALSE, 5 );
    $vbox->pack_start( $sw,    TRUE,  TRUE,  0 );
    $dialog->get_content_area()->add($vbox);
    $dialog->show_all;

    # get response
    my $response = $dialog->run;
    my $name;
    if ( $response eq 'ok' ) {
        my $iter = $treeselection->get_selected;
        if ( defined $iter ) {

            # we want data at the model's columns where the iter is pointing
            $name = $model->get( $iter, 0 );
        }
    }
    $dialog->destroy;

    return $name;
}

sub dialog_choose_db {
    my $self = shift;

    # Create the dialog and init two buttons
    my $dialog = Gtk3::Dialog->new(
        "Choose DB...",
        $self->win, [qw{modal destroy-with-parent}],
        'gtk-ok'    => 'ok',
        'gtk-close' => 'close',
    );

    my $label = Gtk3::Label->new("Choose an AlignDB database:");
    $label->set_alignment( 0, 0.5 );
    $label->set_justify('left');

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");

    my @dbs = grep {/vs/}
        map { s/dbi\:mysql\://i; $_ }
        DBI->data_sources( "mysql",
        { host => $server, port => $port, user => $username, password => $password } );

    my $model = Gtk3::ListStore->new('Glib::String');
    for (@dbs) {
        my $iter = $model->append;
        $model->set( $iter, 0 => $_, );
    }

    my $treeview = Gtk3::TreeView->new_with_model($model);

    # Add columns
    $treeview->insert_column_with_attributes( 0, "DB name", Gtk3::CellRendererText->new,
        text => 0 );
    $treeview->get_column(0)->set_sort_column_id(0);

    # get the Gtk3::TreeSelection of $treeview
    my $treeselection = $treeview->get_selection;
    $treeselection->set_mode('single');

    # add TreeView to a scrolledwindow
    my $sw = Gtk3::ScrolledWindow->new;
    $sw->set_size_request( 360, 240 );
    $sw->set_policy( 'automatic', 'automatic' );
    $sw->add($treeview);

    # Create a vbox
    my $vbox = Gtk3::Box->new( 'vertical', 5 );
    $vbox->pack_start( $label, FALSE, FALSE, 5 );
    $vbox->pack_start( $sw,    TRUE,  TRUE,  0 );
    $dialog->get_content_area()->add($vbox);
    $dialog->show_all;

    # get response
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

        my ( $target_name, $query_name ) = $obj->get_names;

        $self->append_text("choose db succeeded\n");
        return {
            target_name => $target_name,
            query_name  => $query_name,
            db_name     => $db_name,
        };
    }
    else {
        $self->append_text("choose db closed\n");
        return;
    }
}

sub dialog_db_meta {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");

    my $db_name = $self->get_value("entry_db_name");

    # Create the dialog and init two buttons
    my $dialog = Gtk3::Dialog->new(
        "Steps done",
        $self->win, [qw{modal destroy-with-parent}],
        'gtk-ok'    => 'ok',
        'gtk-close' => 'close',
    );

    my $model = Gtk3::ListStore->new( 'Glib::Int', 'Glib::String', 'Glib::String' );

    # retrieve data for $model
    {
        # get dbh
        my $obj = AlignDB->new(
            mysql  => "$db_name:$server",
            user   => $username,
            passwd => $password,
        );
        my $dbh = $obj->dbh;

        # query table meta
        my $sql = q{
            SELECT meta_value
            FROM meta
            WHERE meta_key LIKE "a_%"
               OR meta_key LIKE "c_%"
        };
        my $sth = $dbh->prepare($sql);
        $sth->execute;

        # get meta info
        my $serial;
        while ( my ($operation) = $sth->fetchrow_array ) {
            $serial++;
            my ($time) = $sth->fetchrow_array;

            # The iter is a pointer in the treestore. We use to add data.
            my $iter = $model->append;
            $model->set(
                $iter,
                0 => $serial,
                1 => $operation,
                2 => $time,
            );
        }
    }

    my $treeview = Gtk3::TreeView->new_with_model($model);

    # Add columns
    $treeview->insert_column_with_attributes( 0, "serial", Gtk3::CellRendererText->new, text => 0 );
    $treeview->insert_column_with_attributes( 1, "operation", Gtk3::CellRendererText->new,
        text => 1 );
    $treeview->insert_column_with_attributes( 2, "time", Gtk3::CellRendererText->new, text => 2 );
    $treeview->get_column(0)->set_sort_column_id(0);
    $treeview->get_column(1)->set_sort_column_id(1);
    $treeview->get_column(2)->set_sort_column_id(2);

    # get the Gtk3::TreeSelection of $treeview
    my $treeselection = $treeview->get_selection;
    $treeselection->set_mode('single');

    # add TreeView to a scrolledwindow
    my $sw = Gtk3::ScrolledWindow->new;
    $sw->set_size_request( 360, 240 );
    $sw->set_policy( 'automatic', 'automatic' );
    $sw->add($treeview);

    # Create a vbox
    my $vbox = Gtk3::Box->new( 'vertical', 5 );
    $vbox->pack_start( $sw, TRUE, TRUE, 0 );
    $dialog->get_content_area()->add($vbox);
    $dialog->show_all;

    # get response
    my $response = $dialog->run;
    my ( $operation, $time );
    if ( $response eq 'ok' ) {
        my $iter = $treeselection->get_selected;
        if ( defined $iter ) {

            # we want data at the model's columns where the iter is pointing
            $operation = $model->get( $iter, 1 );
            $time      = $model->get( $iter, 2 );
        }
    }
    $dialog->destroy;

    return ( $operation, $time );
}

#----------------------------#
# button callback events
#----------------------------#
sub on_button_db_meta_clicked {
    my $self = shift;

    my ( $operation, $time ) = $self->dialog_db_meta;

    if ($operation) {
        $self->append_text("operation:\t");
        $self->append_text( "$operation\n", "bold", "blue" );
        $self->append_text("time:\t\t");
        $self->append_text( "$time\n", "bold", "blue" );
    }
    return;
}

sub on_button_load_target_clicked {
    my $self = shift;

    my $name = $self->dialog_taxon;

    if ( defined $name ) {
        $self->set_value( "entry_target_name", $name );
    }

    return;
}

sub on_button_load_query_clicked {
    my $self = shift;

    my $name = $self->dialog_taxon;

    if ( defined $name ) {
        $self->set_value( "entry_query_name", $name );
    }

    return;
}

sub on_button_auto_db_name_clicked {
    my $self = shift;

    my $active = $self->get_value("notebook_generate");
    my @pages = ( "Gen. axt", "Gen. fas", "Join dbs", "GC", "Feature" );

    $self->append_text("Active page is [$pages[$active]]\n");

    if ( $active == 0 ) {
        my $target_name = $self->get_value("entry_target_name");
        my $query_name  = $self->get_value("entry_query_name");
        my $db_name     = "$target_name" . "vs" . "$query_name";

        $self->set_value( "entry_db_name", $db_name );
        $self->append_text("db_name set to ");
        $self->append_text( "$db_name\n", "bold", "blue" );
    }
    elsif ( $active == 1 ) {
        my $dir_align_fas = $self->get_value("entry_dir_align_fas");
        my $db_name       = path($dir_align_fas)->basename;

        $self->set_value( "entry_db_name", $db_name );
        $self->append_text("db_name set to ");
        $self->append_text( "$db_name\n", "bold", "blue" );
    }
    elsif ( $active == 2 ) {
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
        $goal_db = $name_of{$first} . 'vs' . $name_of{$second} . 'ref' . $name_of{$outgroup};

        $self->set_value( "entry_db_name", $goal_db );
        $self->append_text("db_name set to ");
        $self->append_text( "$goal_db\n", "bold", "blue" );
    }
    else {
        $self->append_text("Can't set db_name on page [$pages[$active]] of Generate.\n");
    }

    return;
}

sub on_button_choose_db_clicked {
    my $self = shift;

    my $result = $self->dialog_choose_db;
    return unless $result;

    $self->set_value( "entry_db_name", $result->{db_name} );

    return;
}

sub on_button_open_dir_align_axt_clicked {
    my $self = shift;

    my $dia = Gtk3::FileChooserDialog->new(
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
                    $self->set_value( "entry_dir_align_axt", $dir );
                }
            }
            $dia->destroy;
            return FALSE;
        }
    );

    return;
}

sub on_button_open_dir_align_fas_clicked {
    my $self = shift;

    my $dia = Gtk3::FileChooserDialog->new(
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
                    $self->set_value( "entry_dir_align_fas", $dir );
                }
            }
            $dia->destroy;
            return FALSE;
        }
    );

    return;
}

sub on_button_auto_stat_file_common_clicked {
    my $self = shift;

    my $db_name = $self->get_value("entry_db_name");
    my $outfile
        = path( $FindBin::RealBin, '..', 'stat', "$db_name.common.xlsx" )->absolute->stringify;
    $self->set_value( "entry_stat_file_common", $outfile );

    return;
}

sub on_button_auto_stat_file_gc_clicked {
    my $self = shift;

    my $db_name = $self->get_value("entry_db_name");
    my $outfile = path( $FindBin::RealBin, '..', 'stat', "$db_name.gc.xlsx" )->absolute->stringify;
    $self->set_value( "entry_stat_file_gc", $outfile );

    return;
}

sub on_button_auto_stat_file_multi_clicked {
    my $self = shift;

    my $db_name = $self->get_value("entry_db_name");
    my $outfile
        = path( $FindBin::RealBin, '..', 'stat', "$db_name.multi.xlsx" )->absolute->stringify;
    $self->set_value( "entry_stat_file_multi", $outfile );

    return;
}

sub on_button_choose_first_db_clicked {
    my $self = shift;

    my $result = $self->dialog_choose_db;
    return unless $result;

    $self->set_value( "entry_first_db", $result->{db_name} );
    return;
}

sub on_button_choose_second_db_clicked {
    my $self = shift;

    my $result = $self->dialog_choose_db;
    return unless $result;

    $self->set_value( "entry_second_db", $result->{db_name} );
    return;
}

sub on_button_clear_clicked {
    my $self = shift;

    my $text = $self->text;
    $text->delete( $text->get_start_iter, $text->get_end_iter );

    return;
}

#----------------------------------------------------------#
# run *.pl events
#----------------------------------------------------------#
sub on_button_test_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");

    my $dbh = AlignDB->new(
        mysql  => "mysql:$server",
        user   => $username,
        passwd => $password,
    )->dbh;
    my @row_ary = $dbh->selectrow_array("SELECT 123456789");
    print YAML::Syck::Dump @row_ary;

    my $success = grep {/123456789/} @row_ary;

    if ($success) {
        $self->append_text("Test connection successful\n");
    }
    else {
        $self->append_text("Test connection failed\n");
    }

    return;
}

sub on_button_init_aligndb_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $cmd
        = "perl $FindBin::RealBin/../init/init_alignDB.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_gen_aligndb_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $length_threshold = $self->get_value("entry_length_threshold");

    my $target_name = $self->get_value("entry_target_name");
    my $query_name  = $self->get_value("entry_query_name");
    my $dir_align   = $self->get_value("entry_dir_align_axt");

    my $cmd
        = "perl $FindBin::RealBin/../init/gen_alignDB.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -t $target_name"
        . " -q $query_name"
        . " -dir $dir_align"
        . " -lt $length_threshold"
        . " --parallel $parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_gen_aligndb_fas_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $length_threshold = $self->get_value("entry_length_threshold");
    my $parallel         = $self->get_value("entry_parallel");

    my $dir_align = $self->get_value("entry_dir_align_fas");

    my $outgroup = $self->get_value("checkbutton_fas_outgroup");

    my $cmd
        = "perl $FindBin::RealBin/../init/gen_alignDB_fas.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " --db $db_name"
        . " --dir $dir_align"
        . " --length $length_threshold"
        . " --parallel $parallel"
        . ( $outgroup ? " --outgroup" : "" );

    $self->exec_cmd($cmd);
    return;
}

sub on_button_insert_isw_axt_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $cmd
        = "perl $FindBin::RealBin/../init/insert_isw.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_insert_isw_fas_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $outgroup = $self->get_value("checkbutton_fas_outgroup");

    my $cmd
        = "perl $FindBin::RealBin/../init/insert_isw.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel"
        . ( $outgroup ? " --outgroup" : "" );

    $self->exec_cmd($cmd);
    return;
}

sub on_button_insert_gc_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $insert_gc      = $self->get_value("checkbutton_insert_gc");
    my $insert_segment = $self->get_value("checkbutton_insert_segment");

    my $cmd
        = "perl $FindBin::RealBin/../init/insert_gc.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --insert_gc $insert_gc"
        . " --insert_segment $insert_segment"
        . " --parallel $parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_swcv_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $parallel = $self->get_value("entry_parallel");

    my $cmd
        = "perl $FindBin::RealBin/../init/update_sw_cv.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_feature_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $ensembl = $self->get_value("entry_ensembl");

    my $parallel = $self->get_value("entry_parallel");

    my $cmd
        = "perl $FindBin::RealBin/../init/update_feature.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -e $ensembl"
        . " --parallel $parallel";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_upd_slippage_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $cmd
        = "perl $FindBin::RealBin/../init/update_indel_slippage.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_stat_common_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $output = $self->get_value("entry_stat_file_common");

    my $cmd
        = "perl $FindBin::RealBin/../stat/common_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $output"
        . " -r all --index --chart";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_stat_multi_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $output = $self->get_value("entry_stat_file_multi");

    my $cmd
        = "perl $FindBin::RealBin/../stat/multi_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $output"
        . " -r all --index --chart";

    $self->exec_cmd($cmd);
    return;
}

sub on_button_stat_gc_clicked {
    my $self = shift;

    my $server   = $self->get_value("entry_server");
    my $port     = $self->get_value("entry_port");
    my $username = $self->get_value("entry_username");
    my $password = $self->get_value("entry_password");
    my $db_name  = $self->get_value("entry_db_name");

    my $output = $self->get_value("entry_stat_file_gc");

    my $cmd
        = "perl $FindBin::RealBin/../stat/gc_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $output"
        . " -r all --index --chart";

    $self->exec_cmd($cmd);
    return;
}

package main;

AlignDB::GUI->new;

exit;
