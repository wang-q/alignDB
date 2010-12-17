package AlignDB::GUI::Frame;
use strict;
use warnings;

use Wx qw(:id :filedialog :dirdialog);
use Cwd;
use File::Spec;

use FindBin;
use lib "$FindBin::Bin/../lib";
require AlignDB;
require AlignDB::GUI::TaxonDialog;

#----------------------------#
# callback events
#----------------------------#
sub about_dialog {
    my $self = shift;
    my $info = Wx::AboutDialogInfo->new;

    $info->SetName('alignDB GUI2');
    $info->SetVersion('0.2.0');
    $info->SetDescription("A GUI shell for alignDB scripts.");
    $info->SetCopyright("(C) 2004-2007 Wang Qiang");
    $info->SetLicense( "This program is free software;\n"
            . "you can redistribute it and/or modify\n"
            . "it under the same terms as Perl itself.\n" );
    $info->SetWebSite( 'http://where.can.I.find.you/',
        'We have no website yet' );
    $info->AddDeveloper('Wang Qiang <wangqiang1997@gmail.com>');
    $info->AddDeveloper('Zhu Liucun');
    $info->AddArtist('Wang Qiang <wangqiang1997@gmail.com>');
    $info->AddArtist('Looks ugly? Help me!');
    $info->AddDocWriter('No documnets! Help me!');

    Wx::AboutBox($info);
    return;
}

sub event_list_process {
    my $self  = shift;
    my $event = shift;

    use YAML qw{Dump};
    
    my $processes = $self->processes;
    return unless ref $processes eq 'ARRAY';
    for my $proc (@$processes) {
        print Dump {
            start => $proc->start_time,
            end => $proc->end_time,
            alive => $proc->alive,
            pid => $proc->pid,
            cmd => $proc->{cmd},
        };
    }
    
    return;
}

sub event_test_connect {
    my $self  = shift;
    my $event = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $cmd
        = "mysql"
        . " -h$server"
        . " -P$port"
        . " -u$username"
        . " -p$password"
        . q{ -e "SELECT 123456789"};
    my $success = grep {/123456789/} split "\n", `$cmd`;

    if ($success) {
        print "Test connection successful\n\n";
    }
    else {
        print "Test connection failed\n\n";
    }

    return;
}

sub event_load_target {
    my $self  = shift;
    my $event = shift;

    my $dialog = AlignDB::GUI::TaxonDialog->new;
    $dialog->build_window("$FindBin::Bin/normal_taxon.csv");

    if ( $dialog->ShowModal == wxID_OK ) {
        $self->set_value( "target_taxon_id", $dialog->id );
        $self->set_value( "target_name",     $dialog->name );
    }

    $dialog->Destroy;
    return;
}

sub event_load_query {
    my $self  = shift;
    my $event = shift;

    my $dialog = AlignDB::GUI::TaxonDialog->new;
    $dialog->build_window("$FindBin::Bin/normal_taxon.csv");

    if ( $dialog->ShowModal == wxID_OK ) {
        $self->set_value( "query_taxon_id", $dialog->id );
        $self->set_value( "query_name",     $dialog->name );
    }

    $dialog->Destroy;
    return;
}

sub choose_db {
    my $self  = shift;

    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $cmd
        = "mysql"
        . " -h$server"
        . " -P$port"
        . " -u$username"
        . " -p$password"
        . ' -e "SHOW DATABASES"';
    my @dbs = grep {/vs/} split "\n", `$cmd`;

    my $dialog
        = Wx::SingleChoiceDialog->new( $self, "Choose an alignDB database",
        "Choose DB", [@dbs], );
        
    my $icon = Wx::Icon->new;
    $icon->CopyFromBitmap( AlignDB::GUI::Bitmap->get_bitmap_db );
    $dialog->SetIcon($icon);

    if ( $dialog->ShowModal == wxID_OK ) {
        my $selected_db = $dialog->GetStringSelection;

        my $obj = AlignDB->new(
            mysql  => "$selected_db:$server",
            user   => $username,
            passwd => $password,
        );

        my ( $target_taxon_id, $query_taxon_id ) = $obj->get_taxon_ids;
        my ( $target_name,     $query_name )     = $obj->get_names;

        return {
            target_taxon_id => $target_taxon_id,
            query_taxon_id  => $query_taxon_id,
            target_name     => $target_name,
            query_name      => $query_name,
            db_name         => $selected_db,
        };
    }

    return;
}

sub event_choose_db {
    my $self  = shift;
    my $event = shift;

    my $result = $self->choose_db;
    
    return unless $result;

    $self->set_value( "target_taxon_id", $result->{target_taxon_id} );
    $self->set_value( "query_taxon_id",  $result->{query_taxon_id} );
    $self->set_value( "target_name",     $result->{target_name} );
    $self->set_value( "query_name",      $result->{query_name} );
    $self->set_value( "db_name",         $result->{db_name} );

    return;
}

sub event_choose_first_db {
    my $self  = shift;
    my $event = shift;

    my $result = $self->choose_db;
    $self->set_value( "first_db", $result->{db_name} );

    return;
}

sub event_choose_second_db {
    my $self  = shift;
    my $event = shift;

    my $result = $self->choose_db;
    $self->set_value( "second_db", $result->{db_name} );
    
    return;
}

sub event_choose_goal_db {
    my $self  = shift;
    my $event = shift;

    my $result = $self->choose_db;
    $self->set_value( "goal_db", $result->{db_name} );
    
    return;
}

sub event_auto_db_name {
    my $self  = shift;
    my $event = shift;

    my $target_name = $self->get_value("target_name");
    my $query_name  = $self->get_value("query_name");
    $self->set_value( "db_name", "$target_name" . "vs" . "$query_name" );
    return;
}

sub event_auto_goal_db_name {
    my $self  = shift;
    my $event = shift;
    
    my $first_db = $self->get_value("first_db");
    my $second_db = $self->get_value("second_db");
    my $goal_db;
    
    my $first = $self->get_value("first");
    my $second = $self->get_value("second");
    my $outgroup = $self->get_value("outgroup");
    
    my $server   = $self->get_value("server");
    my $port     = $self->get_value("port");
    my $username = $self->get_value("username");
    my $password = $self->get_value("password");

    my $first_obj = AlignDB->new(
        mysql  => "$first_db:$server",
        user   => $username,
        passwd => $password,
    );
    my ( $first_target_name,     $first_query_name )     = $first_obj->get_names;

    my $second_obj = AlignDB->new(
        mysql  => "$second_db:$server",
        user   => $username,
        passwd => $password,
    );
    my ( $second_target_name,     $second_query_name )     = $second_obj->get_names;
    
    my %name_of = (
        '0target' => $first_target_name,
        '0query' => $first_query_name,
        '1target' => $second_target_name,
        '1query' => $second_query_name,
    );
    $goal_db = $name_of{$first} . 'vs' . $name_of{$second} . 'ref' . $name_of{$outgroup};
    
    $self->set_value( "goal_db", $goal_db);
    return;
}

sub event_choose_aim_db {
    my $self  = shift;
    my $event = shift;

    my $result = $self->choose_db;
    $self->set_value( "aim_db", $result->{db_name} );

    return;
}

sub event_choose_ref_db {
    my $self  = shift;
    my $event = shift;

    my $result = $self->choose_db;
    $self->set_value( "ref_db", $result->{db_name} );

    return;
}

sub event_open_axt_dir {
    my $self  = shift;
    my $event = shift;

    my $dialog
        = Wx::DirDialog->new( $self, "Select a dir",
        $self->previous_directory || cwd,
        );

    if ( $dialog->ShowModal == wxID_OK ) {
        my $dir = $dialog->GetPath;
        $self->set_value( "axt_dir", $dir );
        $self->previous_directory($dir);
    }

    $dialog->Destroy;
    return;
}

sub event_open_sql_file {
    my $self  = shift;
    my $event = shift;

    my $dialog = Wx::FileDialog->new(
        $self,
        "Select a file",
        $self->previous_directory || cwd,
        '',
        (   join '|',
            'SQL files (*.sql)|*.sql',
            'Text files (*.txt)|*.txt',
            'All files (*.*)|*.*'
        ),
        wxFD_OPEN | wxFD_MULTIPLE
    );

    if ( $dialog->ShowModal == wxID_OK ) {
        my ($file) = $dialog->GetPaths;    # May return multiple files
        $self->set_value( "sql_file", $file );
        $self->previous_directory( $dialog->GetDirectory );
    }

    $dialog->Destroy;
    return;
}

sub event_open_output_file {
    my $self  = shift;
    my $event = shift;

    my @filters = (
        'CSV files (*.csv)|*.csv',
        'Neat files (*.neat)|*.neat',
        'Table files (*.table)|*.table',
        'Box files (*.box)|*.box',
        'HTML files (*.html)|*.html',
    );
    my $dialog = Wx::FileDialog->new(
        $self,
        "Select a file",
        $self->previous_directory || cwd,
        '',
        ( join '|', @filters ),
        wxFD_OPEN | wxFD_MULTIPLE
    );

    if ( $dialog->ShowModal == wxID_OK ) {
        my ($file) = $dialog->GetPaths;    # May return multiple files
        $self->set_value( "query_output", $file );
        $self->previous_directory( $dialog->GetDirectory );

        # Returns the index into the list of filters
        # Transform it to query_sql.pl types
        my $filter_index = $dialog->GetFilterIndex;
        my $filter       = $filters[$filter_index];
        $filter =~ s/^(\w+).+/$1/;
        $filter = lc $filter;
        $self->query_output_type($filter);
    }

    $dialog->Destroy;
    return;
}

sub event_auto_common_stat_file {
    my $self  = shift;
    my $event = shift;

    my $db_name = $self->get_value("db_name");
    my $outfile = "$FindBin::Bin/../stat/$db_name.common.xls";
    $outfile = File::Spec->rel2abs($outfile);
    $self->set_value( "common_stat_file", $outfile );

    return;
}

sub event_auto_gc_stat_file {
    my $self  = shift;
    my $event = shift;

    my $db_name = $self->get_value("db_name");
    my $outfile = "$FindBin::Bin/../stat/$db_name.gc.xls";
    $outfile = File::Spec->rel2abs($outfile);
    $self->set_value( "gc_stat_file", $outfile );

    return;
}

sub event_auto_gene_stat_file {
    my $self  = shift;
    my $event = shift;

    my $db_name = $self->get_value("db_name");
    my $outfile = "$FindBin::Bin/../stat/$db_name.gene.xls";
    $outfile = File::Spec->rel2abs($outfile);
    $self->set_value( "gene_stat_file", $outfile );

    return;
}

sub event_auto_three_stat_file {
    my $self  = shift;
    my $event = shift;

    my $goal_db = $self->get_value("goal_db");
    my $outfile = "$FindBin::Bin/../stat/$goal_db.three.xls";
    $outfile = File::Spec->rel2abs($outfile);
    $self->set_value( "three_stat_file", $outfile );

    return;
}

1;
