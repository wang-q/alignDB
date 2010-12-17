package AlignDB::GUI::Frame;
use strict;
use warnings;

use Wx qw(
    :sizer :textctrl :staticline :id :font :color :combobox :toolbar
);
use Wx::Event qw(
    EVT_BUTTON EVT_CLOSE EVT_MENU EVT_COLLAPSIBLEPANE_CHANGED EVT_CHOICE
    EVT_TOOL_ENTER
);
use base qw(Wx::Frame);
use base qw(Class::Accessor);

use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
require AlignDB::GUI::Bitmap;
require AlignDB::GUI::Frame_cmd;
require AlignDB::GUI::Frame_config;
require AlignDB::GUI::Frame_event;

# these attributes are all wxWidgets object, which should be accessed by
# set_value and get_value
__PACKAGE__->mk_accessors(
    qw(
        server port username password
        target_taxon_id target_name query_taxon_id query_name
        db_name ensembl
        first_db second_db goal_db first second outgroup
        axt_dir insert_dG axt_threshold parallel
        insert_gc insert_segment
        insert_genesw insert_exonsw insert_codingsw gene_ensembl
        process_align process_indel process_isw process_snp process_window
        aim_db ref_db
        sql_file query_sql query_output
        )
);

# stat related attributes
__PACKAGE__->mk_accessors(
    qw(
        common_stat_file common_run common_threshold
        common_jc common_time_stamp common_add_index_sheet
        gc_stat_file gc_run gc_threshold
        gc_jc gc_time_stamp gc_add_index_sheet
        gene_stat_file gene_run
        gene_jc gene_time_stamp gene_add_index_sheet
        three_stat_file three_run three_threshold
        three_jc three_time_stamp three_add_index_sheet
        )
);

# these attributes should be accessed directly
__PACKAGE__->mk_accessors(
    qw(
        previous_directory query_output_type
        processes
        )
);

sub new {
    my $ref  = shift;
    my $self = $ref->SUPER::new(
        undef,             # parent window
        -1,                # ID -1 means any
        'alignDB GUI2',    # title
        [ -1, -1 ],        # default position
        [ -1, -1 ],        # size
    );

    $self->build_icon;
    $self->build_menu;
    $self->build_toolbar;
    $self->build_window;
    $self->read_config;

    EVT_CLOSE( $self, \&event_close_window );

    return $self;
}

sub build_icon {
    my $self = shift;

    my $icon = Wx::Icon->new;
    $icon->CopyFromBitmap( AlignDB::GUI::Bitmap->get_bitmap_main );
    $self->SetIcon($icon);

    return;
}

sub build_menu {
    my $self = shift;

    my $file_menu = Wx::Menu->new;
    my $help_menu = Wx::Menu->new;

    # Using these special contants will automatically use the doc/view
    # The EVT_MENU sets events called from clicking on menu items
    EVT_MENU( $self, $file_menu->Append( wxID_DEFAULT, "&Default" ),
        \&read_config );
    EVT_MENU( $self, $file_menu->Append( wxID_REFRESH, "&Process" ),
        \&event_list_process );
    $file_menu->AppendSeparator;
    EVT_MENU( $self, $file_menu->Append( wxID_EXIT, "E&xit" ),
        \&event_close_window );

    $help_menu->AppendSeparator;
    EVT_MENU( $self, $help_menu->Append( wxID_ABOUT, "&About" ),
        \&about_dialog );

    my $menu_bar = Wx::MenuBar->new;
    $menu_bar->Append( $file_menu, "&File" );
    $menu_bar->Append( $help_menu, "&Help" );

    $self->SetMenuBar($menu_bar);

    return;
}

sub build_toolbar {
    my $self = shift;

    my $style = wxTB_HORIZONTAL | wxTB_FLAT | wxTB_DOCKABLE | wxNO_BORDER;

    my $toolbar = $self->CreateToolBar( $style, -1 );
    $toolbar->SetMargins( 4, 4 );
    $toolbar->SetToolBitmapSize( [ 20, 20 ] );

    $toolbar->AddTool( wxID_DEFAULT, '',
        AlignDB::GUI::Bitmap->get_bitmap_default, 'Default'
    );
    $toolbar->AddTool( wxID_REFRESH, '',
        AlignDB::GUI::Bitmap->get_bitmap_process, 'Process'
    );

    {
        $toolbar->AddSeparator;
        my $static_text = Wx::StaticText->new( $toolbar, -1, "parallel: " );
        $toolbar->AddControl($static_text);

        my $text_ctrl
            = Wx::TextCtrl->new( $toolbar, -1, '', [ -1, -1 ], [ 30, -1 ],
            0 );
        $toolbar->AddControl($text_ctrl);
        $self->{parallel} = $text_ctrl;
    }

    $toolbar->AddSeparator;
    $toolbar->AddTool( wxID_ABOUT, '', AlignDB::GUI::Bitmap->get_bitmap_info,
        'About' );

    $toolbar->Realize;
    $self->SetToolBar($toolbar);

    return;
}

#----------------------------#
# pane: Database
#----------------------------#
sub build_database_pane {
    my $self         = shift;
    my $parent_panel = shift;
    my $parent_sizer = shift;

    my ( $pane, $nb )
        = $self->add_collapsible_nb( $parent_panel, $parent_sizer,
        'Database' );

    {    # DB Server
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "DB Server" );

        $sizer->AddStretchSpacer;

        {    # use GridBagSizer
            my $gb_sizer = Wx::GridBagSizer->new;
            $gb_sizer->AddGrowableCol(1);
            $sizer->Add( $gb_sizer, 0, wxGROW, 0 );

            # server and port
            $self->add_gb_static_text( $panel, $gb_sizer, "server:",
                [ 0, 0 ] );
            $self->add_gb_text_ctrl( $panel, $gb_sizer, "server", [ 0, 1 ] );

            $self->add_gb_static_text( $panel, $gb_sizer, "port:", [ 0, 2 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "port",
                [ 0,  3 ],
                [ 1,  1 ],
                [ 50, -1 ]
            );

            # username and password
            $self->add_gb_static_text( $panel, $gb_sizer, "username:",
                [ 1, 0 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "username",
                [ 1, 1 ],
                [ 1, 3 ]
            );

            $self->add_gb_static_text( $panel, $gb_sizer, "password:",
                [ 2, 0 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "password",
                [ 2, 1 ],
                [ 1, 3 ]
            );

        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Test",
                \&event_test_connect );
        }
    }

    {        # Initiate
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Initiate" );

        $sizer->AddStretchSpacer;

        {    # use GridBagSizer
            my $gb_sizer = Wx::GridBagSizer->new;
            $gb_sizer->AddGrowableCol(3);
            $sizer->Add( $gb_sizer, 0, wxGROW, 0 );

            # target
            $self->add_gb_static_text( $panel, $gb_sizer, "target id:",
                [ 0, 0 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "target_taxon_id",
                [ 0,  1 ],
                [ 1,  1 ],
                [ 80, -1 ]
            );

            $self->add_gb_static_text( $panel, $gb_sizer, "name:", [ 0, 2 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "target_name",
                [ 0, 3 ],
                [ 1, 2 ]
            );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_load,
                \&event_load_target, [ 0, 5 ],
            );

            # query
            $self->add_gb_static_text( $panel, $gb_sizer, "query id:",
                [ 1, 0 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "query_taxon_id",
                [ 1,  1 ],
                [ 1,  1 ],
                [ 80, -1 ]
            );

            $self->add_gb_static_text( $panel, $gb_sizer, "name:", [ 1, 2 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "query_name",
                [ 1, 3 ],
                [ 1, 2 ]
            );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_load,
                \&event_load_query, [ 1, 5 ],
            );

            # db name
            $self->add_gb_static_text( $panel, $gb_sizer, "db name:",
                [ 2, 0 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "db_name",
                [ 2, 1 ],
                [ 1, 3 ]
            );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_auto,
                \&event_auto_db_name, [ 2, 4 ],
            );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_db,
                \&event_choose_db, [ 2, 5 ],
            );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Init. alignDB",
                \&event_init_alignDB );
        }
    }

    return ( $pane, $nb );
}

#----------------------------#
# pane: Generate
#----------------------------#
sub build_generate_pane {
    my $self         = shift;
    my $parent_panel = shift;
    my $parent_sizer = shift;

    my ( $pane, $nb )
        = $self->add_collapsible_nb( $parent_panel, $parent_sizer,
        'Generate' );

    {    # generate
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Generate" );

        $sizer->AddStretchSpacer;

        {    # axt dir
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, ".axt dir:" );
            $self->add_text_ctrl( $panel, $boxsizer, "axt_dir", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_bitmap_button( $panel, $boxsizer,
                AlignDB::GUI::Bitmap->get_bitmap_open,
                \&event_open_axt_dir );
        }

        {    # dG, threshold and parallel
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "insert_dG", "DeltaG" );

            $boxsizer->AddSpacer(10);

            $self->add_static_text( $panel, $boxsizer, "threshold:" );
            $self->add_text_ctrl( $panel, $boxsizer, "axt_threshold", 1 );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Gen. alignDB",
                \&event_gen_alignDB );
        }
    }

    {        # insert GC
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Insert GC" );

        $sizer->AddStretchSpacer;

        {    # insert_gc and insert_segment
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "insert_gc", "GC" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "insert_segment",
                "Segment" );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Ins. GC",
                \&event_insert_gc );
        }
    }

    {        # insert gene
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Insert Gene" );

        $sizer->AddStretchSpacer;

        {    # insert_genesw, insert_exonsw and insert_codingsw
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "insert_genesw",
                "GeneSW" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "insert_exonsw",
                "ExonSW" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "insert_codingsw",
                "CodingSW" );
        }

        {    # ensembl
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "ensembl:" );
            $self->add_text_ctrl( $panel, $boxsizer, "gene_ensembl" );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Ins. Gene",
                \&event_insert_gene );
        }
    }

    return ( $pane, $nb );
}

#----------------------------#
# pane: Update
#----------------------------#
sub build_update_pane {
    my $self         = shift;
    my $parent_panel = shift;
    my $parent_sizer = shift;

    my ( $pane, $nb )
        = $self->add_collapsible_nb( $parent_panel, $parent_sizer, 'Update' );

    {    # update feature
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Feature" );

        $sizer->AddStretchSpacer;

        {    # insert_genesw, insert_exonsw and insert_codingsw
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "process_align",
                "align" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "process_indel",
                "indel" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "process_isw", "isw" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "process_snp", "snp" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "process_window",
                "window" );
        }

        {    # ensembl db
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "ensembl:" );
            $self->add_text_ctrl( $panel, $boxsizer, "ensembl" );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Upd. feature",
                \&event_update_feature );
        }
    }

    {        # update indel slippage, isw indel id and segment extreme
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Misc" );

        $sizer->AddStretchSpacer;

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            # indel slippage
            $boxsizer->AddStretchSpacer(1);
            $self->add_button(
                $panel, $boxsizer,
                "Upd. slippage",
                \&event_update_indel_slippage
            );

            # isw indel id
            $boxsizer->AddSpacer(10);
            $self->add_button( $panel, $boxsizer, "Upd. isw",
                \&event_update_isw_indel_id );

            # segment extreme
            $boxsizer->AddSpacer(10);
            $self->add_button( $panel, $boxsizer, "Upd. segment",
                \&event_update_segment );
        }
    }

    return ( $pane, $nb );
}

#----------------------------#
# pane: Stat
#----------------------------#
sub build_stat_pane {
    my $self         = shift;
    my $parent_panel = shift;
    my $parent_sizer = shift;

    my ( $pane, $nb )
        = $self->add_collapsible_nb( $parent_panel, $parent_sizer, 'Stat' );

    {    # common
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Common" );

        $sizer->AddStretchSpacer;

        {    # run and threshold
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "run:" );
            $self->add_text_ctrl( $panel, $boxsizer, "common_run", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_static_text( $panel, $boxsizer, "threshold:" );
            $self->add_text_ctrl( $panel, $boxsizer, "common_threshold" );
        }

        {    # stat file
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "stat file:" );
            $self->add_text_ctrl( $panel, $boxsizer, "common_stat_file", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_bitmap_button(
                $panel, $boxsizer,
                AlignDB::GUI::Bitmap->get_bitmap_auto,
                \&event_auto_common_stat_file
            );
        }

        {    # jc, time_stamp and add_index_sheet
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "common_jc",
                "JC correction" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "common_time_stamp",
                "time stamp" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "common_add_index_sheet",
                "add index sheet" );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Common stat",
                \&event_common_stat );

            $boxsizer->AddSpacer(10);
            $self->add_button( $panel, $boxsizer, "Common chart",
                \&event_common_chart );
        }
    }

    {    # gc
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "GC" );

        $sizer->AddStretchSpacer;

        {    # run and threshold
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "run:" );
            $self->add_text_ctrl( $panel, $boxsizer, "gc_run", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_static_text( $panel, $boxsizer, "threshold:" );
            $self->add_text_ctrl( $panel, $boxsizer, "gc_threshold" );
        }

        {    # stat file
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "stat file:" );
            $self->add_text_ctrl( $panel, $boxsizer, "gc_stat_file", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_bitmap_button( $panel, $boxsizer,
                AlignDB::GUI::Bitmap->get_bitmap_auto,
                \&event_auto_gc_stat_file );
        }

        {    # jc, time_stamp and add_index_sheet
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "gc_jc",
                "JC correction" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "gc_time_stamp",
                "time stamp" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "gc_add_index_sheet",
                "add index sheet" );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "GC stat",
                \&event_gc_stat );

            $boxsizer->AddSpacer(10);
            $self->add_button( $panel, $boxsizer, "GC chart",
                \&event_gc_chart );
        }
    }

    {    # gene
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Gene" );

        $sizer->AddStretchSpacer;

        {    # run and threshold
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "run:" );
            $self->add_text_ctrl( $panel, $boxsizer, "gene_run", 1 );

            $boxsizer->AddSpacer(10);
        }

        {    # stat file
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "stat file:" );
            $self->add_text_ctrl( $panel, $boxsizer, "gene_stat_file", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_bitmap_button( $panel, $boxsizer,
                AlignDB::GUI::Bitmap->get_bitmap_auto,
                \&event_auto_gene_stat_file );
        }

        {    # jc, time_stamp and add_index_sheet
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "gene_jc",
                "JC correction" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "gene_time_stamp",
                "time stamp" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "gene_add_index_sheet",
                "add index sheet" );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Gene stat",
                \&event_gene_stat );

            $boxsizer->AddSpacer(10);
            $self->add_button( $panel, $boxsizer, "Gene chart",
                \&event_gene_chart );
        }
    }

    return ( $pane, $nb );
}

#----------------------------#
# pane: Util
#----------------------------#
sub build_util_pane {
    my $self         = shift;
    my $parent_panel = shift;
    my $parent_sizer = shift;

    my ( $pane, $nb )
        = $self->add_collapsible_nb( $parent_panel, $parent_sizer, 'Util' );

    {    # apply sql
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Apply SQL" );

        $sizer->AddStretchSpacer;

        {    # sql_file
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "SQL file:" );
            $self->add_text_ctrl( $panel, $boxsizer, "sql_file", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_bitmap_button( $panel, $boxsizer,
                AlignDB::GUI::Bitmap->get_bitmap_open,
                \&event_open_sql_file );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Apply SQL",
                \&event_apply_sql );
        }
    }

    {        # query sql
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Query SQL" );

        {    # sql statements
            $self->add_text_ctrl( $panel, $sizer, "query_sql", 1,
                wxTE_MULTILINE, [ -1, 100 ] );
            $self->query_sql->AppendText("SQL here ...");
        }

        {    # output
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "output:" );
            $self->add_text_ctrl( $panel, $boxsizer, "query_output", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_bitmap_button( $panel, $boxsizer,
                AlignDB::GUI::Bitmap->get_bitmap_open,
                \&event_open_output_file );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Query SQL",
                \&event_query_sql );
        }
    }

    return ( $pane, $nb );
}

#----------------------------#
# pane: Three
#----------------------------#
sub build_three_pane {
    my $self         = shift;
    my $parent_panel = shift;
    my $parent_sizer = shift;

    my ( $pane, $nb )
        = $self->add_collapsible_nb( $parent_panel, $parent_sizer,
        'Three-way' );

    {    # Three-way
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Three-way" );

        $sizer->AddStretchSpacer;

        {    # use GridBagSizer
            my $gb_sizer = Wx::GridBagSizer->new;
            $gb_sizer->AddGrowableCol(1);
            $sizer->Add( $gb_sizer, 0, wxGROW, 0 );

            # first db
            $self->add_gb_static_text( $panel, $gb_sizer, "first db:",
                [ 0, 0 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "first_db",
                [ 0, 1 ],
                [ 1, 2 ]
            );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_db,
                \&event_choose_first_db, [ 0, 3 ],
            );

            # second db
            $self->add_gb_static_text( $panel, $gb_sizer, "second db:",
                [ 1, 0 ] );
            $self->add_gb_text_ctrl(
                $panel, $gb_sizer, "second_db",
                [ 1, 1 ],
                [ 1, 2 ]
            );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_db,
                \&event_choose_second_db, [ 1, 3 ],
            );

            # goal db
            $self->add_gb_static_text( $panel, $gb_sizer, "goal db:",
                [ 2, 0 ] );
            $self->add_gb_text_ctrl( $panel, $gb_sizer, "goal_db", [ 2, 1 ] );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_auto,
                \&event_auto_goal_db_name, [ 2, 2 ],
            );
            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_db,
                \&event_choose_goal_db, [ 2, 3 ],
            );
        }

        {    # first, second and outgroup
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "T/Q/R:" );
            $boxsizer->AddSpacer(10);

            my $choices = [ '0target', '0query', '1target', '1query' ];

            # choose first
            $self->add_choice( $panel, $boxsizer, "first", $choices, 0 );
            $boxsizer->AddSpacer(10);

            # choose second
            $self->add_choice( $panel, $boxsizer, "second", $choices, 1 );
            $boxsizer->AddSpacer(10);

            # choose outgroup
            $self->add_choice( $panel, $boxsizer, "outgroup", $choices, 3 );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Ref outgroup",
                \&event_ref_outgroup );
        }
    }

    {        # update isw dxr and cpg
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Misc" );

        $sizer->AddStretchSpacer;

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);
            $boxsizer->AddStretchSpacer(1);

            # isw dxr
            $self->add_button( $panel, $boxsizer, "Upd. isw Dxr",
                \&event_update_isw_dxr );

            # cpg
            $boxsizer->AddSpacer(10);
            $self->add_button( $panel, $boxsizer, "Upd. CpG",
                \&event_update_snp_cpg );
        }
    }

    {    # update indel occured
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Occured" );

        $sizer->AddStretchSpacer;

        {    # use GridBagSizer
            my $gb_sizer = Wx::GridBagSizer->new;
            $gb_sizer->AddGrowableCol(1);
            $sizer->Add( $gb_sizer, 0, wxGROW, 0 );

            # aim db
            $self->add_gb_static_text( $panel, $gb_sizer, "aim db:",
                [ 0, 0 ] );
            $self->add_gb_text_ctrl( $panel, $gb_sizer, "aim_db", [ 0, 1 ] );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_db,
                \&event_choose_aim_db, [ 0, 2 ],
            );

            # ref db
            $self->add_gb_static_text( $panel, $gb_sizer, "ref db:",
                [ 1, 0 ] );
            $self->add_gb_text_ctrl( $panel, $gb_sizer, "ref_db", [ 1, 1 ] );

            $self->add_gb_bitmap_button(
                $panel, $gb_sizer, AlignDB::GUI::Bitmap->get_bitmap_db,
                \&event_choose_ref_db, [ 1, 2 ],
            );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Upd. occured",
                \&event_update_indel_occured );
        }
    }

    {        # three
        my ( $panel, $sizer ) = $self->add_nb_page( $nb, "Stat" );

        $sizer->AddStretchSpacer;

        {    # run and threshold
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "run:" );
            $self->add_text_ctrl( $panel, $boxsizer, "three_run", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_static_text( $panel, $boxsizer, "threshold:" );
            $self->add_text_ctrl( $panel, $boxsizer, "three_threshold" );
        }

        {    # stat file
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_static_text( $panel, $boxsizer, "stat file:" );
            $self->add_text_ctrl( $panel, $boxsizer, "three_stat_file", 1 );

            $boxsizer->AddSpacer(10);

            $self->add_bitmap_button(
                $panel, $boxsizer,
                AlignDB::GUI::Bitmap->get_bitmap_auto,
                \&event_auto_three_stat_file
            );
        }

        {    # jc, time_stamp and add_index_sheet
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $self->add_check_box( $panel, $boxsizer, "three_jc",
                "JC correction" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "three_time_stamp",
                "time stamp" );
            $boxsizer->AddSpacer(10);
            $self->add_check_box( $panel, $boxsizer, "three_add_index_sheet",
                "add index sheet" );
        }

        {    # execute button
            $self->add_static_line( $panel, $sizer );
            my $boxsizer = $self->add_boxsizer_h($sizer);

            $boxsizer->AddStretchSpacer(1);
            $self->add_button( $panel, $boxsizer, "Three stat",
                \&event_three_stat );

            $boxsizer->AddSpacer(10);
            $self->add_button( $panel, $boxsizer, "Three chart",
                \&event_three_chart );
        }
    }

    return ( $pane, $nb );
}

sub build_window {
    my $self = shift;

    # init sizer
    # Naming $self and $main_sizer as $top_panel and $top_sizer is a
    #   historical issue, but this works well and make $top_panel is more
    #   clean than $self
    my $main_sizer = Wx::BoxSizer->new(wxVERTICAL);
    $self->SetSizer($main_sizer);
    my $top_panel = $self;
    my $top_sizer = $main_sizer;

    # build collapsible notebooks
    my ( $database_pane, $database_nb )
        = $self->build_database_pane( $top_panel, $top_sizer );
    my ( $generate_pane, $generate_nb )
        = $self->build_generate_pane( $top_panel, $top_sizer );
    my ( $update_pane, $update_nb )
        = $self->build_update_pane( $top_panel, $top_sizer );
    my ( $stat_pane, $stat_nb )
        = $self->build_stat_pane( $top_panel, $top_sizer );
    my ( $util_pane, $util_nb )
        = $self->build_util_pane( $top_panel, $top_sizer );
    my ( $three_pane, $three_nb )
        = $self->build_three_pane( $top_panel, $top_sizer );

    # set default styles of every panes and notebooks
    $database_pane->Expand;
    $generate_pane->Expand;
    $update_pane->Collapse;
    $three_pane->Collapse;
    $stat_pane->Collapse;
    $util_pane->Collapse;
    $database_nb->ChangeSelection(1);    # select Initiate tab by default

    # automatic layout, size the window optimally and set its minimal size
    $top_panel->SetAutoLayout(1);
    $top_panel->SetSizerAndFit($top_sizer);
    $top_sizer->SetSizeHints($top_panel);

    return;
}

#----------------------------#
# setter and getter of object attrs which are all textctrl or checkbox
#----------------------------#
sub set_value {
    my $self  = shift;
    my $attr  = shift;
    my $value = shift;

    $self->{$attr}->SetValue($value) if defined $value;

    return;
}

sub get_value {
    my $self = shift;
    my $attr = shift;

    if ( ref $self->{$attr} eq "Wx::Choice" ) {
        return $self->{$attr}->GetStringSelection;
    }
    return $self->{$attr}->GetValue;
}

#----------------------------#
# methods creating unified widgets
#----------------------------#
sub add_static_text {
    my $self   = shift;
    my $parent = shift;
    my $sizer  = shift;
    my $text   = shift;

    my $static_text = Wx::StaticText->new( $parent, -1, $text );
    my $flag = wxALIGN_CENTER_VERTICAL | wxALIGN_LEFT | wxALL;
    $sizer->Add( $static_text, 0, $flag, 2 );

    return;
}

sub add_gb_static_text {
    my $self        = shift;
    my $parent      = shift;
    my $sizer       = shift;
    my $text        = shift;
    my $gb_position = shift;
    my $gb_span     = shift || [ 1, 1 ];

    my $static_text = Wx::StaticText->new( $parent, -1, $text );
    my $flag = wxALIGN_CENTER_VERTICAL | wxALIGN_RIGHT | wxALL;
    $sizer->AddWindow(
        $static_text,
        Wx::GBPosition->new( @{$gb_position} ),
        Wx::GBSpan->new( @{$gb_span} ),
        $flag, 2,
    );

    return;
}

sub add_text_ctrl {
    my $self           = shift;
    my $parent         = shift;
    my $sizer          = shift;
    my $object_attr    = shift;
    my $stretch_factor = shift || 0;
    my $text_style     = shift || 0;
    my $size           = shift || [ -1, -1 ];

    my $text_ctrl = Wx::TextCtrl->new( $parent, -1, '', [ -1, -1 ], $size,
        $text_style );
    my $flag = wxGROW | wxALIGN_CENTER_VERTICAL | wxALL;
    $sizer->Add( $text_ctrl, $stretch_factor, $flag, 2 );

    my $font = Wx::Font->new( 10, wxSWISS, wxNORMAL, wxNORMAL, 0, "Arial" );
    $text_ctrl->SetFont($font);

    $self->{$object_attr} = $text_ctrl;

    return;
}

sub add_gb_text_ctrl {
    my $self        = shift;
    my $parent      = shift;
    my $sizer       = shift;
    my $object_attr = shift;
    my $gb_position = shift;
    my $gb_span     = shift || [ 1, 1 ];
    my $size        = shift || [ -1, -1 ];

    my $text_ctrl = Wx::TextCtrl->new( $parent, -1, '', [ -1, -1 ], $size );
    my $flag = wxGROW | wxALIGN_CENTER_VERTICAL | wxALL;

    $sizer->AddWindow(
        $text_ctrl,
        Wx::GBPosition->new( @{$gb_position} ),
        Wx::GBSpan->new( @{$gb_span} ),
        $flag, 2,
    );

    my $font = Wx::Font->new( 10, wxSWISS, wxNORMAL, wxNORMAL, 0, "Arial" );
    $text_ctrl->SetFont($font);

    $self->{$object_attr} = $text_ctrl;

    return;
}

sub add_button {
    my $self   = shift;
    my $parent = shift;
    my $sizer  = shift;
    my $text   = shift;
    my $event  = shift;
    my $stretch_factor = shift || 0;

    my $button = Wx::Button->new( $parent, -1, $text );
    my $flag = wxALIGN_CENTER_VERTICAL | wxALL;
    $sizer->Add( $button, $stretch_factor, $flag, 0 );

    EVT_BUTTON( $self, $button, $event );

    return;
}

sub add_bitmap_button {
    my $self   = shift;
    my $parent = shift;
    my $sizer  = shift;
    my $bitmap = shift;
    my $event  = shift;

    my $button = Wx::BitmapButton->new( $parent, -1, $bitmap );
    my $flag = wxALIGN_CENTER_VERTICAL | wxALL;
    $sizer->Add( $button, 0, $flag, 0 );

    EVT_BUTTON( $self, $button, $event );

    return;
}

sub add_gb_bitmap_button {
    my $self        = shift;
    my $parent      = shift;
    my $sizer       = shift;
    my $bitmap      = shift;
    my $event       = shift;
    my $gb_position = shift;
    my $gb_span     = shift || [ 1, 1 ];

    my $button = Wx::BitmapButton->new( $parent, -1, $bitmap );
    my $flag = wxALIGN_CENTER_VERTICAL | wxALL;
    $sizer->AddWindow(
        $button,
        Wx::GBPosition->new( @{$gb_position} ),
        Wx::GBSpan->new( @{$gb_span} ),
        $flag, 2,
    );

    EVT_BUTTON( $self, $button, $event );

    return;
}

sub add_check_box {
    my $self        = shift;
    my $parent      = shift;
    my $sizer       = shift;
    my $object_attr = shift;
    my $text        = shift || $object_attr;

    my $check_box = Wx::CheckBox->new( $parent, -1, $text, );
    my $flag = wxGROW | wxALIGN_CENTER_VERTICAL | wxALL;
    $sizer->Add( $check_box, 0, $flag, 2 );

    $self->{$object_attr} = $check_box;

    return;
}

sub add_choice {
    my $self        = shift;
    my $parent      = shift;
    my $sizer       = shift;
    my $object_attr = shift;
    my $choices     = shift || [$object_attr];
    my $selection   = shift || 0;

    my $choice
        = Wx::Choice->new( $parent, -1, [ -1, -1 ], [ 70, -1 ], $choices, );
    $choice->SetSelection($selection);
    my $flag = wxGROW | wxALIGN_CENTER_VERTICAL | wxALL;
    $sizer->Add( $choice, 0, $flag, 2 );

    $self->{$object_attr} = $choice;

    return;
}

sub add_static_line {
    my $self   = shift;
    my $parent = shift;
    my $sizer  = shift;

    my $static_line = Wx::StaticLine->new(
        $parent, -1,
        [ -1, -1 ],
        [ -1, -1 ],
        wxLI_HORIZONTAL
    );
    $sizer->AddWindow( $static_line, 0,
        wxGROW | wxALIGN_CENTER_VERTICAL | wxALL, 5 );

    return;
}

sub add_boxsizer_h {
    my $self         = shift;
    my $parent_sizer = shift;

    my $sizer = Wx::BoxSizer->new(wxHORIZONTAL);
    $parent_sizer->Add( $sizer, 0, wxGROW, 0 );

    return $sizer;
}

sub add_boxsizer_v {
    my $self         = shift;
    my $parent_sizer = shift;

    my $sizer = Wx::BoxSizer->new(wxVERTICAL);
    $parent_sizer->Add( $sizer, 0, wxGROW, 0 );

    return $sizer;
}

sub add_static_panel {
    my $self         = shift;
    my $parent       = shift;
    my $parent_sizer = shift;
    my $title        = shift;

    my $sub_panel = Wx::Panel->new( $parent, -1 );
    $parent_sizer->Add( $sub_panel, 0, wxGROW, 0 );

    my $sub_sizer = Wx::StaticBoxSizer->new(
        Wx::StaticBox->new( $sub_panel, -1, $title ), wxVERTICAL );

    $sub_panel->SetSizer($sub_sizer);

    return ( $sub_panel, $sub_sizer );
}

sub add_collapsible_nb {
    my $self         = shift;
    my $parent       = shift;
    my $parent_sizer = shift;
    my $title        = shift;

    # create a CollapsiblePane that will contain the following Notebook
    my $sub_pane = Wx::CollapsiblePane->new( $parent, -1, $title );
    $parent_sizer->Add( $sub_pane, 0, wxGROW, 0 );
    my $sub_window = $sub_pane->GetPane;
    my $sub_sizer  = Wx::BoxSizer->new(wxVERTICAL);
    my $sub_nb     = Wx::Notebook->new( $sub_window, -1 );
    $sub_sizer->Add( $sub_nb, 0, wxGROW | wxALL, 0 );

    $sub_window->SetSizer($sub_sizer);
    $sub_sizer->SetSizeHints($sub_window);

    EVT_COLLAPSIBLEPANE_CHANGED( $parent, $sub_pane, \&event_pane_changed );

    return ( $sub_pane, $sub_nb );
}

sub add_nb_page {
    my $self     = shift;
    my $notebook = shift;
    my $title    = shift;

    my $sub_panel = Wx::Panel->new( $notebook, -1 );
    $notebook->AddPage( $sub_panel, $title );

    my $sub_sizer = Wx::BoxSizer->new(wxVERTICAL);
    $sub_panel->SetSizer($sub_sizer);

    return ( $sub_panel, $sub_sizer );
}

#----------------------------#
# top panel callback events
#----------------------------#
sub event_pane_changed {
    my ( $self, $event ) = @_;

    $self->GetSizer->Layout;
    return;
}

sub event_close_window {
    my ( $self, $event ) = @_;

    $self->Destroy;
    return;
}

1;

