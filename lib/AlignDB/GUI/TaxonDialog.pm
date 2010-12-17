package AlignDB::GUI::TaxonDialog;
use strict;
use warnings;

use Wx qw(:sizer :listctrl :id :staticline wxSUNKEN_BORDER);
use Wx::Event qw(EVT_LIST_ITEM_SELECTED EVT_LIST_COL_CLICK );

use Parse::CSV;
use Scalar::Util qw(looks_like_number);

use base qw(Wx::Dialog);
use base qw(Class::Accessor);

use FindBin;
use lib "$FindBin::Bin/../lib";
require AlignDB::GUI::Bitmap;

__PACKAGE__->mk_accessors(qw(id name listctrl sort_order item_data));

sub new {
    my $ref = shift;

    my $self = $ref->SUPER::new(
        undef,              # parent window
        -1,                 # ID -1 means any
        'Load taxon...',    # title
        [ -1, -1 ],         # default position
        [ -1, -1 ],         # default size
    );
    
    my $icon = Wx::Icon->new;
    $icon->CopyFromBitmap( AlignDB::GUI::Bitmap->get_bitmap_load );
    $self->SetIcon($icon);
    
    $self->{sort_order} = 'asc';
    return $self;
}

sub build_window {
    my $self     = shift;
    my $filename = shift;

    my $top_sizer = Wx::BoxSizer->new(wxVERTICAL);

    # add static text
    $top_sizer->Add( Wx::StaticText->new( $self, -1, "Choose a taxon" ),
        0, wxALIGN_LEFT | wxALL, 5 );

    # build taxon listctrl
    {
        my $listctrl = Wx::ListCtrl->new(
            $self, -1,
            [ -1,  -1 ],
            [ 360, 200 ],
            wxLC_REPORT | wxLC_SINGLE_SEL | wxSUNKEN_BORDER
        );
        $self->listctrl($listctrl);
        $top_sizer->Add( $listctrl, 0, wxGROW | wxALIGN_CENTER | wxALL, 5 );

        # use this in sort
        my $item_data = {};
        $self->{item_data} = $item_data;

        # create columns
        my $item_col = Wx::ListItem->new();
        $item_col->SetText("id");
        $listctrl->InsertColumn( 0, $item_col );
        $item_col->SetText("name");
        $listctrl->InsertColumn( 1, $item_col );
        $item_col->SetText("genus");
        $listctrl->InsertColumn( 2, $item_col );
        $item_col->SetText("species");
        $listctrl->InsertColumn( 3, $item_col );

        # read out normal taxons from filename
        # and then put them into listctrl
        my $csv = Parse::CSV->new( file => $filename );
        $csv->fetch;    # bypass title

        my $i = 0;
        while ( my $record = $csv->fetch ) {
            last unless @$record;
            my $idx = $listctrl->InsertStringItem( $i, $record->[0] );
            $listctrl->SetItemData( $idx, $i );
            $item_data->{$i} = { 0 => $record->[0] };
            for ( 1 .. 3 ) {
                $listctrl->SetItem( $idx, $_, $record->[$_] );
                $item_data->{$i}{$_} = $record->[$_];
            }
            $i++;
        }

        # auto adjust column width
        for ( 0 .. 3 ) {
            $listctrl->SetColumnWidth( $_, wxLIST_AUTOSIZE );
        }
    }

    my $static_line = Wx::StaticLine->new(
        $self, -1,
        [ -1, -1 ],
        [ -1, -1 ],
        wxLI_HORIZONTAL
    );
    $top_sizer->AddWindow( $static_line, 0,
        wxGROW | wxALIGN_CENTER_VERTICAL | wxALL, 5 );

    ## A StdDialogButtonSizer containing OK and Cancel buttons
    my $ok_cancal_box = Wx::StdDialogButtonSizer->new;
    $top_sizer->Add( $ok_cancal_box, 0, wxALIGN_CENTER_HORIZONTAL | wxALL,
        5 );

    # the OK button
    my $ok_button = Wx::Button->new( $self, wxID_OK, "&OK" );
    $ok_cancal_box->AddButton($ok_button);

    # The Cancel button
    my $cancel_button = Wx::Button->new( $self, wxID_CANCEL, "&Cancel" );
    $ok_cancal_box->AddButton($cancel_button);

    # Rearranges the buttons and applies proper spacing between buttons to
    # make them match the platform or toolkit's interface guidelines.
    $ok_cancal_box->Realize;

    # tell we want automatic layout
    $self->SetAutoLayout(1);
    $self->SetSizer($top_sizer);

    # size the window optimally and set its minimal size
    $top_sizer->Fit($self);
    $top_sizer->SetSizeHints($self);

    EVT_LIST_ITEM_SELECTED( $self, $self->listctrl, \&OnSelected );
    EVT_LIST_COL_CLICK( $self, $self->listctrl, \&OnColClick );

    return;
}

sub OnSelected {
    my ( $self, $event ) = @_;

    my $listctrl = $self->listctrl;
    $self->{id}   = $listctrl->GetItem( $event->GetIndex, 0 )->GetText;
    $self->{name} = $listctrl->GetItem( $event->GetIndex, 1 )->GetText;

    return;
}

sub OnColClick {
    my ( $self, $event ) = @_;

    my $listctrl   = $self->listctrl;
    my $sort_order = $self->sort_order;

    my $sorter = sub {
        my $item0 = $self->item_data->{ $_[0] }{ $event->GetColumn };
        my $item1 = $self->item_data->{ $_[1] }{ $event->GetColumn };

        if ( $sort_order eq 'asc' ) {
            $self->sort_order('desc');
        }
        else {
            $self->sort_order('asc');
            ( $item0, $item1 ) = ( $item1, $item0 );
        }

        if ( looks_like_number($item0) and looks_like_number($item1) ) {
            return $item0 < $item1;
        }
        else {
            return $item0 lt $item1;
        }
    };

    $listctrl->SortItems($sorter);

    return;
}

1;
