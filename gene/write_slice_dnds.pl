#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Multi;
use AlignDB::Position;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

my $outdir;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=s'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'output=s'   => \$outdir,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
my $obj = AlignDB::Multi->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
my $dbh = $obj->dbh;
my $pos_obj = AlignDB::Position->new( dbh => $dbh );

my $group = {
    group_dn0 => {},
    group_ds0 => {},
};

my $taxon_id;
{
    my $sth = $dbh->prepare(
        q{
        SELECT 
                c.taxon_id
        FROM
            sequence s
            INNER JOIN target t
                ON s.seq_id = t.seq_id
            INNER JOIN chromosome c
                ON s.chr_id = c.chr_id
        }
    );
    $sth->execute;
    ($taxon_id) = $sth->fetchrow_array;
}

{    # dn = 0, ds > 0
    my $sth = $dbh->prepare(
        q{
        SELECT 
                    c.chr_name,   w.align_id,   g.gene_tl_runlist
            FROM
                gene g INNER JOIN window w
                    ON g.window_id = w.window_id
                INNER JOIN sequence s
                    ON s.align_id = w.align_id
                INNER JOIN target t
                    ON s.seq_id = t.seq_id
                INNER JOIN chromosome c
                    ON s.chr_id = c.chr_id
            WHERE
                1 = 1 
                AND gene_nsy = 0 
                AND gene_syn > 0 
                AND w.window_pi > 0.001
        }
    );
    $sth->execute;

    while ( my ( $chr_name, $align_id, $runlist ) = $sth->fetchrow_array ) {
        my $set = to_target_chr( $align_id, $runlist );

        if ( !defined $group->{group_dn0}{$chr_name} ) {
            $group->{group_dn0}{$chr_name} = $set;
        }
        else {
            $group->{group_dn0}{$chr_name}->add($set);
        }
    }
}

{    # dn > 0, ds = 0

    my $sth = $dbh->prepare(
        q{
        SELECT 
                    c.chr_name,   w.align_id,   g.gene_tl_runlist
            FROM
                gene g INNER JOIN window w
                    ON g.window_id = w.window_id
                INNER JOIN sequence s
                    ON s.align_id = w.align_id
                INNER JOIN target t
                    ON s.seq_id = t.seq_id
                INNER JOIN chromosome c
                    ON s.chr_id = c.chr_id
            WHERE
                1 = 1 
                AND gene_nsy > 0 
                AND gene_syn > 0 
                AND w.window_pi > 0.001
        }
    );
    $sth->execute;
    while ( my ( $chr_name, $align_id, $runlist ) = $sth->fetchrow_array ) {
        my $set = to_target_chr( $align_id, $runlist );
        if ( !defined $group->{group_ds0}{$chr_name} ) {
            $group->{group_ds0}{$chr_name} = $set;
        }
        else {
            $group->{group_ds0}{$chr_name}->add($set);
        }
    }
}

{    # dnds groups

    my @dnds_levels
        = ( [ 1, 0, 0.3 ], [ 2, 0.3, 0.6 ], [ 3, 0.6, 1 ], [ 4, 1, 9999 ], );

    for my $level (@dnds_levels) {
        my ( $order, $low_border, $high_border ) = @{$level};
        my $group_name = "group_$order-$low_border-$high_border";

        my $sth = $dbh->prepare(
            q{
            SELECT 
                        c.chr_name,   w.align_id,   g.gene_tl_runlist
                FROM
                    gene g INNER JOIN window w
                        ON g.window_id = w.window_id
                    INNER JOIN sequence s
                        ON s.align_id = w.align_id
                    INNER JOIN target t
                        ON s.seq_id = t.seq_id
                    INNER JOIN chromosome c
                        ON s.chr_id = c.chr_id
                WHERE
                    1 = 1 
                    AND gene_nsy > 0 
                    AND gene_syn > 0 
                    AND w.window_pi > 0.001
                    AND gene_nsy / gene_syn >= ?
                    AND gene_nsy / gene_syn <= ?
            }
        );
        $sth->execute( $low_border, $high_border );

        while ( my ( $chr_name, $align_id, $runlist ) = $sth->fetchrow_array )
        {
            my $set = to_target_chr( $align_id, $runlist );
            if ( !defined $group->{$group_name}{$chr_name} ) {
                $group->{$group_name}{$chr_name} = $set;
            }
            else {
                $group->{$group_name}{$chr_name}->add($set);
            }
        }
    }
}

for my $g ( keys %{$group} ) {
    for my $chr ( keys %{ $group->{$g} } ) {
        $group->{$g}{$chr} = $group->{$g}{$chr}->pad(1000);
        $group->{$g}{$chr} = $group->{$g}{$chr}->diff("-9999-0");
        $group->{$g}{$chr} = $group->{$g}{$chr}->runlist;
    }
}

for my $g ( keys %{$group} ) {
    DumpFile( "$g.yml", $group->{$g} );
}

sub to_target_chr {
    my $align_id = shift;
    my $runlist  = shift;

    my $set     = AlignDB::IntSpan->new($runlist);
    my $chr_set = AlignDB::IntSpan->new;
    for my $span ( $set->spans ) {
        my ( $lower, $upper ) = @{$span};

        # align coordinates to target & query chromosome coordinates
        my $chr_lower = $pos_obj->at_target_chr( $align_id, $lower );
        my $chr_upper = $pos_obj->at_target_chr( $align_id, $upper );
        $chr_set->add_range( $chr_lower, $chr_upper );
    }

    return $chr_set;

}


__END__

write yml

