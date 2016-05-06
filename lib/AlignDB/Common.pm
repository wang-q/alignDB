# Static functions for AlignDB. Export nothing.
package AlignDB::Common;
use strict;
use warnings;
use autodie;

use Carp;
use List::Util;
use List::MoreUtils::PP;
use YAML::Syck;

use AlignDB::IntSpan;

sub mean {
    @_ = grep { defined $_ } @_;
    return unless @_;
    return $_[0] unless @_ > 1;
    return List::Util::sum(@_) / scalar(@_);
}

sub calc_gc_ratio {
    my $seq_refs = shift;

    my $seq_count = scalar @{$seq_refs};

    my @ratios;
    for my $i ( 0 .. $seq_count - 1 ) {

        # Count all four bases
        my $a_count = $seq_refs->[$i] =~ tr/Aa/Aa/;
        my $g_count = $seq_refs->[$i] =~ tr/Gg/Gg/;
        my $c_count = $seq_refs->[$i] =~ tr/Cc/Cc/;
        my $t_count = $seq_refs->[$i] =~ tr/Tt/Tt/;

        my $four_count = $a_count + $g_count + $c_count + $t_count;
        my $gc_count   = $g_count + $c_count;

        if ( $four_count == 0 ) {
            next;
        }
        else {
            my $gc_ratio = $gc_count / $four_count;
            push @ratios, $gc_ratio;
        }
    }

    return mean(@ratios);
}

sub pair_D {
    my $seq_refs = shift;

    my $seq_count = scalar @{$seq_refs};
    if ( $seq_count != 2 ) {
        Carp::confess "Need two sequences\n";
    }

    my $legnth = length $seq_refs->[0];

    my ( $comparable_bases, $differences ) = (0) x 2;

    for my $pos ( 1 .. $legnth ) {
        my $base0 = substr $seq_refs->[0], $pos - 1, 1;
        my $base1 = substr $seq_refs->[1], $pos - 1, 1;

        if (   $base0 =~ /[atcg]/i
            && $base1 =~ /[atcg]/i )
        {
            $comparable_bases++;
            if ( $base0 ne $base1 ) {
                $differences++;
            }
        }
    }

    if ( $comparable_bases == 0 ) {
        Carp::carp "comparable_bases == 0\n";
        return 0;
    }
    else {
        return $differences / $comparable_bases;
    }
}

sub multi_seq_stat {
    my $seq_refs = shift;

    my $seq_count = scalar @{$seq_refs};
    if ( $seq_count < 2 ) {
        Carp::confess "Need two or more sequences\n";
    }

    my $legnth = length $seq_refs->[0];

    # For every positions, search for polymorphism_site
    my ( $comparable_bases, $identities, $differences, $gaps, $ns, $align_errors, ) = (0) x 6;
    for my $pos ( 1 .. $legnth ) {
        my @bases = ();
        for my $i ( 0 .. $seq_count - 1 ) {
            my $base = substr( $seq_refs->[$i], $pos - 1, 1 );
            push @bases, $base;
        }
        @bases = List::MoreUtils::PP::uniq(@bases);

        if ( List::MoreUtils::PP::all { $_ =~ /[agct]/i } @bases ) {
            $comparable_bases++;
            if ( List::MoreUtils::PP::all { $_ eq $bases[0] } @bases ) {
                $identities++;
            }
            else {
                $differences++;
            }
        }
        elsif ( List::MoreUtils::PP::any { $_ eq '-' } @bases ) {
            $gaps++;
        }
        else {
            $ns++;
        }
    }
    if ( $comparable_bases == 0 ) {
        print YAML::Syck::Dump { seqs => $seq_refs, };
        Carp::carp "number_of_comparable_bases == 0!!\n";
        return [
            $legnth, $comparable_bases, $identities, $differences, $gaps,
            $ns,     $legnth,           undef,       undef,        undef,
        ];
    }

    my @all_Ds;
    for ( my $i = 0; $i < $seq_count; $i++ ) {
        for ( my $j = $i + 1; $j < $seq_count; $j++ ) {
            my $D = pair_D( [ $seq_refs->[$i], $seq_refs->[$j] ] );
            push @all_Ds, $D;
        }
    }

    my $D = mean(@all_Ds);

    my $target_gc = calc_gc_ratio( [ $seq_refs->[0] ] );
    my $average_gc = calc_gc_ratio($seq_refs);

    return [
        $legnth, $comparable_bases, $identities, $differences, $gaps,
        $ns,     $align_errors,     $D,          $target_gc,   $average_gc,
    ];
}

sub find_indel_set {
    my $seq = shift;

    my $length = length $seq;

    my $indel_set    = AlignDB::IntSpan->new;
    my $indel_offset = 0;
    my $indel_start  = 0;
    my $indel_end    = 0;

    for my $pos ( 1 .. $length ) {
        my $base = substr( $seq, $pos - 1, 1 );
        if ( $base eq '-' ) {
            if ( $indel_offset == 0 ) {
                $indel_start = $pos;
            }
            $indel_offset++;
        }
        else {
            if ( $indel_offset != 0 ) {
                $indel_end = $pos - 1;
                $indel_set->add_pair( $indel_start, $indel_end );
            }
            $indel_offset = 0;
        }
    }

    if ( $indel_offset != 0 ) {
        $indel_end = $length;
        $indel_set->add_pair( $indel_start, $indel_end );
    }

    $indel_set = $indel_set->intersect("1-$length");
    return $indel_set;
}


sub decode_header {
    my $header = shift;

    tie my %info, "Tie::IxHash";

    # S288.chrI(+):27070-29557|species=S288C
    my $head_qr = qr{
        (?:(?P<name>[\w_]+)\.)?
            (?P<chr_name>[\w-]+)
            (?:\((?P<chr_strand>.+)\))?
            [\:]                        # spacer
            (?P<chr_start>\d+)
            [\_\-]?                      # spacer
            (?P<chr_end>\d+)?
        }xi;

    $header =~ $head_qr;
    my $chr_name  = $2;
    my $chr_start = $4;
    my $chr_end   = $5;

    if ( defined $chr_name and defined $chr_start ) {
        if ( !defined $chr_end ) {
            $chr_end = $chr_start;
        }
        %info = (
            name       => $1,
            chr_name   => $chr_name,
            chr_strand => $3,
            chr_start  => $chr_start,
            chr_end    => $chr_end,
        );
        if ( defined $info{chr_strand} ) {
            if ( $info{chr_strand} eq '1' ) {
                $info{chr_strand} = '+';
            }
            elsif ( $info{chr_strand} eq '-1' ) {
                $info{chr_strand} = '-';
            }
        }
    }
    else {
        $header =~ /^(\S+)/;
        my $name = $1;
        %info = (
            name       => $name,
            chr_name   => undef,
            chr_strand => undef,
            chr_start  => undef,
            chr_end    => undef,
        );
    }

    # additional keys
    if ( $header =~ /\|(.+)/ ) {
        my @parts = grep {defined} split /;/, $1;
        for my $part (@parts) {
            my ( $key, $value ) = split /=/, $part;
            if ( defined $key and defined $value ) {
                $info{$key} = $value;
            }
        }
    }

    return \%info;
}

1;
