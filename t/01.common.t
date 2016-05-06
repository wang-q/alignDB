use strict;
use warnings;
use autodie;

use Test::More;
use Test::Number::Delta within => 1e-2;

use AlignDB::Common;

{
    print "#calc_gc_ratio\n";

    my @data = (
        [ "ATAA",            0 ],
        [ "AtaA",            0 ],
        [ "CCGC",            1 ],
        [ "CcGc",            1 ],
        [ "TAGggATaaC",      0.4 ],
        [ "GCaN--NN--NNNaC", 0.6 ],
    );

    for my $i ( 0 .. $#data ) {
        my ( $ori, $expected ) = @{ $data[$i] };
        my $result = AlignDB::Common::calc_gc_ratio( [$ori] );
        print "original: $ori\n";
        is( $result, $expected, "calc_gc_ratio_$i" );
    }
}

{
    print "#multi_seq_stat\n";

    #$seq_legnth,            $number_of_comparable_bases,
    #$number_of_identities,  $number_of_differences,
    #$number_of_gaps,        $number_of_n,
    #$number_of_align_error, $pi,
    #$first_seq_gc,          $average_gc,
    my @data = (

        #AAAATTTTGG
        #AAAATTTTTG
        [ [qw{ AAAATTTTGG AAAATTTTTG }], [ 10, 10, 9, 1, 0, 0, 0, 0.1, 0.2, 0.15 ], ],

        #TTAGCCGCTGAGAAGC
        #GTAGCCGCTGA-AGGC
        [   [qw{ TTAGCCGCTGAGAAGC GTAGCCGCTGA-AGGC }],
            [ 16, 15, 13, 2, 1, 0, 0, 0.1333, 0.5625, 0.6146 ],
        ],

        #GATTATCATCACCCCAGCCACATA
        #GATTTT--TCACTCCATTCGCATA
        [   [qw{ GATTATCATCACCCCAGCCACATW GATTTT--TCACTCCATTCGCATA }],
            [ 24, 21, 16, 5, 2, 1, 0, 0.2381, 0.4783, 0.4209 ],
        ],

    );

    for my $i ( 0 .. $#data ) {
        my ( $seq_pair_ref, $except_ref ) = @{ $data[$i] };

        my $result_multi_ref = AlignDB::Common::multi_seq_stat($seq_pair_ref);
        Test::Number::Delta::delta_ok( $result_multi_ref, $except_ref, "stat $i" );
    }
}

{
    print "#find_indel_set\n";

    my @data = (
        [ "ATAA",            "-" ],
        [ "CcGc",            "-" ],
        [ "TAGggATaaC",      "-" ],
        [ "C-Gc",            "2" ],
        [ "C--c",            "2-3" ],
        [ "---c",            "1-3" ],
        [ "C---",            "2-4" ],
        [ "GCaN--NN--NNNaC", "5-6,9-10" ],
    );

    for my $i ( 0 .. $#data ) {
        my ( $ori, $expected ) = @{ $data[$i] };
        my $result = AlignDB::Common::find_indel_set($ori);
        print "original: $ori\n";
        is( $result->runlist, $expected, "calc_gc_ratio_$i" );
    }
}

done_testing();
