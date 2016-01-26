#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use AlignDB::SQL;
use AlignDB::SQL::Library;

#----------------------------------------------------------#
# SQL
#----------------------------------------------------------#

# Object headers in sql_library are named under the following rules:
#   TYEP-NAME-BINDINGs
# e.g.: common-distance-0
#       multi-distance-0

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

ld_stat_factory.pl - LD stats for alignDB

=head1 SYNOPSIS

    perl ld_stat_factory.pl [options]
      Options:
        --help      -?          brief help message
        --lib           STR     Path to sql.lib
        --verbose               verbose mode
        
=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'lib=s' => \( my $lib_file = "$FindBin::RealBin/sql.lib" ),
    'verbose' => \my $verbose,
) or HelpMessage(1);

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
unlink $lib_file if -e $lib_file;
my $sql_file = AlignDB::SQL::Library->new( lib => $lib_file );

sub ns { return AlignDB::SQL->new; }

#----------------------------------------------------------#
# stat_factory.pl SQL
#----------------------------------------------------------#

#SELECT
#  isw.isw_distance distance,
#  AVG(isw.isw_pi) AVG_pi,
#  STD(isw.isw_pi) STD_pi,
#  AVG(isw.isw_average_gc) AVG_gc,
#  STD(isw.isw_average_gc) STD_gc,
#  AVG(isw.isw_cv) AVG_cv,
#  STD(isw.isw_cv) STD_cv,
#  COUNT(*) COUNT
#FROM isw
#GROUP BY
#  isw_isw_distance
{
    my $name = 'common-d1_pi_gc_cv-0';

    my $sql = ns();
    $sql->add_select( 'isw.isw_distance',        'distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',         'AVG_pi' );
    $sql->add_select( 'STD(isw.isw_pi)',         'STD_pi' );
    $sql->add_select( 'AVG(isw.isw_average_gc)', 'AVG_gc' );
    $sql->add_select( 'STD(isw.isw_average_gc)', 'STD_gc' );
    $sql->add_select( 'AVG(isw.isw_cv)',         'AVG_cv' );
    $sql->add_select( 'STD(isw.isw_cv)',         'STD_cv' );
    $sql->add_select( 'COUNT(*)',                'COUNT' );
    $sql->from( ['isw'] );
    $sql->group( { column => 'isw.isw_distance' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT
#  isw.isw_distance distance,
#  COUNT(*) COUNT,
#  SUM(isw.isw_length) SUM_length
#FROM isw
#GROUP BY
#  isw.isw_distance
{
    my $sql = ns();
    $sql->add_select( 'isw.isw_distance', 'distance' );
    $sql->add_select( 'COUNT(*)',         'COUNT' );

    # for group last portion
    $sql->add_select( 'SUM(isw_length)', 'SUM_length' );

    $sql->from( ['isw'] );
    $sql->group( { column => 'isw.isw_distance' } );

    $sql_file->set( 'common-d1_combine-0', $sql );
    print $sql->as_sql if $verbose;
}

#SELECT
#  isw.isw_distance distance,
#  AVG(isw.isw_pi) AVG_pi,
#  STD(isw.isw_pi) STD_pi,
#  AVG(isw.isw_average_gc) AVG_gc,
#  STD(isw.isw_average_gc) STD_gc,
#  AVG(isw.isw_cv) AVG_cv,
#  STD(isw.isw_cv) STD_cv,
#  COUNT(*) COUNT
#FROM isw
#GROUP BY
#  isw_isw_distance
{
    my $name = 'common-d1_comb_pi_gc_cv-0';

    my $sql = ns();
    $sql->add_select( 'AVG(isw.isw_distance)',   'AVG_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',         'AVG_pi' );
    $sql->add_select( 'STD(isw.isw_pi)',         'STD_pi' );
    $sql->add_select( 'AVG(isw.isw_average_gc)', 'AVG_gc' );
    $sql->add_select( 'STD(isw.isw_average_gc)', 'STD_gc' );
    $sql->add_select( 'AVG(isw.isw_cv)',         'AVG_cv' );
    $sql->add_select( 'STD(isw.isw_cv)',         'STD_cv' );
    $sql->add_select( 'COUNT(*)',                'COUNT' );
    $sql->from( ['isw'] );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT
#  AVG(isw.isw_distance) AVG_distance,
#  AVG(isw.isw_pi) AVG_pi,
#  COUNT(*) COUNT,
#  STD(isw.isw_pi) STD_pi
#FROM isw
{
    my $name = 'common-d1_pi_avg-0';

    my $sql = ns();
    $sql->add_select( 'AVG(isw.isw_distance)', 'AVG_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',       'AVG_pi' );
    $sql->add_select( 'COUNT(*)',              'COUNT' );
    $sql->add_select( 'STD(isw.isw_pi)',       'STD_pi' );
    $sql->from( ['isw'] );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT
#  isw.isw_distance distance,
#  COUNT(*) COUNT
#FROM isw
#WHERE (isw.isw_coding >= ?)
#  AND (isw.isw_coding <= ?)
#GROUP BY
#  isw.isw_distance
{
    my $name = 'common-d1_make_combine_coding-0';

    my $sql = ns();
    $sql->add_select( 'isw.isw_distance', 'distance' );
    $sql->add_select( 'COUNT(*)',         'COUNT' );
    $sql->from( ['isw'] );
    $sql->add_where( 'isw.isw_coding' => { op => '>=', value => '1' } );
    $sql->add_where( 'isw.isw_coding' => { op => '<=', value => '1' } );
    $sql->group( { column => 'isw.isw_distance' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT
#  AVG(isw.isw_distance) AVG_distance,
#  AVG(isw.isw_pi) AVG_pi,
#  STD(isw.isw_pi) STD_pi,
#  AVG(isw.isw_average_gc) AVG_gc,
#  STD(isw.isw_average_gc) STD_gc,
#  AVG(isw.isw_cv) AVG_cv,
#  STD(isw.isw_cv) STD_cv,
#  COUNT(*) COUNT
#FROM isw
#WHERE (isw.isw_coding >= ?)
#  AND (isw.isw_coding <= ?)
{
    my $name = 'common-d1_comb_coding-0';

    my $sql = ns();
    $sql->add_select( 'AVG(isw.isw_distance)',   'AVG_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',         'AVG_pi' );
    $sql->add_select( 'STD(isw.isw_pi)',         'STD_pi' );
    $sql->add_select( 'AVG(isw.isw_average_gc)', 'AVG_gc' );
    $sql->add_select( 'STD(isw.isw_average_gc)', 'STD_gc' );
    $sql->add_select( 'AVG(isw.isw_cv)',         'AVG_cv' );
    $sql->add_select( 'STD(isw.isw_cv)',         'STD_cv' );
    $sql->add_select( 'COUNT(*)',                'COUNT' );
    $sql->from( ['isw'] );

    $sql->add_where( 'isw.isw_coding' => { op => '>=', value => '1' } );
    $sql->add_where( 'isw.isw_coding' => { op => '<=', value => '1' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT
#  isw.isw_distance distance,
#  COUNT(*) COUNT
#FROM isw
#  INNER JOIN indel ON
#    isw.isw_indel_id = indel.indel_id
#WHERE (indel.indel_slippage >= ?)
#  AND (indel.indel_slippage <= ?)
#GROUP BY
#  isw.isw_distance
{
    my $name = 'common-d1_make_combine_slippage-0';

    my $sql = ns();
    $sql->add_select( 'isw.isw_distance', 'distance' );
    $sql->add_select( 'COUNT(*)',         'COUNT' );

    $sql->add_join(
        isw => {
            type      => 'inner',
            table     => 'indel',
            condition => 'isw.isw_indel_id = indel.indel_id',
        }
    );
    $sql->add_where( 'indel.indel_slippage' => { op => '>=', value => '1' } );
    $sql->add_where( 'indel.indel_slippage' => { op => '<=', value => '1' } );
    $sql->group( { column => 'isw.isw_distance' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT
#  AVG(isw.isw_distance) AVG_distance,
#  AVG(isw.isw_pi) AVG_pi,
#  STD(isw.isw_pi) STD_pi,
#  AVG(isw.isw_average_gc) AVG_gc,
#  STD(isw.isw_average_gc) STD_gc,
#  AVG(isw.isw_cv) AVG_cv,
#  STD(isw.isw_cv) STD_cv,
#  COUNT(*) COUNT
#FROM isw
#  INNER JOIN indel ON
#    isw.isw_indel_id = indel.indel_id
#WHERE (indel.indel_slippage >= ?)
#  AND (indel.indel_slippage <= ?)
{
    my $name = 'common-d1_comb_slippage-0';

    my $sql = ns();
    $sql->add_select( 'AVG(isw.isw_distance)',   'AVG_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',         'AVG_pi' );
    $sql->add_select( 'STD(isw.isw_pi)',         'STD_pi' );
    $sql->add_select( 'AVG(isw.isw_average_gc)', 'AVG_gc' );
    $sql->add_select( 'STD(isw.isw_average_gc)', 'STD_gc' );
    $sql->add_select( 'AVG(isw.isw_cv)',         'AVG_cv' );
    $sql->add_select( 'STD(isw.isw_cv)',         'STD_cv' );
    $sql->add_select( 'COUNT(*)',                'COUNT' );

    $sql->add_join(
        isw => {
            type      => 'inner',
            table     => 'indel',
            condition => 'isw.isw_indel_id = indel.indel_id',
        }
    );
    $sql->add_where( 'indel.indel_slippage' => { op => '>=', value => '1' } );
    $sql->add_where( 'indel.indel_slippage' => { op => '<=', value => '1' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT CONCAT(isw_type, isw_distance) isw_type_distance,
#       AVG(isw_pi) AVG_pi,
#       COUNT(isw_pi) COUNT,
#       STD(isw_pi) STD_pi
#FROM isw
#WHERE isw_type = ?
#AND isw_density BETWEEN ? AND ?
#AND isw_distance <= (? + 1) / 2
#GROUP BY CONCAT(isw_type, isw_distance)
#ORDER BY isw_distance
{
    my $sql = ns();
    $sql->add_select( 'CONCAT(isw_type, isw_distance)', 'isw_type_distance' );
    $sql->add_select( 'AVG(isw_pi)',                    'AVG_pi' );
    $sql->add_select( 'COUNT(*)',                       'COUNT' );
    $sql->add_select( 'STD(isw_pi)',                    'STD_pi' );
    $sql->from( ['isw'] );
    $sql->add_where( 'isw_type'     => 'L' );
    $sql->add_where( 'isw_density'  => { op => '>=', value => '1' } );
    $sql->add_where( 'isw_density'  => { op => '<=', value => '2' } );
    $sql->add_where( 'isw_distance' => \'<= (? + 1) / 2' );

    $sql->group( { column => 'CONCAT(isw_type, isw_distance)' } );
    $sql->order( { column => 'isw_distance' } );

    $sql_file->set( 'common-dd_group-4', $sql );
    print $sql->as_sql if $verbose;
}

#SELECT CONCAT(isw.isw_type, isw.isw_distance) isw_type_distance,
#       AVG(isw_pi) AVG_pi,
#       COUNT(isw_pi) COUNT,
#       STD(isw_pi) STD_pi
#FROM indel INNER JOIN isw ON indel.indel_id = isw.indel_id
#WHERE isw.isw_density > 9
#AND isw.isw_distance <= 5
#AND isw.isw_type = 'R'
##AND indel.indel_length >= ?
##AND indel.indel_length <= ?
#GROUP BY CONCAT(isw.isw_type, isw.isw_distance) DESC
{
    my $sql = ns();
    $sql->add_select( 'CONCAT(isw_type, isw_distance)', 'isw_type_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',                'AVG_pi' );
    $sql->add_select( 'COUNT(*)',                       'COUNT' );
    $sql->add_select( 'STD(isw.isw_pi)',                'STD_pi' );
    $sql->add_where( 'isw.isw_density'  => \'> 9' );
    $sql->add_where( 'isw.isw_distance' => \'<= 5' );

    #$sql->add_where( 'indel.indel_length' => { op => '>=', value => '1' } );
    #$sql->add_where( 'indel.indel_length' => { op => '<=', value => '5' } );

    my $sql_R = $sql->copy;
    $sql_R->add_join(
        indel => [
            {   type      => 'inner',
                table     => 'isw',
                condition => 'indel.indel_id = isw.indel_id',
            },
        ]
    );
    $sql_R->add_where( 'isw.isw_type' => \'= \'R\'' );
    $sql_R->group(
        {   column => 'CONCAT(isw.isw_type, isw.isw_distance)',
            desc   => 'DESC'
        }
    );

    my $sql_L = $sql->copy;
    $sql_L->add_join(
        indel => [
            {   type      => 'inner',
                table     => 'isw',
                condition => 'indel.indel_id = isw.prev_indel_id',
            },
        ]
    );
    $sql_L->add_where( 'isw.isw_type' => \'= \'L\'' );
    $sql_L->group( { column => 'CONCAT(isw.isw_type, isw.isw_distance)' } );

    $sql_file->set( 'common-indel_size_r-0', $sql_R );
    $sql_file->set( 'common-indel_size_l-0', $sql_L );
    print $sql_R->as_sql if $verbose;
    print $sql_L->as_sql if $verbose;
}

#SELECT
#  CONCAT(isw_type, isw_distance) isw_type_distance,
#  AVG(isw.isw_pi) AVG_pi,
#  COUNT(*) COUNT,
#  STD(isw.isw_pi) STD_pi
#FROM indel
#  INNER JOIN isw ON
#    indel.indel_id = isw.indel_id
#WHERE (isw.isw_density > 9)
#  AND (isw.isw_distance <= 5)
#  AND (isw.isw_type = 'R')
#GROUP BY
#  CONCAT(isw.isw_type, isw.isw_distance) DESC
{
    my $sql = ns();
    $sql->add_select( 'CONCAT(isw_type, isw_distance)', 'isw_type_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',                'AVG_pi' );
    $sql->add_select( 'COUNT(*)',                       'COUNT' );
    $sql->add_select( 'STD(isw.isw_pi)',                'STD_pi' );
    $sql->add_where( 'isw.isw_density'  => \'> 9' );
    $sql->add_where( 'isw.isw_distance' => \'<= 5' );

    my $sql_R = $sql->copy;
    $sql_R->add_join(
        indel => [
            {   type      => 'inner',
                table     => 'isw',
                condition => 'indel.indel_id = isw.indel_id',
            },
        ]
    );
    $sql_R->add_where( 'isw.isw_type' => \'= \'R\'' );
    $sql_R->group(
        {   column => 'CONCAT(isw.isw_type, isw.isw_distance)',
            desc   => 'DESC'
        }
    );

    my $sql_L = $sql->copy;
    $sql_L->add_join(
        indel => [
            {   type      => 'inner',
                table     => 'isw',
                condition => 'indel.indel_id = isw.prev_indel_id',
            },
        ]
    );
    $sql_L->add_where( 'isw.isw_type' => \'= \'L\'' );
    $sql_L->group( { column => 'CONCAT(isw.isw_type, isw.isw_distance)' } );

    $sql_file->set( 'common-indel_feature_r-0', $sql_R );
    $sql_file->set( 'common-indel_feature_l-0', $sql_L );
    print $sql_R->as_sql if $verbose;
    print $sql_L->as_sql if $verbose;
}

#SELECT  indel_length,
#        COUNT(*) indel_number,
#        AVG(indel_gc) AVG_gc,
#        SUM(indel_length) indel_sum
#FROM indel
#GROUP BY indel_length
{
    my $sql = ns();
    $sql->add_select('indel_length');
    $sql->add_select( 'COUNT(*)',          'indel_number' );
    $sql->add_select( 'AVG(indel_gc)',     'AVG_indel_gc' );
    $sql->add_select( 'SUM(indel_length)', 'indel_sum' );
    $sql->from( ['indel'] );
    $sql->group( { column => 'indel_length' } );

    $sql_file->set( 'common-indel_length-0', $sql );
    print $sql->as_sql if $verbose;
}

#SELECT
#  isw.isw_distance distance,
#  AVG(isw.isw_pi) AVG_pi,
#  COUNT(*) COUNT,
#  STD(isw.isw_pi) STD_pi
#FROM indel
#  INNER JOIN isw ON
#    isw.indel_id = indel.indel_id
#  INNER JOIN align ON
#    align.align_id = indel.align_id
{
    my $sql = ns();
    $sql->add_select( 'isw.isw_distance', 'distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',  'AVG_pi' );
    $sql->add_select( 'COUNT(*)',         'COUNT' );
    $sql->add_select( 'STD(isw.isw_pi)',  'STD_pi' );

    $sql->add_join(
        indel => [
            {   type      => 'inner',
                table     => 'isw',
                condition => 'isw.indel_id = indel.indel_id',
            },
            {   type      => 'inner',
                table     => 'align',
                condition => 'align.align_id = indel.align_id',
            },
        ]
    );

    $sql->group( { column => 'isw.isw_distance' } );

    $sql_file->set( 'common-align-0', $sql );
    print $sql->as_sql if $verbose;
}

#----------------------------------------------------------#
# multi_stat_factory.pl SQL
#----------------------------------------------------------#

#SELECT isw_distance distance,
#       AVG(isw_pi) AVG_D,
#       AVG(isw_d_indel) AVG_Di,
#       AVG(isw_d_noindel) AVG_Dni,
#       AVG(isw_d_bii)/2 `AVG_Dbii/2`,
#       AVG(isw_d_bnn)/2 `AVG_Dbnn/2`,
#       AVG(isw_d_complex) AVG_Dc,
#       AVG(isw_d_indel)/AVG(isw_d_noindel) `Di/Dn`,
#       COUNT(*) COUNT
#FROM isw s, indel i
#WHERE s.isw_indel_id = i.indel_id
#AND i.indel_slippage = 0
#GROUP BY isw_distance
{
    my $sql = ns();
    $sql->add_select( 'isw_distance',                          'distance' );
    $sql->add_select( 'AVG(isw_pi)',                           'AVG_D' );
    $sql->add_select( 'AVG(isw_d_indel)',                      'AVG_Di' );
    $sql->add_select( 'AVG(isw_d_noindel)',                    'AVG_Dn' );
    $sql->add_select( 'AVG(isw_d_bii)/2',                      '`AVG_Dbii/2`' );
    $sql->add_select( 'AVG(isw_d_bnn)/2',                      '`AVG_Dbnn/2`' );
    $sql->add_select( 'AVG(isw_d_complex)',                    'AVG_Dc' );
    $sql->add_select( 'AVG(isw_d_indel) / AVG(isw_d_noindel)', '`Di/Dn`' );
    $sql->add_select( 'COUNT(*)',                              'COUNT' );

    $sql->add_join(
        isw => {
            type      => 'inner',
            table     => 'indel',
            condition => 'isw.isw_indel_id = indel.indel_id',
        }
    );
    $sql->add_where( 'indel.indel_slippage' => \'= 0' );

    $sql->group( { column => 'isw_distance' } );

    $sql_file->set( 'multi-distance-0', $sql );
    print $sql->as_sql if $verbose;
}

{
    my $sql = ns();
    $sql->add_select( 'isw_distance',                            'distance' );
    $sql->add_select( 'AVG(isw_pi)',                             'AVG_D' );
    $sql->add_select( 'AVG(isw_d_indel2)',                       'AVG_Di2' );
    $sql->add_select( 'AVG(isw_d_noindel2)',                     'AVG_Dn2' );
    $sql->add_select( 'AVG(isw_d_bii2)/2',                       '`AVG_Dbii2/2`' );
    $sql->add_select( 'AVG(isw_d_bnn2)/2',                       '`AVG_Dbnn2/2`' );
    $sql->add_select( 'AVG(isw_d_complex2)',                     'AVG_Dc2' );
    $sql->add_select( 'AVG(isw_d_indel2) / AVG(isw_d_noindel2)', '`Di2/Dn2`' );
    $sql->add_select( 'COUNT(*)',                                'COUNT' );

    $sql->add_join(
        isw => {
            type      => 'inner',
            table     => 'indel',
            condition => 'isw.isw_indel_id = indel.indel_id',
        }
    );
    $sql->add_where( 'indel.indel_slippage' => \'= 0' );

    $sql->group( { column => 'isw_distance' } );

    $sql_file->set( 'multi-distance2-0', $sql );
    print $sql->as_sql if $verbose;
}

{
    my $sql = ns();
    $sql->add_select( 'isw_distance',                            'distance' );
    $sql->add_select( 'AVG(isw_pi)',                             'AVG_D' );
    $sql->add_select( 'AVG(isw_d_indel3)',                       'AVG_Di3' );
    $sql->add_select( 'AVG(isw_d_noindel3)',                     'AVG_Dn3' );
    $sql->add_select( 'AVG(isw_d_bii3)/2',                       '`AVG_Dbii3/2`' );
    $sql->add_select( 'AVG(isw_d_bnn3)/2',                       '`AVG_Dbnn3/2`' );
    $sql->add_select( 'AVG(isw_d_complex3)',                     'AVG_Dc3' );
    $sql->add_select( 'AVG(isw_d_indel3) / AVG(isw_d_noindel3)', '`Di3/Dn3`' );
    $sql->add_select( 'COUNT(*)',                                'COUNT' );

    $sql->add_join(
        isw => {
            type      => 'inner',
            table     => 'indel',
            condition => 'isw.isw_indel_id = indel.indel_id',
        }
    );
    $sql->add_where( 'indel.indel_slippage' => \'= 0' );

    $sql->group( { column => 'isw_distance' } );

    $sql_file->set( 'multi-distance3-0', $sql );
    print $sql->as_sql if $verbose;
}

#SELECT  indel_length,
#        COUNT(*) indel_number,
#        AVG(indel_gc) AVG_gc,
#        SUM(indel_length) indel_sum
#FROM indel
#WHERE indel.indel_slippage = 0
#GROUP BY indel_length
{
    my $sql = ns();
    $sql->add_select( 'indel_length',      'indel_length' );
    $sql->add_select( 'COUNT(*)',          'indel_number' );
    $sql->add_select( 'AVG(indel_gc)',     'AVG_gc' );
    $sql->add_select( 'SUM(indel_length)', 'indel_sum' );
    $sql->from( ['indel'] );
    $sql->add_where( 'indel.indel_slippage' => \'= 0' );
    $sql->group( { column => 'indel_length' } );

    $sql_file->set( 'multi-indel_length-0', $sql );
    print $sql->as_sql if $verbose;
}

#----------------------------------------------------------#
# gene dnds
#----------------------------------------------------------#

#SELECT
#    AVG(i.isw_distance) AVG_distance,
#    AVG(i.isw_pi) AVG_pi,
#    AVG(i.isw_syn) AVG_d_syn,
#    AVG(i.isw_nsy) AVG_d_nsy,
#    AVG(i.isw_stop) AVG_d_stop,
#    COUNT(*) COUNT,
#    AVG(i.isw_nsy) / AVG(i.isw_syn) `dn/ds`
#FROM
#    isw i
{
    my $name = 'dnds-d1_comb_dn_ds-0';

    my $sql = ns();
    $sql->add_select( 'AVG(isw.isw_distance)',               'AVG_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',                     'AVG_pi' );
    $sql->add_select( 'AVG(isw.isw_syn)',                    'AVG_d_syn' );
    $sql->add_select( 'AVG(isw.isw_nsy)',                    'AVG_d_nsy' );
    $sql->add_select( 'AVG(isw.isw_stop)',                   'AVG_d_stop' );
    $sql->add_select( 'COUNT(*)',                            'COUNT' );
    $sql->add_select( 'AVG(isw.isw_nsy) / AVG(isw.isw_syn)', '`dn/ds`' );
    $sql->from( ['isw'] );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

{
    my $name = 'dnds-d1_dn_ds-0';

    my $sql = ns();
    $sql->add_select( 'isw.isw_distance',                    'isw_distance' );
    $sql->add_select( 'AVG(isw.isw_pi)',                     'AVG_pi' );
    $sql->add_select( 'AVG(isw.isw_syn)',                    'AVG_d_syn' );
    $sql->add_select( 'AVG(isw.isw_nsy)',                    'AVG_d_nsy' );
    $sql->add_select( 'AVG(isw.isw_stop)',                   'AVG_d_stop' );
    $sql->add_select( 'COUNT(*)',                            'COUNT' );
    $sql->add_select( 'AVG(isw.isw_nsy) / AVG(isw.isw_syn)', '`dn/ds`' );
    $sql->from( ['isw'] );
    $sql->group( { column => 'isw_distance' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#----------------------------------------------------------#
# ld_stat_factory.pl SQL
#----------------------------------------------------------#

#SELECT
#  isw.isw_distance distance,
#  AVG(snp.snp_r) AVG_r,
#  AVG(POWER(snp.snp_r, 2)) AVG_r2,
#  AVG(snp.snp_dprime) AVG_Dprime,
#  AVG(ABS(snp.snp_dprime)) AVG_Dprime_abs,
#  COUNT(*) COUNT
#FROM isw
#  INNER JOIN indel ON
#    isw.isw_indel_id = indel.indel_id
#  INNER JOIN snp ON
#    isw.isw_id = snp.isw_id
#WHERE (indel.indel_occured != 'unknown')
#  AND (snp.snp_occured != 'unknown')
#GROUP BY
#  isw.isw_distance
{
    my $name = 'ld-indel_ld-0';

    my $sql = ns();
    $sql->add_select( 'isw.isw_distance',         'distance' );
    $sql->add_select( 'AVG(snp.snp_r)',           'AVG_r' );
    $sql->add_select( 'AVG(POWER(snp.snp_r, 2))', 'AVG_r2' );
    $sql->add_select( 'AVG(snp.snp_dprime)',      'AVG_Dprime' );
    $sql->add_select( 'AVG(ABS(snp.snp_dprime))', 'AVG_Dprime_abs' );
    $sql->add_select( 'COUNT(*)',                 'COUNT' );

    $sql->add_join(
        isw => [
            {   type      => 'inner',
                table     => 'indel',
                condition => 'isw.isw_indel_id = indel.indel_id',
            },
            {   type      => 'inner',
                table     => 'snp',
                condition => 'isw.isw_id = snp.isw_id',
            },
        ]
    );
    $sql->add_where( 'indel.indel_occured' => \"!= 'unknown'" );
    $sql->add_where( 'snp.snp_occured'     => \"!= 'unknown'" );
    $sql->group( { column => 'isw.isw_distance' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

#SELECT
#  isw.isw_distance distance,
#  AVG(POWER(snp.snp_r, 2)) AVG_r2,
#  AVG(snp.snp_r2_s) AVG_r2_s,
#  AVG(ABS(snp.snp_dprime)) AVG_Dprime_abs,
#  AVG(snp.snp_dprime_abs_s) AVG_Dprime_abs_s,
#  AVG(snp.snp_r2_i) AVG_r2_i,
#  AVG(snp.snp_r2_ni) AVG_r2_ni,
#  AVG(snp.snp_dprime_abs_i) AVG_Dprime_abs_i,
#  AVG(snp.snp_dprime_abs_ni) AVG_Dprime_abs_ni,
#  COUNT(*) COUNT
#FROM isw
#  INNER JOIN indel ON
#    isw.isw_indel_id = indel.indel_id
#  INNER JOIN snp ON
#    isw.isw_id = snp.isw_id
#WHERE (indel.indel_occured != 'unknown')
#  AND (snp.snp_occured != 'unknown')
#GROUP BY
#  isw.isw_distance
{
    my $name = 'ld-snps_ld-0';

    my $sql = ns();
    $sql->add_select( 'isw.isw_distance',           'distance' );
    $sql->add_select( 'AVG(POWER(snp.snp_r, 2))',   'AVG_r2' );
    $sql->add_select( 'AVG(snp.snp_r2_s)',          'AVG_r2_s' );
    $sql->add_select( 'AVG(ABS(snp.snp_dprime))',   'AVG_Dprime_abs' );
    $sql->add_select( 'AVG(snp.snp_dprime_abs_s)',  'AVG_Dprime_abs_s' );
    $sql->add_select( 'AVG(snp.snp_r2_i)',          'AVG_r2_i' );
    $sql->add_select( 'AVG(snp.snp_r2_ni)',         'AVG_r2_ni' );
    $sql->add_select( 'AVG(snp.snp_dprime_abs_i)',  'AVG_Dprime_abs_i' );
    $sql->add_select( 'AVG(snp.snp_dprime_abs_ni)', 'AVG_Dprime_abs_ni' );
    $sql->add_select( 'COUNT(*)',                   'COUNT' );

    $sql->add_join(
        isw => [
            {   type      => 'inner',
                table     => 'indel',
                condition => 'isw.isw_indel_id = indel.indel_id',
            },
            {   type      => 'inner',
                table     => 'snp',
                condition => 'isw.isw_id = snp.isw_id',
            },
        ]
    );
    $sql->add_where( 'indel.indel_occured' => \"!= 'unknown'" );
    $sql->add_where( 'snp.snp_occured'     => \"!= 'unknown'" );
    $sql->group( { column => 'isw.isw_distance' } );

    $sql_file->set( $name, $sql );

    print "\n[$name]\n";
    print $sql->as_sql if $verbose;
}

END {
    $sql_file->write;
}
