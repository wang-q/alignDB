#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Scalar::Util qw(looks_like_number);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Update isw_dxr of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh();

#----------------------------#
# Add columns
#----------------------------#
# Add columns to isw_extra
{
    # Dir
    $obj->create_column("isw_extra", "isw_feature4", "DOUBLE");
    # Dnr
    $obj->create_column("isw_extra", "isw_feature5", "DOUBLE");
    # Dtr
    $obj->create_column("isw_extra", "isw_feature6", "DOUBLE");
    # Dqr
    $obj->create_column("isw_extra", "isw_feature7", "DOUBLE");
    # Dtotal
    $obj->create_column("isw_extra", "isw_feature8", "DOUBLE");

    print "Table isw_extra altered\n";
}

#----------------------------#
# Prefill tables if necessarily
#----------------------------#
# check rows in table isw_extra
{
    my $sql_query = qq{
        SELECT COUNT(isw_extra_id)
        FROM isw_extra
    };
    my $sth = $dbh->prepare($sql_query);
    $sth->execute();
    my ($count) = $sth->fetchrow_array;

    unless ($count) {
        $sql_query = qq{
            INSERT INTO isw_extra (isw_id)
            SELECT isw.isw_id
            FROM isw
        };
        $sth = $dbh->prepare($sql_query);
        $sth->execute();
    }
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
{

    # alignments
    my $align_query = q{
        SELECT align_id
        FROM align 
    };
    my $align_sth = $dbh->prepare($align_query);

    # sequence
    my $seq_query = q{
        SELECT t.target_seq, q.query_seq, r.ref_seq
        FROM align a, target t, query q, reference r
        WHERE a.align_id = ?
        AND a.align_id = q.align_id
        AND a.align_id = t.align_id
        AND a.align_id = r.align_id
    };
    my $seq_sth = $dbh->prepare($seq_query);

    # select all indels in this alignment
    my $indel_query = q{
        SELECT indel_id, indel_occured
        FROM indel
        WHERE align_id = ?
    };
    my $indel_sth = $dbh->prepare($indel_query);

    # select all isw in this indel
    my $isw_L_query = q{
        SELECT isw_id, isw_start, isw_end, isw_length
        FROM isw
        WHERE indel_id = ?
        AND isw_type = "L"
    };
    my $isw_L_sth = $dbh->prepare($isw_L_query);

    my $isw_R_query = q{
        SELECT isw_id, isw_start, isw_end, isw_length
        FROM isw
        WHERE foregoing_indel_id = ?
        AND isw_type = "R"
    };
    my $isw_R_sth = $dbh->prepare($isw_R_query);

    # update isw table in the new feature column
    my $isw_extra = q{
        UPDATE isw_extra
        SET isw_feature4 = ?,
            isw_feature5 = ?,
            isw_feature6 = ?,
            isw_feature7 = ?,
            isw_feature8 = ?
        WHERE isw_id = ?
    };
    my $isw_extra_sth = $dbh->prepare($isw_extra);

    $align_sth->execute();

    # for each align
    while ( my @row = $align_sth->fetchrow_array ) {
        my ($align_id) = @row;
        print "Processing align_id $align_id\n";

        $seq_sth->execute($align_id);
        my ( $target_seq, $query_seq, $ref_seq ) = $seq_sth->fetchrow_array;

        $indel_sth->execute($align_id);
        while ( my @row = $indel_sth->fetchrow_array ) {
            my ( $indel_id, $indel_occured ) = @row;

            my ( $seq_ref, $seq_indel, $seq_noindel );
            $seq_ref = \$ref_seq;

            if ( $indel_occured eq "T" ) {
                $seq_indel   = \$target_seq;
                $seq_noindel = \$query_seq;
            }
            elsif ( $indel_occured eq "Q" ) {
                $seq_indel   = \$query_seq;
                $seq_noindel = \$target_seq;
            }
            else {
                next;
            }

            $isw_L_sth->execute($indel_id);
            while ( my @row = $isw_L_sth->fetchrow_array ) {
                my ( $isw_id, $isw_start, $isw_end, $isw_length ) = @row;

                # dir & dnr
                my $slice_ref
                    = substr( $$seq_ref, $isw_start - 1, $isw_length );
                my $slice_indel
                    = substr( $$seq_indel, $isw_start - 1, $isw_length );
                my $slice_noindel
                    = substr( $$seq_noindel, $isw_start - 1, $isw_length );

                my $result_indel = &pair_seq_stat( $slice_ref, $slice_indel );
                my $dir = $result_indel->[7];
                my $result_noindel
                    = &pair_seq_stat( $slice_ref, $slice_noindel );
                my $dnr = $result_noindel->[7];

                # dtr & dqr
                my $slice_target
                    = substr( $target_seq, $isw_start - 1, $isw_length );
                my $slice_query
                    = substr( $query_seq, $isw_start - 1, $isw_length );

                my $result_target
                    = &pair_seq_stat( $slice_ref, $slice_target );
                my $dtr          = $result_target->[7];
                my $result_query = &pair_seq_stat( $slice_ref, $slice_query );
                my $dqr          = $result_query->[7];

                # dtotal
                my $result_tq = &pair_seq_stat( $slice_target, $slice_query );
                my $dtq = $result_tq->[7];
                my @each_ds
                    = grep { looks_like_number($_) } ( $dtq, $dtr, $dqr );
                my $valid_d_number = scalar @each_ds;
                my $dtotal =
                    $valid_d_number < 2
                    ? 'NULL'
                    : sum(@each_ds) / $valid_d_number;

                #print Dump {
                #    dir    => $dir,
                #    dnr    => $dnr,
                #    dtr    => $dtr,
                #    dqr    => $dqr,
                #    dtotal => $dtotal,
                #    ds     => \@each_ds,
                #};
                $isw_extra_sth->execute( $dir, $dnr, $dtr, $dqr, $dtotal,
                    $isw_id );
            }

            $isw_R_sth->execute($indel_id);
            while ( my @row = $isw_R_sth->fetchrow_array ) {
                my ( $isw_id, $isw_start, $isw_end, $isw_length ) = @row;

                # dir & dnr
                my $slice_ref
                    = substr( $$seq_ref, $isw_start - 1, $isw_length );
                my $slice_indel
                    = substr( $$seq_indel, $isw_start - 1, $isw_length );
                my $slice_noindel
                    = substr( $$seq_noindel, $isw_start - 1, $isw_length );

                my $result_indel = &pair_seq_stat( $slice_ref, $slice_indel );
                my $dir = $result_indel->[7];
                my $result_noindel
                    = &pair_seq_stat( $slice_ref, $slice_noindel );
                my $dnr = $result_noindel->[7];

                # dtr & dqr
                my $slice_target
                    = substr( $target_seq, $isw_start - 1, $isw_length );
                my $slice_query
                    = substr( $query_seq, $isw_start - 1, $isw_length );

                my $result_target
                    = &pair_seq_stat( $slice_ref, $slice_target );
                my $dtr          = $result_target->[7];
                my $result_query = &pair_seq_stat( $slice_ref, $slice_query );
                my $dqr          = $result_query->[7];

                # dtotal
                my $result_tq = &pair_seq_stat( $slice_target, $slice_query );
                my $dtq = $result_tq->[7];
                my @each_ds
                    = grep { looks_like_number($_) } ( $dtq, $dtr, $dqr );
                my $valid_d_number = scalar @each_ds;
                my $dtotal =
                    $valid_d_number < 2
                    ? 'NULL'
                    : sum(@each_ds) / $valid_d_number;

                $isw_extra_sth->execute( $dir, $dnr, $dtr, $dqr, $dtotal,
                    $isw_id );
            }

        }
    }

    $isw_extra_sth->finish;
    $isw_R_sth->finish;
    $isw_L_sth->finish;
    $indel_sth->finish;
    $align_sth->finish;
}

$stopwatch->end_message();

__END__

=head1 NAME

    update_isw_dxr.pl - Add additional Dir(Dindel_ref), Dnr(Dnoindel_ref)
                            Dtr(Dtarget_ref) and Dqr(Dquery_ref)
                            to table isw_extra
                        Add Dtr and Dqr to table window
    
=head1 SYNOPSIS

    update_isw_dxr.pl [options]
        Options:
            --help              brief help message
            --man               full documentation
            --server            MySQL server IP/Domain name
            --db                database name
            --username          username
            --password          password

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut

