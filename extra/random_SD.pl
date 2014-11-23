#!/usr/bin/perl
use strict;
use warnings;

use DBI;

use Getopt::Long;
use Pod::Usage;

my $help        = 0;
my $man         = 0;
my $server      = 'localhost';
my $port        = 3306;
my $db_file     = 'random_db';
my $username    = 'alignDB';
my $password    = 'alignDB';
my $out_file    = 'random_SD';

GetOptions(
    'help|?'        => \$help,
    'server=s'      => \$server,
    'Port=i'        => \$port,
    'db=s'          => \$db_file,
    'username=s'    => \$username,
    'password=s'    => \$password,
    'output=s'      => \$out_file
)
or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

$|++;

open LIST, "< $db_file";
chomp (my @data = <LIST>);
close LIST;

open OUT, "> $out_file";
open OUT2, '>temp';

foreach my $window ((-1..5)){
    my $count = 0;
    my @pi;
    foreach my $db (@data){
        print "$db\n";
        
        my $dbh = DBI->connect("dbi:mysql:$db:$server",$username,$password);
        
        my $fetch_sql = $dbh->prepare(
            'SELECT isw_pi
            FROM isw i
            WHERE isw_distance = ?'
        );
        $fetch_sql->execute($window);
        
        my $db_cnt = 0;
        
        while(my ($pi) = $fetch_sql->fetchrow_array){
            push @pi,$pi;
            $count++;
            $db_cnt++;
        }
        
        print OUT2 "$db,$window,$db_cnt\n";
    }
    
    my $avg = &average(@pi);
    my $std = &std(@pi);
    
    print OUT "$window,$avg,$std,$count\n";
}
close OUT2;
close OUT;

sub std{
    my @number = @_;
    my @square;
    
    my $avg = &average(@number);
    my $cnt = scalar @number;
    
    foreach (@number){
        push @square, ($_-$avg)*($_-$avg);
    }
    
    return sqrt(&sum(@square)/($cnt-1));
}

sub average{
    my @number      = @_;
    my $total_cnt   = scalar @number;
    
    return &sum(@number)/$total_cnt;
}

sub sum{
    my @number  = @_;
    my $sum     = shift @number;
    
    foreach (@number){
        $sum+=$_;
    }
    
    return $sum;
}