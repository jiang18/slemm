# sim_ped.pl

use strict;
use warnings;

@ARGV == 1 or die "One argument needed: sample size\n";
my ($np) = @ARGV;
my $st = 1000000;

open OUT,">pedigree.file";
for(1..$np){
        my $id = $_+$st;
        print OUT "F     $id        -1       -2 20000101  0.0    0\n";
}
close OUT;

open OUT,">genotype.data0";
print OUT "idnum  chip\n\n";
for(1..$np){
        my $id = $_+$st;
        print  OUT "    $id    1\n";
}
close OUT;
