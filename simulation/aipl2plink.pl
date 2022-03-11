# aipl2plink.pl

use strict;
use warnings;

@ARGV == 1 or die "One argument needed: plink file prefix\n";
my ($plink) = @ARGV;

open IN,"genotypes.true";
open OUT,">$plink.ped";
while(<IN>)
{
        chomp;
        my @c = split /\s+/;
        my @g = split //,$c[-1];
        print OUT "0 $c[1] 0 0 2 0";
        
        $c[-1] =~ s/2/ 22/g;
	$c[-1] =~ s/1/ 12/g;
	$c[-1] =~ s/0/ 11/g;
        print OUT $c[-1],"\n";
}
close IN;
close OUT;

open IN,"chromosome.data";
open OUT,">$plink.map";
$_=<IN>;
while(<IN>)
{
        chomp;
        my @c = split /\s+/;
        print OUT "$c[1]\t$c[0]\t0\t$c[4]\n";
}
close IN;
close OUT;
