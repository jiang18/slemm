# sim_snp_effects.pl

use strict;
use warnings;

my @snps;
my %freq;
open IN,"true.frequency";
while(<IN>)
{
        chomp;
        my @c = split /\s+/;
        push @snps, $c[0];
        $freq{$c[0]} = $c[2];
}
close IN;

my $PI = 2 * atan2 1, 0;
my @nums = map {
    sqrt(-2 * log rand) * cos(2 * $PI * rand)
} 1..@snps;

@ARGV == 1 or die "One argument needed: output filename\n";
my ($out) = @ARGV;

open OUT,">$out";
print OUT "SNP,Chr,Pos,Allele1,Allele2,MAF,HWE_Pval,Group,Weight,Effect\n";

for(0..$#snps) {
        my $f = $freq{$snps[$_]};
        my $e = $nums[$_]/sqrt(2*$f*(1-$f));
        print OUT "$snps[$_],0,0,1,2,$f,0,NULL,1,$e\n";
}
close OUT;
