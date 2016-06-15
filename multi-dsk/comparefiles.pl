# compare the counts from two datasets
open(file,"$ARGV[0]");
while($l=<file>)
{
	chomp($l);
	($kmer,$c)=split(/ /,$l);
	#$hash{substr($kmer,-8)} ++;
	$hash{$kmer} = $c;
}
close file;

open(file,$ARGV[1]);
while($l=<file>)
{
	chomp($l);
	($kmer,$c)=split(/ /,$l);
#	if($hash{$kmer}==$c) { print "$kmer $hash{$kmer} $c\n";}
#	if(!exists($hash{$kmer})) { print "$k\n"; }
	#$hashm{substr($kmer,-8)}++;
	$hashm{$kmer} =$c;
}
for $k (keys %hash)
{
	if($hash{$k} !=$hashm{$k}) {
		print "$k $hash{$k} $hashm{$k}\n";
	}
} 
close file;
