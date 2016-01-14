open(file,$ARGV[0]);
%hash = ();
while($l=<file>)
{
	chomp($l);
	($kmer,$v)=split(/ /,$l);
	$temp = revcomplement($kmer);
	if(($kmer cmp $temp)<0) 
	{
		$hash{$kmer} = $hash{$kmer} + $v;
#		print "$kmer $temp $hash{$kmer} $v\n";
	}else
	{
		$hash{$temp} = $hash{$temp} + $v;
#		print "$kmer $temp $hash{$temp} $v\n";
	}
}
close file;


for (keys %hash) 
{
#	print "$_ ".revcomplement($_)." $hash{$_}\n";
	print "$_ $hash{$_}\n";
}
sub revcomplement
{
	my($str) = $_[0];
	$str = scalar reverse $str;
	$rev_str="";
	for (my($i)=0;$i<length($str);$i++) {
		if(substr($str,$i,1) eq "A") {$rev_str = $rev_str."T";} 
		elsif(substr($str,$i,1) eq "T") {$rev_str = $rev_str."A";}
		elsif(substr($str,$i,1) eq "C") {$rev_str = $rev_str."G";}
		else {$rev_str = $rev_str."C";}
	}
	return $rev_str;
}

