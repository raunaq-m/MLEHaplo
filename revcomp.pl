# compute reverse complement of reads and print them
open(file,$ARGV[0]);
$read="";
$l=<file>;
$count=0;
while($l=<file>)
{
	chomp($l);
	if(substr($l,0,1) ne ">")
	{
		$read=$read.$l;
	}else
	{
		$count++;
		print ">$count"."_rev\n";
		print revcomplement($read)."\n";
		$read="";
	}
}
$count++;
print ">$count"."_rev\n";
print revcomplement($read)."\n";

sub revcomplement
{
	my($str) = $_[0];
	$str = scalar reverse uc($str);
	$rev_str="";
	for (my($i)=0;$i<length($str);$i++) {
		if(substr($str,$i,1) eq "A") {$rev_str = $rev_str."T";} 
		elsif(substr($str,$i,1) eq "T") {$rev_str = $rev_str."A";}
		elsif(substr($str,$i,1) eq "C") {$rev_str = $rev_str."G";}
		elsif(substr($str,$i,1) eq "G") {$rev_str = $rev_str."C";}
		else {$rev_str = $rev_str."N";}
	}
	return $rev_str;
}

	
