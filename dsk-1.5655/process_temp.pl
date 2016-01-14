open(file,$ARGV[0]); #open the temp file and process the k-mers visited
$count=0;
while($l=<file>)
{
	chomp($l);
	@vs = split(/ /,$l);
	#if(scalar(@vs)==3)
	if(scalar(@vs)==1)
	{
		$hash{$vs[0]}++;
#		$pass{$vs[2]}++;	
	}
	$count++;
	#if($count%1000000==0) { print $count." ".time."\n"; }
}
close file;

for $y ( sort {$a<=>$b} keys %hash )
{
	print "$y $hash{$y}\n";

}
#for $y( sort {$a<=>$b}  keys %pass )
#{
#	print "$y $pass{$y}\n";
#}
