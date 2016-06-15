
open(file,$ARGV[0]);
$count=1;
$read="";
$l=<file>;
while($l=<file>)
{
	chomp($l);
	if(substr($l,0,1) ne ">")
	{
		$read=$read.$l;
	}
	else 
	{
		print ">$count"."_".length($read)."\n$read\n";
		$count++;
		$read="";
	}
}
close file;

#print "Reads in the file are $count\n";
