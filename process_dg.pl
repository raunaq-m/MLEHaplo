open(file,$ARGV[0]); #file that contains dg.fasta, the cover 
$l=<file>;
while(substr($l,0,1) ne ">")
{
	chomp($l);
	$l =<file>;
}
print "$l";
while($l=<file>)
{
	chomp($l);
	if (!($l =~ /factor/ | $l =~ /^1:/) )
	{
		print "$l\n";
	}
}
close file;
