# Remove short contigs from the list of contigs
# USAGE perl remove_short_contigs.pl <Fasta_file> <size_to_remove> 
open(file,$ARGV[0]);
while($l=<file>)
{
	chomp($l);
	$seq =<file>;
	chomp($seq);
	#($number,$length) = split(/__/,$l);
	if(length($seq)  > $ARGV[1])
	{
		print "$l\n$seq\n";
	}
}
close file;
