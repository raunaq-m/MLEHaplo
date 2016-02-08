#Open an alignment file and trim it from given to given bases
#USAGE perl trim_alignment.pl <filename_fasta> <start_pos> <end_pos>

open(file,$ARGV[0]); #Open fasta file

$srt = $ARGV[1]; $end=$ARGV[2];
while($l=<file>)
{
	chomp($l);
	$seq = <file>;
	chomp($seq);
#	print "$l\n";
#	print substr($seq,$srt-1,$end-$srt+1)."\n";
	$sub_seq = substr($seq,$srt-1,$end-$srt+1);
	push @{$read_seq{$sub_seq}} , $l;
}
close file;
for $k(keys %read_seq)
{
	$total = scalar(@{$read_seq{$k}});
	$id = "";
	foreach (@{$read_seq{$k}})
	{
		$id=$id." ".substr($_,1);
	}	
	print ">$id\n$k\n";
} 

