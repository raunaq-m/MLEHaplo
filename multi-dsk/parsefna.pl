#This file will read through the FNA file and find out the frequencies of l-tuples present in all the reads/sequences in the fasta file.
# Usage: perl parsefasta.pl <filename> <LengthofTuple>
open(file,@ARGV[0]);  #fasta file
$L= @ARGV[1]; #k-mer length 
$line = <file>;
$count=0;
$read = "";
%tupleh = ();
%freq = ();
$sum =0;
open(file1,">$ARGV[2]"); # Write file
while($line=<file>){
	chomp($line);
	if(!(substr($line,0,1) eq ">")){
#		$read = join("",$read,substr($line,0,-1));
		$read = join("",$read,$line);
	}else{
	$read=uc($read);
	if(length($read)-$L >=0) {
		for($i=0;$i<length($read)-$L+1;$i++){
			$tmp = substr($read,$i,$L);
			if(index($tmp,'N')==-1) {
				$count++;
				$tupleh{$tmp}++;  #contains counts of tuples from all the reads 
			}
#			print $tmp."\n";
		}
## CHANGE FOR SINGLE READ Print out the tupleh variable right now to obtain counts for the given read and initialize tupleh to () again
	}    
	$read ="";
	}
}
if(length($read)-$L >=0)
{
	for($i=0;$i<length($read)-$L+1;$i++)
	{
		$tmp = substr($read,$i,$L);
		if(index($tmp,'N')==-1)
		{
			$count++;
			$tupleh{$tmp}++;
		}
	}
	$read="";
}

for $k(sort {$tupleh{$a}<=>$tupleh{$b}} keys %tupleh){
	print file1 "$k $tupleh{$k}\n";
#	$sum = $sum + scalar($tupleh{$k});
	$freq{$tupleh{$k}}++; #counts the number of l-tuples that have a particular count
}
#for $k(sort {$a<=>$b} keys %freq){
#	print "$k ".(scalar($freq{$k})*scalar($k))."\n";
#	print "$k $freq{$k}\n";
#}
#print $sum."\n";
close file1;

