#extract the sequences of the maximum likelihood set of paths
use Getopt::Long;
$fastafile ="";
$likelihoodfile = "";

GetOptions("f=s"=>\$fastafile,"l=s"=>\$likelihoodfile);

main();

sub main
{
	loadlikelihoodfile();
	#print "Number of elements is $set_val\n";
	extract_mle_sequences();
}
sub loadlikelihoodfile
{
	open(file,$likelihoodfile);
	$max = -9**9**9;
	while($l=<file>)
	{
		chomp($l);
		@vs = split(/ /,$l);
		if($vs[1] > $max)
		{
			$max = $vs[1];
			$set_val = $vs[0];
			# Number of elements in the MLE set
		}
		$id{$vs[2]} = $vs[0];
		#print "$vs[0] $vs[1] $vs[2] $id{$vs[2]} \n";
	}
	close file;
}

sub extract_mle_sequences
{
	open(file,$fastafile);
	while($l=<file>)
	{
		chomp($l);
		$temp = substr($l,1);
		$seq = <file>;
		chomp($seq);
		($temp_id,$length) = split(/__/,$temp);
		if($id{$temp_id} > $set_val)
		{
			#print "ID $temp_id $set_val \n";		
			delete $id{$temp_id};
		}
	}
	close file;
	open(file,$fastafile);
	while($l=<file>)
	{
		chomp($l);
		$temp = substr($l,1);
		$seq = <file>;
		chomp($seq);
		($temp_id,$length) = split(/__/,$temp);
		if(exists($id{$temp_id}))
		{
			print "$l\n$seq\n";
		}
	}
	close file;	
}