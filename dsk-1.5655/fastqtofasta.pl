use Bio::SeqIO;

$seqio=Bio::SeqIO->new(-file=>"$ARGV[0]",-format=>"fastq");
$filetostore = substr($ARGV[0],0,index($ARGV[0],"fastq"))."fasta";
open(file,">$filetostore");

while($obj = $seqio->next_seq){
	$id = $obj->id;
	$seq = $obj->seq;
	print file ">$id\n";
	print file $seq."\n";
}
print "Fasta file written in $filetostore\n";
# The file converts the given fastq file to fasta file 

