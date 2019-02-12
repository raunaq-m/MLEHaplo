# use lib "/Users/raunaq/perl5/lib/perl5/";
use Bio::SeqIO;

$seqio=Bio::SeqIO->new(-file=>"$ARGV[0]",-format=>"fastq");
$filetostore = $ARGV[0];
$in = index($filetostore,"fq");
if($in==-1) {
	$filetostore =~ s/fastq/fasta/;
}else
{
	$filetostore =~ s/fq/fasta/;
}
print $filetostore."\n";
open(file,">$filetostore");

while($obj = $seqio->next_seq){
	$id = $obj->id;
	$seq = $obj->seq;
	print file ">$id\n";
	print file $seq."\n";
}
print "Fasta file written in $filetostore\n";
# The file converts the given fastq file to fasta file 

