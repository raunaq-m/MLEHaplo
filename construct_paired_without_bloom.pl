# Store information from paired reads by picking up reads in pairs from a read file 

use Getopt::Long;
use Class::Struct;
use Bio::SeqIO;
#use Bloom::Filter;

########################################################################
my($fastafile) = ""; #(paired reads present in the file)
my($kmerfile) = "";
my($thresh) = 0;
my($twofiles) = ''; #If there are two files corresponding to paired reads
my($f1) = "";
my($f2) = "";
$wrfile = $kmerfile."PS.pk.txt";
$store_threshold = 1;
GetOptions( "fasta:s" =>\$fastafile,
	    "kmerfile=s" => \$kmerfile,
	    "thresh:i" => \$thresh,
	    "paired" => \$twofiles,
	    "file1:s" => \$f1, "file2:s" => \$f2,
		"st:i"=>\$store_threshold,
		"wr:s"=>\$wrfile)
	or die("Error in command line arguements \n");
die("USAGE:perl construct_paired_without_bloom.pl -file1 file1.fastq -file2 file2.fastq -paired -kmerfile file1.kvalue -thresh \"number\" -wr \"outputPairedSetfile\"\nOutput is a file that contains pairs of k-mers on a line and the number of times such pair is observed:\nkmer1 kmer2 count\n") if $f1 eq "";
die("USAGE:perl construct_paired_without_bloom.pl -file1 file1.fastq -file2 file2.fastq -paired -kmerfile file1.kvalue -thresh \"number\" -wr \"outputPairedSetfile\"\nOutput is a file that contains pairs of k-mers on a line and the number of times such pair is observed:\nkmer1 kmer2 count\n") if $f2 eq "";
die("USAGE:perl construct_paired_without_bloom.pl -file1 file1.fastq -file2 file2.fastq -paired -kmerfile file1.kvalue -thresh \"number\" -wr \"outputPairedSetfile\"\nOutput is a file that contains pairs of k-mers on a line and the number of times such pair is observed:\nkmer1 kmer2 count\n") if $kmerfile eq "";

#######################################################################
%kmerhash = (); $K=0;
%pairedinfo = ();
loadkmerfile($thresh);
$key_total = scalar(keys %kmerhash);

open(logfile,">log_file.txt");
#$filter = Bloom::Filter->new(error_rate => 0.01, capacity => $key_total*200);

print "$key_total is the total no. of keys in the kmerhash \n";


if (!$twofiles) {
	print logfile "In fasta file \n";
	pairedinfo_single($fastafile);
}else 
{
	print logfile "In paired file \n";
	pairedinfo_twofiles($f1,$f2);
}
print_paired_kmers();

sub pairedinfo_twofiles
{
	my($file1) = $_[0];	my($file2) = $_[1];
	my($num_reads) = 0;
	$seq1 = Bio::SeqIO->new(-file=>"$file1",-format=>"fastq");
	$seq2 = Bio::SeqIO->new(-file=>"$file2",-format=>"fastq");
	print "$seq1 $seq2 files \n";
	while($obj1=$seq1->next_seq)
	{
		$obj2=$seq2->next_seq;
		$r1 = $obj1->seq;
		$r2 = $obj2->seq;
		add_paired_reads($r1,$r2);
		$num_reads++;
		if($num_reads %1000 ==0 ) { 
			print "Processed $num_reads in ".time." \n" ; 
		}
	}
}

sub pairedinfo_single
{
	my($fas) = $_[0];
	my($num_reads) = 0;
	open(file,$fas);
	while($l1=<file>) {
		chomp($l1);
		$r1 = <file>; chomp($r1);
		$l2 = <file>; chomp($l2);
		$r2 = <file>; chomp($r2);
#		print length($r1)." ".length($r2)."\n";
		add_paired_reads($r1,$r2);
		$num_reads++;
		if($num_reads %1000 ==0 ) { 
			print "Processed $num_reads in ".time." \n" ; 
		}
	}
	close(file);
}
sub add_paired_reads
{
	my($rr1) = $_[0]; my($rr2) = $_[1]; my($len) = length($rr1);
	my($revrr1) = revcomplement($rr1); my($revrr2) = revcomplement($rr2);
#	my($rev_1) = scalar reverse $revrr1; my($rev_2) = scalar reverse $revrr2; 
	#Adding paired k-mers from the read that are at a fixed distance to each other. 
	#For example if the reads are 150 bps each , rr1 being forward strand and rr2 being reverse strand, 
	#then add the kmer pair (1st from rr1, 1st from revcomplement(rr2) ) and so on 
	for(my($i)=0; $i < $len - $K+1; $i++)
	{
		for(my($j) = $i+5; $j< $len - $K+1; $j++ ) 
		{
			my($p1a) = substr($rr1,$i,$K);
			my($p1b) = substr($rr1,$j,$K);
			my($p2a) = substr($revrr1,$len -$K - $i,$K);
			my($p2b) = substr($revrr1,$len -$K - $j,$K);
			if(!($p1a =~ /N/ || $p1b =~ /N/ || $p2b =~ /N/ || $p2a =~ /N/))
			{	
				pairedinfopair($p1a,$p2a,$p1b,$p2b);
			}
			my($p1a) = substr($rr2,$i,$K);
			my($p1b) = substr($rr2,$j,$K);
			my($p2a) = substr($revrr2,$len -$K - $i,$K);
			my($p2b) = substr($revrr2,$len -$K - $j,$K);
			if(!($p1a =~ /N/ || $p1b =~ /N/ || $p2b =~ /N/ || $p2a =~ /N/))
			{
				pairedinfopair($p1a,$p2a,$p1b,$p2b);
			}
		}
	}
	for(my($i)=0; $i<$len-$K+1; $i++) {
		my($p1a) = substr($rr1,$i,$K);
		my($p2a) = substr($revrr1,$len - $K - $i,$K);
#		my($j) = $i ; 
		for(my($j)=0;$j<$len-$K+1;$j++) {
			my($p1b) = substr($revrr2,$j,$K);
			my($p2b) = substr($rr2,$len - $K - $j,$K);
			if(!($p1a =~ /N/ || $p1b =~ /N/ || $p2b =~ /N/ || $p2a =~ /N/))
			{
				pairedinfopair($p1a,$p2a,$p1b,$p2b);
			}
		}
	}
}
sub pairedinfopair
{
	my($p1a) = $_[0];
	my($p2a) = $_[1];
	my($p1b) = $_[2];
	my($p2b) = $_[3];
	if( (exists($kmerhash{$p1a}) | exists($kmerhash{$p2a})) && (exists($kmerhash{$p1b}) | exists($kmerhash{$p2b}) ) ) 
	{
		# Add kmer pair to the paired hash only if the bloom filter says they are present. 
		# Bloom filter version is very slow!!! Removing bloom filter implementation for simulation datasets 
		#if($filter->check($p1a.$p1b) || $filter->check($p2b.$p2a) )
		#{
			$pairedinfo{$p1a}{$p1b} ++;
			$pairedinfo{$p2b}{$p2a} ++;
		#}else 
		#{
		#	$filter->add($p1a.$p1b);
		#	$filter->add($p2b.$p2a);
		#}
	}
}

sub print_paired_kmers
{
	open(file,">$wrfile");
	$count = 0;
	for $p1(keys %pairedinfo)
	{
		$count += scalar(keys %{$pairedinfo{$p1}});
		for $p2(keys %{$pairedinfo{$p1}}) 
		{
			if($pairedinfo{$p1}{$p2} > $store_threshold) 
			{
				print file "$p1 $p2 $pairedinfo{$p1}{$p2}\n";
			}
		}
	}
	close file;
	print logfile "Total pairs stored is $count\n";
}

sub revcomplement
{
	my($str) = $_[0];
	$str = scalar reverse $str;
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
sub loadkmerfile
{
	#only load the kmers that are above a certain threshold 
	my($t) = $_[0];
	my($iter) = 0;
	open(file,$kmerfile);
	while($l=<file>)
	{
		chomp($l);
		($kmer,$c) = split(/ /,$l);
		if($c>$t) 
		{
			$kmerhash{$kmer} = $c;
		}
	}
	$K = length($kmer);
	close file;
}
