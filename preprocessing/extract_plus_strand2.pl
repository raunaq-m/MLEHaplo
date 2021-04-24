open(file,$ARGV[0]); #sam
open(file1, $ARGV[1]); #karect_whole.fa
open(file2,">$ARGV[2]");#plus_strand.fa
open(file3,">$ARGV[3]"); # new_paired_reads.fa
my %strands;
my %sequences;
while($l=<file>)
{
	chomp($l);
	@txt = split(/\t/,$l);
	@seq = split (/\//,$txt[0]);
	$strands{$seq[0]}{int($txt[1])} = $seq[1];	
}
close file;

while($line=<file1>)
{
	chomp($line);
	if(substr($line,0,1) eq ">"){
		$h=$line;
	}else{
		$sequences{$h} = $line;	
	}
}
close file1;

my $pairedno;
my $singleno;
my $rcno;
foreach $k (keys %strands){
	if(exists($strands{$k}{0})){
		$header= ">".$k."/".$strands{$k}{0};
		$s = $sequences{$header};
		print file2 "$header\n$s\n";		
	}elsif(exists($strands{$k}{16})){
		$header= ">".$k."/".$strands{$k}{16};
		$s = $sequences{$header};
		$s2 = revcomplement($s);
		$rcno++;
		print file2 "$header\n$s2\n";
	}else{
		print "0/16 not exist\n";
	}
	
	if(exists($strands{$k}{0}) & exists($strands{$k}{16})){
		$pairedno++;
		$header1= ">".$k."/".$strands{$k}{0};
		$s1 = $sequences{$header1};
		$header2= ">".$k."/".$strands{$k}{16};
		$s2 = $sequences{$header2};
		print file3 "$header1\n$s1\n$header2\n$s2\n";
	}else{
		$singleno++;
		if(exists($strands{$k}{0})){
			$header1= ">".$k."/".$strands{$k}{0};
			$s1 = $sequences{$header1};
			$no = int(2/int($strands{$k}{0}));
			$header2 = ">".$k."/".$no;
			$s2 = revcomplement($s1);
			print file3 "$header1\n$s1\n$header2\n$s2\n";
		}elsif(exists($strands{$k}{16})){
			$header2= ">".$k."/".$strands{$k}{16};
			$s2 = $sequences{$header2};
			$no = int(2/int($strands{$k}{16}));
			$header1 = ">".$k."/".$no;
			$s1 = revcomplement($s2);
			print file3 "$header1\n$s1\n$header2\n$s2\n";
		}else{
			print "error\n";
		}
	}
}
print "paired: $pairedno\n single: $singleno\n rc: $rcno\n";

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