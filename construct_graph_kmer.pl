# Construct graph from kmer+1 file directly 
# perl construct_graph_kmer.pl -k1 <kmer+1countsfile> -k <kmerfile> -t <threshold> -w <graphwritefile>

use Getopt::Long;
use Graph;

$kmer_1file = "";
$kmerfile = "";
$graphwritefile = "";
$threshold = 0;
GetOptions("k1=s",\$kmer_1file,"k=s",\$kmerfile,"t:i",\$threshold,"w=s",\$graphwritefile) or die("USAGE:perl construct_graph_kmer.pl -k1 <kmer+1countsfile> -k <kmerfile> -t <threshold> -w <graphwritefile>\n");

$kmerlength = 0;
%kmer_counts = ();
$g = Graph->new(directed=>1);
main();

sub main
{
	# Create graph data structure
	# Remove graph file if it exists previously (For now)
	unlink $graphwritefile or print "Could not remove $graphwritefile\n";
	# load k-mers from the kmer coutns file 
	load_kmers_from_kmercounts();
	construct_graph_from_k_1_mer_counts();
}

sub load_kmers_from_kmercounts
{
	open(file,$kmerfile) or die("Error opening $kmerfile kmer counts file\n");
	print "Loading kmers from the kmer counts file $kmerfile\n";
	while($l=<file>)
	{
		chomp($l);
		($kmer,$c) = split(/ /,$l);
		$kmer_counts{$kmer} = $c;
	}
	$kmerlength = length($kmer);
	close file;
	print "Kmer counts loaded \n";
}

sub construct_graph_from_k_1_mer_counts
{
	open(file,$kmer_1file) or die ("Could not load k+1-mer file $kmer_1file\n");
	print "Creating graph \n";
	open(wrfile,">$graphwritefile");
	while($l=<file>)
	{
		chomp($l);
		($kmer,$c) = split(/ /,$l);
		$u = substr($kmer,0,$kmerlength);
		$v = substr($kmer,1);
		#Only proceed if (k+1)-mer count is greater than the threshold
		$rev_kmer = revcomplement($kmer);
		$u_r = substr($rev_kmer,1);
		$v_r = substr($rev_kmer,0,$kmerlength); 
		if(($kmer_counts{$u}>$threshold | $kmer_counts{$u_r} > $threshold) && ($kmer_counts{$v_r} > $threhsold | $kmer_counts{$v} > $threshold))
		{
			if(!$g->has_edge($u,$v)) 
			{
				print wrfile "$u $v\n";
			}
			if(!$g->has_edge($v_r,$u_r))
			{
				print wrfile "$v_r $u_r\n";
			}
			$g->add_edge($u,$v);
			$g->add_edge($v_r,$u_r);
		}
	}
	close file;
	close wrfile;
	print "Graph written to $graphwritefile\n";
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
