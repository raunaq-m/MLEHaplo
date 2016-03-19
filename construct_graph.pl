# Construct the graph from kmer file
# perl construct_graph.pl <fasta_file> <kmercountsfile> <threshold_forgraph> <graphwritefile>
#
use lib "/i3c/alisa/rom5161/rh6/perl/share/perl5";
use Graph;

if(scalar(@ARGV)<5)
{
	print "Error in Input Arguments\n";
	print "USAGE:\nperl construct_graph.pl \"fastafile\" \"kmerfile\" \"threshold\" \"graphfile\" \"s\"\n";
	exit(0);
}

$g = Graph->new(directed=>1);
$mode = $ARGV[4];
$kmer_threshold = $ARGV[2];
$write_file = $ARGV[3];

unlink $write_file or print "Could not remove $write_file \n";

#$no_of_vertices=0; $no_of_edges=0; %vertices =();
if ($mode eq "s")
{
	$fasta_file = $ARGV[0];
	$kmer_count_file = $ARGV[1];
	count_kmers_from_kmercounts($kmer_count_file);
	construct_graph_from_read_file($fasta_file);
}
else 
{
	#### DEPRICATED CODE, Convert the paried fasta file into a single file for analysis. 
	$fasta_file = $ARGV[0];
	$kmer_count_file = $ARGV[1];
	$fasta_file2 = $ARGV[5];
	$kmer_count_file2 = $ARGV[6];
	count_kmers_from_kmercounts($kmer_count_file);
	count_kmers_from_kmercounts($kmer_count_file2);
	print "Creating graph from firt paired end file \n";
	construct_graph_from_read_file($fasta_file);
	print "Creating graph from second paired end file \n";
	construct_graph_from_read_file($fasta_file2);
}

$kmerlength = 0;
%kmer_counts = ();

#compute the kmer counts 


#construct graph from read file 


print "Number of vertices ".$g->vertices." number of edges ".$g->edges."\n";
#$no_of_vertices = scalar(keys %vertices);
#print "Number of vertices $no_of_vertices number of edges $no_of_edges \n";
sub count_kmers_from_kmercounts
{
	print "Counting kmer counts from the kmer file\n";
	my($file) = $_[0];
	open(file_count,$file);
	while($l=<file_count>)
	{
		chomp($l);
		($kmer,$c) = split(/ /,$l);
		$kmer_counts{$kmer} = $c;
	}
	$kmerlength = length($kmer);
	close file_count;
	print "Kmer counts loaded\n";
}

sub construct_graph_from_read_file
{
	my($f) = $_[0];
	print "Creating graph \n";
	open(file,$f) or print "file not open $f \n";
	open(wr_file,">>$write_file");
	$l=<file>;
	$read="";
	$read_count =0;
	while($l=<file>) 
	{	#Open fasta file and construct the read graph
		chomp($l);
		if(substr($l,0,1) ne ">") 
		{
			$read=$read.$l;
		}
		else 
		{
	#		print "Adding $read \n";
			addkmers_from_read($read);
			$read_count++;
			if(($read_count %1000)==0)
			{
				print "processed $read_count reads in ".time."\n";
			}
			$read = "";
		}

	}
	addkmers_from_read($read);
	print "Graph done \n";
	close wr_file;
	close file;
}

sub addkmers_from_read
{
	my($r) = $_[0];
	$r = uc($r);
	my($rev_read) = revcomplement($r);
	my($len) = length($r);
	for(my($i)=0; $i<$len-$kmerlength; $i++)
	{
		$u = substr($r,$i,$kmerlength);
		$v = substr($r,$i+1,$kmerlength);
		if(!($u =~ /N/ || $v =~/N/))
		{
			$u_r = substr($rev_read,$len-$kmerlength-$i,$kmerlength);
			$v_r = substr($rev_read,$len-$kmerlength-$i-1,$kmerlength);
			if( ($kmer_counts{$u} > $kmer_threshold |$kmer_counts{$u_r} > $kmer_threshold) && ($kmer_counts{$v} > $kmer_threshold |$kmer_counts{$v_r} > $kmer_threshold) )
			{
				if(!$g->has_edge($u,$v)) { print wr_file "$u $v\n"; }
				if(!$g->has_edge($v_r,$u_r)) { print wr_file "$v_r $u_r\n"; } # reverse complement kmer 
				$g->add_edge($u,$v);
				$g->add_edge($v_r,$u_r);
			}
		}
	}
}
#sub addkmers_from_read
#{	#the lexicographically lower k-mer's count is used for computation of the graph
#	my($r) = $_[0]; # get 
#	for(my($i)=0;$i<length($r)-$kmerlength; $i++) 
#	{
#		$u = substr($r,$i,$kmerlength);
#		$u_= $u;
#		$u_r = revcomplement($u);
#		if ( ($u cmp $u_r)>0 ) { $u_ = $u_r; }
#		$v = substr($r,$i+1,$kmerlength);
#		$v_= $v;
#		$v_r = revcomplement($v);
#		if ( ($v cmp $v_r)>0 ) { $v_ = $v_r; }
#		#print $kmer_counts{convert_to_number($u_)}." ".$kmer_counts{convert_to_number($v_)}."\n";
#		if ( $kmer_counts{$u_} > $kmer_threshold && $kmer_counts{$v_} > $kmer_threshold ) 
#		{
#			if(!$g->has_edge($u,$v))
#			{ print wr_file "$u $v\n";  } #$no_of_edges++;  $vertices{$u}=1; $vertices{$v}=1;}
#			if(!$g->has_edge($v_r,$u_r)) 
#			{ print wr_file "$v_r $u_r\n"; } #$no_of_edges++; $vertices{$v_r}=1; $vertices{$u_r}=1;}
#			$g->add_edge($u,$v);
#			$g->add_edge($v_r,$u_r);
#			
#		}
#	}
#}
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
