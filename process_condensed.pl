# Compute the average coverage of the condensed nodes in the graph 
# Remove condensed nodes from the graph which consist of the tips of the graph
# USAGE: perl process_condensed.pl -f condensedgraph -k kmerfile -t thresholdForTips -ct CoverageThreshold -rem
#  RemoveTipsOnly -ck CheckFlag  
use Getopt::Long;
$condgraph_file = "";
$kmerfile = "";
$threshold = 0;
$coverage_threshold = 35;
$check = '';
$remove_tips= '';
$writeoutput ="output.graph";
GetOptions("f=s",\$condgraph_file,
			"k=s",\$kmerfile,
			"t:i",\$threshold,
			"ct:i",\$coverage_threshold,
			"rem",\$remove_tips, # remove all nodes below threshold of coverage if this parameter is mentioned
			"ck",\$check,
			"o:s",\$writeoutput)
			or die("Error in Input Arguements");
die("perl process_condensed.pl -f condensedgraph -k kmerfile -t thresholdForTips -ct CoverageThreshold -rem") if $condgraph_file eq "";
# Variables that will be in use
%cond_ver_average = ();
%cond_ver = ();
%kmerhash = ();
%cond_ver_suffix = ();
%cond_ver_prefix = ();
%vertex = ();
%cond_graphout = ();
%cond_graphin = ();


main();

sub main
{
	loadkmerfile($threshold);
	computeaveragecoverage();
	#printcondgraphaverage();
	loadcondgraph();
	if(!$remove_tips)
	{
		remove_low_coverage_tips($coverage_threshold,1); # Set second parameter to 1 if you only want to remove tips/ or anything else if you want to remove all nodes with coverage below converage_threshold
	}
	else
	{
		remove_low_coverage_tips($coverage_threshold,10); # Set second parameter to 10 if you want to remove all nodes with coverage below converage_threshold
	}
	if($check)
	{
		printremovednodesstats();
	}
	else
	{
		update_cond_graph();
	}
}

sub computeaveragecoverage
{
	open(file,$condgraph_file)  or die("Cannot load $condgraph_file\n");
	while($l=<file>)
	{
		chomp($l);
		($n1,$n2) = split(/ /,$l);
		if($n2 =~ /[AGCT]/)
		{
			$cond_ver{$n1} = $n2;
			$nokmers = length($n2) - $K +1 ; $count_total = 0;
			for(my($i)=0;$i< $nokmers;$i++)
			{
				$count_total += $kmerhash{substr($n2,$i,$K)} ;
			}
			$average = $count_total/$nokmers;
			$cond_ver_average{$n1} = $average;
			# print "Vtx $n1, Seq $n2\n";
		}
	}
	close file;
}

sub loadkmerfile
{
	#only load the kmers that are above in the depth file  
	my($t) = $_[0];
	open(file,$kmerfile) or die("Cannot load $kmerfile\n");
	while($l=<file>)
	{
		chomp($l);
		($kmer,$c) = split(/ /,$l);
		if($c>$t) 
		# if(exists($depth{$kmer}))
		{
			$kmerhash{$kmer} = $c;
			$rev_kmer = revcomplement($kmer); # ADD the reverse complement directly to the kmer hash 6/18/2015
			if (! exists($kmerhash{$rev_kmer}) ) 
			{
				$kmerhash{$rev_kmer} = $c;
			}
		}
	}
	$totalkmers = scalar(keys %kmerhash);
	$K = length($kmer);
	# print "#Total kmers from $kmerfile is $totalkmers and $K is k length\n";
	close file;
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
		else	{$rev_str = $rev_str."N"; } 
		
	}
	return $rev_str;
}

sub printcondgraphaverage
{
	for $k (sort {$cond_ver_average{$a} <=> $cond_ver_average{$b}} keys %cond_ver_average)
	{
		print "$k $cond_ver_average{$k} ".length($cond_ver{$k})."\n";
	}
}

sub loadcondgraph
{
	# Load condensed graph from the text file
	open(file,$condgraph_file)  or die("Cannot load $condgraph_file\n");
	while($l=<file>)
	{
		chomp($l);
		($n1,$n2) = split(/ /,$l);
		if($n2 =~ /[AGCT]/)
		{
			$cond_ver{$n1} = $n2;
			$cond_ver_suffix{substr($n2,-$K)} = $n1;
			$cond_ver_prefix{substr($n2,0,$K)} = $n1;
			for(my($i)=0;$i<length($n2)-$K+1;$i++)
			{
				$vertex{substr($n2,$i,$K)} = $n1;
			}
			# print "Vtx $n1, Seq $n2\n";
		}
		else
		{
			# print "Vtx1 $n1, Vtx2 $n2\n";
			$cond_graphout{$n1}{$n2} = 1;
			$cond_graphin{$n2}{$n1} = 1;
			$source{$n2} = 0;
			$sink{$n1} = 0;
		}
	}
	close file;
	for $vv(keys %cond_ver) 
	{
		if(!exists($source{$vv}))
		{
			push @cond_source_nodes, $vv;
		}
		if(!exists($sink{$vv}))
		{
			push @cond_sink_nodes , $vv;
		}
	}
	undef %source;
	undef %sink;
	$total_ver = scalar(keys %cond_ver);
}

sub remove_low_coverage_tips
{
	# This function will remove nodes that are sink/source nodes
	# and have a small average coverage 
	$coverage_thresh = $_[0];
	$tips_only = $_[1];
	$remove_node = 0;
	$number_source = 0;
	$number_sink = 0;
	%mark_node = ();
	for $k( sort {$cond_ver_average{$a} <=> $cond_ver_average{$b}} keys %cond_ver_average) 
	{
		# Remove a condensed vertex if it is a source/sink and has average coverage less than threshold
		if($tips_only ==1)
		{
			if ( (scalar(%{$cond_graphin{$k}})==0 || scalar(%{$cond_graphout{$k}})==0) && $cond_ver_average{$k} < $coverage_thresh )
			{
				$remove_node++;
				$mark_node{$k} = 1;
				delete $cond_ver{$k};
			}
		}
		else
		{
			if( $cond_ver_average{$k} < $coverage_thresh )
			{
				$remove_node++;
				$mark_node{$k} = 1;
				delete $cond_ver{$k};
			}
		}	
	}
}
sub update_cond_graph
{
	open(wrfile,">$writeoutput") or die ("Cannot write to file $writeoutput\n");

	for $k ( sort {$a <=> $b} keys %cond_ver)
	{
		print wrfile "$k $cond_ver{$k}\n";
	}
	
	for $k ( keys %cond_graphout)
	{
		if(exists($cond_ver{$k}))
		{
			for $y(keys %{$cond_graphout{$k}})
			{
				if(exists($cond_ver{$y}))
				{
					print wrfile "$k $y\n";
				}
			}
		}
	}
	close wrfile;
	print "Written updated condensed graph to $writeoutput\n";
}

sub printremovednodesstats
{
	for $k(@cond_source_nodes)
	{
		if(exists($mark_node{$k}))
		{
			$number_source++;
		}else
		{
			print "Source $k $cond_ver_average{$k}\n";
		}
	}
	for $k(@cond_sink_nodes)
	{
		if(exists($mark_node{$k}))
		{
			$number_sink++;
		}else
		{
			print "Sink $k $cond_ver_average{$k}\n";
		}
	}
	print "The total number of condensed nodes is ".scalar(keys %cond_ver).", compared to $total_ver before \n";
	print "The number of source and sink nodes are ".scalar(@cond_source_nodes)." ".scalar(@cond_sink_nodes)."\n";
	print "The number of nodes that have to be removed is $remove_node, Number sink $number_sink, Number Source $number_source \n";
}
