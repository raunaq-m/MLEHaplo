########## Algorithm for cover of all nodes using the condensed graph
# usage : perl dg_cover_condensed.pl 
use Bio::SeqIO;
use List::Util qw(shuffle max min);
use POSIX qw(ceil);
use Getopt::Long;
#############VARIABLES########################
$condgraph_file = "";
$kmerfile = "";
$nodedepthfile = "";
$thresh = 0;
$paired_kmers_file = "";
$factor = 0;
$numpaths = '';
$depth_of_end = 0;
$score = -10;
$IS_factor = 500;

###############GET OPTIONS####################
GetOptions( "cg=s",\$condgraph_file, 
			"nd=s",\$nodedepthfile,
			"k=s",\$kmerfile,
			"t:i",\$thresh,
			"paired:s",\$paired_kmers_file,
			"fact:f",\$factor,
			"IS:i",\$IS_factor)
			or die("Error in input arguments\n");
###############################################

############### HASHES AND PARAMETERS #########
%kmerhash = ();
%depth = ();
%vertex = ();
%start = ();
%finish = ();
%pred = ();
%cond_graphin = ();
%cond_graphout = ();
%paired_kmers = ();
@cond_source_nodes = ();
@cond_sink_nodes = ();
$K = 0;

################################################
main();

sub main
{
	loadkmerfile($thresh);
	loadcondgraph();
	$nokmers = scalar(keys %vertex); #nokmers is the number of vertices in the graph
	print "Parameters : graph $condgraph_file kmer $kmerfile thresh $thresh paired $paired_kmers_file factor $factor IS $IS_factor depth $nodedepthfile\n";
	print "Graph loaded with no of vertices $nokmers \n";
	print "K value of the kmer is $K \n";
	print "Number of source nodes is ".scalar(@cond_source_nodes)."\n";
	print "Number of sink nodes is ".scalar(@cond_sink_nodes)."\n";
	loaddepthfile();
	loadpaired();
	$num_comp_pairs = load_paired_feature_select_updated() if $paired_kmers_file ne "";
	$num_comp_pairs = generate_compatible_set() if $paired_kmers_file ne "";
	print "Total paired k-mers selected is ".$num_comp_pairs."\n";
	condgraph_number_of_paths_wrapper();
	_path_all_nodes_backward();
	path_all_nodes(); #print_constraints();
	_redundant_constraints();
	print_constraints();
	print "$factor is factor\n";
}

sub loaddepthfile
{
	open(file,$nodedepthfile) or die("Cannot load $nodedepthfile\n");
	while($l=<file>)
	{
		chomp($l);
		($kmer,$depth_temp,$root) = split(/ /,$l);
		if(exists($vertex{$kmer}))
		{
			$depth{$kmer} = $depth_temp;
			$root{$kmer} = $root;
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
}

sub loadpaired 
{
	open(filepaired,$paired_kmers_file)  or die("Cannot load $paired_kmers_file\n");
	while($l=<filepaired>) 
	{
		chomp($l);
		my ($lk1,$rk2,$cc) = split(/ /,$l);
		#if( exists($kmerhash{$lk1}) || exists($kmerhash{revcomplement($lk1)}) ) && ( exists($kmerhash{$rk2}) || exists($kmerhash{revcomplement($rk2)}) )
		# Only store k-mer pairs if they are present in the vertex set 
		if(exists($vertex{$lk1}) && exists($vertex{$rk2})) 
		{
			$paired_kmers{$lk1}{$rk2} = $cc;
			if ($depth{$rk2} >= $depth{$lk1})
			{
				$paired_depth{$lk1}{$rk2} = $depth{$rk2} - $depth{$lk1} ;
			}
		}
	}
	close filepaired;
	$avg = 0; $count=0;
	for $k1(keys %paired_depth) {
		for $k2(keys %{$paired_depth{$k1}}) { $avg = $avg + $paired_depth{$k1}{$k2} ; $count++;}
	}
	$size_insert = $avg/$count;
	print "#Avg insert size is $avg/$count =  ".$avg/$count."\n";
	undef %paired_depth;
}

sub load_paired_feature_select_updated
{
	my(%nodes_bubbles) = ();
	my(%add_to_bubbles) = ();
	$depth{"source"} = 50000; $depth{"sink"} = 50000;
	for $k1( sort {$depth{$a}<=>$depth{$b}} keys %cond_ver_suffix )
	{
		# if ( scalar(keys %{$cond_graphin{ $cond_ver_suffix{$k1} } } ) == 1 && scalar(keys %{$cond_graphout{ $cond_ver_suffix{$k1} } } ) == 1 ) 
		# { 
			# $nodes_bubbles{$vertex{$k1}} = 1 ;
		# }
		push @{$add_to_bubbles{$depth{substr($cond_ver{$cond_ver_suffix{$k1}},0,$K)}}}, $vertex{$k1};
	}
	for $k(keys %add_to_bubbles)
	{
		if(scalar(@{$add_to_bubbles{$k}}) >1)
		{
			foreach (@{$add_to_bubbles{$k}}) { $nodes_bubbles{$_} =1; }
		}
	}
#	for (keys %nodes_bubbles) { print "Bubble: $_ \n"; }
	# SECOND: Add sink and source nodes also to the bubble set, does this mean all the edges are accounted in the above two steps?? 
	foreach (@cond_source_nodes) { $nodes_bubbles{$_} = 1; }
	foreach (@cond_sink_nodes)   { $nodes_bubbles{$_} = 1; }
	print "Number of vertices ".scalar(keys %cond_ver).", number of seleccted bubbles ".scalar(keys %nodes_bubbles)."\n";
	# THIRD: Find pairs in the paired-kmers set, where both of the k-mers are present in one of the bubbles
#	%selective_pairs = ();
	%bubble_sets = ();
	my ($sel_no) = 0;
	for $n1( keys %paired_kmers ) 
	{	
		for $n2( keys %{$paired_kmers{$n1}} ) 
		{
			if (exists($nodes_bubbles{$vertex{$n1}}) && exists($nodes_bubbles{$vertex{$n2}} ) ) 
			{
		#		$selective_pairs{$n1}{$n2} = $paired_kmers{$n1}{$n2}; # support for two k-mers occurring together
				$sel_no ++;
				$bubble_sets{$vertex{$n1}}{$vertex{$n2}} += $paired_kmers{$n1}{$n2}; 
			}
		}
	}
	# FOURTH : Add neighbors to the compatible sets as default 
	for $n1 ( keys %cond_graphout) 
	{
		for $n2 ( keys %{$cond_graphout{$n1}} )
		{
			if(!exists($bubble_sets{$n1}{$n2}) ) {
				$bubble_sets{$n1}{$n2} = 1 ; 
			}
		}
	}
	undef %nodes_bubbles;
	my($comp_no) = 0; 
	foreach $n1 (keys %bubble_sets) {
		$comp_no += scalar(keys %{$bubble_sets{$n1}});
	}
	print_compatible_sets("bubble.txt");
	print "Pairs used are $sel_no ";
	return $comp_no;	
}

sub generate_compatible_set
{
	#Generate the compatible set file, add all pairs that were observed in reads 
	%compatible_set = ();
	my ($sel_no) = 0;
	for $n1( keys %paired_kmers ) 
	{	
		for $n2( keys %{$paired_kmers{$n1}} ) 
		{
			if($depth{$n2} > $depth{$n1})
			{
				
				$compatible_set{$vertex{$n1}}{$vertex{$n2}} += $paired_kmers{$n1}{$n2}; 
			}
			else
			{
				
				$compatible_set{$vertex{$n2}}{$vertex{$n1}} +=$paired_kmers{$n1}{$n2};
			}
			$sel_no++;
		}
	}
	my($comp_no) = 0; 
	foreach $n1 (keys %compatible_set) {
		$comp_no += scalar(keys %{$compatible_set{$n1}});
	}
	print_compatible_sets("comp.txt");
	print "Pairs used are $sel_no ";
	return $comp_no;
}

sub print_compatible_sets 
{
	my($suffix) = shift;
 	my($fname) = $condgraph_file;
	$fname =~ s/graph/$suffix/;
	open(fname,">$fname");
	for $k1 ( sort { $depth{$a}<=>$depth{$b} } keys %cond_ver_suffix) 
	{
#		print "Vertex $cond_ver_suffix{$k1} $cond_ver{$cond_ver_suffix{$k1}} $paths_vertex{$cond_ver_suffix{$k1}} Compatible set :";
		print fname "$cond_ver_suffix{$k1}:".$depth{substr($cond_ver{$vertex{$k1}},0,$K)}." ";
		my($last) = $k1;
		if ($suffix =~ /comp/)
		{
			foreach $k2(sort {$depth{substr($cond_ver{$a},0,$K)}<=>$depth{substr($cond_ver{$b},0,$K)}} keys %{$compatible_set{$cond_ver_suffix{$k1}}}) 
			{ 
				print fname "$k2:".$depth{substr($cond_ver{$k2},0,$K)}.":".$compatible_set{$cond_ver_suffix{$k1}}{$k2}." "; 
				$last = $k2
			}
		}
		else 
		{
			foreach $k2(sort {$depth{substr($cond_ver{$a},0,$K)}<=>$depth{substr($cond_ver{$b},0,$K)}} keys %{$bubble_sets{$cond_ver_suffix{$k1}}}) 
			{ 
				print fname "$k2:".$depth{substr($cond_ver{$k2},0,$K)}.":".$bubble_sets{$cond_ver_suffix{$k1}}{$k2}." "; 
				$last = $k2
			}
		}
		$max_depth{$cond_ver_suffix{$k1}} = $depth{substr($cond_ver{$last},0,$K)} - $depth{substr($cond_ver{$vertex{$k1}},0,$K)};
		#print fname "$max_depth{$cond_ver_suffix{$k1}} ";
#		print "Edge set :";
#		foreach (keys %{$cond_graphout{$cond_ver_suffix{$k1}}}) { print "$_ "; }
		print fname "\n";
	}
	close fname;
}

sub condgraph_number_of_paths_wrapper
{
	print "Number of condensed vertices is                                               ".scalar(keys %cond_ver)."\n";
	for $source( @cond_source_nodes) {
		$cond_graphout{"source"}{$source} = 1;
		$cond_graphin{$source}{"source"} = 1;
	}
	for $k(@cond_sink_nodes) {
		$cond_graphout{$k}{"sink"} = 1;
		$cond_graphin{"sink"}{$k} = 1;
	}
	$cond_ver{"source"} = $vno+1;
	$cond_ver{"sink"} = $vno+2;
	undef %paths_vertex;
	for $k ( keys %cond_ver) { $paths_vertex{$k} = -1; }
	$number_of_paths = find_number_paths("source","sink",\%cond_graphout);
	print "Total number of paths 							$number_of_paths\n";	
}

sub find_number_paths
{
	my($first) = shift;
	my($second) = shift;
	my($arhash) = shift;
#	my($hashchoice) = shift;
	if($first eq $second ) 
	{
		$paths_vertex{$first} = 1;
		return 1;
	}else
	{
		if($paths_vertex{$first} == -1) {
			my(@arr) = keys %{$arhash->{$first}};
			if(scalar(@arr) !=0) {
				$paths_vertex{$first} = 0;
				foreach (@arr) {
				#	if($hashchoice==2) {
					 $paths_vertex{$first} += find_number_paths($_,$second,\%cond_graphout);  
				#	else { $paths_vertex{$first} += find_number_paths($_,$second,\%graphout,1); }
				}
				return $paths_vertex{$first};
			}else {
				return 0;
			}
		}else { return $paths_vertex{$first}; }
	}
}

sub _path_all_nodes_backward
{
	_compute_membership_paired_nodes(); # Probability based on just node membership
	initialize_two_nodes();
	#Marked is to make sure that sink nodes are processed only once in the dataset
	my(%marked) = ();
	for (keys %cond_ver)
	{
		$marked{$_} = 0;
	}
	$marked{"source"} = 1; 
	$marked{"sink"} = 1;
	for $sink( @cond_sink_nodes )
	{
		path_two_nodes_dg_backward_cover($sink);
		print "sink $sink $depth{substr($cond_ver{$sink},0,$K)}\n";
		foreach (@{$node_paths_back{$sink}}) { print $_."\t";}
		print "\n";		
		foreach (@{$node_paths_like_back{$sink}}) { print $_."\t";}
		print "\n";
		$marked{$sink} = 1;
		#last;
	}
	for $ss ( @cond_source_nodes )
	{
		$marked{$ss} = 0; #Visit nodes that are both source and sink twice
	}
	for $node ( sort {$depth{substr($cond_ver{$b},0,$K)} <=> $depth{substr($cond_ver{$a},0,$K)}} keys %cond_ver ) 
	{
		if($marked{$node} ==0) 
		{
			path_two_nodes_dg_backward_cover($node);
			print $node."back depth".$depth{substr($cond_ver{$node},0,$K)}."\nNeighbors";
			foreach (keys %{$cond_graphout{$node}}) { print "$_ ";} print "\n";
			foreach (@{$node_paths_back{$node}}) { print $_."\t";}
			print "\n";		
			foreach (@{$node_paths_like_back{$node}}) { print $_."\t";}
			print "\n";
			$marked{$node} = 1;
		}
	}
}


sub _compute_membership_paired_nodes
{
	for $k(keys %compatible_set)
	{
		# print $k." compatible\t";
		# my($total) = scalar(keys %{$compatible_set{$k}});
		my($total) = 0;
		foreach (keys %{$compatible_set{$k}})
		{
			$total +=$compatible_set{$k}{$_};
		}
		foreach (keys %{$compatible_set{$k}})
		{
			$compatible_set_mem{$k}{$_} = $compatible_set{$k}{$_}/$total; # Modified 4/2/2015
			# print $compatible_set_mem{$k}{$_}."\t";
		}
		# print "\n";
	}
}


sub initialize_two_nodes
{
	#Initializes each node to have an empty path set associated with it, along with 0 probability to the path
	for $k(keys %cond_ver)
	{
		# if($k eq "source") 
		# {
			# print $k."\n";
		# }	
		push (@{$node_paths{$k}} , (""));
		push (@{$node_paths_like{$k}} , (0));
		push (@{$node_paths_back{$k}} , (""));
		push (@{$node_paths_like_back{$k}} , (0));
		#foreach (@{$node_paths_like{$k}}) { print $_." "; }
		#print "\n";
	}
	for $k( keys %compatible_set)
	{
		for $y( keys %{$compatible_set{$k}})
		{
			$assigned{$k}{$y} = -1;
		}
	}
}

sub path_two_nodes_dg_backward_cover
{
	#initialization assumed that every node has at most $factor top paths saved 
	my($start) = shift; 
	my($epsilon) = 0.000001;
	my(@opt_sub_paths) = ();
	my(@prob_sub_paths) = ();
	#Create a list of possible paths from previous node to current node
	# print "dg_ $start ";
	for $k( keys %{$cond_graphout{$start}} ) 
	{
		# print "$k ";
		push (@opt_sub_paths, @{$node_paths_back{$k}});
		push (@prob_sub_paths, @{$node_paths_like_back{$k}});
	}
	
	for(my($u)=0;$u<=$#opt_sub_paths; $u++)
	{
		@nodes = split(/_/,$opt_sub_paths[$u]);
		my($temp_like) = 0;
		my($path_length) = 0;
		my($cc) = scalar(@nodes);
		if($cc==1) 
		{
			$temp_like = $prob_sub_paths[$u];
		}
		else
		{
			$temp_like = $prob_sub_paths[$u]*($cc-1)*($cc)/2;
		}
		print $temp_like." is the like \n";		
		for(my($iter)=0;$iter<$cc; $iter++)
		{
			if($depth{substr($cond_ver{$nodes[$iter]},0,$K)} - $depth{substr($cond_ver{$start},0,$K)} < $IS_factor)
			{
				$path_length += length($cond_ver{$nodes[$iter]})-$K;
				if(exists($compatible_set_mem{$start}{$nodes[$iter]}))
				{
					if(exists($bubble_sets{$start}{$nodes[$iter]}))
					{
						# bubble sets indicate two nodes that are present in a bubble
						# $temp_like += 2;
						$temp_like +=2*$compatible_set_mem{$start}{$nodes[$iter]};
					}
					else
					{
						# $temp_like +=1;
						$temp_like += 1*$compatible_set_mem{$start}{$nodes[$iter]};
					}
				}
				else
				{
					$temp_like -= ($cc);
				}
			}
			elsif(exists($compatible_set_mem{$start}{$nodes[$iter]}))
			{
				if(exists($bubble_sets{$start}{$nodes[$iter]}))
				{
					# bubble sets indicate two nodes that are present in a bubble
					# $temp_like += 1;
					$temp_like += 1*$compatible_set_mem{$start}{$nodes[$iter]};
				}
				else
				{
					# $temp_like +=0.5;
					$temp_like +=0.5*$compatible_set_mem{$start}{$nodes[$iter]};
				}
			}else
			{
				$temp_like += 0.01;
				# $temp_like += 0.0001;
			}
			# if(exists($compatible_set_mem{$start}{$nodes[$iter]}))
			# {
				# # print "$start $nodes[$iter] yes\n";
				# $temp_like += 1;
			# }
			# elsif($depth{substr($cond_ver{$nodes[$iter]},0,$K)} - $depth{substr($cond_ver{$start},0,$K)} > $max_depth{$start})
			# {
				# # print "$start $nodes[$iter] \n";
				# $temp_like += 1;
			# }
			# else
			# {
				# # print "penalty $start $nodes[$iter] \n";
				# $temp_like -= ($cc); #Penalize the path as it contains un-supported pairs. 
			# }
		}
		if($cc>1)
		{
			$prob_sub_paths[$u]  = 2*$temp_like/($cc*($cc+1));
		}else
		{
			$prob_sub_paths[$u] = $temp_like;
		}
		# Adjust $prob_sub_paths to 1 if > 1
		# if($prob_sub_paths[$u] > 1) 
		# {	
			# $prob_sub_paths[$u] = 1;
		# }
		$path_length += $K;
		print "Like after evaluation $prob_sub_paths[$u] $temp_like $path_length \n";
	}
	my($last_max) = -9**9**9; my($pos) = -100;
	#Find top 5 values of probabilities now 
	@{$node_paths_back{$start}} = ();
	@{$node_paths_like_back{$start}} = ();
	@unsorted = (0..$#prob_sub_paths);
	my @sorted_paths_indices = sort {$prob_sub_paths[$b] <=> $prob_sub_paths[$a]} 0..$#prob_sub_paths;
	my @sorted_paths_like = @prob_sub_paths[ @sorted_paths_indices];
	$no_elements  = min($factor, scalar(@prob_sub_paths));
	$total_constraints = (keys %{$compatible_set{$start}}) -1;
	$total_count = 0;
	my(%constraints_per_path) = ();
	my($i)=0;my($iterate)=0;
	#for (my($i) = 0; $i< $no_elements; $i ++)
	# Adding cutoff for path scores for each node
	$initial_val = $sorted_paths_like[0];
	$cutoff = $initial_val*0.25;
	while($i<$no_elements && $iterate <=$#sorted_paths_indices)
	{
		# print "$i Next val $sorted_paths_like[$i] $sorted_paths_indices[$i]\n";
		$pos = $sorted_paths_indices[$iterate];
		$val = $sorted_paths_like[$iterate];
		if ($val >= $cutoff)
		{
			@added_nodes = split(/_/, $opt_sub_paths[$pos]);
			$path_count = 0;
			foreach (@added_nodes) 
			{
				if($assigned{$start}{$_} == -1)
				{
					$assigned{$start}{$_} = 1;
					$path_count ++;
				}
			}
			$constraints_per_path{$pos} = $path_count;
			$total_count += $path_count;
			#print "count $path_count\t$total_count/$total_constraints\n";
			if( (scalar(@added_nodes) < 2 && $path_count ==0) || ($path_count != 0) || ($val > 0) )
			{
				push @{$node_paths_back{$start}} , $start."_".$opt_sub_paths[$pos];
				if($initial_val !=0)
				{
					$val = $val/$initial_val;
				}
				push @{$node_paths_like_back{$start}}, $val;
				$i++;	
			}
			# if($path_count != 0)
			# {
				# push @{$node_paths_back{$start}} , $start."_".$opt_sub_paths[$pos];
				# push @{$node_paths_like_back{$start}}, $val;
				# $i++;
			# }elsif($val > 0)
			# {
				# push @{$node_paths_back{$start}} , $start."_".$opt_sub_paths[$pos];
				# push @{$node_paths_like_back{$start}}, $val;
				# $i++;
			# }
		}
		$iterate++; 
	}
	if (scalar( @{$node_paths_back{$start}} ) == 0  ) 
	{
		push @{$node_paths_back{$start}}, $start."_";
		push @{$node_paths_like_back{$start}}, 0; 
	}
	if ($total_count < $total_constraints)
	{
		print "Add more paths for $start node\n";
	}
}

sub path_all_nodes
{
	#_compute_probability_paired_nodes(); # Probability based on the paired reads supported
	_compute_membership_paired_nodes(); # Probability based on just node membership
	initialize_two_nodes();
	#initialize from the source nodes
	%marked = (); # To make sure source nodes are not processed twice
	for (keys %cond_ver)
	{
		$marked{$_} = 0;
	}
	$marked{"source"} = 1; 
	$marked{"sink"} = 1;
	for $source( @cond_source_nodes )
	{
		path_two_nodes_dg_forward($source);
		# print "source $source $depth{substr($cond_ver{$node},0,$K)}\n";
		# foreach (@{$node_paths{$source}}) { print $_."\t";}
		# print "\n";		
		# foreach (@{$node_paths_like{$source}}) { print $_."\t";}
		# print "\n";
		$marked{$source} = 1;
		#last;
	}
	for $node ( sort {$depth{substr($cond_ver{$a},0,$K)} <=> $depth{substr($cond_ver{$b},0,$K)}} keys %cond_ver ) 
	{
		if($marked{$node} ==0) 
		{
			path_two_nodes_dg_forward($node);
			# print $node."node depth".$depth{substr($cond_ver{$node},0,$K)}."\nNeighbors";
			# foreach (keys %{$cond_graphin{$node}}) { print "$_ ";} print "\n";
			# foreach (@{$node_paths{$node}}) { print $_."\t";}
			# print "\n";		
			# foreach (@{$node_paths_like{$node}}) { print $_."\t";}
			# print "\n";
			$marked{$node} = 1;
		}
	}
	# Do the same analysis on backward direction
	#_path_all_nodes_backward();
	#_redundant_constraints();
}

sub path_two_nodes_dg_forward
{
	#initialization assumed that every node has at most 5 top paths saved 
	my($start) = shift; 
	my($epsilon) = 0.000001;
	my(@opt_sub_paths) = ();
	my(@prob_sub_paths) = ();
	#Create a list of possible paths from previous node to current node
	for $k( keys %{$cond_graphin{$start}} ) 
	{
		# print "$k ";
		push (@opt_sub_paths, @{$node_paths{$k}});
		push (@prob_sub_paths, @{$node_paths_like{$k}});
	}
	# print "Number of paths is ".scalar(@opt_sub_paths)."\n";
	for(my($u)=0;$u<=$#opt_sub_paths; $u++)
	{
		@nodes = split(/_/,$opt_sub_paths[$u]);
		my($temp_like) = 0;
		my($path_length) = 0;
		my($cc) = scalar(@nodes);
		if($cc==1) 
		{
			$temp_like = $prob_sub_paths[$u];
		}
		else
		{
			$temp_like = $prob_sub_paths[$u]*($cc-1)*($cc)/2;
		}
		# print $temp_like." is the like \n";
		
		for(my($iter)=1;$iter<$cc; $iter++)
		{
			if($depth{substr($cond_ver{$start},0,$K)} - $depth{substr($cond_ver{$nodes[$iter]},0,$K)} < $IS_factor)
			{
				$path_length += length($cond_ver{$nodes[$iter]}) -$K;
				if(exists($compatible_set_mem{$nodes[$iter]}{$start}))
				{
					if(exists($bubble_sets{$nodes[$iter]}{$start}))
					{
						# $temp_like +=2;
						$temp_like += 2*$compatible_set_mem{$nodes[$iter]}{$start};
					}
					else
					{
						# $temp_like +=1;
						$temp_like += 1*$compatible_set_mem{$nodes[$iter]}{$start};
					}
				}
				else
				{
					$temp_like -=$cc;
				}
			}
			elsif(exists($compatible_set_mem{$nodes[$iter]}{$start}))
			{
				if(exists($bubble_sets{$nodes[$iter]}{$start}))
				{
					 # $temp_like +=1;
					 $temp_like += 1*$compatible_set_mem{$nodes[$iter]}{$start};
				}
				else
				{
					# $temp_like +=0.5;
					$temp_like += 0.5*$compatible_set_mem{$nodes[$iter]}{$start};
				}
			}
			else
			{
				$temp_like += 0.01;
				# $temp_like += 0.0001;
			}
			# if(exists($compatible_set_mem{$nodes[$iter]}{$start}))
			# {
				# #$temp_like += log($compatible_set_mem{$nodes[$iter]}{$start});
				# $temp_like +=1;
			# }
			# elsif($depth{substr($cond_ver{$start},0,$K)} - $depth{substr($cond_ver{$nodes[$iter]},0,$K)} > $max_depth{$nodes[$iter]})
			# {
				# # print "$start $nodes[$iter] \n";
				# #Penalty less as two nodes outside the range of size_insert
				# $temp_like += 1; 
			# }
			# else
			# {
				# # print "penalty $start $nodes[$iter] \n";
				# #Within size range but support not present, then penalize 
				# $temp_like -= ($cc-1);
			# }
		}
		if($cc>1)
		{
			$prob_sub_paths[$u]  = 2*$temp_like/($cc*($cc+1));
		}else
		{
			$prob_sub_paths[$u] = $temp_like;
		}
		if($prob_sub_paths[$u] > 1)
		{
			$prob_sub_paths[$u] = 1;
		}
		# print "Like after evaluation $prob_sub_paths[$u] \n";
	}
	my($last_max) = -9**9**9; my($pos) = -100;
	#Find top 5 values of probabilities now 
	@{$node_paths{$start}} = ();
	@{$node_paths_like{$start}} = ();
	@unsorted = (0..$#prob_sub_paths);
	my @sorted_paths_indices = sort {$prob_sub_paths[$b] <=> $prob_sub_paths[$a]} 0..$#prob_sub_paths;
	my @sorted_paths_like = @prob_sub_paths[ @sorted_paths_indices];
	$no_elements  = min($factor, scalar(@prob_sub_paths));
	$initial_val = $sorted_paths_like[0];
	$cutoff = $initial_val*0.25;
	for (my($i) = 0; $i< $no_elements; $i ++)
	{
		# print "$i Next val $sorted_paths_like[$i] $sorted_paths_indices[$i]\n";
		$pos = $sorted_paths_indices[$i];
		$val = $sorted_paths_like[$i];
		if($val >=$cutoff)
		{
			if($initial_val !=0)
			{
				$val = $val/$initial_val;
			}
			push @{$node_paths{$start}} , $opt_sub_paths[$pos]."_".$start;
			push @{$node_paths_like{$start}}, $val;
		}
		# $last_max = $val; $last_val = $pos ; 
	}
}

sub _redundant_constraints
{
	$all_constraints = "";
	$all_constraints_like = "";
	$back_constraints = "";
		for $node( sort {$depth{substr($cond_ver{$a},0,$K)} <=> $depth{substr($cond_ver{$b},0,$K)}} @cond_source_nodes )
	{
		# print "constraints $node ";
		for (my($i) = 0; $i < scalar(@{$node_paths_back{$node}}); $i++)
		{
			${$node_paths_back{$node}}[$i] =~ s/^_//;
			${$node_paths_back{$node}}[$i] =~ s/_$//;
			# print ${$node_paths_back{$node}}[$i]."\t";
			${$node_paths_back{$node}}[$i] = "#".${$node_paths_back{$node}}[$i]."#";
			if($all_constraints !~ /${$node_paths_back{$node}}[$i]/ )
			{
				$all_constraints = $all_constraints.substr(${$node_paths_back{$node}}[$i],0,-1);
				$all_constraints_like = $all_constraints_like."#".${$node_paths_like_back{$node}}[$i];
			}
		}
		# print "\n";
		# print "CC ".$all_constraints."\n";
	}
	for $node( sort {$depth{substr($cond_ver{$b},0,$K)} <=> $depth{substr($cond_ver{$a},0,$K)}} keys %cond_ver)
	{
		$str1 = "_".$node."_";
		$str2 = "#".$node."#";
		$str3 = "_".$node."#";
		$str4 = "#".$node."_";
		if(!($all_constraints =~ /$str1/ || $all_constraints =~ /$str2/ || $all_constraints =~ /$str3/ ||$all_constraints =~ /$str4/))
		{
			print "$node is absent\n";
			my($num_to_pick1) = min(10,scalar(@{$node_paths{$node}}));
			my($num_to_pick2) = min(10,scalar(@{$node_paths_back{$node}}));
			# my($num_to_pick1) = scalar(@{$node_paths{$node}});
			# my($num_to_pick2) = scalar(@{$node_paths_back{$node}});			
			for(my($i)=0; $i< $num_to_pick1; $i++)
			{
				$path1 = ${$node_paths{$node}}[$i];
				for(my($j)=0; $j < $num_to_pick2 ; $j++)
				{
					$path2_temp = ${$node_paths_back{$node}}[$j];
					$path2 = substr($path2_temp,index($path2_temp,"_"));
					if($path1 eq $path2) 
					{
						$new_path = $path1;
					}
					else
					{
						$new_path = $path1.$path2;
					}
					$new_path =~ s/^_//;
					$new_path =~ s/_$//;
					$vals = split(/_/,$new_path);
					if($all_constraints !~ /$new_path/  ||$vals ==1) 
					{
						$all_constraints = $all_constraints."#".$new_path;
						$path_like = ${$node_paths_like{$node}}[$i] + ${$node_paths_like_back{$node}}[$j];
						$all_constraints_like = $all_constraints_like."#".$path_like;
					}
				}
			}
		}
	}
	$count = split(/#/,$all_constraints);
	# $count_back = split(/#/,$back_constraints);
	print $all_constraints."\n".$all_constraints_like."\n";
	# print $back_constraints."\n".$back_constraints_like."\n";
	print "Total count of constratins is $count and $count_back\n";
}

sub print_constraints
{
	@all_paths = split(/#/,$all_constraints);
	@paths_like = split(/#/,$all_constraints_like);
	$iter = 1;
	for(my($c)=0;$c<scalar(@all_paths);$c++)
	{
		if($paths_like[$c] >= $score )
		{
			@nodes = split(/_/,$all_paths[$c]); #Break the paths
			$count = 1;
			if(scalar(@nodes)>=1) 
			{
				for(my($i)=0;$i<scalar(@nodes); $i++)
				{
					print "$count:$nodes[$i] $depth{substr($cond_ver{$nodes[$i]},-$K)} ";
					$count++;
				}
				print "\n";
				$length = $depth{substr($cond_ver{$nodes[$#nodes]},-$K)} + $K;
				#print ">".$iter."_".$nodes[$#nodes]."_".$length."\n";
				#unshift @nodes, "start";
				# if($nodes[0] eq "")
				{
					push @nodes, "sink";
				}
				# else
				{
					unshift @nodes, "start";
				}
				print ">$paths_like[$c]";
				print_path_in_array_condensed(@nodes);
				$iter++;
			}
		}
	}
}

sub print_path_in_array_condensed {
	my(@arr) = @_;
	$pathcounter++;
#	if($pathcounter %10000 ==0) { $duration = time - $initial_time; $initial_time= time ; print " Written $pathcounter paths in $duration seconds \n"; }
	my($str)  = $cond_ver{$arr[1]};
#	my($str) = $arr[0];
	for(my($i)=2; $i<scalar(@arr)-1;$i++) 
	{
		$str = $str.substr($cond_ver{$arr[$i]},$K-1);
	#	$str = $str." ".$arr[$i];
	#	print allpath substr($arr[$i],$K-1);
	}
	print ">$pathcounter"."_"."$cond_ver{$arr[0]}_".length($str)."\n";
	print $str;
	print "\n"; 
}