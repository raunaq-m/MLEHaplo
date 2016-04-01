#### Algorithm for cover of all nodes up till the given node in computation


# usage: perl dominant_path.pl <edge_file> <kmer_counts_file>

use lib "/i3c/alisa/rom5161/rh6/perl/lib/perl5";
use Bio::SeqIO;
use List::Util qw(shuffle max min);
use POSIX qw(ceil);
use Getopt::Long;
 
$graph_file = "";
$kmerfile = "";
#if( defined($ARGV[2])) { $thresh = $ARGV[2];} else { $thresh = 2; }
$thresh = 0;
$wr_file="";
$paired_kmers_file = ""; # change this later only temporary 
$factor = 0; ## USE FACTOR TO SELECT $FACTOR TOP PATHS AT A NODE
$numpaths=''; 
$specific_node='';
$node1 = -1;
$depth_of_specific = 0;
$depth_of_end = 0;
$score = -10;
$start_preference='';
$start_order_file = "";
# Print the condensed graph nodes
$IS_factor = 500; # InsertSize - ReadLength + Kmer (350-150+60) 
#$IS_factor = 291; # InsertSize -Kmer + 1  when largest k-mer pair is used
#$file_base_seq = $ARGV[4];
#######################################################################
GetOptions ( "graph=s", \$graph_file, "kmer=s",\$kmerfile, "thresh:i",\$thresh,
		"paired:s",\$paired_kmers_file,"fact:f",\$factor, "sp" =>\$specific_node, "node:i", \$node1 , "d1:i" => \$depth_of_specific, "d2:i" => \$depth_of_end, "IS:i",\$IS_factor, 
		"start_order:s",\$start_order_file,"so" =>\$start_preference)
		or die("Error in inputting parameters \n");


#######################################################################
%kmerhash = (); # Contains the kmers to be used in the graph nodes
%graphout = (); %graphin = ();# Double hash containing the graph structure
%vertex = (); # Hash of the vertices in the graph 
%pred = ();%start=();%finish=(); %depth = ();
@source_nodes = (); # Hash for ruling out source/sink nodes
@sink_nodes = ();
$K=0; # Kmer length on which graph is based
print "Threshold for kmer filtering is $thresh\n"; 
loadkmerfile($thresh); # load the kmer counts data above the threshold
loadgraph();  # Create the graph of the reads
$nokmers = scalar(keys %vertex); #nokmers is the number of vertices in the graph
print "Parameters : graph $graph_file kmer $kmerfile thresh $thresh paired $paired_kmers_file factor $factor IS $IS_factor starts $start_order_file\n";
print "Graph loaded with no of vertices $nokmers \n";
print "K value of the kmer is $K \n";
print "Number of source nodes is ".scalar(@source_nodes)."\n";
print "Number of sink nodes is ".scalar(@sink_nodes)."\n";
print "Computing the maxflow paths in the graph \n";

# Set the order of visits in depth first search and print out paths that follow the depth first search order
depth_order();
#depthfirst();
#####################################################################
loadpaired() if $paired_kmers_file ne "";
condense_graph_wrapper();

if($specific)
{
	$num_comp_pairs = generate_compatible_set() if $paired_kmers_file ne "";
	print "Total paired k-mers selected is ".$num_comp_pairs."\n";
	$pathcounter = 0;
	$paths_traversed = 0;
	condgraph_number_of_paths_wrapper();
	print "Computing paths in range less than $depth_of_specific\n";
	_path_all_nodes_backward_specific($depth_of_specific,$depth_of_end);
	_redundant_constraints();
	print_constraints();
	print "$factor is factor\n";
}
else 
{
	$num_comp_pairs = load_paired_feature_select_updated() if $paired_kmers_file ne "";
	$num_comp_pairs = generate_compatible_set() if $paired_kmers_file ne "";
	print "Total paired k-mers selected is ".$num_comp_pairs."\n";
	$pathcounter = 0;
	$paths_traversed = 0;
	condgraph_number_of_paths_wrapper();
	_path_all_nodes_backward();
	path_all_nodes(); #print_constraints();
	_redundant_constraints();
	print_constraints();
	print "$factor is factor\n";

	
	$initial_time  = time;
	@tmp = ();
}


sub loadgraph
{
	open(file,$graph_file);
	print "Opened $graph_file\n";
	$edges=0;
	while($l=<file>)
	{
		chomp($l);
		($u,$v) = split(/ /,$l);
		if((exists($kmerhash{$u})||exists($kmerhash{revcomplement($u)})) && (exists($kmerhash{$v})||exists($kmerhash{revcomplement($v)})) )
		{
			$graphout{$u}{$v} = 1;
			$graphin{$v}{$u} = 1;
			$vertex{$u} = 1;
			$vertex{$v} = 1;
			$source{$v} = 0;
			$sink{$u} = 0;
			$edges++;
		}
#		else { 
#			print "$u $kmerhash{$u} $v $kmerhash{$v} \n";
#		}
	}
	close file;
	for $k(keys %vertex) 
	{
		if(!exists($source{$k})) { push @source_nodes, $k; }
		if(!exists($sink{$k})) { push @sink_nodes, $k; }
	}	
	print "Number of edges read is $edges\n";
	undef %sink; # forget %sink ever existed
	undef %source; # forget %source ever existed
}



sub loadkmerfile
{
	#only load the kmers that are above a certain threshold 
	my($t) = $_[0];
	open(file,$kmerfile);
	while($l=<file>)
	{
		chomp($l);
		($kmer,$c) = split(/ /,$l);
		if($c>$t) 
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
	print "Total kmers from $kmerfile is $totalkmers and $K is k length\n";
	close file;
}

sub loadpaired 
{
	open(filepaired,$paired_kmers_file);
	while($l=<filepaired>) 
	{
		chomp($l);
		my ($lk1,$rk2,$cc) = split(/ /,$l);
		#if( exists($kmerhash{$lk1}) || exists($kmerhash{revcomplement($lk1)}) ) && ( exists($kmerhash{$rk2}) || exists($kmerhash{revcomplement($rk2)}) )
		if(exists($kmerhash{$lk1}) && exists($kmerhash{$rk2})) 
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
	print "Avg insert size is $avg/$count =  ".$avg/$count."\n";
	undef %paired_depth;
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


sub depth_order
{
 #function to find the highest depth nodes in the graph
         $time=0;
         $start_time=0;
         $counter=0;
         $time_=time;
         for $kv(keys %vertex) {
                 $discovered{$kv}=0;
                 $depth{$kv} =0;
         }
	 $depth_file = $graph_file;
	 $depth_file =~ s/graph/depth/;
         open(wr_fi,">$depth_file");
		 if ($start_preference) 
		 {
			@start_node_order = ();
			open(strt_order, $start_order_file);
			while($l=<strt_order>)
			{
				chomp($l);
				if(exists($vertex{$l}))
				{
					push @start_node_order, $l;
				}
			}
			close strt_order;
			for $k( @start_node_order)
			{
				$counter++;
				$pred{$k} = undef;
				$depth{$k} = 0;
				$root{$k} = $k;
				colourneighbours($k,$k);
			}
		 }
		 else
		 {
			 #DEFAULT : 
			 # Visiting the source nodes in the order of their counts ( higher counts first )
			 for $k( sort {$kmerhash{$b} <=> $kmerhash{$a}} @source_nodes) {
					$counter++;
					$pred{$k}=undef;
					$depth{$k}=0;
					$root{$k}=$k;
					colourneighbours($k,$k);   
	 #              if(($counter%1000)==0) 
	 #				{
	 #                        $duration=time-$time_;
	 #                        print "Finished $counter/".scalar(@source_nodes)." source nodes duration $duration seconds\n";
	 #              }
			 }
		 }
        for $k(keys %vertex )
         {
                 $visited_order{$k} = 0;
         }
         $order=0;$counter=0;
         for $k( @source_nodes ) { $start_order{$k} = 0 ; } # Define the start order of the nodes
	  print wr_fi "#".scalar(keys %vertex)." nodes and depth nodes ".scalar(keys %depth)."\n";
         for $k( sort {$depth{$b} <=> $depth{$a}} keys %depth ) # order nodes in decreasing order of depth 
         {
		$counter++;
                 my($temp)= $root{$k};
		 if($start_depth{$temp} < $depth{$k} ) { $start_depth{$temp} = $depth{$k}; }
                 if ($visited_order{$k}==0) {
                 if($start_order{$temp} ==0) { $order++; $start_order{$temp} = $order; }
                 #print wr_fi "$k $depth{$k} $temp $start_order{$temp}\n"; 
                 } 
                 print wr_fi "$k $depth{$k} $temp $start_order{$temp} \n";
		 #else {
                 #print wr_fi "\n";
                 #}
#                 if(($counter%1000)==0) {
 #                        $duration=time-$time_;
  #                       print "Finished $counter/".scalar(keys %depth)."nodes duration $duration seconds\n";
 #               }
         }
        print wr_fi "#Starting order for graph start nodes \n";
         for $k( sort {$start_order{$a} <=> $start_order{$b} } keys %start_order) { print wr_fi "$k $start_order{$k} start depth $start_depth{$k}\n"; }
	 $node_depth_file = $graph_file;
	 $node_depth_file =~ s/graph/nodedepth/; 
	 open(depth_of_nodes,">$node_depth_file");
	for $k(sort {$depth{$a}<=>$depth{$b}} keys %depth) {
		print depth_of_nodes "$k $depth{$k} $root{$k}\n";
	}
}


sub colourneighbours(){
	my($v) = $_[0];
	my($head) = $_[1];
	$discovered{$v} = 1;
	$time++;
	$start{$v} = $time;
	foreach $u(sort {$kmerhash{$b} <=> $kmerhash{$a}} keys %{$graphout{$v}}) {
		if($discovered{$u}==0) {
			$pred{$u} = $v;
			$root{$u} = $head;
			$depth{$u} = $depth{$v}+1;
			colourneighbours($u,$head);
		}
	}
	$discovered{$v}=2;
	$time++;
	$finish{$v} = $time;
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


sub load_paired_feature_select_updated
{
	# Load only those paired k-mers that span across two bubbles, as the rest of them are uninformative in regards to resolving the bubbles 
	# Use the depth of the paired k-mers to inform which sets of paired k-mers will be informative 
	# Traverse the nodes that form a bubble, and fill up their compatible sets by traversing other nodes that are within
	# insert size distance 
	#
	# FIRST: obtain a list of nodes that form bubbles: Bubbles are those that have only one incoming and one outgoing edge 
	# Assuming the following variables are already defined: %depth, %cond_graphout %cond_graphin, %cond_ver_suffix
	# %cond_ver_prefix, %cond_ver
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
 	my($fname) = $graph_file;
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


sub condense_graph_wrapper
{
	%graphin = ();%cond_ver = ();
	create_ingraph();
	$cond_file = $graph_file;
	$cond_file =~ s/graph/cond.graph/;
        open(wr_fi,">$cond_file");
	for (keys %vertex) {
		$condensed{$_} = 0;
	}
	for $k( @source_nodes) 
	{
		push @queue, $k; #every source node will start a new condensed node
	}
	$vno = 0;
	while(scalar(@queue)!=0) {
		my($n1) = shift(@queue);
	#	print "Opening queue $n1\n";
		if($condensed{$n1} == 0) {
#			print "Condensing node $n1\n";
			my($tt) = condense($n1);
			print wr_fi "$vno $tt\n";
			$cond_ver{$vno} = $tt;
			$cond_ver_suffix{substr($tt,-$K) } = $vno;
			$cond_ver_prefix{substr($tt,0,$K)} = $vno;
			for (my($i)=0;$i<length($tt)-$K+1;$i++) { #Assign each k-mer in a condense node to its vertex ID
				$vertex{substr($tt,$i,$K)} = $vno; 
			}
			$vno++;
			$condensed{$n1} = 1;
		}
	}
	%source = (); %sink = ();
	for $uu(keys %cond_ver_suffix)
	{
		for $vv(keys %{$graphout{$uu}}) {
			$cond_graphout{$cond_ver_suffix{$uu}}{$cond_ver_prefix{$vv}} = 1;
			$cond_graphin{$cond_ver_prefix{$vv}}{$cond_ver_suffix{$uu}} = 1;
			$source{$cond_ver_prefix{$vv}} = 0; $sink{$cond_ver_suffix{$uu}} = 0;
		}
	}
	for $vv(keys %cond_ver) { if (!exists($source{$vv})) { push @cond_source_nodes, $vv; }	}
	for $vv(keys %cond_ver) { if (!exists($sink{$vv}) )  { push @cond_sink_nodes , $vv; }	}
	undef %source;
	undef %sink;
#	for (@cond_source_nodes ) {print wr_fi "Source $_ $cond_ver{$_}\n";}
#	for (@cond_sink_nodes ) {print wr_fi "Sink $_ $cond_ver{$_}\n";}
#	for (@sink_nodes ) {print wr_fi "Sink $_ \n";}
	print "Number of source and sink nodes ".scalar(@cond_source_nodes)." ".scalar(@cond_sink_nodes)."\n";
# Prin the condensed graph if you need it 
	for $uu(keys %cond_graphout) { for $vv(keys %{$cond_graphout{$uu}} ) { print wr_fi "$uu $vv\n"; } }
	undef %graphin ;
	close wr_fi;
}


sub condense
{
	my($start) = shift;
	my($condensed_start) = $start;
	@next = keys %{$graphout{$start}};
	while(scalar(@next)==1) {
		$nnode = $next[0];
		if(scalar(keys %{$graphin{$nnode}})==1 && scalar(keys %{$graphout{$nnode}})==1) {
			$condensed{$nnode} = 1;
			$condensed_start = $condensed_start.substr($nnode,-1);
		}
		else {
			if(scalar(keys %{$graphin{$nnode}})==1  && scalar(keys %{$graphout{$nnode}})<1)
			{
			#	push @queue, $nnode;
				$condensed{$nnode} = 1;
				$condensed_start = $condensed_start.substr($nnode,-1);
				return $condensed_start;
			} 
			if(scalar(keys %{$graphin{$nnode}})>1) # && scalar(keys %{$graphout{$nnode}})<=1) 
			{ 
				push @queue, $nnode;
				return $condensed_start;
			}
			if(scalar(keys %{$graphout{$nnode}})>1) {
				$condensed{$nnode} =1;
				$condensed_start = $condensed_start.substr($nnode,-1);
			}
		}
		@next = keys %{$graphout{$nnode}};
	}
	foreach (@next) {
		if($_ ne "") {	 push @queue, $_; }
	#			print "Adding queue $nnode\n";}
	} 
	return $condensed_start;
}


sub create_ingraph
{
	for $p1(keys %graphout)
	{
		for $p2(keys %{$graphout{$p1}}) 
		{
			$graphin{$p2}{$p1} = 1;
		}
	}
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
	#BLOCKED for dynamic testing 
	$number_of_paths = find_number_paths("source","sink",\%cond_graphout);
	print "Total number of paths 							$number_of_paths\n";
#	$file_all_paths = $graph_file;
#	$file_all_paths =~ s/graph/allc.fasta/;
#	open(allpath,">$file_all_paths");
#	close allpath;
}


sub path_all_nodes
{
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
		$marked{$source} = 1;
	}
	for $node ( sort {$depth{substr($cond_ver{$a},0,$K)} <=> $depth{substr($cond_ver{$b},0,$K)}} keys %cond_ver ) 
	{
		if($marked{$node} ==0) 
		{
			path_two_nodes_dg_forward($node);
			$marked{$node} = 1;
		}
	}
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
						$temp_like +=2;
						# $temp_like += 2*$compatible_set_mem{$nodes[$iter]}{$start};
					}
					else
					{
						$temp_like +=1;
						# $temp_like += 1*$compatible_set_mem{$nodes[$iter]}{$start};
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
					 $temp_like +=1;
					 # $temp_like += 1*$compatible_set_mem{$nodes[$iter]}{$start};
				}
				else
				{
					$temp_like +=0.5;
					# $temp_like += 0.5*$compatible_set_mem{$nodes[$iter]}{$start};
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

sub _compute_membership_paired_nodes
{
	for $k(keys %compatible_set)
	{
		# print $k." compatible\t";
		my($total) = scalar(keys %{$compatible_set{$k}});
		foreach (keys %{$compatible_set{$k}})
		{
			$compatible_set_mem{$k}{$_} = $compatible_set{$k}{$_}/$total; # Modified 4/2/2015
			# print $compatible_set_mem{$k}{$_}."\t";
		}
		# print "\n";
	}
}

sub _redundant_constraints
{
	$all_constraints = "";
	$all_constraints_like = "";
	$back_constraints = "";
	# for $node ( sort {$depth{substr($cond_ver{$b},0,$K)} <=> $depth{substr($cond_ver{$a},0,$K)}} @cond_sink_nodes )
	# {
		# $i = 0;
		# for(my($i)=0; $i< scalar(@{$node_paths{$node}}); $i++)
		# {
			# ${$node_paths{$node}}[$i] =~ s/^_//;
			# ${$node_paths{$node}}[$i] =~ s/_$//;
			# if( $all_constraints !~ /${$node_paths{$node}}[$i]/ )
			# {
				# $all_constraints = $all_constraints."#".${$node_paths{$node}}[$i];
				# $all_constraints_like = $all_constraints_like."#".${$node_paths_like{$node}}[$i];
			# }
		# }
	# }
	
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
			# for(my($i)=0; $i<scalar(@{$node_paths_back{$node}}); $i++)
			# {
				# $back_constraints = $back_constraints."#".$
				# $back_constraints_like = $back_constraints_like."#".;
			# }
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
						$temp_like += 2;
						# $temp_like +=2*$compatible_set_mem{$start}{$nodes[$iter]};
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
					$temp_like += 1;
					# $temp_like += 1*$compatible_set_mem{$start}{$nodes[$iter]};
				}
				else
				{
					$temp_like +=0.5;
					# $temp_like +=0.5*$compatible_set_mem{$start}{$nodes[$iter]};
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


