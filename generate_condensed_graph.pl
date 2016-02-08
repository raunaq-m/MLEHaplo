# Generate a condensed graph 
use Bio::SeqIO;
use List::Util qw(shuffle max min);
use POSIX qw(ceil);
use Getopt::Long;

# File name variables
$graph_file = "";
$kmerfile = "";
$thresh = 0;
$wr_file="";
#######################################################################
GetOptions ( "graph=s", \$graph_file, "kmer=s",\$kmerfile, "thresh:i",\$thresh)
		or die("Error in inputting parameters \n");

#######################################################################

%kmerhash = (); # Contains the kmers to be used in the graph nodes
%graphout = (); %graphin = ();# Double hash containing the graph structure
%vertex = (); # Hash of the vertices in the graph 
%pred = ();%start=();%finish=(); %depth = ();
@source_nodes = (); # Hash for ruling out source/sink nodes
@sink_nodes = ();
$start_preference = 0;
$K=0; # Kmer length on which graph is based 
loadkmerfile($thresh); # load the kmer counts data above the threshold
loadgraph();  # Create the graph of the reads
$nokmers = scalar(keys %vertex); #nokmers is the number of vertices in the graph
print "Parameters : graph $graph_file kmer $kmerfile thresh $thresh\n";
print "Graph loaded with no of vertices $nokmers \n";
print "K value of the kmer is $K \n";
print "Number of source nodes is ".scalar(@source_nodes)."\n";
print "Number of sink nodes is ".scalar(@sink_nodes)."\n";
#####################################################################

##################################################################### 
# Set the order of visits in depth first search and print out paths that follow the depth first search order
depth_order();
#depthfirst();
#####################################################################
condense_graph_wrapper();
#####################################################################


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
			# $rev_kmer = revcomplement($kmer); # ADD the reverse complement directly to the kmer hash 6/18/2015
			# if (!exists($kmerhash{$rev_kmer}) ) 
			# {
				# $kmerhash{$rev_kmer} = $c;
			# }
		}
	}
	$totalkmers = scalar(keys %kmerhash);
	$K = length($kmer);
	print "Total kmers from $kmerfile is $totalkmers and $K is k length\n";
	close file;
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
		# if((exists($kmerhash{$u})|| exists($kmerhash{revcomplement($u)})) && (exists($kmerhash{$v})||exists($kmerhash{revcomplement($v)})) )
		if(exists($kmerhash{$u}) && exists($kmerhash{$v}))
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
#			print "$kmerhash{$u} $kmerhash{$v} \n";
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
#	node_stats();
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
