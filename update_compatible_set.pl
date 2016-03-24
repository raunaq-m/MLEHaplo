# Update the compset file and remove the vertex pairs with low counts. 
# vertex pair to remove is the one which together have low counts, remove that pair from the set
# USAGE: perl update_compatible_set.pl -cs compatible_set -w writeoutputfile -t thresholdForPairedCounts -wg UpdatedGraph
use Getopt::Long;

$comp_set = "";
$threshold = 10;
$updated_compset = "";
$updated_graph = "";
$cond_graph = "";
GetOptions("cs=s",\$comp_set,"w=s",\$updated_compset,"t:i",\$threshold,"wg:s",\$updated_graph) or die("Enter input correctly\n");
$cond_graph = $comp_set;
$graph = $comp_set;
if( $cond_graph =~ /comp.txt/)
{
	$cond_graph =~ s/comp.txt/cond.graph/;
	$graph =~ s/comp.txt/graph/;
	print $cond_graph."\n";;
}
elsif ($cond_graph =~ /bubble.txt/) 
{
	$cond_graph =~ s/bubble.txt/cond.graph/;
	$graph =~ s/bubble.txt/graph/;
	print $cond_graph."\n";
}
else 
{
	$cond_graph = "";
	$graph ="";
}

main();

sub main
{
	load_compatible_set();
	open(temp,$graph) if ($graph ne ""); 
	$line = <temp>; chomp($line); ($k1,$k2) = split(/ /,$line); 
	$Kmersize = length($k1); 
	close temp;
	load_cond_graph() if($cond_graph ne "");
	load_graph() if($graph ne "");
}

sub load_compatible_set
{
	# Load the compatible set first
	open(file,$comp_set) or die("Error opening $comp_set file\n");
	open(wrfile,">$updated_compset");
	while($l=<file>)
	{
		chomp($l);
		@vs = split(/ /,$l);
		%valid_position = ();
		%valid_counts = ();
		($id,$position) = split(/:/,$vs[0]); # vertex and its depth
		for ($i=1;$i<scalar(@vs); $i++) # iterate through all the vertex pairs in the adjacency list 
		{
			($id2,$position2,$count) = split(/:/,$vs[$i]);
			if($count > $threshold)
			{
				$valid_position{$id2} = $position2;
				$valid_counts{$id2} = $count;
			}
		}
		if (scalar(keys %valid_counts)>0)
		{
			print wrfile $vs[0]." ";
			for $k(sort {$valid_position{$a}<=>$valid_position{$b}} keys %valid_position)
			{
				print wrfile "$k:$valid_position{$k}:$valid_counts{$k} ";
				$condvertex{$id} = 1;
				$condvertex{$k} = 1;
			}
			print wrfile "\n";
		}
	}
	close file;
	close wrfile;
}

sub load_cond_graph
{
	open(file, $cond_graph) or die("Condensed graph didn't load");
	$not_existing = 0;$no_vertex = 0;
	while($l=<file>)
	{
		chomp($l);
		if($l =~ /[AGCT]/)
		{
			($vertex,$seq) = split(/ /,$l);
			$cond_ver{$vertex} = $seq;
			if(!exists($condvertex{$vertex}))
			{
				$len = length($seq); 
				for ($i=0;$i<$len-$Kmersize+1; $i++)
				{
					$remove_vertex{substr($seq,$i,$Kmersize)} = 1;
				}
				$no_vertex ++;
			}
		}
		else
		{
			($ver1,$ver2) = split(/ /,$l);
			if(!exists($condvertex{$ver1}) || !exists($condvertex{$ver2}))
			{
				$not_existing++;
			}
		}
	}
	close file;
	print "$not_existing edges were lost\n$no_vertex vertices were lost\n";
}

sub load_graph
{
	open(file,$graph) or die("Graph file $graph did not load\n"); 
	open(wrfile,">$updated_graph");
	
	while($l=<file>)
	{
		chomp($l);
		($k1,$k2) = split(/ /,$l);
		if(!(exists($remove_vertex{$k1}) || exists($remove_vertex{$k2})))
		{
			print wrfile $l."\n";
		}
	}
	close file;
	close wrfile;
}
