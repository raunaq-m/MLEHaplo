# Perl script to convert a condensed graph generated for the Viral Populations
# to a fastg format, so that it can be visualized and observed easily. 
# USAGE : perl convert_condensed_to_fastg.pl -condensed condensedgraphfile -output outputfastgfile
use Getopt::Long;

$condensedgraphfile = "";
$outputfastgfile = "Condensedgraph.fastg";
$kmer_size = 0;

GetOptions("condensed=s",\$condensedgraphfile,"output:s",\$outputfastgfile,"kmer_size=i",\$kmer_size)
	or die("Error in input arguments, Enter a condensed graph file\n");
die("Error in input arguments, Enter a condensed graph file\nperl convert_condensed_to_fastg.pl -condensed condensedgraphfile -output outputfastgfile\n") if $condensedgraphfile eq "";
main();

sub loadcondensedgraph
{
	open(file,$condensedgraphfile) or die("Error in loading file $condensedgraphfile\n");
	while($l=<file>)
	{
		chomp($l);
		($v1,$v2) = split(/ /,$l);
		if($v2 =~ /[AGCT]/)
		{
			# node descriptions 
			$node_seq{$v1}  = $v2;
		}
		else
		{
			# store the condensed graph
			$graphout{$v1}{$v2} = 1 ;
			$source{$v2} = 0;
			$sink{$v1} = 0;
		}
	}
	# Creating list of source nodes
	for $v(keys %node_seq)
	{
		if(!exists($source{$v})) 
		{
			$source_nodes{$v}=1;
		}
	}
	close file;
	undef %source;
	print "The graph has ".scalar(keys %node_seq)." nodes\n";
	print "There are ".scalar(keys %source_nodes)." source nodes\n";
}

sub writefastg
{
	open(file,">$outputfastgfile");
	for $k(sort {$a<=>$b} keys %node_seq)
	{
		$line = ">NODE_$k"."_length_".length($node_seq{$k})."_cov_1.0:";
		# print file ">$k:";
		foreach (keys %{$graphout{$k}})
		{
			# print file "$_,";
			$line = $line."NODE_$_"."_length_".length($node_seq{$_})."_cov_1.0,";
		}
		$line =~ s/,$/;/;
		$line =~ s/:$/;/;
		print file "$line\n";
		if(!exists($source_nodes{$k}))
		{
			print file substr($node_seq{$k},$kmer_size-1)."\n";
		}
		else 
		{
			print file $node_seq{$k}."\n";
		}
	}
}
sub main
{
	loadcondensedgraph();
	writefastg();
}
