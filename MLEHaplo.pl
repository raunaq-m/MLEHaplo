# Wrapper main file for MLEHaplo: A Maximum Likelihood Estimation for viral population reconstruction
# USAGE: perl MLEHaplo.pl -fas <fastafile> -out <mlehaplofile> -k <kmersize> -IS <insert_size> -M <ViPRA M Factor> -GL <genome length>

# The program takes as input fasta file and outputs a MLEHaplo predicted output 
#
use Getopt::Long;
use File::Basename;

$fastafile = "";
$fastqfile = "";
$flag_fasta = 0;
$kmersize = 55;
$base_dir = "";
$single_file = -1;
$M_factor = 15;
$genome_length = 1200;
$insert_size = 400; 
$multiresoutput = "";

## Threshold for k-mer storage in DE BRUIJN GRAPH 
$threshold = 0;
$threshold_paired = 5;
# File to store k-mer counts
@writefiles = ();

# Graph file name & Paired Set file 
$graphwritefile = "";
$pairedsetfile = "";


GetOptions("fas=s",\$fastafile,"fas1:s",\$fasta_paired1,"fas2:s",\$fasta_paired2,"MR:s",\$multiresoutput,"k=i",\$kmersize,"M=i",\$M_factor,"IS=i",\$insert_size,"GL=i",\$genome_length,"out=s",\$outputfile);

main();

sub main
{
	parse_args();
	if($multiresoutput ne "")
	{
		print "K-mer counts obtained from output of MultiRes\nUsing file $multiresoutput for k-mer counts\n";
	}
	# Run the k-mer counting as before.
	run_kmer_counting();
	#Storing the graph
	# Generating the paired-set file
	$graphwritefile = $writefiles[0].".graph";
	$pairedsetfile =  $writefiles[0].".pk.txt";
	if($multiresoutput eq "")
	{
		generate_de_bruijn_graph($writefiles[0],$graphwritefile);
		generate_paired_set($writefiles[0]);
	}
	else
	{
		generate_de_bruijn_graph($multiresoutput,$graphwritefile);
		generate_paired_set($multiresoutput);
	}
		
	# Running ViPRA algorithm 
	# Input files: graph file 
	# 			   kmers file 
	# 			   pairedsetfile
	# $M_factor = 15;
	$ViPRAoutputfile = $writefiles[0].".fact$M_factor".".txt";
	# $insert_size = 400;
	print "Running ViPRA\n";
	if($multiresoutput eq "")
	{
		print "Command: perl dg_cover.pl -graph $graphwritefile -kmer $writefiles[0] -paired $pairedsetfile -fact $M_factor -thresh $threshold -IS $insert_size > $ViPRAoutputfile\n";
		$cmd = `perl dg_cover.pl -graph $graphwritefile -kmer $writefiles[0] -paired $pairedsetfile -fact $M_factor -thresh $threshold -IS $insert_size > $ViPRAoutputfile`;
	}
	else
	{
		print "Command: perl dg_cover.pl -graph $graphwritefile -kmer $multiresoutput -paired $pairedsetfile -fact $M_factor -thresh $threshold -IS $insert_size > $ViPRAoutputfile\n";
		$cmd = `perl dg_cover.pl -graph $graphwritefile -kmer $multiresoutput -paired $pairedsetfile -fact $M_factor -thresh $threshold -IS $insert_size > $ViPRAoutputfile`;
	}
	print "ViPRA Execution Done\n";
	
	# Parsing output files from ViPRA
	$ViPRAfastafile = $writefiles[0].".fact$M_factor".".fasta";
	$ViPRApathsfile = $writefiles[0].".fact$M_factor".".paths.txt";
	print "Processing ViPRA output to two files $ViPRAfastafile and $ViPRApathsfile\n";
	$cmd = `perl process_dg.pl $ViPRAoutputfile > $ViPRAfastafile`;
	$cmd = `perl get_paths_dgcover.pl -f $ViPRAoutputfile -w $ViPRApathsfile`;
	print "ViPRA Output Processing Done\n";	
	
	# Running MLEHaplo maximum likelihood estimation 
	$condensedgraphfile = $writefiles[0].".cond.graph";
	$compatiblesetfile  = $writefiles[0].".comp.txt";
	$likelihoodoutputfile = $writefiles[0].".smxlik.txt";
	# $genome_length = 1200;
	print "Running MLEHAplo with command:";
	print "perl likelihood_singles_wrapper.pl -condgraph $condensedgraphfile -compset $compatiblesetfile -pathsfile $ViPRApathsfile -back -gl $genome_length -slow > $likelihoodoutputfile\n";
	$cmd = `perl likelihood_singles_wrapper.pl -condgraph $condensedgraphfile -compset $compatiblesetfile -pathsfile $ViPRApathsfile -back -gl $genome_length -slow > $likelihoodoutputfile`;
	print "MLEHAplo Done\n";
	
	# Generating final output file 
	$cmd = `perl extract_MLE.pl -f $ViPRAfastafile -l $likelihoodoutputfile > $outputfile`;
	print "Storing Final output in $outputfile\n";
}

sub generate_paired_set
{
	$kmercountsfile = $_[0];
	if ($single_file ==1)
	{
		print "Generating Paired Set file $pairedsetfile\n";
		print "Command: perl construct_paired_without_bloom.pl -fasta $fastafile -kmerfile $kmercountsfile -thresh $threshold_paired -wr $pairedsetfile\n";
		#$cmd = `perl construct_paired_without_bloom.pl -fasta $fastafile -kmerfile $kmercountsfile -thresh $threshold_paired -wr $pairedsetfile`;
		system("perl construct_paired_without_bloom.pl -fasta $fastafile -kmerfile $kmercountsfile -thresh $threshold_paired -wr $pairedsetfile");
		print "Paired set generated \n";
	}
	elsif($single_file ==0) 
	{
		print "Generating Paired Set file from paired fasta files\n";
		print "Command: perl construct_paired_without_bloom.pl -file1 $fasta_paired1 -file2 $fasta_paired2 -paired -kmerfile $kmercountsfile -thresh $threshold_paired -wr $pairedsetfile\n";
		system("perl construct_paired_without_bloom.pl -file1 $fasta_paired1 -file2 $fasta_paired2 -paired -kmerfile $kmercountsfile -thresh $threshold_paired -wr $pairedsetfile");
		print "Paired Set Generated \n";
	}
}

sub generate_de_bruijn_graph
{
	my($kmercountsfile) = $_[0];
	my($output_file) = $_[1];
	# Store the output of the graph in the file $output_file
	print "Generating the De Bruijn graph \n";
	print "Command:\nperl construct_graph_kmer.pl -k1 $writefiles[1] -k $kmercountsfile -t $threshold -w $output_file\n";
	$cmd = `perl construct_graph_kmer.pl -k1 $writefiles[1] -k $kmercountsfile -t $threshold -w $output_file`;
	print "Graph generation done\n";
}

sub run_kmer_counting
{
	if( $single_file ==1)
	{
		print "Running multi-dsk counting on single fasta file\n";
		print "Executing commmand: \"multi-dsk/multi-dsk $fastafile listofkmers.txt -m 8192 -d 10000\"\n";
		system("multi-dsk/multi-dsk $fastafile listofkmers.txt -m 8192 -d 10000");
	}
	elsif( $single_file == 0)
	{
		print "Running multi-dsk for the two fasta file $fasta_paired1 and $fasta_paired2\n";
		open(temp,">paired_files.txt");
		print temp "$fasta_paired1\n";
		print temp "$fasta_paired2\n";
		close temp;
		print "Executing commmand: \"multi-dsk/multi-dsk paired_files.txt listofkmers.txt -m 8192 -d 10000\"";
		system("multi-dsk/multi-dsk paired_files.txt listofkmers.txt -m 8192 -d 10000\n");
	}
	else 
	{
		die("Function: run_kmer_counting: Error in fasta single file or paired end files\n");
	}
	@kmerlists = `ls *solid_kmers_binary*`;
	@writefiles = ();
	for $k (@kmerlists)
	{
		chomp($k);
		print "$k\n";
		$writefile =$k;
		$writefile =~ s/solid_kmers_binary.//;
		print "Writing kmer counts of $k to $writefile\n";
		push @writefiles, $dirname."/".$writefile;
		$cmd = `python multi-dsk/parse_results.py $k > $dirname/$writefile`;
		$cmd = `rm $k`;
	}
	# Clean up temporary files 
	$cmd = `rm write_counts.txt`;
	$cmd = `rm *reads_binary`;
}
sub parse_args
{
	# Parse the arguments entered and generate files need for processing
	# If no fasta file provided, kill the program
	if($fastafile eq "" && $fasta_paired1 eq "")
	{
		$single_file = -1;
		die("Enter at least a fasta file containing the reads or a paired of fasta files\nUSAGE: perl MLEHaplo.pl -fas <fastafile> -out <mlehaplofile> -k <kmersize> -IS <insert_size> -M <ViPRA M Factor> -GL <genome length>\n");
	}
	elsif ($fastafile ne "" )
	{
		$single_file = 1;
		$dirname = dirname($fastafile);
		print "Single fasta file entered...\nAssuming the single file as an interleaved paired-end fasta file\n";
	}
	elsif ($fasta_paired1 ne "" && $fasta_paired2 ne "")
	{
		$single_file = 0;
		$dirname = dirname($fasta_paired1);
		print "Two fasta files entered...\nAssuming them as forward and reverse reads of a paried end files\n";
	}
	else
	{
		$single_file = -1;
		die("Enter at two paired-end fasta files as follows:\nUSAGE: perl MLEHaplo.pl -fas1 <fasta_file> -fas2 <fasta_file> -k <kmer_size> -out <output_file> -IS <insert_size> -M <ViPRA M Factor> -GL <genome length>\n");
	}
	
	# Setting the file for k-mer counting using MultiDsk. Use the k-mer size given by user and (k+1)
	$kmer_graph = $kmersize+1;
	open(file,">listofkmers.txt");
	print file "$kmer_graph\n";
	print file "$kmersize\n";
	close file;
}
