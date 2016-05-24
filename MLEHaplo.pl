# Wrapper main file for MLEHaplo: A Maximum Likelihood Estimation for viral population reconstruction
# USAGE: perl MLEHaplo.pl -fas <fastafile> -out <mlehaplofile> -k <kmersize>
# The program takes as input fasta file and outputs a MLEHaplo predicted output 
#
use Getopt::Long;
use File::Basename;

$fastafile = "";
$fastqfile = "";
$flag_fasta = 0;
$kmersize = 55;
$base_dir = "";

GetOptions("fas=s",\$fastafile,"k=i",\$kmersize,"out=s",\$outputfile);

main();

sub main
{

}


