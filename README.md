# MLEHaplo
Maximum Likelihood Estimation for Viral Populations

VIPRA and MLEHaplo README file

Pre-requisites:

1. multi-dsk: k-mer counting software ( Extension of dsk [http://minia.genouest.org/dsk/] version 1.5655)
  - Counts k-mers for multiple values of k simultaneously. The software doesn't combines the counts of forward and reverse complement k-mers, as is performed in traditional k-mer counting softwares
2. perl with modules Bio::Perl, Getopt::Long, Graph.   
  - BioPerl is available at [http://bioperl.org/](http://bioperl.org/)

# Dockerfile 
A dockerfile containing all the dependencies is now available. Find more documenation on Docker [here](https://docs.docker.com/).

To build a docker image, you would need Docker installed in your HPC or cloud environment. Create an empty folder with just the Dockerfile in it and type
```
docker build --tag mlehaplo:1.0 . 
```

This will create the docker container containing all relevant scripts and libraries for you to run MLEHaplo. You can now start a live container by typing: 

```
docker run -it mlehaplo:1.0
```

# MLEHaplo workflow 

The following is a list of steps/commands you'd have to follow to run MLEHaplo. 

## Preliminaries

#### Step 1: Generate k-mer counts file


**Command:** `multi-k fasta/fastq_file list_of_kvalues -d diskspace_limit -m memory_limit`

Example:
```
./multi-dsk/multi-dsk Example/paired-reads.fasta  Example/listofkmers.txt  -m 8192 -d 10000
```

- **fasta/fastq file**: is the file name of the fasta/fastq file. If there are a list of files, use a text file containing locations of all the files, one file per line.
For example, if there are two files containing paired reads, `file1.fastq` and `file2.fastq`, create a file `list_of_files.txt` containing following:
```
file1.fastq
file2.fastq
```

- **list_of_kvalues** : contains a list of k values for which k-mer counting has to be performed. For example if k-mer counting is desired for k-values 55,45,35,25. Create a file `list_of_kmers.txt` containing these numbers in decreasing order:
```
55
45
35
25
```

- **diskspace_limit** : defines the limit of temporary disk space (in MB) used while performing k-mer counting
- **memory_limit** : defines the memory limit (in MB) for storage of temporary hashes while performing k-mer counting

- Output of multi-dsk is a collection of files with extension `prefix.solid_kmers_binary."kvalue"` in compressed format, which contains counts of k-mers present in the fasta/fastq file.

#### Step 2: De-compress the output of multi-dsk
**Command:**
`parse_results prefix.solid_kmers_binary.kvalue_file  > prefix_file.kvalue`

Example:
```
./multi-dsk/parse_results Example/paired-reads.solid_kmers_binary.60 > paired-reads.60
```
- The file `fasta/fastq_file.kvalue` now contains the k-mer counts for the fasta/fastq file in the format "k-mer count" per line.

#### Step 3: Generate the De Bruijn graph
Generating the graph needs two files and a parameter. This will combine paired files into a single file.

1. `fasta_file` containing all the reads.
2. `kmer_file` generated above.
3. `threshold` value for ignoring erroneous k-mers.

**Command:**
`perl construct_graph.pl fasta_file kmer_file threshold graph_file "s"`

Example:
```
perl construct_graph.pl  Example/paired-reads.fasta paired-reads.60 0 paired-reads.60.graph "s"
```
- **"s"** parameter tells the perl script that there is a single fasta file of reads
- Output is the `graph_file` containing pairs of k-mers that form edges in the De Bruijn graph.

**TODO:** Add "From (k+1)-mer counts file"

#### Step 4: Create the paired set file
Create the paired set using the paired reads. It takes as input the two paired files,
`file1.fasta` and `file2.fasta`, the k-mer counts file
`file1.kvalue`, and a threshold for ignoring erroneous k-mers

**Command:**
`perl construct_paired_without_bloom.pl -file1 file1.fasta -file2 file2.fasta -paired -kmerfile file1.kvalue -thresh number -wr output_paired_set_file`

Example:
```
perl construct_paired_without_bloom.pl -fasta Example/paired-reads.fasta -kmerfile paired-reads.60 -thresh 0 -wr paired-reads.60.pk.txt
```

- Choice of threshold : Dependent on sequencing coverage. Lower threshold includes more erroneous k-mers in the graph, while higher threshold decreases the number of true k-mers and size of the graph.

- Output is a file that contains pairs of k-mers on a line and the number of times such pair is observed:
`kmer1 kmer2 count`


## VIPRA

#### Step A: Running VIPRA
Running the VIPRA algorithm takes inputs generated above and a parameter for the average insert size, threshold parameter and a value for M (factor) which decides the number of paths to generate per vertex

**Command:**
`perl dg_cover.pl -graph graph_file_from_step3 -kmer kmer_file_from_step2 -paired paired_set_from_step4 -fact M -thresh threshold_value -IS insert_size > vipra_output_file`

Example:
```
perl dg_cover.pl -graph paired-reads.60.graph -kmer paired-reads.60 -paired paired-reads.60.pk.txt -fact 15 -thresh 0 -IS 400 > paired-reads.60.fact15.txt
```
- `graph_file_from_step3` - output_file from [**Step 3: Generate the De Bruijn graph**](#step-3-generate-the-de-bruijn-graph)
- `kmer_file_from_step2` - output_file from [**Step 2: De-compress the output of multi-dsk**](#step-2-de-compress-the-output-of-multi-dsk)
- `paired_set_from_step4` - output_file from [**Step 4: Create the paired set file**](#step-4-create-the-paired-set-file)
- `vipra_output_file`: Contains the paths generated from the graph with high paired end supports.
- `prefix.comp.txt` : Contains a sets of paired vertices in the condensed graph that are compatible with each other based on the paired set.
- `prefix.cond.graph` : Contains the condensed version of De Bruijn graph, with non-branching paths condensed to a single vertex.

- Temporary output files:
  1. `prefix.bubble.txt`: Contains a sets of paired bubbles in the condensed graph that are compatible with each other based on the paired set.
  2. `prefix.depth` : Contains the depth first search traversal of the graph.
  3. `prefix.nodedepth` : Temporary file for debugging of code.


#### Step B: Generate fasta file
Extracting fasta file from outputfile

**Command:**
`perl process_dg.pl vipra_output_file > paths_fasta_file`

Example:
```
perl process_dg.pl paired-reads.60.fact15.txt > paired-reads.60.fact15.fasta
```

- Output: `paths_fasta_file`. The paths generated by VIPRA

#### Step C: Generate paths file for maximum likelihood estimation
Extracting just the paths in terms of nodes in the graph

**Command:**
`perl get_paths_dgcover.pl -f vipra_output_file -w paths_write_file`

Example:
```
perl get_paths_dgcover.pl -f paired-reads.60.fact15.txt -w paired-reads.60.fact15.paths.txt
```
- Output: `paths_write_file`


## MLEHaplo
Running MLEHaplo takes as input intermediate files generated by VIPRA [Step A: Running VIPRA](#step-a-running-vipra) and `paths_write_file` generated in [Step C: Generate paths file for maximum likelihood estimation](#step-c-generate-paths-file-for-maximum-likelihood-estimation)

**Command:**
`perl likelihood_singles_wrapper_parallel.pl -condgraph prefix.cond.graph -compset prefix.comp.txt -pathsfile paths_write_file -back -slow -gl approximate_genome_size  > MLE_textfile_output`

**Non Parallel Version Command:**
`perl likelihood_singles_wrapper.pl -condgraph prefix.cond.graph -compset prefix.comp.txt -pathsfile paths_write_file -back -gl approximate_genome_size -slow  > MLE_textfile_output`

Example
```
perl likelihood_singles_wrapper.pl -condgraph paired-reads.60.cond.graph -compset paired-reads.60.comp.txt -pathsfile paired-reads.60.fact15.paths.txt -back -gl 1200 -slow  > paired-reads.60.smxlik.txt
```
Required input files & Parameters

1. `prefix.cond.graph`: Condensed Graph file generated by VIPRA .
2. `prefix.comp.txt`: Compatible Set generated by VIPRA.
3. `paths_write_file`: Paths file from ViPRA Step C.
4. `approximate_genome_size`: Approximate genome Length.

Final viral population generation using `MLE_textfile_output`

**Command:**
`perl extract_MLE.pl -f paths_fasta_file -l MLE_textfile_output  > MLE_Haplo_fasta_OUTPUT`

Example
```
perl extract_MLE.pl -f paired-reads.60.fact15.fasta -l paired-reads.60.smxlik.txt > paired-reads.60.MLE.fasta
```

- `paths_fasta_file`. The paths generated by [VIPRA Step B: Generate fasta file](#step-b-generate-fasta-file)

## ~~ update for the ViPRA-Haplo paper (submitted in 2022) ~~
## Error Correction
Install Karect

karect -correct -inputfile=*_1.fa -inputfile=*_2.fa -celltype=haploid -matchtype=hamming -aggressive=5.0 -numstages=2 -errorrate=0.25 -errorratesec=0.25 


## VSEARCH clustering
centroid-based clustering by similarity. A centroid in a cluster is a contig.

**Command:**
`vsearch --cluster_fast paths_fasta_file --id --centroids cluster_centroid_sequences --consout cluster_consensus_sequences`

Example
```
vsearch --cluster_fast paired-reads.60.fact15.fasta --id 0.995 --centroids paired-reads.60.fact15.centroids.0.995.fa --consout paired-reads.60.fact15.consout.0.995.fa
```

## update MLEHaplo: 
Input: paths reconstructed by ViPRA, or paths after reduced by VSEARCH when taking paths reduced by VSEARCH as input.

**Command:**
`perl dg_subpath.pl paths_write_file cluster_centroid_sequences new_paths_file`

Example
```
perl dg_subpath.pl paired-reads.60.fact15.paths.txt  paired-reads.60.fact15.centroids.0.995.fa paired-reads.60.fact15.0.995.paths.txt 
```
use this new paths file for likelihood_singles_wrapper.pl
