# use strict;
# use warnings;
use List::Util qw(shuffle min );
use Getopt::Long;
use Parallel::Loops;
#require "/i3c/alisa/rom5161/haplotype/MD_virus/likelihood_singles_functions.pl";

my($cond_graph) = "";
my($comp_sets) = "";
my($paths_file) = "";
my($truefile) = "";
my($iter) = 100000; 
my($mode) = '';  
my($back) = ''; #backward elimination mode
my($slow) = '';
my($kmer_size) = 60;
my($genome_length) = 0; 
GetOptions( "condgraph=s" => \$cond_graph,
	    "compset=s" => \$comp_sets,
	    "pathsfile:s" => \$paths_file,
	    "trueset:s" => \$truefile,
	     "iter:i" => \$iter,
	     "kmer:i" => \$kmer_size,
		 "gl=i" =>\$genome_length,
	    "random" => \$mode, "back" => \$back,
		"slow" =>\$slow)  
	or die("Error in input data \n");

load_condensed_graph($cond_graph);
$count = load_paired_compatibles($comp_sets);

pair_complement_nodes();
load_all_paths() if $paths_file ne "";  

$max_processes = 8;
$pl = Parallel::Loops->new($max_processes);
	

%set_paths = ();
#print_compatible_set();
#print "$count loaded ".scalar(keys %cond_graphout)." nodes loaded ".scalar(keys %compatible_set)." compatible set loaded \n";
if($mode)
{	
	$wrfile =~ s/trueset/like_rand.txt/;
	open(wr,">$wrfile"); 
	likelihood_subset_random(14,$iter); 
}elsif($back)
{
	#print "Estimating backward elimination \n";
#	print "number of possible haplotypes is ".scalar(keys %allpaths)."\n";
#	my(@start_set) = keys %allpaths;
#	print "Likelihood all ".likelihood_current_set(@start_set)."\n";
	likelihood_subset_backward_elimination();
}else
{
	$wrfile =~ s/trueset/like_replace.txt/;
	open(wr,">$wrfile"); 
	likelihood_subset(14,$iter)
}

sub likelihood_subset_backward_elimination
{
	# Start by initializing the set of all paths as the start set for analysis 	
	my($current_like) = initialize_hashes();
	#print "$current_like is starting likelihood\n";
	$remove_haps = backward_elimination_fixed_size();
	$old_like = $temp_like;
	while($temp_like != -9**9**9 && scalar(keys %set_paths)>0) 
	{
		$old_like = $temp_like;
		# Remove the haplotype for which remaining set has highest likelihood
		remove_one_d_hashtable($remove_haps);
		#($errorhigh,$errorlow) = likelihood_error_bar(@start_set);
		#print scalar(keys %set_paths)." $temp_like $errorlow $errorhigh $remove_haps $path_pair{$remove_haps}\n";
		print scalar(keys %set_paths)." $temp_like $remove_haps $percent_missed\n";
		$remove_haps = backward_elimination_fixed_size();
	}
#	print "printing start set \n";
#	foreach (sort @start_set) { print "$_ "; } print "\n";
	#print_haplotypes_final(@start_set);
#	foreach (sort @trueset) { print "$_ "; } print "\n";
#	print "$current_like with ".$#current_haps_." haplotypes and $temp_like with best minus one. Removed one is $current_haps_[$imp_var] \n";
}

sub backward_elimination_fixed_size
{
	my($current_like ) = compute_set_likelihood_using_d();
	#print $current_like."\n";
	my($new_like) = 0; 
	my(@current_haps_) = keys %set_paths; 
	# foreach (sort {$a <=> $b } @current_haps_) { print "$_ ";}
	# print "\n";
	my($imp_var) = -10;
	my(%likeli_hash) = ();
	my(%likeli_pcmissed_hash) = ();
	
	#my(%temp_set_reduced) = map{ $_ => 0} @current_haps_;
	$temp_like  = -9**9**9;
	$size = scalar(@current_haps_);
	# print $size." is the size\n";
	$pl->share(\%likeli_hash,\%likeli_pcmissed_hash);
	$pl->foreach(\@current_haps_, sub {
		# Remove an haplotype from the current set and compute likelihood of the remaining set 
		my($k) = $current_haps_[$_];
		remove_one_d_hashtable($k);
		# Change the likelihood function here
		if($slow)
		{
			$likeli_hash{$k} = compute_set_likelihood_using_d();
			$likeli_pcmissed_hash{$k} = $percent_pairs_missed;
			#print "$k $likeli_hash{$k} \n";
		}else
		{
			$likeli_hash{$k} = compute_set_likelihood_using_d_rem($current_like);
			$likeli_pcmissed_hash{$k} = $percent_pairs_missed;
		}
		# print "$new_like when removing $k\n";
		# if($new_like !=0 && $new_like > $temp_like) 
		# {
			# $temp_like = $new_like; $imp_var = $k;
			# $percent_missed = $percent_pairs_missed;
		# }
		# Add the removed haplotype again from the current set to restore the original set 
		add_one_d_hashtable($k);
	});
	# undef(%temp_set_reduced);
	my( @vals )  = sort { $likeli_hash{$b} <=> $likeli_hash{$a} } keys %likeli_hash;
	$temp_like = $likeli_hash{$vals[0]};
	$percent_missed = $likeli_pcmissed_hash{$vals[0]};
	# foreach (sort {$a <=> $b } keys %likeli_hash) { print "$_ ";}
	# print "\n";
	return $vals[0];
}
sub compute_d_hashtable_initial
{
	for $n1( keys %set_paths ) 
	{
		#pick a path in the haplotype set
		my(@one_path) = @{$set_paths{$n1}};
		#print "Operating on $n1 ".scalar(@one_path)." nodes in it \n";
		#my(%temp_set) = ();
		for (my($i)=0;$i<scalar(@one_path);$i++)
		{
			for (my($j)=$i; $j < scalar(@one_path); $j++)
			#for $k(keys %{$compatible_set{$one_path[$i]}})
			{
				if(exists($compatible_set{$one_path[$i]}{$one_path[$j]})) 
				{
					$d_hashtable{$one_path[$i]}{$one_path[$j]}++;
				}
			}
		}
	}
}

sub remove_one_d_hashtable
{
	# update d_hashtable by removing a path provided
	my($remove_path) = shift;
	
	# update set_paths
	#print "Removing path $remove_path\n";
	my(@one_path) = @{$set_paths{$remove_path}};
	%rem_newhashtable = ();
	%rem_oldhashtable = ();
	delete($set_paths{$remove_path});
	for (my($i)=0; $i<scalar(@one_path); $i++)
	{
		for(my($j)=$i; $j<scalar(@one_path); $j++)
		{
			if(exists( $compatible_set{$one_path[$i]}{$one_path[$j]} ) )
			{
				$rem_oldhashtable{$one_path[$i]}{$one_path[$j]} = $d_hashtable{$one_path[$i]}{$one_path[$j]};
				# print "$one_path[$i]:$one_path[$j]:$rem_oldhashtable{$one_path[$i]}{$one_path[$j]} ";
				$d_hashtable{$one_path[$i]}{$one_path[$j]} --;
				# IF d_hashtable goes to zero then just remove it
				$rem_newhashtable{$one_path[$i]}{$one_path[$j]} = $d_hashtable{$one_path[$i]}{$one_path[$j]};
				if($d_hashtable{$one_path[$i]}{$one_path[$j]}==0)
				{
					# print "dd:$one_path[$i]:$one_path[$j]:$d_hashtable{$one_path[$i]}{$one_path[$j]} ";
					delete $d_hashtable{$one_path[$i]}{$one_path[$j]};
					delete $rem_newhashtable{$one_path[$i]}{$one_path[$j]};

				}
			}
		}
	}
	# print "rem_old_new\n";
#	print "Paths remaining remove path ".scalar(keys %set_paths);
}

sub add_one_d_hashtable 
{
	my($add_path) = shift;
	my(@one_path) = @{$allpaths{$add_path}};
	# Add the path and its nodes to the current set of paths 
	@{$set_paths{$add_path}} = @one_path;
	%add_newhashtable = ();
	%add_oldhashtable = ();
	for (my($i)=0; $i<scalar(@one_path); $i++)
	{
		for(my($j)=$i; $j<scalar(@one_path); $j++)
		{
			if(exists( $compatible_set{$one_path[$i]}{$one_path[$j]} ) )
			{
				if(exists $d_hashtable{$one_path[$i]}{$one_path[$j]}) 
				{
					$add_oldhashtable{$one_path[$i]}{$one_path[$j]} = $d_hashtable{$one_path[$i]}{$one_path[$j]};
				}
				$d_hashtable{$one_path[$i]}{$one_path[$j]} ++;
				$add_newhashtable{$one_path[$i]}{$one_path[$j]} = $d_hashtable{$one_path[$i]}{$one_path[$j]};
				# print "$one_path[$i]:$one_path[$j]:$add_newhashtable{$one_path[$i]}{$one_path[$j]} ";
			}
		}
	}
	# print "add_old_new\n";
#	print "Paths remaining add path ".scalar(keys %set_paths);
}

sub compute_set_likelihood_using_d_rem
{
	my($llik) = shift;
	my($start_like) =$llik;
	$num_haps = scalar(keys %set_paths);
	$Scale = $num_haps*$genome_length;
	my(%visit_fwd_rev) = ();
	my($remove_add_count) = 0;
	for $k1 ( keys %rem_newhashtable) 
	{
		for $k2( keys %{$rem_newhashtable{$k1}} ) 
		{
			 if(exists($compatible_set{$k1}{$k2})) 
			 {	$visit_fwd_rev{$k1}{$k2} = -1; }
		}
	}
	# for $k1 ( keys %rem_oldhashtable) 
	# {
		# for $k2( keys %{$rem_oldhashtable{$k1}} ) 
		# {
			 # if(exists($compatible_set{$k1}{$k2})) 
			 # {	$visit_fwd_rev{$k1}{$k2} = -1; }
		# }
	# }	
	# Update $llik start with rem_newhashtable
	for $k1( keys %compatible_set)
	{
		for $k2( keys %{$compatible_set{$k1}})
		{
			if(exists($d_hashtable{$k1}{$k2}))
			{
				if($visit_fwd_rev{$k1}{$k2}==-1 )
				{
					$qij = $rem_newhashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
					$temp1 = $total_count - $compatible_set{$k1}{$k2};
					$temp2 = 1 - $qij;
					#print "$k1 $k2 rem_newhashtable $rem_newhashtable{$k1}{$k2} PS $compatible_set{$k1}{$k2} $temp2 $qij\n";
					$val = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
					#print "$val \n";
					$llik += $val;
					# Substract the old hash part
					$qij = $rem_oldhashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
					$temp1 = $total_count - $compatible_set{$k1}{$k2};
					$temp2 = 1 - $qij; 
					#print "$k1 $k2 rem_oldhashtable $rem_oldhashtable{$k1}{$k2} PS $compatible_set{$k1}{$k2} $temp2 $qij\n";
					$val  = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
					$llik -=$val;
					$visit_fwd_rev{$k1}{$k2} = 1;
					#$visit_fwd_rev{$paired_nodes{$k2}}{$paired_nodes{$k1}} = 1;
				}
			}
			else
			{
				if( exists($rem_oldhashtable{$k1}{$k2}) )
				{
					$qij = $rem_oldhashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
					$temp1 = $total_count - $compatible_set{$k1}{$k2};
					$temp2 = 1 - $qij;
					#print "$k1 $k2 singles rem_oldhashtable $rem_oldhashtable{$k1}{$k2} PS $compatible_set{$k1}{$k2} $temp2 $qij\n";
					$val = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
					#print "$val \n";
					$llik -= $val;
					$visit_fwd_rev{$k1}{$k2} = 1;
					$visit_fwd_rev{$paired_nodes{$k2}}{$paired_nodes{$k1}} = 1;
					$llik +=$epsilon;
					# print "rem$k1:$k2:$length_diff{$k1}{$k2}\t";
					$num_pairs_missed += $compatible_set{$k1}{$k2};
				}
			}
		}
	}
	#print "Removing $remove_add_count pairs \n";	
	# for $k1(keys %rem_oldhashtable)
	# {
		# for $k2(keys %rem_oldhashtable)
		# {
			# if(exists($visit_fwd_rev{$k1}{$k2}) && $visit_fwd_rev{$k1}{$k2}==0 && exists($compatible_set{$k1}{$k2}))
			# {
				# $qij = $rem_oldhashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
				# $temp1 = $total_count - $compatible_set{$k1}{$k2};
				# $temp2 = 1 - $qij;
				# #print "$k1 $k2 singles rem_oldhashtable $rem_oldhashtable{$k1}{$k2} PS $compatible_set{$k1}{$k2} $temp2 $qij\n";
				# $val = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
				# #print "$val \n";
				# $llik -= $val;
				# #if(exists($compatible_set{$k1}{$k2}))
				
				
				
				# $visit_fwd_rev{$k1}{$k2} = 1;
				# $visit_fwd_rev{$paired_nodes{$k2}}{$paired_nodes{$k1}} = 1;
			# }
		# }
	# }
	#print "Starting like $start_like End like $llik\n";
	$percent_pairs_missed = ( $num_pairs_missed/$total_pairs );
	return $llik;
}

sub compute_set_likelihood_using_d
{
	# Compute likelihood using the d_hashtable 
	# d_hashtable contains information about the number of times a vertex pair occurs in current %set_paths
	# num_haps is the input from 
	$num_haps = scalar(keys %set_paths); 
	my($llik) = 0; 
	my($min) = 0; 
	my($max) = -9**9**9; 
	$num_pairs_missed = 0;
	$total_pairs = 0;
	# Penalty term when a vertex pair in the reads is not present in the current %set_paths
	$epsilon = -1e10;
	$Scale = $num_haps*$genome_length;
	# Make sure that one -forward and reverse node pair is only visited once. 
	my(%visit_fwd_rev) = ();
	for $k1 ( keys %compatible_set) 
	{
		for $k2( keys %{$compatible_set{$k1}} ) 
		{
			$total_pairs += $compatible_set{$k1}{$k2}; 
			$visit_fwd_rev{$k1}{$k2} = 0;
		}
	}
	
	#print "Genome length $genome_length\n";
	for $k1( keys %compatible_set ) 
	{
		for $k2(keys %{$compatible_set{$k1}}) 
		{
			if(exists($d_hashtable{$k1}{$k2}) )
			{
				if($visit_fwd_rev{$k1}{$k2}==0) 
				{
					$qij = $d_hashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
					$temp1 = $total_count - $compatible_set{$k1}{$k2};
					$temp2 = 1 - $qij;
					#print "$k1 $k2 d_hashtable $d_hashtable{$k1}{$k2} PS $compatible_set{$k1}{$k2} $temp2 $qij\n";
					$val = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
					#print "$val \n";
					$llik += $val;
					$visit_fwd_rev{$k1}{$k2} = 1;
					#$visit_fwd_rev{$paired_nodes{$k2}}{$paired_nodes{$k1}} = 1;
				}
			}
			else 
			{
				#compatible set contains a haplotype which is not supported by the set of haplotypes given, reject the set of haplotypes 	
				$llik +=$epsilon;
				$num_pairs_missed += $compatible_set{$k1}{$k2};
				# print "beg$k1:$k2:$length_diff{$k1}{$k2}\t";
				#return $min;
			}
		}
	}
	# print "\n";
	#print "Max $max and Min $min $t1 $t2 $compatible_set{$t1}{$t2}\n";
	#print "Likelihood is $llik \n";
	$percent_pairs_missed = ( $num_pairs_missed/$total_pairs );
	#print "Pairs missed $percent_pairs_missed\n";
	return $llik;
}

sub initialize_hashes
{
	# Initialize set_paths, d_hashtable and other hashes 
	# Take as the initial set of paths
	%set_paths = (); 
	foreach (keys %allpaths) {
		push @{$set_paths{$_}} , @{$allpaths{$_}};
	}
	# Initialize d_hashtable
	%d_hashtable = ();
	compute_d_hashtable_initial();
	$initial_like = compute_set_likelihood_using_d();
	return $initial_like;
}



sub load_condensed_graph
{
	my($cond_graph) = $_[0]; # Condensed graph file
	open(file,$cond_graph);
	%cond_ver = ();
	%cond_graphout = ();
	while($l=<file>)
	{
		chomp($l);
		($v1, $v2) = split(/ /,$l);
		if( $v2 =~ /^[A,G,C,T]+$/) 
		{
			$cond_ver{$v1} = $v2;
			$length_ver{$v1} = length($v2);
		}else
		{
			$cond_graphout{$v1}{$v2}=1;
		}
	}
	close file;
}

sub load_paired_compatibles
{
	# Load the node pairs that are mapped in the compatiblle sets ... 
	# The file contains the following information.
	# x1:Depth of x1 compatiblevertex:Depth:SupportCount
	my($pairedfile) = $_[0];
	open(paired_compatible_file,"$pairedfile");
	%compatible_set = ();
	%length_diff = ();
	$total_count = 0;
	while($l=<paired_compatible_file>)
	{
		chomp($l);
		@arr = split(/ /,$l);
		($v1,$v1_depth) = split(/:/,$arr[0]);
		for(my($i)=1;$i<scalar(@arr); $i++) 
		{
			($v2,$v2_depth,$count) = split(/:/,$arr[$i]);
			
			if ($v2_depth >= $v1_depth) 
			{
				$total_count +=$count;
				$compatible_set{$v1}{$v2} = $count;
				$length_diff{$v1}{$v2} = $v2_depth - $v1_depth + $length_ver{$v2} ;
			}
		}
	}
	close paired_compatible_file;
	return $total_count;
}

sub pair_complement_nodes
{
	# take all the nodes and compute their reverse complements, and a paired hash table be created. 
	my(%temp_hash) = ();
	my(%inverse_hash) = ();
	for $k( keys %cond_ver) {
		$temp_hash{$k} = 0;
		$inverse_hash{$cond_ver{$k}} = $k;
	} 
	#Link nodes which are reverse complements of each other. 
	for $k(keys %cond_ver)
	{
		if($temp_hash{$k}==0) {
			my($temp) = revcomplement($cond_ver{$k});
			$paired_nodes{$k} = $inverse_hash{$temp};
			$paired_nodes{$inverse_hash{$temp}} = $k ; 
			$temp_hash{$k} = 1; 
			$temp_hash{$inverse_hash{$temp}} = 1;
		}
	}
	undef %inverse_hash;
	undef %temp_hash;
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
		else {$rev_str = $rev_str."C";}
	}
	return $rev_str;
}

sub load_all_paths
{
	#Loading only paired (reverse and forward strain) paths as they were artificially created in the alogrithm to make sure we don't miss information from a single read 
	%allpaths = ();
	open(file,$paths_file);
	while($l=<file>)
	{
		chomp($l);
		@arr = split(/ /,$l);
		$l = <file>;
		chomp($l);
		($pathno,$t,$length) = split(/_/,substr($l,1));
			for(my($i)=0;$i<scalar(@arr);$i=$i+2){ 
				my($n,$node) = split(/:/,$arr[$i]);
				push @{$allpaths{$pathno}}, $node;
			}
	}
	close file;
	#print "Number of paired paths is ".scalar(keys %allpaths)."\n";
}

sub compute_set_likelihood_using_d_add
{
	my($llik) = shift;
	$num_haps = scalar(keys %set_paths);
	$Scale = $num_haps*$genome_length;
	my(%visit_fwd_rev) = ();
	# only visit things 
	for $k1 ( keys %add_newhashtable) 
	{
		for $k2( keys %{$add_newhashtable{$k1}} ) 
		{
			$visit_fwd_rev{$k1}{$k2} = 0;
		}
	}
	for $k1 ( keys %add_oldhashtable) 
	{
		for $k2( keys %{$add_oldhashtable{$k1}} ) 
		{
			if(!exists($visit_fwd_rev{$k1}{$k2})) { $visit_fwd_rev{$k1}{$k2} = 0; }
		}
	}
	# Update llik value 
	for $k1( keys %add_oldhashtable ) 
	{
		for $k2(keys %{$add_oldhashtable{$k1}}) 
		{
			#if(exists(${$k1}{$k2}) )
			if($visit_fwd_rev{$k1}{$k2}==0) 
			{
				$qij = $add_newhashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
				$temp1 = $total_count - $compatible_set{$k1}{$k2};
				$temp2 = 1 - $qij;
				#print "$k1 $k2 d_hashtable $d_hashtable{$k1}{$k2} PS $compatible_set{$k1}{$k2} $temp2 $qij\n";
				$val = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
				#print "$val \n";
				$llik += $val;
		# Substract the old hash part
				$qij = $add_oldhashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
				$temp1 = $total_count - $compatible_set{$k1}{$k2};
				$temp2 = 1 - $qij; 
				$val  = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
				$llik -=$val;
			
				$visit_fwd_rev{$k1}{$k2} = 1;
				$visit_fwd_rev{$paired_nodes{$k2}}{$paired_nodes{$k1}} = 1;
			#compatible set contains a haplotype which is not supported by the set of haplotypes given, reject the set of haplotypes 	
			#return $min;
			}
		}
	}
	#Something that was not in add_old but became positive in add_new 
	#Action : Add val and substract epsilon if present in compatible_set
	for $k1(keys %add_newhashtable)
	{
		for $k2(keys %add_newhashtable)
		{
			if($visit_fwd_rev{$k1}{$k2}==0)
			{
				$qij = $add_newhashtable{$k1}{$k2}/($Scale-$num_haps*$length_diff{$k1}{$k2}+$num_haps*2);
				$temp1 = $total_count - $compatible_set{$k1}{$k2};
				$temp2 = 1 - $qij;
				#print "$k1 $k2 d_hashtable $d_hashtable{$k1}{$k2} PS $compatible_set{$k1}{$k2} $temp2 $qij\n";
				$val = $temp1*log($temp2) + $compatible_set{$k1}{$k2} * log($qij);
				#print "$val \n";
				$llik += $val;
				if(exists($compatible_set{$k1}{$k2}))
				{
					$llik -=$epsilon;
					print "p$k1:$k2:$length_diff{$k1}{$k2}\t";
					$num_pairs_missed -= $compatible_set{$k1}{$k2};					
				}
				$visit_fwd_rev{$k1}{$k2} = 1;
				$visit_fwd_rev{$paired_nodes{$k2}}{$paired_nodes{$k1}} = 1;
			}
		}
	}
	return $llik;
}
