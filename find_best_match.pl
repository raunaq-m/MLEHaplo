# process the sam file dump to find the best predicted haplotype for a true haplotype
# perl find_best_match.pl <sam file>
open(file,$ARGV[0]);

main();

sub main
{
	%identity_all = ();
	%identity_pos = ();
	%identity_len = ();
	analyze_samfile();
	# Print what is the NM to true haplotypes
	# print_identity_to_true();
	# Count the haplotypes that occur at distance "x" from true haplotypes
	# count_haps_matched_at_distance();
	# Print haps that map to each reference, sorted by either mapping position or by minimum edit distance 
	print_all_haps_mapped_to_true();
	# Print average distance of the closest distance haplotype to a reference ? 
	# print_average_distances();
	
}
sub analyze_samfile
{
	while($l=<file>)
	{
		chomp($l);
		if(substr($l,0,1) ne "@")
		{
			@vs = split(/\t/,$l);
			@msmatch = split(/:/,$vs[11]);
			@tlen = split(/__/,$vs[0]);
			# store the best match for a reference sequence in identity hash,
			# the length of the best match in identity_length
			# the name of the haplotype in id hash
			if(exists($identity{$vs[2]})) 
			{	
				# if($identity{$vs[2]}> $msmatch[2] || $identity_length{$vs[2]} < $tlen[1])
				if($identity{$vs[2]}> $msmatch[2])
				# if ($identity_length{$vs[2]} < $tlen[1] && $identity{$vs[2]} > $msmatch[2])
				{
					$identity{$vs[2]} = $msmatch[2];
					$identity_length{$vs[2]} = $tlen[1];
					$id{$vs[2]} = $vs[0];
				}
				
			}
			else 
			{
				$identity{$vs[2]} = $msmatch[2];
				$identity_length{$vs[2]} = $tlen[1];
				$id{$vs[2]} = $vs[0];
			}
			# Store the edit distance of every hapltype $vs[0] to the reference that it maps to $vs[2] in identity_all
			# Store the position of mapping in the hash identity_pos
			$identity_all{$vs[2]}{$vs[0]} = $msmatch[2];
			$identity_pos{$vs[2]}{$vs[0]} = $vs[3];
			# Length of the haplotype 
			$identity_len{$vs[0]} = $tlen[1];
			# path_identity_score stores the edit distance of haplotype that is at max distance? to a reference sequence
			# path_identity stores the name of the reference which is at path_identity_score edit distance
			if(exists($path_identity_score{$vs[0]}))
			{
				if($path_identity_score{$vs[0]} > $msmatch[2])
				{
					$path_identity_score{$vs[0]}  = $msmatch[2];
					$path_identity{$vs[0]} = $vs[2];
				}
			}
			else
			{
				$path_identity_score{$vs[0]}  = $msmatch[2];
				$path_identity{$vs[0]} = $vs[2];
			}
		}
	}
	close file;
}
# Print distances of true haplotypes to the predicted haplotypes

sub print_identity_to_true
{
	for (keys %identity)
	{
		print "$_ $identity{$_} $identity_length{$_} $id{$_}\n";
	}
}

sub print_all_haps_mapped_to_true
{

	for $k(keys %identity_all)
	{
		print $k."\n";
		# for $y(sort {$a <=> $b} keys %{$identity_all{$k}})
		for $y (sort {$identity_all{$k}{$a} <=> $identity_all{$k}{$b}} keys %{$identity_all{$k}})
		{
			#if($identity_all{$k}{$y} <=10)
			{
				print "$y,$identity_all{$k}{$y},$identity_pos{$k}{$y}\n";
			}
		}
		print "\n";
	}
}

sub print_average_distances
{
	# Print distances of the predicted haplotypes to the true haplotypes
	$count =1 ; 
	$distance = 1;
	$temp="";
	for (sort {$path_identity{$a} cmp $path_identity{$b}} keys %path_identity_score)
	{
		if($temp eq $path_identity{$_}) 
		{
			#print "$path_identity_score{$_}\t";
			$distance += $path_identity_score{$_};
			$count++;
		}
		else
		{
			$last_average = $distance/$count;
			if ($temp ne "")
			{
				print "$last_average\t$count\n$path_identity{$_}\t";
			}
			else 
			{
				print "$path_identity{$_}\t";
			}
			$distance = $path_identity_score{$_};
			$count = 1;
			$temp = $path_identity{$_};
		}
	}
	print $distance/$count;
	print "\t$count\n";

}


sub count_haps_matched_at_distance
{
	# purpose: Count haplotypes that are mapped to all references at distance "x"
	# Use
	# $identity_all{$vs[2]}{$vs[0]}	 $identity_pos{$vs[2]}{$vs[0]}
	# Define hash distance 
	%distance = ();
	for $ref(sort keys %identity_all)
	{
		for $haps(keys %{$identity_all{$ref}})
		{
			# print "$ref $haps $identity_all{$ref}{$haps}\n";
			$distance{$identity_all{$ref}{$haps}}++;
			$length{$identity_len{$haps}}++;
		}
	}
	for $k (sort {$a <=>$b} keys %distance)
	{
		print "$k\t$distance{$k}\n";
	}
	print "####\n";
	for $k(sort {$a <=>$b} keys %length)
	{
		print "$k\t$length{$k}\n";
	}
}
