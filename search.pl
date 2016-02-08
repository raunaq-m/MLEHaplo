open(file,$ARGV[0]); #true haplotype file
%truehash = ();
while($l=<file>)
{
	chomp($l);
	$id = substr($l,1);
	$l=<file>;
	chomp($l);
	$truehash{$id} = $l;
	push @{$unique{$l}} , $id; 
	$found_hash{$id} = 0;
	$found_unique{$l} = 0;
	# also look for reverse complements
	$rev_string = revcomplement($l);
	push @{$unique{$rev_string}} , $id;
}
close file;

open(file,$ARGV[1]); #predicted file 
$count = 0;$total_haps = 0;
while($l=<file>)
{
	chomp($l);
	if( $l =~ /[A,G,C,T][A,G,C,T]/)
	{
		$found=0; @id = ();
		#foreach $k( keys %truehash )
		foreach $k(keys %unique)
		{
			if ($k =~ /$l/)
			{ 
				$found=1; 
				push @id, @{$unique{$k}};
			 	$seq = $k; 
			}
		}
		#if($found ==1 && $found_unique{$seq} ==0 ) { print "$id found in $store\n" ; $count++; $found_hash{$id} = 1; $found_unique{$seq} = 1; } 
		if($found == 1) {
			$count += scalar(@id);
			foreach (@id) { 
				print "$_ ";
				$found_hash{$_} = 1; 
			}
			print "found in $store\n";
		}
	}
	if( substr($l,0,1) eq ">" ) {  $store = $l; $total_haps ++ ; }
}
close file;
print "Count is $count  total is $total_haps\n";

#foreach $k(keys %found_hash)
#{
#	if($found_hash{$k} ==0)
#	{
#		print "$k key not found \n";
#	}
#}
$total_not_found=0;
for $k( keys %unique)
{
	
	$checker=  0;
	foreach (@{$unique{$k}})
	{
		if($found_hash{$_}==0)
		{	
			$checker = 1;
			print "$_ ";
		}
	}
	if($checker ==1) 
	{
		$total_not_found++;
		print "key not found \n";
	}
}

print "$total_not_found keys not found\n";

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
