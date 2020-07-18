open(file,$ARGV[0]); #path
open(file1,$ARGV[1]); #filtered fa
open(file2,">$ARGV[2]"); #filtered path

my %paths;

while($l=<file>)
{
	chomp($l);
	$l =~ s/\v//l;
	if(substr($l,0,1) eq ">"){
		if(exists($paths{$l})){
			print "hash error!\n";
		}else{
			$paths{$l}=$p;
		}
	}else{
		$p=$l;
	}
}
close file;

#foreach $key (keys %paths)
#{
#  print "$key\n";
#}

while($l1=<file1>)
{
	chomp($l1);
	if(substr($l1,0,1) eq ">"){
		@array1 = split /\s/,$l1;
		if(!exists($paths{$array1[0]})){
			print "codes error!\n";
		}else{
			$path1=$paths{$array1[0]};
		}
		print file2 "$path1\n$array1[0]\n";
	}
}
close file1;
close file2;