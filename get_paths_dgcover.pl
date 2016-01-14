# Process dg_cover file to get paths 
use Getopt::Long;

$filename = "";
$writefile = "";
$length_threshold = 0;
GetOptions("f=s"=>\$filename,"w:s"=>\$writefile,"t:s",\$length_threshold) or die ("Options not given");
if ($writefile eq "")
{
	$writefile = $filename;
	$writefile =~ s/fact.*fasta/paths.txt/;
}

main();

sub main
{
	open(file,$filename);
	open(wrfile,">$writefile");
	$l=<file>;
	while(substr($l,0,2) ne "1:")
	{
		$l=<file>;
	}
	chomp ($l);
	$path = $l;
	while($l=<file>)
	{
		chomp($l);
		if(substr($l,0,1) eq ">")
		{
			@vs = split(/__/,$l);
			if($vs[1] > $length_threshold)
			{
				print wrfile $path."\n";
				print wrfile "$l\n";
			}
		}
		if(substr($l,0,2) eq "1:")
		{
			$path = $l;
			#print wrfile "$l\n";
		}
	}
	close file;
	close wrfile;
}