#!/usr/bin/perl

$usage  = "Usage:   perl $0 <fq.list> <num_file> <out_prefix>\n";
$usage .= "Author:  Li Fang (fangli\@grandomics.com)\n";
$usage .= "Version: 0.4.0\n";

die $usage if (@ARGV < 3);

$list = shift (@ARGV);
$num_file = shift(@ARGV);
$out_prefix = shift(@ARGV);
open (LIST, $list) or die $!;


$i = 0;
for ($fidx = 0; $fidx < $num_file; $fidx++) {
	$fq[$fidx] = "";
	$out = "$out_prefix.$fidx.fastq";
	open (OUT, "> $out") or die $!;
	print OUT $fq[$fidx];
	close OUT;
}

while ($file = <LIST>)
{
	chomp $file;
	open (FN, $file) or die $!;
	while ($line1 = <FN>)
	{
		$line2 = <FN>;
		$line3 = <FN>;
		$line4 = <FN>;
		$i++;
		$fidx = int(rand($num_file));
		$fq[$fidx] .= $line1.$line2.$line3.$line4;
		if ($i % 10000 == 0)
		{
			for ($fidx = 0; $fidx < $num_file; $fidx++) 
			{
				$out = "$out_prefix.$fidx.fastq";
				open (OUT, ">> $out") or die $!;
				print OUT $fq[$fidx];
				close OUT;
				$fq[$fidx] = "";
			}
		}
	}
}

for ($fidx = 0; $fidx < $num_file; $fidx++) 
{
	$out = "$out_prefix.$fidx.fastq";
	open (OUT, ">> $out") or die $!;
	print OUT $fq[$fidx];
	close OUT;
	$fq[$fidx] = "";
}

close LIST;
