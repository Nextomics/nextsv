#!/usr/bin/perl

use strict;
my $usage  = "Usage:   perl $0 <in1.bed> <in2.bed> <out.intersect.bed> <out.union.bed> <sv_type>\n";
$usage    .= "Version: 0.4.0\n";
$usage    .= "Contact: Li Fang (fangli\@grandomics.com)\n";

die $usage if (@ARGV < 2);

my $in1  = shift(@ARGV); 
my $in2  = shift(@ARGV);
my $out1 = shift(@ARGV); # output file name for intersect SV calls
my $out2 = shift(@ARGV); # output file name for union SV calls
my $svtype = shift(@ARGV);

open (IN1, $in1);
open (IN2, $in2);
open (OUT1, "> $out1");
open (OUT2, "> $out2");

my $i = 0;
my $line;
my @chr1, my @start1, my @end1, my @info1, my @merged1;
my @chr2, my @start2, my @end2, my @info2, my @merged2;
my $ref, my $alt;
my $max_dist = 1000; 

while ($line = <IN1>)
{
	next if ($line =~ /^#/);
	chomp $line;
	($chr1[$i], $start1[$i], $end1[$i], $ref, $alt, $info1[$i]) = split ("\t", $line);
	$info1[$i] = substr($info1[$i], length("comments:"));
	$merged1[$i] = 0;
	$i++;
}
my $n_sv1 = $i;

my $j = 0;
while ($line = <IN2>)
{
	next if ($line =~ /^#/);
	chomp $line;
	($chr2[$j], $start2[$j], $end2[$j], $ref, $alt, $info2[$j]) = split ("\t", $line);
	$info2[$j] = substr($info2[$j], length("comments:"));
	$merged2[$j] = 0;
	$j++;
}
my $n_sv2 = $j;

my $start_max, my $end_min, my $length_intersect;
my $start_merge, my $end_merge, my $length_merge;
for ($i = 0; $i < $n_sv1; $i++)
{
	for ($j = 0; $j < $n_sv2; $j++)
	{
		next if ($merged1[$i]);
		next if ($merged2[$j]);
		if ($chr1[$i] eq $chr2[$j])
		{
			$start_max = $start1[$i] > $start2[$j] ? $start1[$i] : $start2[$j];
			$end_min  = $end1[$i] < $end2[$j] ? $end1[$i] : $end2[$j];
			$length_intersect = $end_min - $start_max;
			if ($svtype eq "INS" or $svtype eq "ins"){
				next if ($length_intersect < -$max_dist);
			}else{
				next if ($length_intersect <= 0);
				next if ($length_intersect < 0.5 * ($end1[$i] - $start1[$i]));
				next if ($length_intersect < 0.5 * ($end2[$j] - $start2[$j]));
			}
			$start_merge = int(($start1[$i] + $start2[$j])/2);
			$end_merge   = int(($end1[$i] + $end2[$j])/2);
			$length_merge = $end_merge - $start_merge;

			print OUT1 "$chr1[$i]\t$start_merge\t$end_merge\t0\t0\tcomments:merged from {$info1[$i]} and {$info2[$j]}\n";
			print OUT2 "$chr1[$i]\t$start_merge\t$end_merge\t0\t0\tcomments:merged from {$info1[$i]} and {$info2[$j]}\n";
			$merged1[$i] = 1;
			$merged2[$j] = 1;
		}
	}
}

for ($i = 0; $i < $n_sv1; $i++)
{
	next if ($merged1[$i]);
	print OUT2 "$chr1[$i]\t$start1[$i]\t$end1[$i]\tcomments:$info1[$i]\n";
}

for ($i = 0; $i < $n_sv2; $i++)
{
	next if ($merged2[$i]);
	print OUT2 "$chr2[$i]\t$start2[$i]\t$end2[$i]\tcomments:$info2[$i]\n";
}

close IN1;
close IN2;
close OUT1;
close OUT2;
