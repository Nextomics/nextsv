#!/usr/bin/perl
use strict;

my $usage  = "Usage:   perl $0 <in.spots>\n";
$usage    .= "Contact: Li Fang (fangli\@grandomics.com)\n";
$usage    .= "Version: 0.4.0\n";

die $usage if (@ARGV < 1);

my $min_length = 50;     # min sv length for output
my $max_length = 1000000; # max sv length for output

my $in = shift(@ARGV);
my $out1 = "$in.INS.bed"; # insertion calls in bed format
my $out2 = "$in.DEL.bed"; # deletion calls in bed format

my %hash = ("1" => 1, "2" => 1, "3" => 1, "4" => 1, "5" => 1, "6" => 1, "7" => 1, "8" => 1, "9" => 1, "10" => 1, "11" => 1, "12" => 1, "13" => 1, "14" => 1, "15" => 1, "16" => 1, "17" => 1, "18" => 1, "19" => 1, "20" => 1, "21" => 1, "22" => 1, "X" => 1, "Y" => 1, "MT" => 1, "chr1" => 1, "chr2" => 1, "chr3" => 1, "chr4" => 1, "chr5" => 1, "chr6" => 1, "chr7" => 1, "chr8" => 1, "chr9" => 1, "chr10" => 1, "chr11" => 1, "chr12" => 1, "chr13" => 1, "chr14" => 1, "chr15" => 1, "chr16" => 1, "chr17" => 1, "chr18" => 1, "chr19" => 1, "chr20" => 1, "chr21" => 1, "chr22" => 1, "chrX" => 1, "chrY" => 1, "chrMT" => 1); 

my $caller = "PBhoney-Spots";

open (IN, $in) or die $!;
open (OUT1, "> $out1") or die $!;
open (OUT2, "> $out2") or die $!;

my $header = "#Chr\tStart\tEnd\tRef\tAlt\tType\tLength\tSV_caller\tInfo\n";
print OUT1 $header;
print OUT2 $header;

print STDERR "Formating $in ... Settings: min sv length = $min_length, max sv length = $max_length\n";

##### read spots file #####
seek IN, 0, 0;
while (my $line = <IN>)
{
	next if ($line =~ /^#/);
	chomp $line;
	(my $chr, my $start, my $end, my $type, my $length, my $info) = split("\t", $line);
	if ($type eq "DEL"){
		$length = $end - $start;
	}
	
	$start = $start - 1;  # The first base in a chromosome is numbered 0 in bed format
	my $index = index($info, "szCount"); # number of reads support
	my $szcount = substr($info, $index+length("szCount="));
	my @b = split(";", $szcount);
	my $rd_sup = $b[0];
	$rd_sup = int($rd_sup);

	my $index = index($info, "coverage");  # depth at the sv site
	my $coverage = substr($info, $index+length("coverage="));
	@b = split(";", $coverage);
	my $depth = $b[0];
	$depth = int($depth);

###### filters ######

	next if (!exists ($hash{$chr}));                   
	next if ($length > $max_length );                  # larger than max sv length, probably false positives
	next if ($length < $min_length );                  # less than min sv length
	next if ($end - $start > 1000 and $type eq "INS"); # the breakpoint of insertion call is fuzzy
	next if ($rd_sup * 5 < $depth);                    # variant frequency < 20%, probably false positives 

##### output results #####
# the  output format can be annotated by ANNOVAR
	if ($type eq "INS") {
		print OUT1 "$chr\t$start\t$end\t0\t0\tcomments:sv_type=$type;sv_length=$length;sv_caller=$caller;coverage=$depth;numReads=$rd_sup;$info\n";
	}elsif ($type eq "DEL") {
		print OUT2 "$chr\t$start\t$end\t0\t0\tcomments:sv_type=$type;sv_length=$length;sv_caller=$caller;coverage=$depth;numReads=$rd_sup;$info\n";
	}
}
close IN;
close OUT1;
close OUT2;
