#!/usr/bin/perl
use strict;

my $usage = "Usage:   perl $0 <in.sniffles.vcf> <out_prefix>\n";
$usage   .= "Contact: Li Fang (fangli\@grandomics.com)\n";
$usage   .= "Version: 0.4.0\n";
die $usage if (@ARGV < 2);

my $min_length = 50;        # min sv length for DEL, DUP and INV
my $max_length = 10000000;   # max sv length for DEL, DUP and INV
my $min_read_supp = 2;
my $in = shift(@ARGV);
my $out_prefix = shift(@ARGV);
my $out1 = "$out_prefix.INS.bed";
my $out2 = "$out_prefix.DEL.bed";
my $out3 = "$out_prefix.INV.bed";
my $out4 = "$out_prefix.TRA.bedpe";   # not bed format
my $out5 = "$out_prefix.DUP.bed";

open (IN, $in) or die $!;
open (INS, "> $out1") or die $!;
open (DEL, "> $out2") or die $!;
open (INV, "> $out3") or die $!;
open (TRA, "> $out4") or die $!;
open (DUP, "> $out5") or die $!;

print STDERR "Format $in ... Settings: min length = $min_length, max length = $max_length\n";

my $header1 = "#Chr\tStart\tEnd\tRef\tAlt\tType\tLength\tSV_caller\tInfo\n";

print INS  $header1;
print DEL  $header1; 
print INV  $header1; 
print DUP  $header1; 

my %hash = ("1" => 1, "2" => 1, "3" => 1, "4" => 1, "5" => 1, "6" => 1, "7" => 1, "8" => 1, "9" => 1, "10" => 1, "11" => 1, "12" => 1, "13" => 1, "14" => 1, "15" => 1, "16" => 1, "17" => 1, "18" => 1, "19" => 1, "20" => 1, "21" => 1, "22" => 1, "X" => 1, "Y" => 1, "MT" => 1, "chr1" => 1, "chr2" => 1, "chr3" => 1, "chr4" => 1, "chr5" => 1, "chr6" => 1, "chr7" => 1, "chr8" => 1, "chr9" => 1, "chr10" => 1, "chr11" => 1, "chr12" => 1, "chr13" => 1, "chr14" => 1, "chr15" => 1, "chr16" => 1, "chr17" => 1, "chr18" => 1, "chr19" => 1, "chr20" => 1, "chr21" => 1, "chr22" => 1, "chrX" => 1, "chrY" => 1, "chrMT" => 1);

while (my $line = <IN>)
{
	chomp $line;
	(my $chr, my $start, my $id, my $ref, my $type, my $qual, my $filter, my $info, my $format, my $sample) = split("\t", $line);

	my $index = index($info, "CHR2=");
	my $chr2  = substr($info, $index+length("CHR2="));
	my @b = split(";", $chr2);
	my $chr2 = $b[0];

	my $index = index($info, "END=");
	my $end = substr($info, $index+length("END="));
	my @b = split(";", $end);
	my $end = $b[0];

	my $index = index($info, "SVLEN=");
	my $sv_length = substr($info, $index+length("SVLEN="));
	my @b = split(";", $sv_length);
	my $sv_length = $b[0];

	my $index = index($info, "RE=");
	my $rd_supp = substr($info, $index+length("RE="));
	my @b = split(";", $rd_supp);
	my $rd_supp = $b[0];
	next if ($rd_supp < $min_read_supp );

	if ($type =~ /TRA/){
		$type = "TRA";
		if (exists ($hash{$chr}) and exists ($hash{$chr2})){
			my $start1 = $start;
			my $end1 = $start1+1;
			my $start2 = $end;
			my $end2 = $start2+1;
			my $item = "$chr\t$start1\t$end1\t$chr2\t$start2\t$end2\tcomments:sv_type=$type;sv_caller=Sniffles;numReads=$rd_supp;$filter;$info;$format;$sample\n";
			print TRA $item;
		}
		next;
	}

	if ($type =~ /INS/){
		$type = "INS";
	}elsif($type =~ /DEL/){
		$type = "DEL";
	}elsif($type =~ /DUP/){
		$type = "DUP";
	}elsif($type =~ /INV/){
		$type = "INV";
	}

	$start = $start - 1; 

	if ($start > $end){
		print STDERR "Notice: start position ($start) is smaller than end ($end) ! Switch start and end.\n";
		my $temp = $start;
		$start = $end;
		$end = $temp;
	} # switch position if start < end;

	my $item = "$chr\t$start\t$end\t0\t0\tcomments:sv_type=$type;sv_length=$sv_length;sv_caller=Sniffles;numReads=$rd_supp;$filter;$info;$format;$sample\n";

######filter results######

	next if (!exists ($hash{$chr})); 
	next if ($sv_length < $min_length ); # filter sv length
	next if ($sv_length > $max_length ); # filter sv length
	next if ($end - $start > 1000 and $type eq "INS");     # fuzzy breakpoint of INS

######output results######

	if ($type eq "INS") {
		print INS $item;
	}elsif ($type eq "DEL") {
		print DEL $item;
	}elsif ($type eq "DUP") {
		print DUP $item;
	}elsif ($type eq "INV") {
		print INV $item;
	}
}

close IN;
close INS;
close DEL;
close INV;
close TRA;
close DUP;
