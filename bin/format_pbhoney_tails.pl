#!/usr/bin/perl
use strict;

my $usage  = "Usage:   perl $0 <in.tails>\n";
$usage    .= "Contact: Li Fang (fangli\@grandomics.com)\n";
$usage    .= "Version: 0.4.0\n";

die $usage if (@ARGV < 1);

my $min_length = 50;       # min sv length required for DEL and INV
my $max_length = 1000000;  # max sv length required for DEL and INV

my $in  = shift(@ARGV);
my $del = "$in.DEL.bed";
my $inv = "$in.INV.bed";
my $tra = "$in.TRA.vcf";
my $ins = "$in.INS.bed";

open (IN,  $in ) or die $!;
open (INS, "> $ins") or die $!;
open (DEL, "> $del") or die $!;
open (INV, "> $inv") or die $!;
open (TRA, "> $tra") or die $!;

print STDERR "Format $in ... Settings: min sv length = $min_length, max sv length = $max_length\n";

my $header1 = "#Chr\tStart\tEnd\tRef\tAlt\tType\tLength\tSV_caller\tInfo\n";
print INS $header1;
print DEL $header1;
print INV $header1;

my %hash = ("1" => 1, "2" => 1, "3" => 1, "4" => 1, "5" => 1, "6" => 1, "7" => 1, "8" => 1, "9" => 1, "10" => 1, "11" => 1, "12" => 1, "13" => 1, "14" => 1, "15" => 1, "16" => 1, "17" => 1, "18" => 1, "19" => 1, "20" => 1, "21" => 1, "22" => 1, "X" => 1, "Y" => 1, "MT" => 1, "chr1" => 1, "chr2" => 1, "chr3" => 1, "chr4" => 1, "chr5" => 1, "chr6" => 1, "chr7" => 1, "chr8" => 1, "chr9" => 1, "chr10" => 1, "chr11" => 1, "chr12" => 1, "chr13" => 1, "chr14" => 1, "chr15" => 1, "chr16" => 1, "chr17" => 1, "chr18" => 1, "chr19" => 1, "chr20" => 1, "chr21" => 1, "chr22" => 1, "chrX" => 1, "chrY" => 1, "chrMT" => 1); 

my $line = <IN>;
my $bamname;
if ($line =~ /^#Args/)
{
	my $index = index($line, "bam="); 
	my $tmp = substr($line, $index+length("bam="));
	my @a = split("'", $tmp);
	$bamname = $a[1];
}else{
	print STDERR "Warning: bam name not recognized.\n";
	$bamname = "Unknown";
}

my $header2 = "##fileformat=VCFv4.1
##ALT=<ID=TRA,Description=\"Translocation\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">
##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">
##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of reads\">
##INFO=<ID=ZMW,Number=1,Type=Integer,Description=\"Number of ZMWs\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">;
#CHROM\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$bamname\n";

print TRA $header2;

while ($line = <IN>)
{
	chomp $line;
	(my $id, my $chrKey, my $uRef, my $uBreak, my $uMapq, my $dRef, my $dBreak, my $dMapq, my $remainSeq, my $type, my $numReads, my $numZMWs, my $evidence) = split("\t", $line);

	next if (!exists ($hash{$uRef}));
	next if (!exists ($hash{$dRef}));

	if ($type eq "TLOC"){
		print TRA "$uRef\t$uBreak\t$id\t.\t<TRA>\t.\tPASS\tIMPRECISE;SVMETHOD=PBHoney-Tails;CHR2=$dRef;END=$dBreak;RE=$numReads;ZMW=$numZMWs\tGT:DV\t./.\n";
		next;
	}

	$uBreak = $uBreak - 1;  # start from 0 for bed format
	my $sv_length = $dBreak - $uBreak;
	if ($type eq "INS"){
		$sv_length = "N.A.";
	}

	my $bed_item = "$uRef\t$uBreak\t$dBreak\t0\t0\tcomments:sv_type=$type;sv_length=$sv_length;sv_caller=PBhoney-Tails;numReads=$numReads;numZMWs=$numZMWs;uMapq=$uMapq;dMapq=$dMapq;remainSeq=$remainSeq;evidence=$evidence\n";	

	if ($type eq "INS"){
		if ($dBreak - $uBreak <= 1000){
			print INS $bed_item;
		}
	}elsif ($type eq "DEL"){
		if ($sv_length >= $min_length and $sv_length <= $max_length){
			print DEL $bed_item;
		}
	}elsif ($type eq "INV"){
		if ($sv_length >= $min_length and $sv_length <= $max_length){
			print INV $bed_item;
		}
	}
}

close IN;
close INS;
close DEL;
close INV;
close TRA;
