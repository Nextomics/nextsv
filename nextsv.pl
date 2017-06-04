#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

my $usage  = "
Program:  NextSV (Automated SV detection for long-read sequencing)
Version:  0.4.0
Usage:    perl $0 <config>
Contact:  Li Fang (fangli\@grandomics.com)\n\n";

die $usage if (@ARGV < 1);

my $cfg = shift (@ARGV);
my $nextsv_dir = abs_path(dirname($0));

my %arg;
open (CFG, $cfg) or die $!;
while (my $line = <CFG>)
{
	chomp $line;
	next if ($line =~ /^#/);

	($line,) = split ("#", $line);
	next if (!$line);

	$line = &rtrim($line);

	my @a = split ("=", $line);
	if ($a[0]){
		$arg{$a[0]} = $a[1];
	}
}
close CFG;

$arg{blasr}                 = "$nextsv_dir/bin/blasr";
$arg{bwa}                   = "$nextsv_dir/bin/bwa";
$arg{ngmlr}                 = "$nextsv_dir/bin/ngmlr";

$arg{honey}                 = "$nextsv_dir/aligners_and_callers/PBSuite_15.8.24/bin/Honey.py";
$arg{sniffles}              = "$nextsv_dir/bin/sniffles";

$arg{bamstat}               = "$nextsv_dir/bin/bamstat";
$arg{samtools}              = "$nextsv_dir/bin/samtools1.3";
$arg{merge2sv}              = "$nextsv_dir/bin/merge2sv.pl";
$arg{format_pbhoney_tails}  = "$nextsv_dir/bin/format_pbhoney_tails.pl";
$arg{format_sniffles_vcf}   = "$nextsv_dir/bin/format_sniffles_vcf.pl";
$arg{format_pbhoney_spots}  = "$nextsv_dir/bin/format_pbhoney_spots.pl";


if ($arg{enable_PBHoney_Spots} == 1 or $arg{enable_PBHoney_Tails} == 1){
	print "Generating scripts for PBHoney...\n";
	&run_pbhoney;	
}

if ($arg{enable_Sniffles} == 1){
	print "Generating scripts for Sniffles...\n";
	&run_sniffles_bwa;
	&run_sniffles_ngmlr;
}

if ($arg{enable_PBHoney_Spots} and $arg{enable_Sniffles}){
	my $combine_dir = "$arg{out_dir}/combination";
	my $combine_sh = "$combine_dir/combine.sh";

	`mkdir -p $combine_dir`;
	open (SH, "> $combine_sh") or die $!;
	print SH "#!/bin/bash\n\n";
	my @c  = split ("/", $arg{fastq_list});
	my $list_name  = $c[-1];

	print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.spots.DEL.bed $combine_dir/ \n";
	print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.spots.INS.bed $combine_dir/ \n";
	print SH "cp $arg{out_dir}/sniffles/3_vcf/$list_name.sniffles.vcf.DEL.bed $combine_dir/ \n";
	print SH "cp $arg{out_dir}/sniffles/3_vcf/$list_name.sniffles.vcf.INS.bed $combine_dir/ \n";
	print SH "$arg{merge2sv} $combine_dir/$list_name.spots.DEL.bed $combine_dir/$list_name.sniffles.vcf.DEL.bed $combine_dir/$list_name.PBHoneySpots-Sniffles.intersect.DEL.bed $combine_dir/$list_name.PBHoneySpots-Sniffles.union.DEL.bed DEL\n"; 
	print SH "$arg{merge2sv} $combine_dir/$list_name.spots.INS.bed $combine_dir/$list_name.sniffles.vcf.INS.bed $combine_dir/$list_name.PBHoneySpots-Sniffles.intersect.INS.bed $combine_dir/$list_name.PBHoneySpots-Sniffles.union.INS.bed INS\n"; 
	close SH;
}

sub rtrim { 
	my $s = shift; 
	$s =~ s/\s+$//;       
	return $s
};

sub run_pbhoney{

	my $blasr_bam_dir = "$arg{out_dir}/pbhoney/1_blasr_bam"; 
	my $tails_bam_dir = "$arg{out_dir}/pbhoney/2_tails_bam";   
	my $merge_bam_dir = "$arg{out_dir}/pbhoney/3_merge_bam";
	my $result_dir    = "$arg{out_dir}/pbhoney/4_results";
	my $sh_dir        = "$arg{out_dir}/pbhoney/sh";

	`mkdir -p $blasr_bam_dir`;
	`mkdir -p $tails_bam_dir`;
	`mkdir -p $merge_bam_dir`;
	`mkdir -p $result_dir`;
	`mkdir -p $sh_dir`;

	open (LIST, $arg{fastq_list}) or die $!;
	open (QSUB, "> $sh_dir/qsub.sh") or die $!;
	print QSUB "#!/bin/bash\n\n";
	my $i = 0;
	my @split_bam;
	my $prefix;
	while (my $line = <LIST>)
	{
		chomp $line;
		my @a = split("/", $line);
		my $fq = $a[-1];
		my $out = "align.$fq.blasr.sh";

		my $blasr_sam      = "$blasr_bam_dir/$fq.blasr.sam";
		my $blasr_bam      = "$blasr_bam_dir/$fq.blasr.bam";

		my $tails_sam      = "$tails_bam_dir/$fq.tails.sam";
		my $tails_bam      = "$tails_bam_dir/$fq.tails.bam";

		my $tails_sort_bam = "$tails_bam_dir/$fq.tails.sort.bam";

		open (OUT, "> $sh_dir/$out") or die $!;
		print OUT "#!/bin/bash\n\n";
		print OUT "echo \"program started at:\"";
		print OUT "&& date\n\n";
		print OUT "\n########## blasr alignment ##########\n";
		print OUT "$arg{blasr} $line $arg{ref_blasr} -sa $arg{ref_sa_blasr} -nproc $arg{n_thread} -bestn 1 -sam -clipping subread -out $blasr_sam\n\n";

		print OUT "\n########## tail realignment ##########\n";
		print OUT "python $arg{honey} pie --nproc $arg{n_thread} --output $tails_sam $blasr_sam $arg{ref_blasr} \n\n";

		print OUT "\n########## sam to bam ##########\n";
		print OUT "$arg{samtools} view -bS  -@ $arg{n_thread} $blasr_sam > $blasr_bam\n\n";
		print OUT "$arg{samtools} view -bS  -@ $arg{n_thread} $tails_sam > $tails_bam\n\n";
		print OUT "rm $blasr_sam\n";

		print OUT "sleep 1s\n\n";

		print OUT "\n########## sort bam files ##########\n";
		print OUT "$arg{samtools} sort -@ $arg{n_thread} -o $tails_sort_bam $tails_bam\n\n";
		print OUT "$arg{samtools} index $tails_sort_bam\n";
		print OUT "sleep 1s\n\n";
		#print OUT "rm $tails_bam\n\n";
		print OUT "echo \"program finished at:\"";
		print OUT "&& date\n\n";
		close OUT;
		$split_bam[$i] = "$tails_sort_bam";
		$i++;
		print QSUB "qsub -cwd -q all.q -S /bin/bash -pe smp $arg{n_thread}  $out && sleep 1s\n";
	}
	my $n_bam = $i;
	close QSUB;
	close LIST;

	my @b  = split ("/", $arg{fastq_list});
	my $list_name  = $b[-1];
	my $merge_sh   = "$sh_dir/merge_and_call.$list_name.pbhoney.sh\n";
	my $merge_bam  = "$merge_bam_dir/$list_name.bam";

	open (OUT, "> $merge_sh") or die $!;
	print OUT "#!/bin/bash\n\n";
	print OUT "echo \"program started at:\"";
	print OUT "&& date\n\n";
	print OUT "\n########## merge sorted bam ##########\n";

	print OUT "$arg{samtools} merge $merge_bam \\\n";
	for ($i = 0; $i < $n_bam; $i++) {
		print OUT "        $split_bam[$i]\\\n";
	}

	print OUT "\n\n";
	print OUT "$arg{samtools} index $merge_bam\n\n";

	if ($arg{enable_PBHoney_Tails} == 1){
		my $tails_out = "$result_dir/$list_name.tails";
		print OUT "\n########## PBhoney-Tails ##########\n";
		print OUT "python $arg{honey} tails --buffer $arg{buffer} --minBreads $arg{minBreads} --minZMWs $arg{minZMWs} --output $tails_out $merge_bam\n";
		print OUT "sleep 1s\n";
		print OUT "perl $arg{format_pbhoney_tails} $tails_out\n";
	}

	if ($arg{enable_PBHoney_Spots} == 1){
		my $spots_out_prefix = "$result_dir/$list_name";
		if ($arg{consensus} eq "none"){
			$arg{consensus} = "None";
		}
		print OUT "\n########## PBhoney-Spots ##########\n";
		print OUT "python $arg{honey} spots --nproc $arg{n_thread} --reference $arg{ref_blasr} --threshold $arg{threshold} --minErrReads $arg{minErrReads}  --consensus $arg{consensus}  --output $spots_out_prefix $merge_bam\n\n";
		print OUT "sleep 1s\n";
		print OUT "perl $arg{format_pbhoney_spots} $spots_out_prefix.spots\n"; 
	}
	print OUT "$arg{bamstat} $merge_bam > $result_dir/$merge_bam.statistics.txt\n";
	print OUT "echo \"program finished at:\"";
	print OUT "&& date\n\n";

	close OUT;

}

sub run_sniffles_bwa{


	my $raw_bam_dir   = "$arg{out_dir}/sniffles_bwa/1_bwa_bam";
	my $merge_bam_dir = "$arg{out_dir}/sniffles_bwa/2_merge_bam";
	my $vcf_dir       = "$arg{out_dir}/sniffles_bwa/3_vcf";
	my $sh_dir        = "$arg{out_dir}/sniffles_bwa/sh";
	my $qsub          = "$sh_dir/qsub.sh";

	`mkdir -p $raw_bam_dir`;
	`mkdir -p $merge_bam_dir`;
	`mkdir -p $vcf_dir`;
	`mkdir -p $sh_dir`;

	my @a = split("/", $arg{fastq_list});
	my $list_name = $a[-1];

	open (IN, $arg{fastq_list}) or die $!;
	open (QSUB, "> $qsub") or die $!;
	print QSUB "#!/bin/bash\n\n";

	my $i = 0;
	my $line;
	my @bwamem_bam;
	while ($line = <IN>)
	{
		chomp $line;
		@a = split("/", $line);
		my $fq = $a[-1];
		my $out1 = "align.$fq.bwa.sh";

		my $bwa_sam      = "$raw_bam_dir/$fq.bwa.sam";
		my $bwa_sort_bam = "$raw_bam_dir/$fq.bwa.sort.bam";

		open (OUT1, "> $sh_dir/$out1") or die $!;
		print OUT1 "#!/bin/bash\n\n";
		print OUT1 "echo \"program started at:\"";
		print OUT1 "&& date\n\n";
		print OUT1 "$arg{bwa} mem -x pacbio -t $arg{n_thread} -M $arg{ref_bwa} $line  > $bwa_sam\n\n";
		print OUT1 "$arg{samtools} sort -@ $arg{n_thread} -o $bwa_sort_bam $bwa_sam\n\n";
		print OUT1 "$arg{samtools} index $bwa_sort_bam\n\n";
		#print OUT1 "rm $bwa_sam\n\n";
		print OUT1 "echo \"program finished at:\"";
		print OUT1 "&& date\n\n";
		close OUT1;
		$bwamem_bam[$i] = $bwa_sort_bam;
		$i++;

		print QSUB "qsub -cwd -S /bin/bash -q all.q -pe smp $arg{n_thread}  $out1 && sleep 1s\n";
	}
	my $n_bam = $i;
	close IN;
	close QSUB;

	my $merge_sh  = "$sh_dir/merge_and_call.$list_name.sniffles_bwa.sh";
	my $merge_bam = "$merge_bam_dir/$list_name.bwa.merge.sort.bam";
	my $vcf       = "$vcf_dir/$list_name.bwa.sniffles.vcf";


	open (OUT2, "> $merge_sh");
	print OUT2 "#!/bin/bash\n\n";
	print OUT2 "echo \"program started at:\"";
	print OUT2 "&& date\n\n";
	print OUT2 "$arg{samtools} merge $merge_bam \\\n";
	for ($i = 0; $i < $n_bam; $i++) {
		print OUT2  "        $bwamem_bam[$i]\\\n";
	}
	print OUT2 "\n\n";

	print OUT2 "$arg{samtools} index $merge_bam\n";
	print OUT2 "$arg{sniffles} -m $merge_bam --vcf $vcf --min_support $arg{min_support} --max_distance $arg{max_distance} --threads $arg{n_thread}\n";
	print OUT2 "perl $arg{format_sniffles_vcf} $vcf\n";
	print OUT2 "echo \"program finished at:\"";
	print OUT2 "&& date\n\n";
	close OUT2;
}

sub run_sniffles_ngmlr{


	my $raw_bam_dir   = "$arg{out_dir}/sniffles_ngmlr/1_ngmlr_bam";
	my $merge_bam_dir = "$arg{out_dir}/sniffles_ngmlr/2_merge_bam";
	my $vcf_dir       = "$arg{out_dir}/sniffles_ngmlr/3_vcf";
	my $sh_dir        = "$arg{out_dir}/sniffles_ngmlr/sh";
	my $qsub          = "$sh_dir/qsub.sh";

	`mkdir -p $raw_bam_dir`;
	`mkdir -p $merge_bam_dir`;
	`mkdir -p $vcf_dir`;
	`mkdir -p $sh_dir`;

	my @a = split("/", $arg{fastq_list});
	my $list_name = $a[-1];

	open (IN, $arg{fastq_list}) or die $!;
	open (QSUB, "> $qsub") or die $!;
	print QSUB "#!/bin/bash\n\n";

	my $i = 0;
	my $line;
	my @ngmlr_bam;
	while ($line = <IN>)
	{
		chomp $line;
		@a = split("/", $line);
		my $fq = $a[-1];
		my $out1 = "align.$fq.ngmlr.sh";

		my $ngmlr_sam      = "$raw_bam_dir/$fq.ngmlr.sam";
		my $ngmlr_sort_bam = "$raw_bam_dir/$fq.ngmlr.sort.bam";

		open (OUT1, "> $sh_dir/$out1") or die $!;
		print OUT1 "#!/bin/bash\n\n";
		print OUT1 "echo \"program started at:\"";
		print OUT1 "&& date\n\n";
		print OUT1 "$arg{ngmlr}  -t $arg{n_thread} -r $arg{ref_bwa} -q $line -o $ngmlr_sam\n\n";
		print OUT1 "$arg{samtools} sort -@ $arg{n_thread} -o $ngmlr_sort_bam $ngmlr_sam\n\n";
		print OUT1 "$arg{samtools} index $ngmlr_sort_bam\n\n";
		print OUT1 "echo \"program finished at:\"";
		print OUT1 "&& date\n\n";
		close OUT1;
		$ngmlr_bam[$i] = $ngmlr_sort_bam;
		$i++;

		print QSUB "qsub -cwd -S /bin/bash -q all.q -pe smp $arg{n_thread}  $out1 && sleep 1s\n";
	}
	my $n_bam = $i;
	close IN;
	close QSUB;

	my $merge_sh  = "$sh_dir/merge_and_call.$list_name.sniffles_ngmlr.sh";
	my $merge_bam = "$merge_bam_dir/$list_name.ngmlr.merge.sort.bam";
	my $vcf       = "$vcf_dir/$list_name.ngmlr.sniffles.vcf";


	open (OUT2, "> $merge_sh");
	print OUT2 "#!/bin/bash\n\n";
	print OUT2 "echo \"program started at:\"";
	print OUT2 "&& date\n\n";
	print OUT2 "$arg{samtools} merge $merge_bam \\\n";
	for ($i = 0; $i < $n_bam; $i++) {
		print OUT2  "        $ngmlr_bam[$i]\\\n";
	}
	print OUT2 "\n\n";

	print OUT2 "$arg{samtools} index $merge_bam\n";
	print OUT2 "$arg{sniffles} -m $merge_bam --vcf $vcf --min_support $arg{min_support} --max_distance $arg{max_distance} --threads $arg{n_thread}\n";
	print OUT2 "perl $arg{format_sniffles_vcf} $vcf\n";
	print OUT2 "echo \"program finished at:\"";
	print OUT2 "&& date\n\n";
	close OUT2;
}
