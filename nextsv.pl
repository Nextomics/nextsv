#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

my $usage  = "
Program:  NextSV (SV detection from long-read sequencing)
Version:  1.0.0
Usage:    perl $0 <config>
Contact:  Li Fang (fangli\@grandomics.com)\n\n";

die $usage if (@ARGV < 1);

my $config = shift (@ARGV);
my $nextsv_dir = abs_path(dirname($0));

##### analysis of config file #####
my %arg;
open (CFG, $config) or die $!;
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


##### parameters and settings #####

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
	print "Generating scripts for BLASR/PBHoney...\n";
	&run_pbhoney;	
}

if ($arg{enable_bwa_Sniffles} == 1){
	print "Generating scripts for bwa/Sniffles...\n";
	&run_bwa_sniffles;
}

if ($arg{enable_ngmlr_Sniffles} == 1){
	print "Generating scripts for ngmlr/Sniffles...\n";
	&run_ngmlr_sniffles;
}

if ($arg{enable_PBHoney_Spots} + $arg{enable_PBHoney_Tails} + $arg{enable_bwa_Sniffles} + $arg{enable_ngmlr_Sniffles} > 1){
	my $nextsv_res_dir = "$arg{out_dir}/nextsv_results";
	my $nextsv_sh      = "$nextsv_res_dir/nextsv.sh";

	`mkdir -p $nextsv_res_dir`;
	`chmod +x $nextsv_sh`;

	open (SH, "> $nextsv_sh") or die $!;
	print SH "#!/bin/bash\n\n";
	my @c = split ("/", $arg{input_file_list});
	my $list_name = $c[-1];

	if ($arg{enable_PBHoney_Spots} == 1){
		print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.spots.DEL.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.spots.INS.bed $nextsv_res_dir/ \n";
	}

	if ($arg{enable_PBHoney_Tails} == 1){
		print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.tails.DEL.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.tails.INS.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.tails.INV.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/pbhoney/4_results/$list_name.tails.TRA.vcf $nextsv_res_dir/ \n";
	}

	if ($arg{enable_bwa_Sniffles} == 1){
		print SH "cp $arg{out_dir}/bwa_sniffles/3_vcf/$list_name.bwa.sniffles.vcf.DEL.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/bwa_sniffles/3_vcf/$list_name.bwa.sniffles.vcf.INS.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/bwa_sniffles/3_vcf/$list_name.bwa.sniffles.vcf.INV.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/bwa_sniffles/3_vcf/$list_name.bwa.sniffles.vcf.DUP.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/bwa_sniffles/3_vcf/$list_name.bwa.sniffles.TRA.vcf $nextsv_res_dir/ \n";
	}
	
	if ($arg{enable_ngmlr_Sniffles} == 1){
		print SH "cp $arg{out_dir}/ngmlr_sniffles/3_vcf/$list_name.ngmlr.sniffles.vcf.DEL.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/ngmlr_sniffles/3_vcf/$list_name.ngmlr.sniffles.vcf.INS.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/ngmlr_sniffles/3_vcf/$list_name.ngmlr.sniffles.vcf.INV.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/ngmlr_sniffles/3_vcf/$list_name.ngmlr.sniffles.vcf.DUP.bed $nextsv_res_dir/ \n";
		print SH "cp $arg{out_dir}/ngmlr_sniffles/3_vcf/$list_name.ngmlr.sniffles.TRA.vcf $nextsv_res_dir/ \n";
	}

	my $pbhoney_del = "";
	my $pbhoney_ins = "";

	if ($arg{enable_PBHoney_Spots} == 1 and $arg{enable_PBHoney_Tails} == 1){

		$pbhoney_del = "$nextsv_res_dir/$list_name.PBHoney.DEL.bed";
		$pbhoney_ins = "$nextsv_res_dir/$list_name.PBHoney.INS.bed";

		print SH "$arg{merge2sv} $nextsv_res_dir/$list_name.spots.DEL.bed $nextsv_res_dir/$list_name.tails.DEL.bed $nextsv_res_dir/$list_name.PBHoney-Spots-Tails.intersect.DEL.bed $pbhoney_del DEL\n";
		print SH "$arg{merge2sv} $nextsv_res_dir/$list_name.spots.INS.bed $nextsv_res_dir/$list_name.tails.INS.bed $nextsv_res_dir/$list_name.PBHoney-Spots-Tails.intersect.INS.bed $pbhoney_ins INS\n";
		print SH "rm $nextsv_res_dir/$list_name.PBHoney-Spots-Tails.intersect.DEL.bed\n";
		print SH "rm $nextsv_res_dir/$list_name.PBHoney-Spots-Tails.intersect.INS.bed\n";

	}elsif ($arg{enable_PBHoney_Spots} == 1){

		$pbhoney_del = "$nextsv_res_dir/$list_name.PBHoney.DEL.bed";
		$pbhoney_ins = "$nextsv_res_dir/$list_name.PBHoney.INS.bed";

		print SH "ln -s $nextsv_res_dir/$list_name.spots.DEL.bed $pbhoney_del\n";
		print SH "ln -s $nextsv_res_dir/$list_name.spots.INS.bed $pbhoney_ins\n";

	}elsif ($arg{enable_PBHoney_Tails} == 1){

		$pbhoney_del = "$nextsv_res_dir/$list_name.PBHoney.DEL.bed";
		$pbhoney_ins = "$nextsv_res_dir/$list_name.PBHoney.INS.bed";

		print SH "ln -s $nextsv_res_dir/$list_name.tails.DEL.bed $pbhoney_del\n";
		print SH "ln -s $nextsv_res_dir/$list_name.tails.INS.bed $pbhoney_ins\n";
	}

	my $nexsv_sensitive_del = "";
	my $nexsv_sensitive_ins = "";
	my $nexsv_stringent_del = "";
	my $nexsv_stringent_ins = "";

	$nexsv_sensitive_del = "$nextsv_res_dir/$list_name.NextSV.sensitive.DEL.bed";
	$nexsv_sensitive_ins = "$nextsv_res_dir/$list_name.NextSV.sensitive.INS.bed";
	$nexsv_stringent_del = "$nextsv_res_dir/$list_name.NextSV.stringent.DEL.bed";
	$nexsv_stringent_ins = "$nextsv_res_dir/$list_name.NextSV.stringent.INS.bed";

	if ( ($arg{enable_PBHoney_Spots} == 1 or $arg{enable_PBHoney_Tails} == 1) and $arg{enable_ngmlr_Sniffles} == 1 ){

		print SH "$arg{merge2sv} $pbhoney_del $nextsv_res_dir/$list_name.ngmlr.sniffles.vcf.DEL.bed $nexsv_stringent_del  $nexsv_sensitive_del DEL\n";	
		print SH "$arg{merge2sv} $pbhoney_ins $nextsv_res_dir/$list_name.ngmlr.sniffles.vcf.INS.bed $nexsv_stringent_ins  $nexsv_sensitive_ins INS\n";	

	}

	if ( ($arg{enable_PBHoney_Spots} == 1 or $arg{enable_PBHoney_Tails} == 1) and $arg{enable_bwa_Sniffles} == 1 ){

		print SH "$arg{merge2sv} $pbhoney_del $nextsv_res_dir/$list_name.bwa.sniffles.vcf.DEL.bed $nexsv_stringent_del  $nexsv_sensitive_del DEL\n";	
		print SH "$arg{merge2sv} $pbhoney_ins $nextsv_res_dir/$list_name.bwa.sniffles.vcf.INS.bed $nexsv_stringent_ins  $nexsv_sensitive_ins INS\n";	

	}

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

	open (LIST, $arg{input_file_list}) or die $!;
	open (QSUB, "> $sh_dir/qsub.sh") or die $!;
	print QSUB "#!/bin/bash\n\n";
	my $i = 0;
	my @split_bam;
	my $prefix;
	my $fq_index = 0;
	while (my $line = <LIST>)
	{
		$fq_index += 1;
		chomp $line;
		my @a = split("/", $line);
		my $fq = $a[-1];
		my $out = "blasr_align.$fq.$fq_index.sh";

		my $blasr_sam      = "$blasr_bam_dir/$fq.$fq_index.blasr.sam";
		my $blasr_bam      = "$blasr_bam_dir/$fq.$fq_index.blasr.bam";

		my $tails_sam      = "$tails_bam_dir/$fq.$fq_index.tails.sam";
		my $tails_bam      = "$tails_bam_dir/$fq.$fq_index.tails.bam";

		my $tails_sort_bam = "$tails_bam_dir/$fq.$fq_index.tails.sort.bam";

		open (OUT, "> $sh_dir/$out") or die $!;
		print OUT "#!/bin/bash\n\n";
		print OUT "echo \"program started at:\"";
		print OUT "&& date\n\n";
		print OUT "\n########## blasr alignment ##########\n";
		print OUT "$arg{blasr} $line $arg{ref_blasr} -sa $arg{ref_sa_blasr} -nproc $arg{n_thread} -bestn 1 -sam -clipping subread -out $blasr_sam\n\n";

		print OUT "\n########## tail realignment ##########\n";
		print OUT "python $arg{honey} pie --nproc $arg{n_thread} --output $tails_sam $blasr_sam $arg{ref_blasr} \n\n";

		print OUT "\n########## sam to bam ##########\n";
		print OUT "$arg{samtools} view -bS -@ $arg{n_thread} $blasr_sam > $blasr_bam\n\n";
		print OUT "$arg{samtools} view -bS -@ $arg{n_thread} $tails_sam > $tails_bam\n\n";

		print OUT "sleep 1s\n\n";

		print OUT "\n########## sort bam files ##########\n";
		print OUT "$arg{samtools} sort -@ $arg{n_thread} -o $tails_sort_bam $tails_bam\n\n";
		print OUT "$arg{samtools} index $tails_sort_bam\n";
		print OUT "sleep 1s\n\n";
		print OUT "echo \"program finished at:\"";
		print OUT "&& date\n\n";
		close OUT;
		$split_bam[$i] = "$tails_sort_bam";
		$i++;
		print QSUB "qsub -V -cwd -S /bin/bash -pe smp $arg{n_thread} $out && sleep 1s\n";
	}
	my $n_bam = $i;
	close QSUB;
	close LIST;

	my @b = split ("/", $arg{input_file_list});
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
		print OUT "python $arg{honey} spots --nproc $arg{n_thread} --reference $arg{ref_blasr} --threshold $arg{threshold} --minErrReads $arg{minErrReads} --consensus $arg{consensus} --output $spots_out_prefix $merge_bam\n\n";
		print OUT "sleep 1s\n";
		print OUT "perl $arg{format_pbhoney_spots} $spots_out_prefix.spots\n"; 
	}
	print OUT "$arg{bamstat} $merge_bam > $result_dir/$merge_bam.statistics.txt\n";
	print OUT "echo \"program finished at:\"";
	print OUT "&& date\n\n";

	close OUT;

}

sub run_bwa_sniffles{

	my $raw_bam_dir   = "$arg{out_dir}/bwa_sniffles/1_bwa_bam";
	my $merge_bam_dir = "$arg{out_dir}/bwa_sniffles/2_merge_bam";
	my $vcf_dir       = "$arg{out_dir}/bwa_sniffles/3_vcf";
	my $sh_dir        = "$arg{out_dir}/bwa_sniffles/sh";
	my $qsub          = "$sh_dir/qsub.sh";

	`mkdir -p $raw_bam_dir`;
	`mkdir -p $merge_bam_dir`;
	`mkdir -p $vcf_dir`;
	`mkdir -p $sh_dir`;

	my @a = split("/", $arg{input_file_list});
	my $list_name = $a[-1];

	open (IN, $arg{input_file_list}) or die $!;
	open (QSUB, "> $qsub") or die $!;
	print QSUB "#!/bin/bash\n\n";

	my $i = 0;
	my $line;
	my @bwamem_bam;
	my $fq_index = 0;
	while ($line = <IN>)
	{
		chomp $line;
		@a = split("/", $line);
		my $fq = $a[-1];
		my $out1 = "bwa_align.$fq.$fq_index.sh";

		my $bwa_sam      = "$raw_bam_dir/$fq.$fq_index.bwa.sam";
		my $bwa_sort_bam = "$raw_bam_dir/$fq.$fq_index.bwa.sort.bam";

		open (OUT1, "> $sh_dir/$out1") or die $!;
		print OUT1 "#!/bin/bash\n\n";
		print OUT1 "echo \"program started at:\"";
		print OUT1 "&& date\n\n";
		print OUT1 "$arg{bwa} mem -x pacbio -t $arg{n_thread} -M $arg{ref_bwa} $line > $bwa_sam\n\n";
		print OUT1 "$arg{samtools} sort -@ $arg{n_thread} -o $bwa_sort_bam $bwa_sam\n\n";
		print OUT1 "$arg{samtools} index $bwa_sort_bam\n\n";
		#print OUT1 "rm $bwa_sam\n\n";
		print OUT1 "echo \"program finished at:\"";
		print OUT1 "&& date\n\n";
		close OUT1;
		$bwamem_bam[$i] = $bwa_sort_bam;
		$i++;

		print QSUB "qsub -V -cwd -S /bin/bash -pe smp $arg{n_thread}  $out1 && sleep 1s\n";
	}
	my $n_bam = $i;
	close IN;
	close QSUB;

	my $merge_sh  = "$sh_dir/merge_and_call.$list_name.bwa_sniffles.sh";
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

sub run_ngmlr_sniffles{


	my $raw_bam_dir   = "$arg{out_dir}/ngmlr_sniffles/1_ngmlr_bam";
	my $merge_bam_dir = "$arg{out_dir}/ngmlr_sniffles/2_merge_bam";
	my $vcf_dir       = "$arg{out_dir}/ngmlr_sniffles/3_vcf";
	my $sh_dir        = "$arg{out_dir}/ngmlr_sniffles/sh";
	my $qsub          = "$sh_dir/qsub.sh";

	`mkdir -p $raw_bam_dir`;
	`mkdir -p $merge_bam_dir`;
	`mkdir -p $vcf_dir`;
	`mkdir -p $sh_dir`;

	my @a = split("/", $arg{input_file_list});
	my $list_name = $a[-1];

	open (IN, $arg{input_file_list}) or die $!;
	open (QSUB, "> $qsub") or die $!;
	print QSUB "#!/bin/bash\n\n";

	my $i = 0;
	my $line;
	my @ngmlr_bam;
	my $fq_index = 0;
	while ($line = <IN>)
	{
		$fq_index += 1;
		chomp $line;
		@a = split("/", $line);
		my $fq = $a[-1];
		my $out1 = "ngmlr_align.$fq.$fq_index.sh";

		my $ngmlr_sam      = "$raw_bam_dir/$fq.$fq_index.ngmlr.sam";
		my $ngmlr_sort_bam = "$raw_bam_dir/$fq.$fq_index.ngmlr.sort.bam";

		open (OUT1, "> $sh_dir/$out1") or die $!;
		print OUT1 "#!/bin/bash\n\n";
		print OUT1 "echo \"program started at:\"";
		print OUT1 "&& date\n\n";
		print OUT1 "$arg{ngmlr} -t $arg{n_thread} -r $arg{ref_bwa} -q $line -o $ngmlr_sam\n\n";
		print OUT1 "$arg{samtools} sort -@ $arg{n_thread} -o $ngmlr_sort_bam $ngmlr_sam\n\n";
		print OUT1 "$arg{samtools} index $ngmlr_sort_bam\n\n";
		print OUT1 "echo \"program finished at:\"";
		print OUT1 "&& date\n\n";
		close OUT1;
		$ngmlr_bam[$i] = $ngmlr_sort_bam;
		$i++;

		print QSUB "qsub -V -cwd -S /bin/bash -pe smp $arg{n_thread} $out1 && sleep 1s\n";
	}
	my $n_bam = $i;
	close IN;
	close QSUB;

	my $merge_sh  = "$sh_dir/merge_and_call.$list_name.ngmlr_sniffles.sh";
	my $merge_bam = "$merge_bam_dir/$list_name.ngmlr.merge.sort.bam";
	my $vcf       = "$vcf_dir/$list_name.ngmlr.sniffles.vcf";


	open (OUT2, "> $merge_sh");
	print OUT2 "#!/bin/bash\n\n";
	print OUT2 "echo \"program started at:\"";
	print OUT2 "&& date\n\n";
	print OUT2 "$arg{samtools} merge $merge_bam \\\n";
	for ($i = 0; $i < $n_bam; $i++) {
		print OUT2 "        $ngmlr_bam[$i]\\\n";
	}
	print OUT2 "\n\n";

	print OUT2 "$arg{samtools} index $merge_bam\n";
	print OUT2 "$arg{sniffles} -m $merge_bam --vcf $vcf --min_support $arg{min_support} --max_distance $arg{max_distance} --threads $arg{n_thread}\n";
	print OUT2 "perl $arg{format_sniffles_vcf} $vcf\n";
	print OUT2 "echo \"program finished at:\"";
	print OUT2 "&& date\n\n";
	close OUT2;
}
