#!/usr/bin/env python3

import os
import sys
from datetime import datetime
import argparse
import glob
import subprocess

TimeFormat = '%m/%d/%Y %H:%M:%S'

class Setting:
    def __init__(self, nextsv_root_dir):

        self.root_dir                  = os.path.abspath(nextsv_root_dir)
        self.longreadqc                = os.path.join(self.root_dir, 'bin/longreadqc')
        self.pigz                      = os.path.join(self.root_dir, 'bin/pigz')
        self.seqtk                     = os.path.join(self.root_dir, 'bin/seqtk')
        self.minimap2                  = os.path.join(self.root_dir, 'bin/minimap2')
        self.ngmlr                     = os.path.join(self.root_dir, 'bin/ngmlr')
        self.meryl                     = os.path.join(self.root_dir, 'bin/meryl')
        self.winnowmap                 = os.path.join(self.root_dir, 'bin/winnowmap')
        self.check_bam_and_remove_file = os.path.join(self.root_dir, 'bin/check_bam_and_remove_file.py')

        ## required arguments
        self.in_dir      = None
        self.out_dir     = None
        self.sample_name = None
        self.ref_fasta   = None
        self.platform    = None
        self.aligners    = None
        
        ## optional arguments
        self.conda_env = ''
        self.samtools  = ''
        self.sniffles  = ''
        self.cuteSV    = ''
        self.threads   = 4

        self.input_fastq_list = []
        self.input_fasta_list = []
        self.clean_input_fastq = None
        
        self.aligner_list = []
        ## shell scripts
        self.get_clean_reads_sh_file = None

        # paths
        self.clean_reads_dir = None
        self.bam_dir         = None
        self.sv_calls_dir    = None

    def clean_input_prefix(self):
        return os.path.join(self.clean_reads_dir, self.sample_name)
    def aligned_sam_file(self, aligner_name):
        return os.path.join(self.bam_dir, f'{self.sample_name}.{aligner_name}.sam')
    def aligned_bam_file(self, aligner_name):
        return os.path.join(self.bam_dir, f'{self.sample_name}.{aligner_name}.bam')
    def sorted_bam_file(self, aligner_name):
        return os.path.join(self.bam_dir, f'{self.sample_name}.{aligner_name}.sorted.bam')
    def aligner_shell_file(self, aligner_name):
        return os.path.join(self.bam_dir, f'run_{aligner_name}.{self.sample_name}.sh')
    def sv_detection_shell_file(self, aligner_name, svcaller_name):
        return os.path.join(self.sv_calls_dir, f'run_{svcaller_name}_for_{aligner_name}.{self.sample_name}.sh')
    def sv_vcf_file(self, aligner_name, svcaller_name):
        return os.path.join(self.sv_calls_dir, f'{self.sample_name}.{aligner_name}.{svcaller_name}.vcf')
    
def main():

    program = 'nextsv3.py'    
    parser = argparse.ArgumentParser(prog = program, description='nextsv3: an automated pipeline for structrual variation detection from long-read sequencing. \nContact: Li Fang(fangli2718@gmail.com)')

    # required
    parser.add_argument('-i', '--in_dir',      required = True, metavar = 'path/to/input_dir',   type = str, help = '(required) path to input FASTQ/FASTA directory (all fastq/fasta files in this directory will be used as input files)')
    parser.add_argument('-o', '--out_dir',     required = True, metavar = 'path/to/output_dir',  type = str, help = '(required) path to output fastq directory')
    parser.add_argument('-s', '--sample_name', required = True, metavar = 'sample_name',         type = str, help = '(required) a unique name or id for the input sample')
    parser.add_argument('-r', '--ref_fasta',   required = True, metavar = 'ref.fasta',           type = str, help = '(required) path to reference genome sequence in FASTA format')
    parser.add_argument('-p', '--platform',    required = True, metavar = 'sequencing_platform', type = str, help = '(required) sequencing platform. Three valid values: ont/clr/hifi ont: oxford nanopore; clr: PacBio CLR; hifi: PacBio HiFi/CCS')
    
    # optional
    parser.add_argument('-a', '--aligners',    required = True, metavar = 'aligners_to_use', type = str, default = 'minimap2', help = '(optional) which aligner(s) to use. Three supported aligners: minimap2, ngmlr, winnowmap. Use "+" to combine multiple aligners. Examples: minimap2+ngmlr, winnowmap+minimap2, minimap2+ngmlr+winnowmap, minimap2, ngmlr, winnowmap')
    
    parser.add_argument('-t', '--threads',     required = False, metavar = 'INT',   type = int, default = 4,  help = '(optional) number of threads (default: 4)')
    parser.add_argument('-e', '--conda_env',   required = False, metavar = 'conda_env',  type = str, default = '', help = '(optional) conda environment name (default: NULL)')

    parser.add_argument('--samtools', required = False, metavar = 'path/to/samtools',  type = str, default = 'samtools', help = '(optional) path to samtools (default: using environment default)')
    parser.add_argument('--sniffles', required = False, metavar = 'path/to/sniffles',  type = str, default = 'sniffles', help = '(optional) path to sniffles (default: using environment default)')
    parser.add_argument('--cuteSV', required = False, metavar = 'path/to/cuteSV',  type = str, default = 'cuteSV', help = '(optional) path to cuteSV (default: using environment default)')
    parser.add_argument('--minimap2', required = False, metavar = 'path/to/minimap2',  type = str, default = '', help = 'this parameter is deprecated as minimap2 is supplied in NextSV now')
  
    ### Version
    parser.add_argument('-v', '--version', action='version', version= f'NextSV-3.1.0')
 
    if len(sys.argv) < 2 or sys.argv[1] in ['help', 'h', '-help', 'usage']:
        input_args = parser.parse_args(['--help'])
    else:
        input_args = parser.parse_args()

    myprint(f'Input arguments are: {sys.argv}\n')
    generate_scripts(input_args)

    return

def generate_aligner_list(settings:Setting):
    
    valid_aligners = ['minimap2', 'ngmlr', 'winnowmap']
    settings.aligner_list = settings.aligners.split('+')
    for aligner in settings.aligner_list:
        if aligner not in valid_aligners:
            myprint(f'ERROR! unknown aligners: {settings.aligners}. Supported aligners: minimap2, ngmlr, winnowmap')
            sys.exit(1)


    
def generate_scripts(input_args):

    nextsv3 = os.path.abspath(__file__)
    nextsv_root_dir = os.path.split(nextsv3)[0]
    
    settings                      = Setting(nextsv_root_dir)
    settings.in_dir               = os.path.abspath(input_args.in_dir)
    settings.out_dir              = os.path.abspath(input_args.out_dir)
    settings.sample_name          = input_args.sample_name
    settings.ref_fasta            = os.path.abspath(input_args.ref_fasta)
    settings.platform             = input_args.platform.lower()
    settings.conda_env            = input_args.conda_env
    settings.samtools             = input_args.samtools
    settings.sniffles             = input_args.sniffles
    settings.cuteSV               = input_args.cuteSV
    settings.threads              = input_args.threads
    settings.aligners             = input_args.aligners.strip()

    settings.clean_reads_dir      = os.path.join(settings.out_dir, '1_clean_reads')
    settings.bam_dir              = os.path.join(settings.out_dir, '2_aligned_bam')
    settings.sv_calls_dir         = os.path.join(settings.out_dir, '3_SV_calls')

    if settings.platform not in ['ont', 'hifi', 'clr']:
        myprint(f'ERROR: unknown platform: {settings.platform}')
        myprint(f'--platform should be one of the following: ont, hifi, clr (not case sensitive)')
        sys.exit(1)

    generate_aligner_list(settings)
    
    myprint('NOTICE: checking external tools')
    get_full_path_of_tools(settings)

    myprint('NOTICE: making output directories')
    os.makedirs(settings.out_dir, exist_ok=True)
    os.makedirs(settings.clean_reads_dir, exist_ok=True)
    os.makedirs(settings.bam_dir, exist_ok=True)
    os.makedirs(settings.sv_calls_dir, exist_ok=True)

    myprint('NOTICE: generating shell scripts')

    list_input_files(settings)
    clean_input_files(settings)
    print('###### Output information ######')
    myprint(f'NOTICE: clean reads will be here: {settings.clean_reads_dir}')

    if 'minimap2' in settings.aligner_list:
        minimap2_align(settings)
    if 'ngmlr' in settings.aligner_list:
        ngmlr_align(settings)
    if 'winnowmap' in settings.aligner_list:
        winnowmap_align(settings)
        
    myprint(f'NOTICE: aligned bam files will be here: {settings.bam_dir}')

    for aligner in settings.aligner_list:
        sniffles_detection(settings, aligner)
        cuteSV_detection(settings, aligner)
    
    myprint(f'NOTICE: SV calls will be here: {settings.sv_calls_dir}')

    work_sh_file = os.path.join(settings.out_dir, 'work.sh')
    work_sh_f = open(work_sh_file, 'w')
    work_sh_f.write('#!/bin/bash\n\n')
    if settings.conda_env != '':
        work_sh_f.write('eval \"$(conda shell.bash hook)\"\n')
        work_sh_f.write(f'conda activate {settings.conda_env}\n\n')

    work_sh_f.write(f'sh {settings.get_clean_reads_sh_file}\n\n')

    for aligner in settings.aligner_list:
        aligner_shell_file = settings.aligner_shell_file(aligner)
        work_sh_f.write(f'sh {aligner_shell_file}\n\n')

        for sv_caller in ['sniffles', 'cuteSV']:
            sv_detection_shell_file = settings.sv_detection_shell_file(aligner, sv_caller)
            work_sh_f.write(f'sh {sv_detection_shell_file}\n\n')
    
    work_sh_f.close()
    
    print('\n################################')
    myprint(f'NOTICE: Please run the following shell script: {work_sh_file}\n\n')

    return

def get_full_path_of_tools(settings:Setting):

    tool_name = settings.minimap2
    if os.path.exists(tool_name) == False or os.path.isfile(tool_name) == False:
        report_tool_not_found_check_installation(tool_name)
        sys.exit(1)
        
    myprint(f'NOTICE: minimap2 is here: {settings.minimap2}')
    
    tool_name = settings.ngmlr
    if os.path.exists(tool_name) == False or os.path.isfile(tool_name) == False:
        report_tool_not_found_check_installation(tool_name)
        sys.exit(1)

    myprint(f'NOTICE: ngmlr is here: {settings.ngmlr}')
    
    if settings.samtools == 'samtools':
        tool_name = 'samtools'
        result = subprocess.run(['which', tool_name], stdout=subprocess.PIPE)
        if result.returncode == 0:
            settings.samtools = result.stdout.decode().strip()
        else:
            report_tool_not_found_in_env(tool_name)
            sys.exit(1)
    else:
        settings.samtools = os.path.abspath(settings.samtools)
    
    myprint(f'NOTICE: samtools is here: {settings.samtools}')
    
    if settings.sniffles == 'sniffles':
        tool_name = 'sniffles'
        result = subprocess.run(['which', tool_name], stdout=subprocess.PIPE)
        if result.returncode == 0:
            settings.sniffles = result.stdout.decode().strip()
        else:
            report_tool_not_found_in_env(tool_name)
            sys.exit(1)
    else:
        settings.sniffles = os.path.abspath(settings.sniffles)

    myprint(f'NOTICE: sniffles is here: {settings.sniffles}')
    
    if settings.cuteSV == 'cuteSV':
        tool_name = 'cuteSV'
        result = subprocess.run(['which', tool_name], stdout=subprocess.PIPE)
        if result.returncode == 0:
            settings.cuteSV = result.stdout.decode().strip()
        else:
            report_tool_not_found_in_env(tool_name)
            sys.exit(1)
    else:
        settings.cuteSV = os.path.abspath(settings.cuteSV)
    
    myprint(f'NOTICE: cuteSV is here: {settings.cuteSV}')
    
    
    return

def report_tool_not_found_in_env(tool_name):

    myprint(f'ERROR: {tool_name} is not found in the environment! Please check if you have activated the conda environment. If the error still exists, please supply --{tool_name} with the full path')

    return

def report_tool_not_found_check_installation(tool_name):
    myprint(f'ERROR: {tool_name} is not found! Please check if your installation is correct.')
    return

def list_input_files(settings:Setting):

    input_dir = os.path.join(settings.in_dir)

    file_list = glob.glob(f'{input_dir}/*')

    for file in file_list:
        if file[-6:].lower() == '.fasta':
            settings.input_fasta_list.append(file)
        elif file[-9:].lower() == '.fasta.gz':
            settings.input_fasta_list.append(file)
        elif file[-3:].lower() == '.fa':
            settings.input_fasta_list.append(file)
        elif file[-6:].lower() == '.fa.gz':
            settings.input_fasta_list.append(file)
        elif file[-6:].lower() == '.fastq':
            settings.input_fastq_list.append(file)
        elif file[-9:].lower() == '.fastq.gz':
            settings.input_fastq_list.append(file)
        elif file[-3:].lower() == '.fq':
            settings.input_fastq_list.append(file)
        elif file[-6:].lower() == '.fq.gz':
            settings.input_fastq_list.append(file)

    myprint(f'NOTICE: {len(settings.input_fastq_list)} FASTQ files were found in input directory {input_dir}')
    myprint(f'NOTICE: {len(settings.input_fasta_list)} FASTA files were found in input directory {input_dir}')

    if len(settings.input_fastq_list) == 0 and len(settings.input_fasta_list) == 0:
        myprint(f'ERROR: FASTQ/FASTA files were NOT found in input directory: {settings.in_dir}')
        sys.exit(1)
    
    return

def clean_input_files(settings:Setting):

    settings.get_clean_reads_sh_file = os.path.join(settings.clean_reads_dir, 'get_clean_reads.sh')
    cmd  = '#!/bin/bash\n\n'

    fake_fastq_list = []
    fake_fastq_dir = None
    # convert FASTA to FASTQ
    if len(settings.input_fasta_list) > 0:
        # convert fasta to fastq as NGMLR has bugs for input fasta files
        fake_fastq_dir = os.path.join(settings.clean_reads_dir, 'FA_to_FQ')
        cmd += f'rm -rf {fake_fastq_dir}\n\n'
        cmd += f'mkdir -p {fake_fastq_dir}\n\n'
        for fa_file in settings.input_fasta_list:
            fa_filename = os.path.split(fa_file)[1]
            fake_fastq = os.path.join(fake_fastq_dir, f'{fa_filename}.converted.fastq.gz')
            cmd += f'{settings.seqtk} seq -F . {fa_file} | {settings.pigz} - > {fake_fastq}\n\n'
            fake_fastq_list.append(fake_fastq)

    # clean FASTQ
    input_fastq_list_file = os.path.join(settings.clean_reads_dir, 'raw_input_fastq.list')
    input_fastq_list_f    = open(input_fastq_list_file, 'w')
    for fastq in settings.input_fastq_list:
        input_fastq_list_f.write(f'{fastq}\n')

    for fastq in fake_fastq_list:
        input_fastq_list_f.write(f'{fastq}\n')
    input_fastq_list_f.close()

    settings.clean_input_fastq = settings.clean_input_prefix() + '.clean.fastq'
    cmd += f'rm -f {settings.clean_input_fastq}*\n\n'

    cmd += f'{settings.longreadqc} filterfq --input_list_file {input_fastq_list_file} -p {settings.clean_input_prefix()} -n 1 \n\n'
    cmd += f'{settings.pigz} --processes {settings.threads} {settings.clean_input_fastq}\n\n'
    cmd += f'{settings.longreadqc} fq -i {settings.clean_input_fastq}.gz -d {settings.clean_reads_dir} -p {settings.sample_name}\n\n'
    
    if fake_fastq_dir:
        cmd += f'rm -rf {fake_fastq_dir}\n\n'

    settings.clean_input_fastq = settings.clean_input_fastq + '.gz' # settings.clean_input_fastq will be used in other places

    get_clean_reads_sh_f = open(settings.get_clean_reads_sh_file, 'w')
    get_clean_reads_sh_f.write(cmd)
    get_clean_reads_sh_f.close()

    return

def cuteSV_detection(settings:Setting, aligner_name):

    svcaller_name = 'cuteSV'
    sv_detection_shell_file = settings.sv_detection_shell_file(aligner_name, svcaller_name)
    input_bam = settings.sorted_bam_file(aligner_name)
    out_vcf = settings.sv_vcf_file(aligner_name, svcaller_name)
    cuteSV_workdir = os.path.join(settings.sv_calls_dir, 'cuteSV_temp')

    if settings.platform == 'ont':
        platform_arguments = ' --max_cluster_bias_INS 100  --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100  --diff_ratio_merging_DEL 0.3 '
    elif settings.platform == 'hifi':
        platform_arguments = ' --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 '
    elif settings.platform == 'clr':
        platform_arguments = ' --max_cluster_bias_INS 100  --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200  --diff_ratio_merging_DEL 0.5 '
    else:
        myprint(f'ERROR: unknown platform: {settings.platform}')
        sys.exit(1)

    cmd  = '#!/bin/bash\n\n'
    cmd += f'mkdir -p {cuteSV_workdir} \n'
    cmd += f'{settings.cuteSV} {platform_arguments} --report_readid --genotype --min_support 5 --threads {settings.threads} --sample {settings.sample_name} {input_bam} {settings.ref_fasta} {out_vcf} {cuteSV_workdir} \n'
    cmd += f'rm -r {cuteSV_workdir} \n'

    sh_fp = open(sv_detection_shell_file, 'w')
    sh_fp.write(cmd)
    sh_fp.close()

def sniffles_detection(settings:Setting, aligner_name):

    svcaller_name = 'sniffles'
    sv_detection_shell_file = settings.sv_detection_shell_file(aligner_name, svcaller_name)
    input_bam = settings.sorted_bam_file(aligner_name)
    out_vcf = settings.sv_vcf_file(aligner_name, svcaller_name)

    cmd  = '#!/bin/bash\n\n'
    cmd += f'{settings.sniffles} --output-rnames --allow-overwrite --input {input_bam} --vcf {out_vcf} --reference {settings.ref_fasta} --threads {settings.threads} \n'

    sh_fp = open(sv_detection_shell_file, 'w')
    sh_fp.write(cmd)
    sh_fp.close()

    return

def ngmlr_align_for1input(settings:Setting, input_file, aligned_sam_file, sorted_bam_file):
    
    if settings.platform == 'ont':
        platform_arguments = ' -x ont '
    elif settings.platform == 'hifi' or 'cls':
        platform_arguments = ' -x pacbio '
    else:
        myprint(f'ERROR: unknown platform: {settings.platform}')
        sys.exit(1)
        
    cmd = f'{settings.ngmlr} -t {settings.threads} -r {settings.ref_fasta} -q {input_file} {platform_arguments} -o {aligned_sam_file}\n\n'
    cmd += f'{settings.samtools} sort -@ {settings.threads} -o {sorted_bam_file} {aligned_sam_file}\n\n'
    cmd += f'{settings.samtools} index -@ {settings.threads} {sorted_bam_file}\n\n'
    cmd += f'{settings.check_bam_and_remove_file} {sorted_bam_file} {aligned_sam_file} {settings.samtools}\n\n'
    
    return cmd

def ngmlr_align(settings:Setting):
    
    aligner_name = 'ngmlr'
    aligned_sam_file   = settings.aligned_sam_file(aligner_name)
    sorted_bam_file    = settings.sorted_bam_file(aligner_name)
    aligner_shell_file = settings.aligner_shell_file(aligner_name)
        
    cmd = '#!/bin/bash\n\n'
    
    cmd += ngmlr_align_for1input(settings, settings.clean_input_fastq, aligned_sam_file, sorted_bam_file)

    sh_fp = open(aligner_shell_file, 'w')
    sh_fp.write(cmd)
    sh_fp.close()
    
    return


def minimap2_align(settings:Setting):

    aligner_name = 'minimap2'
    aligned_bam_file   = settings.aligned_bam_file(aligner_name)
    sorted_bam_file    = settings.sorted_bam_file(aligner_name)
    aligner_shell_file = settings.aligner_shell_file(aligner_name)

    if settings.platform == 'ont':
        platform_arguments = '-x map-ont'
    elif settings.platform == 'hifi':
        platform_arguments = '-x map-hifi'
    elif settings.platform == 'clr':
        platform_arguments = '-x map-pb'
    else:
        myprint(f'ERROR: unknown platform: {settings.platform}')
        sys.exit(1)
    
    cmd = '#!/bin/bash\n\n'
    cmd += f'{settings.minimap2} --MD -t {settings.threads} -a {platform_arguments} -N 10 {settings.ref_fasta} {settings.clean_input_fastq} | {settings.samtools} view -@ 2 -bS - > {aligned_bam_file}\n\n'
    cmd += f'{settings.samtools} sort -@ {settings.threads} -o {sorted_bam_file} {aligned_bam_file}\n\n'
    cmd += f'{settings.samtools} index -@ {settings.threads} {sorted_bam_file} \n\n'
    cmd += f'{settings.check_bam_and_remove_file} {sorted_bam_file} {aligned_bam_file} {settings.samtools}\n\n'

    sh_fp = open(aligner_shell_file, 'w')
    sh_fp.write(cmd)
    sh_fp.close()

    return

def winnowmap_align(settings:Setting):
    
    aligner_name = 'winnowmap'
    aligned_bam_file   = settings.aligned_bam_file(aligner_name)
    sorted_bam_file    = settings.sorted_bam_file(aligner_name)
    aligner_shell_file = settings.aligner_shell_file(aligner_name)

    if settings.platform == 'ont':
        platform_arguments = '-x map-ont'
    elif settings.platform == 'hifi':
        platform_arguments = '-x map-pb'
    elif settings.platform == 'clr':
        platform_arguments = '-x map-pb-clr'
    else:
        myprint(f'ERROR: unknown platform: {settings.platform}')
        sys.exit(1)
    
    ref_dir             = os.path.split(settings.ref_fasta)[0]
    merylDB_dir         = os.path.join(settings.ref_fasta, '.merylDB') 
    repetitive_k15_file = os.path.join(settings.ref_fasta, '.repetitive_k15.txt')
    
    repetitive_k15_file_is_valid = False
    if os.path.exists(repetitive_k15_file) and os.path.getsize(repetitive_k15_file) > 0:
        repetitive_k15_file_is_valid = True
        
    if repetitive_k15_file_is_valid == False and os.access(ref_dir, os.W_OK) == False:
        myprint(f'ERROR! Failed to create a folder for storing meryl output, which is needed by winnowmap alignment. Please make sure you have written permission in the reference FASTA folder: \n{ref_dir}\nIf you are unable to get written permission in this folder, you can create a new reference FASTA folder in your own space, soft link the FASTA files inside the new folder, and supply NextSV with the path to the new FASTA folder.')
        sys.exit(1)
    
    cmd = '#!/bin/bash\n\n'
    if repetitive_k15_file_is_valid == False:
        cmd += f'{settings.meryl} count k=15 output {merylDB_dir} {settings.ref_fasta}\n\n'
        cmd += f'{settings.meryl} print greater-than distinct=0.9998 {merylDB_dir} > {repetitive_k15_file}\n\n'
    
    cmd += f'{settings.winnowmap} --MD -t {settings.threads} -a {platform_arguments} {settings.ref_fasta} {settings.clean_input_fastq} | {settings.samtools} view -@ 2 -bS - > {aligned_bam_file}\n\n'
    cmd += f'{settings.samtools} sort -@ {settings.threads} -o {sorted_bam_file} {aligned_bam_file}\n\n'
    cmd += f'{settings.samtools} index -@ {settings.threads} {sorted_bam_file} \n\n'
    cmd += f'{settings.check_bam_and_remove_file} {sorted_bam_file} {aligned_bam_file} {settings.samtools}\n\n'

    sh_fp = open(aligner_shell_file, 'w')
    sh_fp.write(cmd)
    sh_fp.close()

    return

def myprint(string):
    print ('[' + datetime.now().strftime(TimeFormat) + '] ' + string )
    return 

if __name__ == '__main__':
    main()
