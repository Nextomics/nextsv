#!/usr/bin/env python

import os
import sys
import random
import time
from datetime import datetime
import subprocess

tab  = '\t'
endl = '\n'

usage = '''Program:  nextsv2 (Pipeline for SV detection from nanopore sequencing)
Usage:    python nextsv2.py <config> <sample_name> <input_fastq_folder>
Contact:  Li Fang (fangli@grandomics.com) '''

class Setting:
    def __init__(self, nextsv_root_dir):

        ## job submission settings ##

        self.n_thread = 8
        self.sample_name = None
        self.job_submission_command = None
        self.wait_time = 30

        ## input and output settings ##
        self.root_dir = os.path.abspath(nextsv_root_dir)
        self.bin_dir  = os.path.join(self.root_dir, 'bin')
        self.input_file_list = None
        self.input_list = list()
        self.input_file_format = None
        self.out_dir = None
        self.total_input_file_size = 0


        self.clean_input_file = ''
        self.clean_input_dir = ''
        self.clean_input_prefix = ''

        self.qc_dir = ''
        self.ngmlr_bam_dir = ''
        self.minimap2_bam_dir = ''

        self.sniffles_calls_dir = ''

        ## aligner and SV caller settings ##
        self.enable_ngmlr_sniffles = 1
        self.enable_minimap2_sniffles = 1

        self.ngmlr_bam = ''
        self.minimap2_bam = ''
        self.ngmlr_cram = ''
        self.minimap2_cram = ''

        self.enable_ngmlr_aligner = 0
        self.enable_minimap2_aligner = 0
        self.enable_sniffles_caller = 0

        self.sniffles_min_support = 2
        self.sniffles_max_distance = 600
        self.nextsv_out_dir = None

        ## path ##
        self.ref_fasta = None
        self.ngmlr = os.path.join(self.root_dir, 'bin/ngmlr')
        self.minimap2 = os.path.join(self.root_dir, 'bin/minimap2')
        self.longreadqc = os.path.join(self.root_dir, 'bin/longreadqc')
        self.samtools = 'samtools'
        self.pigz = os.path.join(self.root_dir, 'bin/pigz')
        self.sniffles = 'sniffles'
        self.runtimekey = None 
        self.format_sniffles_vcf = os.path.join(self.root_dir, 'bin/format_sniffles_vcf.pl')
        self.merge2sv = os.path.join(self.root_dir, 'bin/merge2sv.pl') 
        self.ngmlr_vcf_file = None
        self.samtools_version = None
        self.bad_file_list = list()
        self.check_bam_and_remove_file = os.path.join(self.root_dir, 'bin/check_bam_and_remove_file.py')
        self.remove_raw_input_fastq = os.path.join(self.root_dir, 'bin/remove_raw_input_fastq.py')

        self.maf_convert = None


class Task:

    def __init__(self, task_id = None, sh_file = None, out_file = None, status = None, dependent_taskid_list = None):
        self.task_id = task_id
        self.sh_file = sh_file
        self.out_file = out_file
        self.status = status
        self.dependent_taskid_list = dependent_taskid_list


all_task_list = list()
task_id = 0

def main():

    if len(sys.argv) < 4: 
        print (usage)
        sys.exit()

    myprint('program started')

    config_file = sys.argv[1]
    sample_name = sys.argv[2]
    input_dir  = os.path.abspath(sys.argv[3])

    nextsv_root_dir = os.path.split(__file__)[0] 
    settings = Setting(nextsv_root_dir)
    settings.sample_name = sample_name


    myprint('reading config file')
    parse_config_file(config_file, settings)

    myprint('creating output directories')
    myprint('output directory is: %s' % settings.out_dir)
    create_output_dirs(settings)

    myprint('listing input files')
    read_input_files(input_dir, settings) 

    random.seed()
    settings.runtimekey = random.randint(10000000000, 99999999999)

    myprint('checking samtools version')
    check_samtools_version(settings)


    work_sh_file = os.path.join(settings.out_dir, 'work.sh')
    work_sh_fp = open(work_sh_file, 'w')
    work_sh_fp.write('#!/bin/bash\n\n')

    settings.clean_input_prefix = os.path.join(settings.clean_input_dir, '%s' % (settings.sample_name) )

    if settings.input_file_format == 'fasta':
        settings.clean_input_file =  settings.clean_input_prefix + '.clean.fasta.gz'
    else:
        settings.clean_input_file =  settings.clean_input_prefix + '.clean.fastq.gz'

    clean_fq_sh_file = get_clean_input_files(settings)  # input file must be in fasta or fastq format

    work_sh_fp.write('time sh %s\n' % clean_fq_sh_file)

    # minimap2 #
    if settings.enable_minimap2_aligner:
        myprint('generating tasks for Minimap2')
        minimap2_sh_file = generate_tasks_minimap2(settings)
        work_sh_fp.write('time sh %s\n' % minimap2_sh_file)

    if settings.enable_minimap2_sniffles:
        myprint('generating tasks for Minimap2-Sniffles')
        minimap2_sniffles_sh_file = generate_tasks_sniffles(settings, settings.minimap2_bam, 'minimap2')
        work_sh_fp.write('time sh %s\n' % minimap2_sniffles_sh_file)

    if settings.enable_minimap2_aligner:
        myprint('bam to cram for Minimap2 bam')
        input_bam        = settings.minimap2_bam
        output_dir       = settings.minimap2_bam_dir
        settings.minimap2_cram = os.path.join (output_dir, '%s.minimap2.sorted.cram' % settings.sample_name)
        output_cram      = settings.minimap2_cram
        minimap2_bam2cram_sh_file = os.path.join(output_dir, 'minimap2_bam2cram.sh')
        generate_tasks_bam2cram (settings, input_bam, output_cram, minimap2_bam2cram_sh_file)
        work_sh_fp.write ('time sh %s\n' % minimap2_bam2cram_sh_file)

    # ngmlr #
    if settings.enable_ngmlr_aligner:
        myprint('generating tasks for NGMLR')
        ngmlr_sh_file = generate_tasks_ngmlr(settings)
        work_sh_fp.write('time sh %s\n' % ngmlr_sh_file)

    if settings.enable_ngmlr_sniffles:
        myprint('generating tasks for NGMLR-Sniffles')
        ngmlr_sniffles_sh_file = generate_tasks_sniffles(settings, settings.ngmlr_bam, 'ngmlr')
        work_sh_fp.write('time sh %s\n' % ngmlr_sniffles_sh_file)

    if settings.enable_ngmlr_aligner:
        myprint('bam to cram for NGMLR bam')
        input_bam        = settings.ngmlr_bam
        output_dir       = settings.ngmlr_bam_dir
        settings.ngmlr_cram = os.path.join (output_dir, '%s.ngmlr.sorted.cram' % settings.sample_name)
        output_cram      = settings.ngmlr_cram
        ngmlr_bam2cram_sh_file = os.path.join(output_dir, 'ngmlr_bam2cram.sh')
        generate_tasks_bam2cram (settings, input_bam, output_cram, ngmlr_bam2cram_sh_file)
        work_sh_fp.write ('time sh %s\n' % ngmlr_bam2cram_sh_file)


    ## bam to cram ##

    work_sh_fp.close()

    return

def generate_tasks_bam2cram (settings, input_bam, output_cram, sh_file):


    sh_fp = open(sh_file, 'w')
    sh_fp.write('#!/bin/bash\n\n')
    cmd = '%s view --reference %s -@ %d --output-fmt CRAM -hS %s > %s' % ( settings.samtools, settings.ref_fasta, settings.n_thread, input_bam, output_cram)
    sh_fp.write(cmd + endl)
    cmd = '%s index %s' % ( settings.samtools, output_cram)
    sh_fp.write(cmd + endl)
    cmd = 'rm %s' % (input_bam)
    sh_fp.write(cmd + endl)
    cmd = 'rm %s.bai' % (input_bam)
    sh_fp.write(cmd + endl)
    sh_fp.close()

    return

def get_clean_input_files (settings):

    clean_fq_sh_file = os.path.join(settings.clean_input_dir, 'get_clean_fastq.sh')

    clean_fq_sh_fp = open(clean_fq_sh_file, 'w')
    clean_fq_sh_fp.write('#!/bin/bash\n\n')
    cmd = '%s filterfq -l %s -p %s -n 1' % (settings.longreadqc, settings.input_file_list, settings.clean_input_prefix)
    clean_fq_sh_fp.write(cmd + endl)
    cmd = '%s fq -i %s -d %s -p %s' % (settings.longreadqc, settings.clean_input_prefix + '.clean.fastq', settings.qc_dir, settings.sample_name)
    clean_fq_sh_fp.write(cmd + endl)

    clean_fq_sh_fp.write(endl)
    cmd = '%s -p %d %s\n' % (settings.pigz, settings.n_thread, settings.clean_input_prefix + '.clean.fastq') 
    clean_fq_sh_fp.write(cmd)
    clean_fq_sh_fp.close()
    
    return clean_fq_sh_file

def bam2fastq(settings, input_bam_file, out_fastq_file):

    cmd = '%s fastq %s > %s' % (settings.samtools, input_bam_file, out_fastq_file)
    return cmd


def wait_outputfile(settings, output_file_list):

    while 1:
        time.sleep(settings.wait_time)
        output_file_exist = True
        for output_file in output_file_list:
            if os.path.exists(output_file) == False:
                output_file_exist = False
                break
        if output_file_exist == True:
            break

    return



def check_samtools_version(settings):

    results = subprocess.check_output([settings.samtools, "sort"]).decode()
    lines   = results.split(endl)
    target_str = " -f "
    settings.samtools_version = 'new'
    for line in lines:
        if line.find(target_str) >= 0:
            settings.samtools_version = 'old'
            break 
    return



def generate_tasks_sniffles(settings, input_bam, aligner_name):

    out_dir = settings.sniffles_calls_dir
    out_vcf = os.path.join(out_dir, '%s.%s.sniffles.vcf' % (settings.sample_name, aligner_name)) 
    cmd = 'time %s -m %s --vcf %s --min_support %d --max_distance %d --threads %d --num_reads_report -1 --genotype --cluster --report_seq ' % (settings.sniffles, input_bam, out_vcf, settings.sniffles_min_support, settings.sniffles_max_distance, settings.n_thread)

    sh_file = os.path.join(out_dir, 'sniffles.%s.%s.sh' %  (aligner_name, settings.sample_name))

    log_file = sh_file + '.log'
    sh_fp = open(sh_file, 'w')
    sh_fp.write('#!/bin/bash\n\n')
    sh_fp.write(cmd + endl)
    sh_fp.write('echo %d > %s\n' % (settings.runtimekey, log_file))
    sh_fp.close()

    return sh_file



def generate_tasks_minimap2(settings):

    out_dir = settings.minimap2_bam_dir
    out_prefix = os.path.join(out_dir, '%s' % settings.sample_name)

    align_bam_file  = out_prefix + '.minimap2.bam' 
    sorted_bam_file = out_prefix + '.minimap2.sorted.bam'
    bam_index_file = sorted_bam_file + '.bai'

    settings.minimap2_bam = sorted_bam_file

    sh_file = os.path.join(out_dir, 'minimap2.%s.sh' % settings.sample_name)
    input_file = settings.clean_input_file

    align_cmd = 'time %s --cs --MD -t %d -ax map-ont -N 8 %s %s | %s view -@ 2 -bS - > %s\n' % (settings.minimap2, settings.n_thread, settings.ref_fasta, input_file, settings.samtools, align_bam_file)

    sort_cmd = samtools_sort_cmd(settings, align_bam_file, sorted_bam_file, settings.n_thread) + endl

    index_cmd = samtools_index_cmd(settings.samtools, sorted_bam_file) + endl
    
    checkbam_cmd = '%s %s %s %s\n' %  (settings.check_bam_and_remove_file, sorted_bam_file, align_bam_file, settings.samtools)  


    log_file = sh_file + '.log'

    sh_fp = open(sh_file, 'w')
    sh_fp.write('#!/bin/bash\n')
    sh_fp.write(align_cmd)
    sh_fp.write(sort_cmd)
    sh_fp.write(index_cmd)
    sh_fp.write(checkbam_cmd)
    sh_fp.write('echo %d > %s\n' % (settings.runtimekey, log_file) )
    sh_fp.close()

    return sh_file

def generate_tasks_ngmlr (settings):


    out_dir = settings.ngmlr_bam_dir
    out_prefix = os.path.join(out_dir, '%s' % settings.sample_name)

    align_bam_file  = out_prefix + '.ngmlr.bam' 
    sorted_bam_file = out_prefix + '.ngmlr.sorted.bam'
    bam_index_file = sorted_bam_file + '.bai'

    settings.ngmlr_bam = sorted_bam_file

    sh_file = os.path.join(out_dir, 'ngmlr.%s.sh' % settings.sample_name)
    input_file = settings.clean_input_file

    align_cmd = 'time %s -t %d -r %s -q %s | %s view -bS - > %s\n' % (settings.ngmlr, settings.n_thread, settings.ref_fasta, input_file, settings.samtools, align_bam_file)
    sort_cmd = samtools_sort_cmd(settings, align_bam_file, sorted_bam_file, settings.n_thread) + endl
    index_cmd = samtools_index_cmd(settings.samtools, sorted_bam_file) + endl
    checkbam_cmd = '%s %s %s %s\n' %  (settings.check_bam_and_remove_file, sorted_bam_file, align_bam_file, settings.samtools)  


    log_file = sh_file + '.log'

    sh_fp = open(sh_file, 'w')
    sh_fp.write('#!/bin/bash\n')
    sh_fp.write(align_cmd)
    sh_fp.write(sort_cmd)
    sh_fp.write(index_cmd)
    sh_fp.write(checkbam_cmd)
    sh_fp.write('echo %d > %s\n' % (settings.runtimekey, log_file))
    sh_fp.close()

    return sh_file

def samtools_merge_bam_cmd(samtools, sort_bam_list, out_bam):
    cmd = '%s merge %s \\\n' % (samtools, out_bam)
    for sort_bam in sort_bam_list: 
        cmd += '    %s\\\n' % sort_bam
    cmd += '\n\n'
    return cmd
    
def samtools_index_cmd(samtools, sort_bam):
    cmd = 'time %s index %s' % (samtools, sort_bam)
    return  cmd

def samtools_sort_cmd(settings, input_bam, output_bam, n_thread):

    if settings.samtools_version  == 'new':
        cmd = 'time %s sort -m 2G -@ %d -o %s %s' % (settings.samtools, n_thread, output_bam, input_bam) 
    else:
        cmd = 'time %s sort -m 2G -@ %d -f %s %s' % (settings.samtools, n_thread, input_bam, output_bam)

    return cmd

def get_file_prefix(input_file):
    file_name = os.path.split(input_file)[1]
    prefix = os.path.splitext(file_name)[0]
    return prefix

def read_input_files(input_dir, settings):

    total_file_size = 0 
    input_dir = os.path.abspath(input_dir)
    file_list = os.listdir(input_dir) 
    input_list = list()
    ext1 = ''
    ext2 = ''
    for f in file_list:
        split_f = f.split('.')
        ext1 = split_f[-1].lower()
        if len(split_f) >= 2: ext2 = split_f[-2].lower()
        if ext1 == 'fastq' or ext1 == 'fq' or (ext1 == 'gz' and (ext2 == 'fastq' or ext2 == 'fq') )  :   
            in_file = os.path.join(input_dir, f)
            total_file_size += os.path.getsize(in_file)
            input_list.append(in_file)
    
    settings.total_input_file_size = total_file_size


    if len(input_list) == 0:
        myprint('ERROR! no fastq files were found in the input folder: %s' % input_dir)
        sys.exit()
    else:
        myprint('%d fastq files were found in the input folder: %s' % (len(input_list), input_dir))

    settings.input_list = input_list

    list_file = os.path.join(settings.out_dir, 'input_files.list')
    list_fp = open(list_file, 'w')
    for f in input_list:
        list_fp.write('%s\n' % f)

    list_fp.close()
    
    settings.input_file_list = list_file

    return 

def create_output_dirs(settings):

    settings.clean_input_dir = os.path.join(settings.out_dir, 'clean_reads')
    settings.nextsv_out_dir = os.path.join(settings.out_dir, 'nextsv_results')
    settings.ngmlr_bam_dir = os.path.join(settings.out_dir, 'ngmlr_bam')
    settings.minimap2_bam_dir = os.path.join(settings.out_dir, 'minimap2_bam')
    settings.sniffles_calls_dir = os.path.join(settings.out_dir, 'sniffles_calls')
    settings.qc_dir = os.path.join(settings.out_dir, 'longread_qc')
    

    os.system('mkdir -p %s' % settings.clean_input_dir)
    os.system('mkdir -p %s' % settings.nextsv_out_dir)

    if settings.enable_ngmlr_aligner:
        os.system('mkdir -p %s' % settings.ngmlr_bam_dir)

    if settings.enable_minimap2_aligner:
        os.system('mkdir -p %s' % settings.minimap2_bam_dir)

    if settings.enable_sniffles_caller:
        os.system('mkdir -p %s' % settings.sniffles_calls_dir)



    return        

def parse_config_file(config_file, settings):
    
    config_file_fp = open(config_file, 'r')
    while 1: 
        line = config_file_fp.readline()
        if not line: break
        line = line.strip()
        if '#' in line: line = line.split('#')[0]
        if len(line) == 0: continue
        if '=' not in line: continue
        line = line.split('=') 
        if len(line) < 2: continue
        key = line[0].strip().lower()
        value = '='.join(line[1:])
        value = value.strip()
        if (not key) or (not value): continue


        if key == 'out_dir':
            settings.out_dir = os.path.join(os.path.abspath(value), settings.sample_name)
        elif key == 'n_thread' and int(value) > 0:
            settings.n_thread = int(value)
        elif key == 'mode':
            settings.mode = value
        elif key == 'job_submission_command':
            settings.job_submission_command = value.strip('[').strip(']')



        elif key == 'enable_ngmlr_sniffles':
            settings.enable_ngmlr_sniffles = int(value)
        elif key == 'enable_minimap2_sniffles':
            settings.enable_minimap2_sniffles = int(value)


        elif key == 'samtools':
            settings.samtools = os.path.abspath(value)
        elif key == 'pigz':
            settings.pigz   = os.path.abspath(value)
        elif key == 'sniffles':
            settings.sniffles = os.path.abspath(value)
        elif key == 'ref_fasta':
            settings.ref_fasta = os.path.abspath(value)

        elif key == 'sniffles_min_support':
            settings.sniffles_min_support = int(value)    
        elif key == 'sniffles_max_distance':
            settings.sniffles_max_distance = int(value)
        else:
            print ('unknown parameter: %s' % key)


    config_file_fp.close()
        
    if settings.enable_ngmlr_sniffles:
        settings.enable_ngmlr_aligner = 1
        settings.enable_sniffles_caller = 1

    if settings.enable_minimap2_sniffles:
        settings.enable_minimap2_aligner = 1
        settings.enable_sniffles_caller = 1
    

    if settings.sample_name == None:
        myprint('ERROR! sample name is not specified')
        sys.exit()

    if settings.out_dir == None:
        myprint('ERROR! output directory was not specified') 
        sys.exit()



    if check_file_specified_and_existed('ref_fasta (reference fasta file)', settings.ref_fasta) == False: sys.exit()

    if check_file_specified_and_existed('samtools', settings.samtools) == False: sys.exit()

    if check_file_specified_and_existed('longreadqc', settings.longreadqc) == False: sys.exit()

    if settings.enable_minimap2_aligner and check_file_specified_and_existed('minimap2', settings.minimap2) == False: sys.exit()

    if settings.enable_ngmlr_aligner and check_file_specified_and_existed('ngmlr', settings.ngmlr) == False: sys.exit()

    if settings.enable_sniffles_caller and check_file_specified_and_existed('sniffles', settings.sniffles) == False: sys.exit()

    return settings 


TimeFormat = '%m/%d/%Y %H:%M:%S'

def check_file_specified_and_existed (file_name, file_path):

    if file_path == None:
        myprint('ERROR! path to %s was not specified' % file_name) 
        return False
    if os.path.exists(file_path) == False:
        myprint('ERROR! %s was not found in the path: %s' % (file_name, file_path))
        return False

    return True

def myprint(string):

    print ('[' + datetime.now().strftime(TimeFormat) + '] ' + string )

    return 

if __name__ == '__main__':
    main()
