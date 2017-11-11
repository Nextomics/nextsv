#!/usr/bin/env python

import os
import sys
import random
import time

tab  = '\t'
endl = '\n'

usage = '''Program:  NextSV (SV detection from long-read sequencing)
Version:  2.0.0 
Usage:    python nextsv.py <config>
Contact:  Li Fang (fangli@grandomics.com) '''

class Setting:
    def __init__(self, nextsv_root_dir):

        self.root_dir = os.path.abspath(nextsv_root_dir)
        self.bin_dir = os.path.join(self.root_dir, 'bin')
        self.input_file_list = None
        self.input_list = list()
        self.input_file_format = None
        self.out_dir = None
        self.n_thread = 1
        self.sample_name = None
        self.job_submission_command = None
        self.wait_time = 30

        self.enable_PBHoney_Spots = 1
        self.spots_threshold = 3
        self.spots_minErrReads = 2
        self.spots_consensus = 'None'

        self.enable_PBHoney_Tails = 1
        self.tails_minBreads = 2
        self.tails_minZMWs = 2
        self.tails_buffer = 600
        self.blasr_pbhoney_dir = None

        self.enable_bwa_Sniffles = 0
        self.enable_ngmlr_Sniffles = 1
        self.sniffles_min_support = 2
        self.sniffles_max_distance = 600
        self.bwa_sniffles_dir = None
        self.ngmlr_sniffles_dir = None
        self.nextsv_out_dir = None

        self.bash5_minLength = 500
        self.bash5_minReadScore = 0.75
        self.extract_fastq_dir = None

        self.ref_fasta = None
        self.ref_blasr = None
        self.ref_sa_blasr = None
        self.ngmlr = os.path.join(self.root_dir, 'bin/ngmlr')
        self.samtools = None
        self.sniffles = os.path.join(self.root_dir, 'bin/sniffles')
        self.bash5tools = None
        self.runtimekey = None 
        self.pbhoney = os.path.join(self.root_dir, 'aligners_and_callers/PBSuite_15.8.24/bin/Honey.py')
        self.blasr = os.path.join(self.root_dir, 'bin/blasr')
        self.format_pbhoney_spots = os.path.join(self.root_dir, 'bin/format_pbhoney_spots.pl')
        self.format_pbhoney_tails = os.path.join(self.root_dir, 'bin/format_pbhoney_tails.pl')
        self.format_sniffles_vcf = os.path.join(self.root_dir, 'bin/format_sniffles_vcf.pl')
        self.merge2sv = os.path.join(self.root_dir, 'bin/merge2sv.pl') 
        self.spots_file = None
        self.tails_file = None
        self.bwa_vcf_file = None
        self.ngmlr_vcf_file = None
        self.samtools_version = None


class Task:

    def __init__(self, task_id = None, cmd = None, sh_file = None, out_file = None, status = None, dependent_taskid_list = None):
        self.task_id = task_id
        self.cmd = cmd
        self.sh_file = sh_file
        self.out_file = out_file
        self.status = status
        self.dependent_taskid_list = dependent_taskid_list


all_task_list = list()
task_id = 0

def main():
    if len(sys.argv) < 2: 
        print usage
        sys.exit()

    global all_task_list
    global task_id 

    config_file = sys.argv[1]
    settings = parse_config_file(config_file)
    random.seed()
    settings.runtimekey = random.randint(10000000000, 99999999999)

    creat_output_dirs(settings)
    check_samtools_version(settings)
    settings.input_list = read_input_files(settings.input_file_list) 
    check_input_file_format(settings)

    extract_fastq_from_rawdata(settings)

    ngmlr_sniffles_tasks = None
    bwa_sniffles_tasks = None
    blasr_pbhoney_tasks = None

    if settings.enable_ngmlr_Sniffles:
        ngmlr_sniffles_tasks = generate_tasks_ngmlr_sniffles(settings)

    if settings.enable_bwa_Sniffles:
        bwa_sniffles_tasks = generate_tasks_bwa_sniffles(settings)

    if settings.enable_PBHoney_Spots or settings.enable_PBHoney_Tails:
        blasr_pbhoney_tasks = generate_tasks_blasr_pbhoney(settings)

    run_alignment_and_svcalling(settings)

    merging_results(settings)

def extract_fastq_from_rawdata(settings):

    settings.extract_fastq_dir = os.path.join(settings.out_dir, 'fastq')
    os.system('mkdir -p ' + settings.extract_fastq_dir)

    if settings.input_file_format == 'hdf5':
        extract_fastq_from_hdf5(settings)
    elif settings.input_file_format == 'bam':
        extract_fastq_from_bam(settings)
    return

def extract_fastq_from_bam(settings):
    extract_fastq_sh_file = os.path.join(settings.extract_fastq_dir, 'extract_fastq.sh')
    extract_done_file = os.path.join(settings.extract_fastq_dir, 'extract_fastq.%s.finished' % settings.runtimekey)
    extract_fastq_sh_fp = open(extract_fastq_sh_file, 'w') 
    extract_fastq_sh_fp.write('#!/bin/bash\n\n')
    bam_list = settings.input_list
    fastq_list = list()
    for input_bam_file in bam_list:
        bam_prefix = get_file_prefix(input_bam_file)
        out_fastq_file = os.path.join(settings.extract_fastq_dir, bam_prefix + '.fastq')
        fastq_list.append(out_fastq_file)
        cmd = bam2fastq(settings, input_bam_file, out_fastq_file)
        extract_fastq_sh_fp.write(cmd + endl)
         
    extract_fastq_sh_fp.write('touch %s' % extract_done_file + endl)
    extract_fastq_sh_fp.close() 

    all_fastq_file_exists = True
    for fastq_file in fastq_list:
        if os.path.exists(fastq_file) == False: 
            all_fastq_file_exists = False
            break
    if all_fastq_file_exists: 
        settings.input_list = fastq_list
        myprint ('extracted FASTQ files existed, skipped extracting FASTQ files from bam files')
        return

    extract_fastq_task = Task(0, cmd, extract_fastq_sh_file, extract_done_file, 'UNK', list()) 
    submit_task(settings, extract_fastq_task)
    wait_outputfile(settings, fastq_list + [extract_done_file])
    settings.input_list = fastq_list

    return

def bam2fastq(settings, input_bam_file, out_fastq_file):

    cmd = '%s fastq %s > %s' % (settings.samtools, input_bam_file, out_fastq_file)
    return cmd

def extract_fastq_from_hdf5(settings):
    extract_fastq_sh_file = os.path.join(settings.extract_fastq_dir, 'extract_fastq.sh')
    extract_done_file = os.path.join(settings.extract_fastq_dir, 'extract_fastq.%s.finished' % settings.runtimekey)
    extract_fastq_sh_fp = open(extract_fastq_sh_file, 'w') 
    extract_fastq_sh_fp.write('#!/bin/bash\n\n')
    fastq_list = list() 
    hdf5_file_list = settings.input_list

    for hdf5_file in hdf5_file_list:
        hdf5_prefix = get_file_prefix(hdf5_file)
        hdf5_prefix = hdf5_prefix.rstrip('.h5')
        hdf5_prefix = hdf5_prefix.rstrip('.hdf5')
        hdf5_prefix = hdf5_prefix.rstrip('bas')
        hdf5_prefix = hdf5_prefix.rstrip('bax')
        hdf5_prefix = hdf5_prefix.rstrip('.')
        out_fastq_prefix = os.path.join(settings.extract_fastq_dir, hdf5_prefix)
        out_fastq = out_fastq_prefix + '.fastq'
        fastq_list.append(out_fastq)
        cmd = hdf5tofastq(settings, hdf5_file, out_fastq_prefix)
        extract_fastq_sh_fp.write(cmd + endl)

    extract_fastq_sh_fp.write('touch %s' % extract_done_file + endl)
    extract_fastq_sh_fp.close() 

    all_fastq_file_exists = True
    for fastq_file in fastq_list:
        if os.path.exists(fastq_file) == False: 
            all_fastq_file_exists = False
            break
    if all_fastq_file_exists: 
        settings.input_list = fastq_list
        myprint ('extracted FASTQ files existed, skipped extracting FASTQ files from hdf files')
        return

    extract_fastq_task = Task(0, cmd, extract_fastq_sh_file, extract_done_file, 'UNK', list()) 
    submit_task(settings, extract_fastq_task)
    wait_outputfile(settings, fastq_list + [extract_done_file])
    settings.input_list = fastq_list

    return
    

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


def hdf5tofastq(settings, hdf5_file, out_fastq_prefix):
    
    if settings.bash5tools == None or os.path.exists(settings.bash5tools) == False:
        myprint('ERROR! bash5tools.py not specified or does not exist!')
        sys.exit()

    cmd = 'python %s --readType subreads --outType fastq --minLength %d --minReadScore %f --outFilePrefix %s %s ' % (settings.bash5tools, settings.bash5_minLength, settings.bash5_minReadScore, out_fastq_prefix, hdf5_file) 
    return cmd

def check_input_file_format(settings):

    input_file = settings.input_list[0]
    input_file_ext = os.path.splitext(input_file)[1].lower()
    if input_file_ext == '.hdf5' or input_file_ext == '.h5': 
        settings.input_file_format = 'hdf5'
    elif input_file_ext == '.bam':
        settings.input_file_format = 'bam'
    elif input_file_ext == '.sam':
        settings.input_file_format = 'sam'
    elif input_file_ext == '.fasta' or input_file_ext == '.fa':
        settings.input_file_format = 'fasta'
    elif input_file_ext == '.fastq' or input_file_ext == '.fq':
        settings.input_file_format = 'fastq'
    else:
        settings.input_file_format = 'UNK'
        myprint ('ERROR! input files can only be fastq, fasta, or hdf5')
        sys.exit()

    return
    
def check_samtools_version(settings):

    tmp_file = os.path.join(settings.out_dir, '.samtools_tmp')
    tmp_fp = open(tmp_file, 'w')
    tmp_fp.close()

    cmd = '%s sort &>> %s ' % (settings.samtools, tmp_file) 
    os.system(cmd)
    tmp_fp = open(tmp_file, 'r')
    lines = list(tmp_fp)
    tmp_fp.close()
    target_str = " -f "
    settings.samtools_version = 'new'
    for line in lines:
        if line.find(target_str) >= 0:
            settings.samtools_version = 'old'
            break 
    os.system ('rm %s' % tmp_file)
    return

def merging_results(settings):
    merge_sh = os.path.join(settings.nextsv_out_dir, 'merge_calls.sh')
    merge_sh_fp = open(merge_sh, 'w')
    merge_sh_fp.write('#!/bin/bash\n\n')

    ngmlr_sni_bed_file = dict()
    bwa_sni_bed_file = dict()
    spots_bed_file = dict()
    tails_bed_file = dict()
    honey_bed_file = dict()
    sni_bed_file = dict()
    nextsv_sensitive_file = dict()
    nextsv_stringent_file = dict()
    sv_type_list = ['DEL', 'INS']
    if settings.enable_ngmlr_Sniffles:
        input_vcf  = settings.ngmlr_vcf_file
        prefix     = os.path.split(input_vcf)[1]
        out_prefix = os.path.join(settings.nextsv_out_dir, prefix)
        cmd = formatting_sniffles_vcf(settings, input_vcf, out_prefix)
        merge_sh_fp.write(cmd + endl)
        for sv_type in sv_type_list:
            ngmlr_sni_bed_file[sv_type] = out_prefix + '.%s.bed' % sv_type

    if settings.enable_bwa_Sniffles:
        input_vcf  = settings.bwa_vcf_file
        prefix     = os.path.split(input_vcf)[1]
        out_prefix = os.path.join(settings.nextsv_out_dir, prefix)
        cmd = formatting_sniffles_vcf(settings, input_vcf, out_prefix)
        merge_sh_fp.write(cmd + endl)
        for sv_type in sv_type_list:
            bwa_sni_bed_file[sv_type] = out_prefix + '.%s.bed' % sv_type

    if settings.enable_PBHoney_Spots: 
        input_file = settings.spots_file
        prefix     = os.path.split(input_vcf)[1]
        out_prefix = os.path.join(settings.nextsv_out_dir, prefix)
        cmd = formatting_spots_file(settings, input_file, out_prefix)
        merge_sh_fp.write(cmd + endl)
        for sv_type in sv_type_list:
            spots_bed_file[sv_type] = out_prefix + '.%s.bed' % sv_type

    if settings.enable_PBHoney_Tails:
        input_file = settings.tails_file
        prefix     = os.path.split(input_vcf)[1]
        out_prefix = os.path.join(settings.nextsv_out_dir, prefix)
        cmd = formatting_tails_file(settings, input_file, out_prefix)
        merge_sh_fp.write(cmd + endl)
        for sv_type in sv_type_list:
            tails_bed_file[sv_type] = out_prefix + '.%s.bed' % sv_type

    temp_file = dict()

    for sv_type in sv_type_list:
        nextsv_sensitive_file[sv_type] = os.path.join(settings.nextsv_out_dir, settings.sample_name + '.nextsv_sensitive.%s.bed' % sv_type)
        nextsv_stringent_file[sv_type] = os.path.join(settings.nextsv_out_dir, settings.sample_name + '.nextsv_stringent.%s.bed' % sv_type)

    if settings.enable_PBHoney_Spots and settings.enable_PBHoney_Tails and settings.enable_ngmlr_Sniffles:
        sni_bed_file = ngmlr_sni_bed_file 
        for sv_type in sv_type_list:
            honey_bed_file[sv_type] = os.path.join(settings.nextsv_out_dir, settings.sample_name + '.pbhoney.%s.bed' % sv_type) 
            temp_file[sv_type] = os.path.join(settings.nextsv_out_dir, settings.sample_name + '.tmp.%s.bed' % sv_type)
            cmd = merge2sv_file(settings, spots_bed_file[sv_type], tails_bed_file[sv_type], temp_file[sv_type], honey_bed_file[sv_type], sv_type) 
            merge_sh_fp.write(cmd + endl)
            merge_sh_fp.write('rm %s\n' % temp_file[sv_type])
            cmd = merge2sv_file(settings, sni_bed_file[sv_type], honey_bed_file[sv_type], nextsv_stringent_file[sv_type], nextsv_sensitive_file[sv_type], sv_type) 
            merge_sh_fp.write(cmd + endl)

    elif settings.enable_PBHoney_Spots and settings.enable_PBHoney_Tails and settings.enable_bwa_Sniffles:
        sni_bed_file = bwa_sni_bed_file 
        for sv_type in sv_type_list:
            honey_bed_file[sv_type] = os.path.join(settings.nextsv_out_dir, settings.sample_name + '.pbhoney.%s.bed' % sv_type) 
            temp_file[sv_type] = os.path.join(settings.nextsv_out_dir, settings.sample_name + '.tmp.%s.bed' % sv_type)
            cmd = merge2sv_file(settings, spots_bed_file[sv_type], tails_bed_file[sv_type], temp_file[sv_type], honey_bed_file[sv_type], sv_type) 
            merge_sh_fp.write(cmd + endl)
            merge_sh_fp.write('rm %s\n' % temp_file[sv_type])
            cmd = merge2sv_file(settings, sni_bed_file[sv_type], honey_bed_file[sv_type], nextsv_stringent_file[sv_type], nextsv_sensitive_file[sv_type], sv_type) 
            merge_sh_fp.write(cmd + endl)

    elif settings.enable_PBHoney_Spots and settings.enable_ngmlr_Sniffles:
        sni_bed_file = ngmlr_sni_bed_file 
        honey_bed_file = spots_bed_file
        for sv_type in sv_type_list:
            cmd = merge2sv_file(settings, sni_bed_file[sv_type], honey_bed_file[sv_type], nextsv_stringent_file[sv_type], nextsv_sensitive_file[sv_type], sv_type) 
            merge_sh_fp.write(cmd + endl)

    elif settings.enable_PBHoney_Spots and settings.enable_bwa_Sniffles:
        sni_bed_file = bwa_sni_bed_file 
        honey_bed_file = spots_bed_file
        for sv_type in sv_type_list:
            cmd = merge2sv_file(settings, sni_bed_file[sv_type], honey_bed_file[sv_type], nextsv_stringent_file[sv_type], nextsv_sensitive_file[sv_type], sv_type) 
            merge_sh_fp.write(cmd + endl)

    elif settings.enable_PBHoney_Tails and settings.enable_ngmlr_Sniffles:
        sni_bed_file = ngmlr_sni_bed_file 
        honey_bed_file = tails_bed_file
        for sv_type in sv_type_list:
            cmd = merge2sv_file(settings, sni_bed_file[sv_type], honey_bed_file[sv_type], nextsv_stringent_file[sv_type], nextsv_sensitive_file[sv_type], sv_type) 
            merge_sh_fp.write(cmd + endl)

    elif settings.enable_PBHoney_Tails and settings.enable_bwa_Sniffles:
        sni_bed_file = bwa_sni_bed_file 
        honey_bed_file = tails_bed_file
        for sv_type in sv_type_list:
            cmd = merge2sv_file(settings, sni_bed_file[sv_type], honey_bed_file[sv_type], nextsv_stringent_file[sv_type], nextsv_sensitive_file[sv_type]) 
            merge_sh_fp.write(cmd + endl)

    merge_sh_fp.close()

    cmd = 'cd %s && ' % settings.nextsv_out_dir
    if settings.mode == 'sge':
        cmd += settings.job_submission_command + ' ' + merge_sh
    else:
        cmd += 'sh ' + merge_sh

    os.system(cmd + endl)

    return

def merge2sv_file(settings, bed1_file, bed2_file, out_intersect_bed_file, out_union_bed_file, sv_type):

    cmd = settings.merge2sv + ' %s %s %s %s %s ' % (bed1_file, bed2_file, out_intersect_bed_file, out_union_bed_file, sv_type)
    return cmd

def formatting_tails_file(settings, input_file, out_prefix):

    cmd = 'perl %s %s %s' % (settings.format_pbhoney_tails, input_file, out_prefix) 
    return cmd

def formatting_spots_file(settings, input_file, out_prefix):

    cmd = 'perl %s %s %s' % (settings.format_pbhoney_spots, input_file, out_prefix) 
    return cmd

def formatting_sniffles_vcf(settings, input_vcf, out_prefix):

    cmd = 'perl %s %s %s' % (settings.format_sniffles_vcf, input_vcf, out_prefix) 
    return cmd

def run_alignment_and_svcalling(settings):

    global all_task_list
    while 1:
        for i in range(0, len(all_task_list)): 
            task = all_task_list[i]
            if ready2submit(task): submit_task(settings, task)

        unfinished_task_cnt = 0
        for i in range(0, len(all_task_list)):
            task = all_task_list[i]
            unfinished_task_cnt += check_task_status(settings, i)

        if unfinished_task_cnt == 0: break
        time.sleep(settings.wait_time) 

    print 'all job finished'
    return 

def ready2submit(task):

    if task.status == 'finished' or task.status == 'submitted': return False
    global all_task_list
    if len(task.dependent_taskid_list) == 0: return True

    for task_id in task.dependent_taskid_list:
        if all_task_list[task_id].status != 'finished': return False

    return True

def check_task_status(settings, task_id):

    global all_task_list
    task = all_task_list[task_id]

    if task.status == 'finished': return 0

    task_status_file = task.out_file + '.finished'
    if not os.path.exists(task_status_file): return 1
    task_status_fp = open(task_status_file, 'r')
    lines = list(task_status_fp)
    runtimekey = int(lines[0].strip())
    if runtimekey == settings.runtimekey:
        all_task_list[task_id].status = 'finished'    
        return 0
    else:
        return 1

def submit_task(settings, task):
    sh_file = task.sh_file
    sh_dir = os.path.split(sh_file)[0]
    submit_cmd = 'cd %s &&' % sh_dir 
    if settings.mode == 'sge':
        submit_cmd += settings.job_submission_command + ' ' + task.sh_file 
    else:
        submit_cmd += 'sh ' + task.sh_file 
    submit_cmd += ' && sleep 1s \n'
    task.status = 'submitted'
    os.system(submit_cmd)
    return

def generate_tasks_blasr_pbhoney(settings):
    global all_task_list
    global task_id 

    out_dir = settings.blasr_pbhoney_dir
    align_bam_dir = os.path.join(out_dir, '1_align_bam')
    merge_bam_dir = os.path.join(out_dir, '2_merge_bam')
    sv_call_dir = os.path.join(out_dir, '3_SV_call')
    sh_dir = os.path.join(out_dir, 'sh')
    os.system('mkdir -p %s' % align_bam_dir)
    os.system('mkdir -p %s' % merge_bam_dir)
    os.system('mkdir -p %s' % sv_call_dir)
    os.system('mkdir -p %s' % sh_dir)
    
    sort_bam_list = list()
    merge_bam_dependent_task_list = list()

    for input_file in settings.input_list:
        align_sh_file = os.path.join(sh_dir, 'align.%d.sh' % (task_id))
        sort_sh_file = os.path.join(sh_dir, 'sort.%d.sh' % (task_id))
        index_sh_file = os.path.join(sh_dir, 'index.%d.sh' % (task_id))

        input_file_prefix = get_file_prefix(input_file)
        blasr_sam_file = os.path.join(align_bam_dir, '%s.%d.blasr.sam' % (input_file_prefix, task_id)) 
        tails_sam_file = os.path.join(align_bam_dir, '%s.%d.tails.sam' % (input_file_prefix, task_id))
        tails_bam_file = os.path.join(align_bam_dir, '%s.%d.tails.bam' % (input_file_prefix, task_id))
        tails_sort_bam_file = os.path.join(align_bam_dir, '%s.%d.blasr.tails.sort.bam' % (input_file_prefix, task_id)) 
        sort_bam_list.append(tails_sort_bam_file)
        bam_index_file = tails_sort_bam_file + '.bai'

        align_cmd = blasr_align_cmd(settings, input_file, blasr_sam_file) + ' && ' + tails_align_cmd(settings, blasr_sam_file, tails_sam_file, tails_bam_file) + endl    
        sort_cmd = samtools_sort_cmd(settings, tails_bam_file, tails_sort_bam_file, settings.n_thread) + endl
        index_cmd = samtools_index_cmd(settings.samtools, tails_sort_bam_file) + endl

        align_task = Task(task_id, align_cmd, align_sh_file, tails_bam_file, 'UNK', list())
        task_id += 1
        sort_task = Task(task_id, sort_cmd, sort_sh_file, tails_sort_bam_file, 'UNK', [task_id-1])
        task_id += 1
        index_task = Task(task_id, index_cmd, index_sh_file, bam_index_file, 'UNK', [task_id-1])
        merge_bam_dependent_task_list.append(task_id)
        task_id += 1

        all_task_list.append(align_task)
        all_task_list.append(sort_task)
        all_task_list.append(index_task)

        generate_shell_file_from_task(settings, align_task)
        generate_shell_file_from_task(settings, sort_task)
        generate_shell_file_from_task(settings, index_task)
        

    merge_bam = os.path.join(merge_bam_dir, 'blasr.%s.merge.bam' % settings.sample_name)
    merge_bam_cmd = samtools_merge_bam_cmd(settings.samtools, sort_bam_list, merge_bam)
    merge_sh = os.path.join(sh_dir, 'blasr.%s.merge_bam.sh' % settings.sample_name)
    merge_task = Task(task_id, merge_bam_cmd, merge_sh, merge_bam , 'UNK', merge_bam_dependent_task_list)
    all_task_list.append(merge_task)
    task_id += 1 

    merge_bam_index_file = merge_bam + '.bai'
    index_merge_bam_cmd = samtools_index_cmd(settings.samtools, merge_bam)
    index_merge_bam_sh = os.path.join(sh_dir, 'blasr.%s.index.sh' % settings.sample_name) 
    index_task = Task(task_id, index_merge_bam_cmd, index_merge_bam_sh, merge_bam_index_file, 'UNK', [task_id-1])
    all_task_list.append(index_task)
    index_task_id = task_id
    task_id += 1 

    generate_shell_file_from_task(settings, merge_task) 
    generate_shell_file_from_task(settings, index_task) 

    if settings.enable_PBHoney_Spots:
        svcall_file = os.path.join(sv_call_dir, '%s.spots' % (settings.sample_name))
        settings.spots_file = svcall_file
        svcall_cmd = spots_cmd(settings, merge_bam, svcall_file) + endl 
        svcall_sh = os.path.join(sh_dir, 'blasr.%s.spots.sh'% settings.sample_name)
        svcall_task = Task(task_id, svcall_cmd, svcall_sh, svcall_file, 'UNK', [index_task_id]) 
        all_task_list.append(svcall_task)
        task_id += 1 
        generate_shell_file_from_task(settings, svcall_task) 

    if settings.enable_PBHoney_Tails:
        svcall_file = os.path.join(sv_call_dir, '%s.tails' % (settings.sample_name))
        settings.tails_file = svcall_file
        svcall_cmd = tails_cmd(settings, merge_bam, svcall_file) + endl 
        svcall_sh = os.path.join(sh_dir, 'blasr.%s.tails.sh'% settings.sample_name)
        svcall_task = Task(task_id, svcall_cmd, svcall_sh, svcall_file, 'UNK', [index_task_id]) 
        all_task_list.append(svcall_task)
        generate_shell_file_from_task(settings, svcall_task) 
        task_id += 1 

    #return (align_tasks, sort_tasks, index_tasks, merge_and_call_tasks)
    return

def spots_cmd(settings, input_bam, out_file):

    out_prefix = os.path.splitext(out_file)[0] 
    spots_n_thread = settings.n_thread / 2
    cmd = 'python %s spots --nproc %d  --reference %s --threshold %d --minErrReads %s --consensus %s --output %s %s' % (settings.pbhoney, spots_n_thread, settings.ref_blasr, settings.spots_threshold, settings.spots_minErrReads, settings.spots_consensus, out_prefix, input_bam)
    return cmd

def tails_cmd(settings, input_bam, out_file):
    
    cmd = 'python %s tails --buffer %d --minBreads %d --minZMWs %d --output %s %s' % (settings.pbhoney, settings.tails_buffer, settings.tails_minBreads, settings.tails_minZMWs, out_file, input_bam)   
    return cmd

def tails_align_cmd(settings, input_file, output_sam, ouput_bam):

    cmd = 'python %s pie --nproc %d --output %s %s %s && %s view -hb -@ %d %s > %s' % (settings.pbhoney, settings.n_thread, output_sam, input_file, settings.ref_blasr, settings.samtools, settings.n_thread, output_sam, ouput_bam)
    return cmd

def blasr_align_cmd(settings, input_file, output_file):

    blasr_n_thread = settings.n_thread / 2
    if blasr_n_thread < 1: blasr_n_thread = 1
    cmd = '%s %s %s -sa %s -nproc %d -bestn 1 -sam -clipping subread -out %s' % (settings.blasr, input_file, settings.ref_blasr, settings.ref_sa_blasr, blasr_n_thread, output_file)
    return cmd

def generate_tasks_bwa_sniffles(settings):

    global all_task_list
    global task_id 

    out_dir = settings.bwa_sniffles_dir
    align_bam_dir = os.path.join(out_dir, '1_align_bam')
    merge_bam_dir = os.path.join(out_dir, '2_merge_bam')
    sv_call_dir = os.path.join(out_dir, '3_SV_call')
    sh_dir = os.path.join(out_dir, 'sh')
    os.system('mkdir -p %s' % align_bam_dir)
    os.system('mkdir -p %s' % merge_bam_dir)
    os.system('mkdir -p %s' % sv_call_dir)
    os.system('mkdir -p %s' % sh_dir)

    #align_tasks = list()
    #sort_tasks = list()
    #index_tasks = list()
    sort_bam_list = list()
    merge_bam_dependent_task_list = list()

    for input_file in settings.input_list:

        align_sh_file = os.path.join(sh_dir, 'align.%d.sh' % (task_id))
        sort_sh_file = os.path.join(sh_dir, 'sort.%d.sh' % (task_id))
        index_sh_file = os.path.join(sh_dir, 'index.%d.sh' % (task_id))

        input_file_prefix = get_file_prefix(input_file)
        out_bam_file = os.path.join(align_bam_dir, '%s.%d.bwa.bam' % (input_file_prefix, task_id)) 
        out_sort_bam_file = os.path.join(align_bam_dir, '%s.%d.bwa.sort.bam' % (input_file_prefix, task_id)) 
        bam_index_file = out_sort_bam_file + '.bai'
        sort_bam_list.append(out_sort_bam_file)

        align_cmd = settings.bwa + ' mem -x pacbio -M -t %d %s %s | %s view -bS - > %s\n' % (settings.n_thread, settings.ref_fasta, input_file, settings.samtools, out_bam_file)
        sort_cmd = samtools_sort_cmd(settings, out_bam_file, out_sort_bam_file, settings.n_thread) + endl
        index_cmd = samtools_index_cmd(settings.samtools, out_sort_bam_file) + endl

        align_task = Task(task_id, align_cmd, align_sh_file, out_bam_file, 'UNK', list())
        task_id += 1
        sort_task = Task(task_id, sort_cmd, sort_sh_file, out_sort_bam_file, 'UNK', [task_id-1])
        task_id += 1
        index_task = Task(task_id, index_cmd, index_sh_file, bam_index_file, 'UNK', [task_id-1])
        merge_bam_dependent_task_list.append(task_id)
        task_id += 1

        all_task_list.append(align_task)
        all_task_list.append(sort_task)
        all_task_list.append(index_task)


        generate_shell_file_from_task(settings, align_task)
        generate_shell_file_from_task(settings, sort_task)
        generate_shell_file_from_task(settings, index_task)


    merge_bam = os.path.join(merge_bam_dir, 'bwa.%s.merge.bam' % settings.sample_name)
    merge_bam_cmd = samtools_merge_bam_cmd(settings.samtools, sort_bam_list, merge_bam)
    merge_sh = os.path.join(sh_dir, 'bwa.%s.merge_bam.sh' % settings.sample_name)
    merge_task = Task(task_id, merge_bam_cmd, merge_sh, merge_bam , 'UNK', merge_bam_dependent_task_list)
    all_task_list.append(merge_task)
    task_id += 1

    merge_bam_index_file = merge_bam + '.bai'
    index_merge_bam_cmd = samtools_index_cmd(settings.samtools, merge_bam)
    index_merge_bam_sh = os.path.join(sh_dir, 'bwa.%s.index.sh' % settings.sample_name) 
    index_task = Task(task_id, index_merge_bam_cmd, index_merge_bam_sh, merge_bam_index_file, 'UNK', [task_id-1])
    all_task_list.append(index_task)
    task_id += 1

    svcall_file = os.path.join(sv_call_dir, '%s.bwa.sniffles.vcf' % (settings.sample_name))
    settings.bwa_vcf_file = svcall_file
    svcall_cmd = sniffles_cmd(settings, merge_bam, svcall_file) 
    svcall_sh = os.path.join(sh_dir, 'bwa.%s.sniffles.sh'% settings.sample_name)
    svcall_task = Task(task_id, svcall_cmd, svcall_sh, svcall_file, 'UNK', [task_id-1]) 
    all_task_list.append(svcall_task)
    task_id += 1

    generate_shell_file_from_task(settings, merge_task) 
    generate_shell_file_from_task(settings, index_task) 
    generate_shell_file_from_task(settings, svcall_task) 

    return


def generate_tasks_ngmlr_sniffles(settings):
    global all_task_list
    global task_id 

    out_dir = settings.ngmlr_sniffles_dir
    align_bam_dir = os.path.join(out_dir, '1_align_bam')
    merge_bam_dir = os.path.join(out_dir, '2_merge_bam')
    sv_call_dir = os.path.join(out_dir, '3_SV_call')
    sh_dir = os.path.join(out_dir, 'sh')

    os.system('mkdir -p %s' % align_bam_dir)
    os.system('mkdir -p %s' % merge_bam_dir)
    os.system('mkdir -p %s' % sv_call_dir)
    os.system('mkdir -p %s' % sh_dir)

    sort_bam_list = list()
    merge_bam_dependent_task_list = list()
    for input_file in settings.input_list:

        align_sh_file = os.path.join(sh_dir, 'align.%d.sh' % (task_id))
        sort_sh_file = os.path.join(sh_dir, 'sort.%d.sh' % (task_id))
        index_sh_file = os.path.join(sh_dir, 'index.%d.sh' % (task_id))

        input_file_prefix = get_file_prefix(input_file)
        out_bam_file = os.path.join(align_bam_dir, '%s.%d.ngmlr.bam' % (input_file_prefix, task_id)) 
        out_sort_bam_file = os.path.join(align_bam_dir, '%s.%d.ngmlr.sort.bam' % (input_file_prefix, task_id)) 
        bam_index_file = out_sort_bam_file + '.bai'
        sort_bam_list.append(out_sort_bam_file)

        align_cmd = settings.ngmlr + ' -t %d -r %s -q %s | %s view -bS - > %s\n' % (settings.n_thread, settings.ref_fasta, input_file, settings.samtools, out_bam_file)
        sort_cmd = samtools_sort_cmd(settings, out_bam_file, out_sort_bam_file, settings.n_thread) + endl
        index_cmd = samtools_index_cmd(settings.samtools, out_sort_bam_file) + endl

        align_task = Task(task_id, align_cmd, align_sh_file, out_bam_file, 'UNK', list())
        task_id += 1
        sort_task = Task(task_id, sort_cmd, sort_sh_file, out_sort_bam_file, 'UNK', [task_id-1])
        task_id += 1
        index_task = Task(task_id, index_cmd, index_sh_file, bam_index_file, 'UNK', [task_id-1])
        merge_bam_dependent_task_list.append(task_id)
        task_id += 1

        all_task_list.append(align_task)
        all_task_list.append(sort_task)
        all_task_list.append(index_task)

        generate_shell_file_from_task(settings, align_task)
        generate_shell_file_from_task(settings, sort_task)
        generate_shell_file_from_task(settings, index_task)


    merge_bam = os.path.join(merge_bam_dir, 'ngmlr.%s.merge.bam' % settings.sample_name)
    merge_bam_cmd = samtools_merge_bam_cmd(settings.samtools, sort_bam_list, merge_bam)
    merge_sh = os.path.join(sh_dir, 'ngmlr.%s.merge_bam.sh' % settings.sample_name)
    merge_task = Task(task_id, merge_bam_cmd, merge_sh, merge_bam , 'UNK', merge_bam_dependent_task_list)
    task_id += 1
    all_task_list.append(merge_task)

    merge_bam_index_file = merge_bam + '.bai'
    index_merge_bam_cmd = samtools_index_cmd(settings.samtools, merge_bam)
    index_merge_bam_sh = os.path.join(sh_dir, 'ngmlr.%s.index.sh' % settings.sample_name) 
    index_task = Task(task_id, index_merge_bam_cmd, index_merge_bam_sh, merge_bam_index_file, 'UNK', [task_id-1])
    task_id += 1
    all_task_list.append(index_task)

    svcall_file = os.path.join(sv_call_dir, '%s.ngmlr.sniffles.vcf' % (settings.sample_name))
    settings.ngmlr_vcf_file = svcall_file
    svcall_cmd = sniffles_cmd(settings, merge_bam, svcall_file) 
    svcall_sh = os.path.join(sh_dir, 'ngmlr.%s.sniffles.sh'% settings.sample_name)
    svcall_task = Task(task_id, svcall_cmd, svcall_sh, svcall_file, 'UNK', [task_id-1]) 
    task_id += 1
    all_task_list.append(svcall_task)

    generate_shell_file_from_task(settings, merge_task) 
    generate_shell_file_from_task(settings, index_task) 
    generate_shell_file_from_task(settings, svcall_task) 

    return

def sniffles_cmd(settings, input_bam, out_vcf): 
    
    cmd = '%s -m %s --vcf %s --min_support %d --max_distance %d --threads %d' % (settings.sniffles, input_bam, out_vcf, settings.sniffles_min_support, settings.sniffles_max_distance, settings.n_thread)
    return cmd

def generate_shell_file_from_task(settings, task):

    sh_file = task.sh_file
    sh_fp = open(sh_file, 'w')
    sh_fp.write('#!/bin/bash\n\n')
    sh_fp.write(task.cmd + endl)
    sh_fp.write('echo %d > %s\n' % (settings.runtimekey, task.out_file + '.finished'))
    sh_fp.close()

    return

def samtools_merge_bam_cmd(samtools, sort_bam_list, out_bam):
    cmd = '%s merge %s \\\n' % (samtools, out_bam)
    for sort_bam in sort_bam_list: 
        cmd += '    %s\\\n' % sort_bam
    cmd += '\n\n'
    return cmd
    
def samtools_index_cmd(samtools, sort_bam):
    cmd = '%s index %s' % (samtools, sort_bam)
    return  cmd

def samtools_sort_cmd(settings, input_bam, output_bam, n_thread):

    if settings.samtools_version  == 'new':
        cmd = '%s sort -@ %d -o %s %s' % (settings.samtools, n_thread, output_bam, input_bam) 
    else:
        cmd = '%s sort -@ %d -f %s %s' % (settings.samtools, n_thread, input_bam, output_bam)

    return cmd

def get_file_prefix(input_file):
    file_name = os.path.split(input_file)[1]
    prefix = os.path.splitext(file_name)[0]
    return prefix

def read_input_files(input_file_list):
    input_list = list()
    input_list_fp = open(input_file_list, 'r')
    while 1:
        line = input_list_fp.readline()
        if not line: break
        line = line.strip()
        if line == '': continue
        input_file = os.path.abspath(line)
        input_list.append(input_file)
    input_list_fp.close()
    return input_list

def creat_output_dirs(settings):

    if settings.enable_PBHoney_Spots or settings.enable_PBHoney_Tails:
        settings.blasr_pbhoney_dir = os.path.join(settings.out_dir, 'blasr_pbhoney')
        os.system('mkdir -p %s' % settings.blasr_pbhoney_dir)

    if settings.enable_bwa_Sniffles:
        settings.bwa_sniffles_dir = os.path.join(settings.out_dir, 'bwa_sniffles')
        os.system('mkdir -p %s' % settings.bwa_sniffles_dir)

    if settings.enable_bwa_Sniffles:
        settings.ngmlr_sniffles_dir = os.path.join(settings.out_dir, 'ngmlr_sniffles')
        os.system('mkdir -p %s' % settings.bwa_sniffles_dir)

    settings.nextsv_out_dir = os.path.join(settings.out_dir, 'nextsv_results')
    os.system('mkdir -p %s' % settings.nextsv_out_dir)

    return        

def parse_config_file(config_file):
    
    nextsv_root_dir = os.path.split(__file__)[0] 
    settings = Setting(nextsv_root_dir)
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
        key = line[0].strip()
        value = '='.join(line[1:])
        value = value.strip()
        if (not key) or (not value): continue
        if key == 'input_file_list': 
            settings.input_file_list = os.path.abspath(value)
        elif key == 'out_dir':
            settings.out_dir = os.path.abspath(value)
        elif key == 'n_thread':
            settings.n_thread = int(value)
        elif key == 'mode':
            settings.mode = value
        elif key == 'job_submission_command':
            settings.job_submission_command = value.strip('[').strip(']')
        elif key == 'sample_name':
            settings.sample_name = value
        elif key == 'ref_blasr':
            settings.ref_blasr = os.path.abspath(value)
        elif key == 'ref_sa_blasr':
            settings.ref_sa_blasr = os.path.abspath(value)
        elif key == 'enable_PBHoney_Spots':
            settings.enable_PBHoney_Spots = int(value)
        elif key == 'spots_threshold':
            settings.spots_threshold = int(value)
        elif key == 'spots_minErrReads':
            settings.spots_minErrReads = int(value)
        elif key == 'spots_consensus':
            settings.spots_consensus = value
        elif key == 'enable_PBHoney_Tails':
            settings.enable_PBHoney_Tails = int(value)
        elif key == 'tails_minBreads':
            settings.tails_minBreads = int(value)
        elif key == 'tails_minZMWs':
            settings.tails_minZMWs = int(value)
        elif key == 'tails_buffer':
            settings.talis_buffer = int(value)
        elif key == 'enable_bwa_Sniffles':
            settings.enable_bwa_Sniffles = int(value)
        elif key == 'enable_ngmlr_Sniffles':
            settings.enable_ngmlr_Sniffles = int(value)
        elif key == 'sniffles_min_support':
            settings.sniffles_min_support = int(value)    
        elif key == 'sniffles_max_distance':
            settings.sniffles_max_distance = int(value)
        elif key == 'bwa':
            settings.bwa = os.path.abspath(value)
        elif key == 'ngmlr':
            settings.ngmlr = os.path.abspath(value)
        elif key == 'sniffles':
            settings.sniffles = os.path.abspath(value)
        elif key == 'samtools':
            settings.samtools = os.path.abspath(value)
        elif key == 'bash5tools':
            settings.bash5tools = os.path.abspath(value)
        elif key == 'minLength':
            settings.bash5_minLength = int(value)
        elif key == 'minReadScore':
            settings.bash5_minReadScore = float(value)
        elif key == 'ref_fasta':
            settings.ref_fasta = value
        elif key == 'ref_blasr':
            settings.ref_blasr = value
        elif key == 'ref_sa_blasr':
            settings.ref_sa_blasr = value
        

    config_file_fp.close()
        

    if settings.input_file_list == None:
        myprint ('ERROR! input file list not specified')
        sys.exit()
    if settings.out_dir == None:
        myprint('ERROR! output directory not specified') 
        sys.exit()
    if settings.samtools == None:
        myprint('ERROR! path to samtools not specified') 
        sys.exit()
    if settings.enable_bwa_Sniffles and settings.bwa == None:
        myprint('ERROR! path to bwa not specified') 
        sys.exit()
    if settings.enable_ngmlr_Sniffles and settings.ngmlr == None:
        myprint('ERROR! path to ngmlr not specified') 
        sys.exit()
    if (settings.enable_bwa_Sniffles or settings.enable_ngmlr_Sniffles) and settings.sniffles == None:
        myprint('ERROR! path to sniffles not specified') 
        sys.exit()
    if settings.enable_bwa_Sniffles and settings.ref_fasta == None:
        myprint('ERROR! ref_fasta (reference fasta file for bwa) not specified')
        sys.exit()
    if (settings.enable_PBHoney_Spots or settings.enable_PBHoney_Tails) and settings.ref_blasr == None:
        myprint('ERROR! ref_blasr (reference fasta file for blasr) not specified')
        sys.exit()
    if (settings.enable_PBHoney_Spots or settings.enable_PBHoney_Tails) and settings.ref_sa_blasr == None:
        myprint('ERROR! ref_sa_blasr (suffix array file for blasr) not specified')
        sys.exit()
    if settings.sample_name == None:
        myprint('ERROR! sample name is not specified')
        sys.exit()

    return settings 

def myprint(string):

    print string
    return 

if __name__ == '__main__':
    main()
