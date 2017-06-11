#!/usr/bin/env python

import os
import sys

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]
arg.reverse()

usage = 'python ' + __file__ + ' <bam.list> <out_bam> <n_thread>'
argc = 3 

samtools = 'time $HOME/bin/samtools'

def main():
    if len(arg) < argc:
        print usage
        sys.exit()
    
    bam_list_file = arg.pop()
    out_bam = arg.pop()
    n_thread = arg.pop()
    bam_list_fp = open(bam_list_file, 'r')
    out_dir = os.path.split(out_bam)[0]
    os.system('mkdir -p ' + out_dir)
    
    bam_list = list()
    bam_header_list = list()
    merge_header = ''
    print 'parse sam headers...'
    while 1:
        line =  bam_list_fp.readline()
        if not line:
            break
        bam = line.strip()
        bam_list.append(bam)
        bam_name = os.path.split(bam)[1]
        bam_header = os.path.join(out_dir, bam_name + '.header')
        bam_header_list.append(bam_header)
        cmd = samtools + ' view -H ' + bam + ' > ' + bam_header
        os.system(cmd)
    
    
    header1_fp = open(bam_header_list[0], 'r')
    header1_list = list(header1_fp)
    for line in header1_list:
        merge_header += line
    header1_fp.close()

    for i in range(1, len(bam_header_list)):
        header_fp = open(bam_header_list[i])
        header_list = list(header_fp)
        for line in header1_list:
            if line[0:3] != '@PG':
                continue
            line = line.strip().split(tab)
            for j in range(1, len(line)):
                if ':' in line[j]:
                    key, value = line[j].split(':')
                    if key == 'ID':
                        value = value + '-' + str(i+1)
                        line[j] = 'ID:' + value
                        break
            PG = tab.join(line) + endl
            merge_header += PG
        header_fp.close()

    merge_header_file = out_bam + '.header.sam'
    merge_header_fp = open(merge_header_file, 'w')
    merge_header_fp.write(merge_header)
    merge_header_fp.close()
    
    print 'concatenating input bam files...'
    cat_bam = out_bam + '.cat.bam'
    cmd = samtools + ' cat -h ' + merge_header_file + ' -o ' + cat_bam 
    for inbam in bam_list:
        cmd += ' ' + inbam
    print cmd + endl
    os.system(cmd)
    print 'sorting and indexing the concatenated bam file...'
    cmd = samtools + ' sort -@ ' + n_thread + ' -o ' + out_bam + ' ' + cat_bam
    print cmd + endl
    os.system(cmd)
    cmd = samtools + ' index ' + out_bam
    print cmd + endl
    os.system(cmd)

    ##### remove temp files #####
    
    cmd  = 'rm ' + merge_header_file + endl 
    cmd += 'rm ' + cat_bam + endl
    for bam_header_file in bam_header_list:
        cmd += 'rm ' + bam_header_file + endl
    print cmd
    os.system(cmd)



if __name__ == '__main__':
    main()
