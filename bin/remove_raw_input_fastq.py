#!/usr/bin/env python

import os
import sys
import subprocess 

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<filtered.fq> <raw_input_list>'
argc  = 2 

def main():
    if len(arg) < argc:
        print (usage)
        sys.exit()

    filtered_fq_file = os.path.abspath(arg.pop(0))
    raw_input_list_file = os.path.abspath(arg.pop(0))

    raw_input_list_fp = open(raw_input_list_file, 'r')
    lines = list(raw_input_list_fp)
    raw_input_list_fp.close() 

    total_file_size = 0
    input_file_list = list()

    for line in lines:
        in_file = line.strip()
        input_file_list.append(in_file)
        total_file_size += os.path.getsize(in_file)

    filtered_fq_file_size = os.path.getsize(filtered_fq_file)

    if filtered_fq_file_size > 0.9 * total_file_size:
        for in_file in input_file_list:
            os.remove(in_file)
        
    return 



if __name__ == '__main__':
    main()
