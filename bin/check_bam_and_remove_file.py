#!/usr/bin/env python

import os
import sys
import subprocess 

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<output.bam> <to_be_removed_file> <samtools>'
argc  = 2 

def main():
    if len(arg) < argc:
        print (usage)
        sys.exit()

    out_bam_file = os.path.abspath(arg.pop(0))
    to_be_removed_file = os.path.abspath(arg.pop(0))
    samtools = os.path.abspath(arg.pop(0))

    ret = subprocess.check_output([samtools, "quickcheck", "-vvv", out_bam_file], stderr=subprocess.STDOUT).decode()
    if 'good EOF' in ret: 
        os.remove(to_be_removed_file)

    return 



if __name__ == '__main__':
    main()
