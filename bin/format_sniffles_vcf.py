#!/usr/bin/env python

import os
import sys

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<in.vcf> <out_prefix> <min_read_supp>'
argc  = 3

class SVCall:

    def __init__(self, attr_list):

        self.chr1, self.pos1, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format, self.genotype_info = attr_list

        self.pos1      = int(self.pos1)

        info_list      = self.info.split(';')
        self.chr2      = ''
        self.pos2      = -1
        self.std1      = -1
        self.std2      = -1
        self.kurt1     = -1000000
        self.kurt2     = -1000000
        self.svtype    = '.'
        self.supp_type = '.'
        self.sv_len    = 0 
        self.reported_sv_len = 0 
        self.strands   = ''
        self.read_supp = 0

        for info in info_list:

            info = info.split('=')
            if len(info) != 2: continue
            key, value = info
            if key == 'CHR2':
                self.chr2 = value
            elif key == 'END':
                self.pos2 = int(value)
            elif key == 'STD_quant_start':
                self.std1 = float(value)
            elif key == 'STD_quant_stop':
                self.std2 = float(value)
            elif key == 'Kurtosis_quant_start':
                self.kurt1 = float(value)
            elif key == 'Kurtosis_quant_stop':
                self.kurt2 = float(value)
            elif key == 'SVTYPE':
                self.svtype = value
            elif key == 'SUPTYPE':
                self.supp_type = value
            elif key == 'SVLEN':
                self.reported_sv_len = int(value)
            elif key == 'STRANDS':
                self.strands = value
            elif key == 'RE':
                self.read_supp = int(value)

        if self.chr1 == self.chr2:
            self.sv_len = abs(self.reported_sv_len)
        else:
            self.sv_len = -1
         
        a = list()
        if '[' in self.alt:
            a = self.alt.split('[')
        elif ']' in self.alt:
            a = self.alt.split(']')

        if len(a) == 3:
            bk2 = a[1]
            self.chr2, self.pos2 = a[1].split(':')
            self.pos2 = int(self.pos2)



    def output(self):

        if self.chr1 == self.chr2:
            outstring = '%s\t%d\t%d\t%s\t%s\t'     %  (self.chr1, self.pos1, self.pos2, self.filter, self.svtype)
        else:
            outstring = '%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t' %  (self.chr1, self.pos1-1, self.pos1, self.chr2, self.pos2-1, self.pos2, self.filter, self.svtype)

        outstring += '%d\t%d\tsupp_type=%s;std1=%.4f;std2=%.4f;kurt1=%.4f;kurt2=%.4f;strand=%s' % (self.sv_len, self.read_supp, self.supp_type, self.std1, self.std2, self.kurt1, self.kurt2, self.strands)

        return outstring

def main():

    if len(arg) < argc:
        print usage
        sys.exit()

    in_vcf_file = arg.pop(0)
    out_prefix  = arg.pop(0)
    min_read_supp = int(arg.pop(0))
      

    primary_chr_name = set()
    for i in range(1, 22):
        primary_chr_name.add(str(i))
        primary_chr_name.add('chr' + str(i))

    primary_chr_name.add('X')
    primary_chr_name.add('Y')
    primary_chr_name.add('M')
    primary_chr_name.add('MT')
    primary_chr_name.add('chrX')
    primary_chr_name.add('chrY')
    primary_chr_name.add('chrM')
    primary_chr_name.add('chrMT')

    max_sv_len = int(2e7)
    in_vcf_fp = open(in_vcf_file, 'r')

    out_ins_file = out_prefix + '.ins.bed'
    out_del_file = out_prefix + '.del.bed'
    out_dup_file = out_prefix + '.dup.bed'
    out_inv_file = out_prefix + '.inv.bed'
    out_tra_file = out_prefix + '.tra.bedpe'
    out_complex_file = out_prefix + '.complex.bed'

    out_ins_fp = open(out_ins_file, 'w')
    out_del_fp = open(out_del_file, 'w')
    out_dup_fp = open(out_dup_file, 'w')
    out_inv_fp = open(out_inv_file, 'w')
    out_tra_fp = open(out_tra_file, 'w')
    out_complex_fp = open(out_complex_file, 'w')

    while 1:
        line = in_vcf_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        sv = SVCall(line)
        if (sv.chr1 not in primary_chr_name) or (sv.chr2 not in primary_chr_name): continue
        
        ### filtering calls ###
        if sv.sv_len > max_sv_len: continue 
        if sv.read_supp < min_read_supp: continue
        if sv.chr1 != sv.chr2:
            out_tra_fp.write(sv.output() + endl)
        elif sv.svtype == 'INS':
            out_ins_fp.write(sv.output() + endl)
        elif sv.svtype == 'DEL':
            out_del_fp.write(sv.output() + endl)
        elif sv.svtype == 'DUP':
            out_dup_fp.write(sv.output() + endl)
        elif sv.svtype == 'INV':
            out_inv_fp.write(sv.output() + endl)
        else:
            out_complex_fp.write(sv.output() + endl)


    in_vcf_fp.close()
    out_ins_fp.close()
    out_del_fp.close()
    out_dup_fp.close()
    out_inv_fp.close()
    out_tra_fp.close()
    out_complex_fp.close()




if __name__ == '__main__':
    main()
