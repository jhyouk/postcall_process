#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# for ID signature
# python2.7

import os,sys,pysam,gzip,collections
from numpy import median


def ID_repeat(input_line,ref_fasta):
    input_split = input_line.split('\t')
    input_chr = input_split[0]
    input_pos = long(input_split[1])
    
    input_ref = input_split[3]
    input_alt = input_split[4]
    
    if len(input_ref) == 1: # for small insertion
        ru = input_alt[1:]
        rc = 0
        
        while rc<5:
            if ref_fasta.fetch(input_chr,input_pos+len(ru)*rc,input_pos+len(ru)*(rc+1)) == ru:
                rc+=1
                #print ref_fasta.fetch(input_chr,input_pos+len(ru)*rc,input_pos+len(ru)*(rc+1))
            else:
                break
        #print ref_fasta.fetch(input_chr,input_pos-1,input_pos) + 'qwe'
        return ['insertion',str(ru),str(rc),'noMH']

    else: # for deletion    
        if len(input_ref) == 2:
            ru = input_ref[1:]
            rc = 0
            
            while rc<5:
                if ref_fasta.fetch(input_chr,input_pos+len(ru)*(rc+1),input_pos+len(ru)*(rc+2)) == ru:
                    rc+=1
                    #print ref_fasta.fetch(input_chr,input_pos+len(ru)*rc,input_pos+len(ru)*(rc+1))
                else:
                    break
            #print ref_fasta.fetch(input_chr,input_pos-2,input_pos) + 'qwe'
            return ['deletion',str(ru),str(rc+1),'noMH']            
            
            
        else:
            #print input_line
            ru = input_ref[1:]
            rc = 0
            #print ru
            #print ref_fasta.fetch(input_chr,input_pos+len(ru)*(rc+1),input_pos+len(ru)*(rc+2))
            
            if ref_fasta.fetch(input_chr,input_pos+len(ru)*(rc+1),input_pos+len(ru)*(rc+2)) == ru:
                while rc<5:
                    if ref_fasta.fetch(input_chr,input_pos+len(ru)*(rc+1),input_pos+len(ru)*(rc+2)) == ru:
                        rc+=1
                        #print ref_fasta.fetch(input_chr,input_pos+len(ru)*rc,input_pos+len(ru)*(rc+1))
                    else:
                        break
                #print ref_fasta.fetch(input_chr,input_pos-2,input_pos) + 'qwe'
                return ['deletion',str(ru),str(rc+1),'noMH']
            else: # microhomology check!
                rt_mh = 0
                lt_mh = 0
                for i in range(0,len(ru)):
                    if ref_fasta.fetch(input_chr,input_pos+len(ru),input_pos+len(ru)+1+i) == input_ref[1:1+1+i]:
                        rt_mh+=1
                    else:
                        break
                for j in range(0,len(ru)):
                    if ref_fasta.fetch(input_chr,input_pos-1-j,input_pos) == input_ref[-1-j:]:
                        lt_mh+=1
                    else:
                        break                                            
                if rt_mh == 0 & lt_mh ==0:
                    return ['deletion',str(ru),str(rc+1),'noMH']
                else:
                    #print input_line
                    return ['deletion',str(ru),'0',str(max(rt_mh,lt_mh))]
        
        
        
    
#start

file_list = [sys.argv[1]]
for input_fn in file_list:
    print input_fn
    #chr_pos_list=[]
    input_file = open(input_fn)
    output_file = file(input_fn.replace(".vcf",".context_anno.vcf"),"w")
    input_line = input_file.readline().strip()
    while input_line[0:2] =='##':
        input_line = input_file.readline().strip()
    output_file.write(input_line+'\tID_signature_type;repeat_unit;repat_count;microhomology\n')
    input_line = input_file.readline().strip()
    
    mm=0
    while input_line:
        ref_fa='/home/users/jhyouk/99_reference/human/GRCh37/human_g1k_v37.fasta'
        r_file=pysam.FastaFile(ref_fa)
        if len(input_line.split('\t')[3]) > 11 or len(input_line.split('\t')[4]) > 11 :
            'blank'
        else:
            ID_list = ID_repeat(input_line,r_file)
            output_file.write(input_line + '\t' + ';'.join(ID_list) + '\n')

        input_line = input_file.readline().strip()


    output_file.close()
    
print 'THE END'

