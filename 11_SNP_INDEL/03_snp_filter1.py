#!/usr/bin/env python
# coding: utf-8

# In[11]:


#Arg1 = vcf_filename

#if you want to extract clonal variants for clonal organoid data, please remove '#' in front of output_file2.write in the last section

#if mean depth of normal > 60, please change cut off of n_depth (normal_depth) in the case of "input_caller == '01'":

import sys

input_fn = sys.argv[1]
#input_fn = '/home/users/jhyouk/06_mm10_SNUH_radiation/31_2_SNP_updated_190315/mm_study4_SI_sham_SO3_snp_union_2.readinfo.readc.rasmy_PanelofNormal.vcf'

input_file = file(input_fn)
#output_file = file(input_fn.replace('.vcf','.filter1.vcf'),'w')
#output_file1 = file(input_fn.replace('_snp_union_2.readinfo.readc.rasmy_PanelofNormal.vcf','_snp_filtered.vcf'),'w')
#output_file2 = file(input_fn.replace('_snp_union_2.readinfo.readc.rasmy_PanelofNormal.vcf','_snp_filtered_mutalisk.vcf'),'w')
output_file3 = file(input_fn.replace('_snp_union_2.readinfo.readc.rasmy_PanelofNormal.vcf','_snp_filtered_subclonal.vcf'),'w')

input_line = input_file.readline().strip()
prev_chr = '0'
while input_line[0:1] == '#':
    input_split = input_line.split('\t')
    #output_file.write(input_line + '\tpairedN_read\tPON\tT_vaf\tfilter1\tblank\n')
    #output_file1.write(input_line + '\tpairedN_read\tPON\tT_vaf\tfilter1\tblank\n')
    #output_file2.write('\t'.join(input_split[0:10])+'\n')
    output_file3.write('\t'.join(input_split[0:10])+'\n')
    
    input_line = input_file.readline().strip()

while input_line:
    input_split = input_line.split('\t')
    
    #if input_split[0] != prev_chr:
    #    print input_split[0]
    #    prev_chr = input_split[0]
    
    #if input_split[1] != '3487989':
        #input_line = input_file.readline().strip()
        #continue
    
    t_var_cor = input_split[29].split(';')[10]
        
    info = '\tNA\tNA\tNA\tU'
    filter1='F'
    
    if t_var_cor =='.' or t_var_cor == '0':
        info = '\tNA\tNA\tNA\tF'
        #output_file.write(input_line + info + '\n')
        input_line = input_file.readline().strip()       
        continue
    else:
        if input_split[30] == 'NA':
            input_pon = 9
            input_pon_numofsample=9
        else:
            input_pon = round(float(input_split[30].split(';')[6])/100,2)
            input_pon_numofsample = int(input_split[30].split(';')[5])

        input_caller = input_split[6]
        
        
        t_var_cor = float(t_var_cor)
        t_ref_cor = float(input_split[29].split(';')[9])
        
        #t_var_max = max(float(input_split[18]),t_var_cor)
        
        
        t_depth = float(input_split[18]) + float(input_split[17])
        n_depth = int(input_split[28].split(';')[0]) + int(input_split[28].split(';')[1])
        
        if float(input_split[18]) > t_var_cor:
            t_var_max = float(input_split[18])
            t_vaf = round(t_var_max/t_depth,2)
        else:
            t_var_max = t_var_cor
            t_vaf = round(float(t_var_max)/(t_var_cor+t_ref_cor),2)
            
        
        n_read = max(int(input_split[28].split(';')[1]),int(input_split[29].split(';')[8]))
        
        #t_vaf = round(float(t_var_max)/(t_var_cor+t_ref_cor),2)
        
        try:
            var_NM = float(input_split[27])
        except:
            var_NM = 99
        try:
            ref_NM = float(input_split[26])
        except:
            ref_NM = 0
        try:
            t_ref_MQ = float(input_split[19])
        except:
            t_ref_MQ = 0
        try:
            t_var_MQ = float(input_split[20])
        except:
            t_var_MQ = 0
        try:
            ref_clip = float(input_split[28].split(';')[3].split(',')[1])
            ref_ins = float(input_split[28].split(';')[4].split(',')[1])
            ref_del = float(input_split[28].split(';')[5].split(',')[1])
        except:
            ref_clip = 100;ref_ins=100;ref_del=100
        try:
            t_var_clip = float(input_split[23].split(';')[3])
            t_var_ins = float(input_split[24].split(';')[3])
            t_var_del = float(input_split[25].split(';')[3])
        except:
            t_var_clip = 0
            t_var_ins = 0
            t_var_del = 0
            
        if n_read >= 2:
            filter1='F'
        elif input_pon >=0.04 and input_pon_numofsample>=3:
            filter1='F'
        elif n_read == 1 and t_var_max <= min((0.3*t_depth),9):
            filter1='F'
        elif n_read == 0 and t_var_max < 3:
            filter1='F'
        elif var_NM > 2 or ref_NM > 2:
            filter1='F'
        elif t_ref_MQ == 0 or t_var_MQ == 0:
            filter1='F'
        elif ref_clip < 10 and t_var_clip > 70:
            filter1='F'
        elif ref_ins < 10 and t_var_ins > 70:
            filter1='F' 
        elif ref_del < 10 and t_var_del > 70:
            filter1='F'
        elif input_caller == '01':
            if input_pon >=0.01:
                filter1='F'
            elif n_depth >=100:
                filter1='F'
            elif ref_clip < 20 and ref_ins < 10 and ref_del < 10:
                filter1='T'
            else:
                filter1='F'
        elif t_var_clip > 100:
            filter1='F'
        #elif input_pon >=0.04:
            #print input_split[0:2]
        else:
            filter1='T'
                
        info = '\t%s\t%s\t%s\t%s' %(n_read,input_pon,t_vaf,filter1)
        #output_file.write(input_line + info + '\t.\n')
        if filter1 == 'T':
            #output_file1.write(input_line + info + '\t.\n')
            if float(t_vaf) >= 0.3:     #when clonal organoid used
                #output_file2.write('\t'.join(input_split[0:10])+'\n')
                'blank'
            else: #when clonal organoid used
                output_file3.write('\t'.join(input_split[0:10])+'\n') #when clonal organoid used
        
    input_line = input_file.readline().strip()

#output_file.close()
#output_file1.close()
#output_file2.close()
output_file3.close()


# In[ ]:





# In[ ]:




