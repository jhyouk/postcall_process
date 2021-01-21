#!/usr/bin/env python
# coding: utf-8

# In[2]:


#Arg1 = vcf_filename

#if mean depth of normal > 60, please change cut off of n_depth (normal_depth) in the case of "input_caller == '01'":

import sys, math

input_fn = sys.argv[1]
#input_fn = 'B3-0_2Gy-3_indel_union_2.readinfo.readc.rasmy.vcf.pon'

input_file = file(input_fn)
output_file = file(input_fn.replace('.vcf','.filter1.vcf'),'w')
output_file1 = file(input_fn.replace('_indel_union_2.readinfo.readc.rasmy_PanelofNormal.vcf','_indel_filtered_strict.vcf'),'w')


input_line = input_file.readline().strip()
prev_chr = '0'

while input_line[0:1] == '#':
    output_file.write(input_line + '\tpairedN_read\tPON\tT_vaf\tfilter1\tread_ref;read_var;PON_ref;PON_var;log10OR\n')
    output_file1.write(input_line + '\tpairedN_read\tPON\tT_vaf\tfilter1\tread_ref;read_var;PON_ref;PON_var;log10OR\n')
    
    input_line = input_file.readline().strip()

while input_line:
    input_split = input_line.split('\t')
    
    #if input_split[0] != prev_chr:
    #    print input_split[0]
    #    prev_chr = input_split[0]
    
    #if input_split[1] != '3318865':
        #input_line = input_file.readline().strip()
        #continue
    
        
    t_var_cor = input_split[29].split(';')[10]
    info = '\tNA\tNA\tNA\tU'
    filter1='F'
    
    if t_var_cor =='.' or t_var_cor == '0':
        info = '\tNA\tNA\tNA\tF'
        output_file.write(input_line + info + '\n')
        input_line = input_file.readline().strip()       
        continue
    else:
        if input_split[30] == 'NA':
            input_pon = 9
            input_pon_numofsample = 9
        else:
            input_pon = round(float(input_split[30].split(';')[4])/100,2)
            input_pon_numofsample = int(input_split[30].split(';')[3])

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
        
        #t_vaf = round(float(t_var_max)/(float(t_var_cor)+t_ref_cor),2)
        
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
        #if input_split[1] == '14492774':
        #    print n_read, t_var_max, t_depth
        
        if n_read >= 2:
            filter1='F'
        elif n_read == 1 and t_var_max <= min((0.3*t_depth),9):
            filter1='F'
        elif n_read == 0 and t_var_max < 3:
            filter1='F'
        elif var_NM > 2 or ref_NM > 2 :
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
            if n_depth >= 100:
                filter1='F'
            elif ref_clip < 20 and ref_ins < 10 and ref_del < 10:
                filter1='T'
            else:
                filter1='F'
        elif t_var_clip > 100:
            filter1='F'
        else:
            filter1='T'
            
        input_PON_ref = 0
        input_PON_var = 0
        if input_split[30].split(';')[0] == 'NA':
            'blank'
        else:
            input_PON_ref = input_PON_ref + int(input_split[30].split(';')[0])
            input_PON_var = input_PON_var + int(input_split[30].split(';')[2])        

        if input_PON_ref == 0:
            input_PON_ref = 9999

        input_OR = (float(input_PON_var)*float(input_split[29].split(';')[9])) / (float(input_PON_ref)*float(input_split[29].split(';')[10]))
                
        info = '\t%s\t%s\t%s\t%s' %(n_read,input_pon,t_vaf,filter1)
        #output_file.write(input_line + info + '\n')
        read_var = int(input_split[29].split(';')[10])
        logOR = float(round(math.log10(input_OR+0.000001),3))
        output_file.write(input_line + info + "\t%s;%s;%s;%s;%s\n" %(input_split[29].split(';')[9],read_var,input_PON_ref,input_PON_var,logOR))

        if filter1 == 'T':        
            if read_var > 2 and logOR < -3:            
                output_file1.write(input_line + info + "\t%s;%s;%s;%s;%s\n" %(input_split[29].split(';')[9],input_split[29].split(';')[10],input_PON_ref,input_PON_var,round(math.log10(input_OR+0.000001),3)))
            
    input_line = input_file.readline().strip()

output_file.close()
output_file1.close()


# In[ ]:





# In[ ]:




