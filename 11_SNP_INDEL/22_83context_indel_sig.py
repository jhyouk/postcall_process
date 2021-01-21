#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#2nd part of indel_signature
#83 context

#'1bp_del'
#'1bp_ins'
#'2bp_del_repeat'
#'2bp_ins'
#'2bp_del_micohomology'

#start
import sys

file_list = [sys.argv[1]]

for input_fn in file_list:   
    print input_fn
    input_file = file(input_fn)
    output_file = file(input_fn.replace(".vcf",".indel_spectrum.vcf"),"w")
    head_list = ['ID','Dose','d1c1', 'd1c2', 'd1c3', 'd1c4', 'd1c5', 'd1c6', 'd1t0', 'd1t1', 'd1t2', 'd1t3', 'd1t4', 'd1t5', 'i1c1', 'i1c2', 'i1c3', 'i1c4', 'i1c5', 'i1c6', 'i1t0', 'i1t1', 'i1t2', 'i1t3', 'i1t4', 'i1t5', 'd22_1', 'd22_2', 'd22_3', 'd22_4', 'd22_5', 'd22_6', 'd23_1', 'd23_2', 'd23_3', 'd23_4', 'd23_5', 'd23_6', 'd24_1', 'd24_2', 'd24_3', 'd24_4', 'd24_5', 'd24_6', 'd25_1', 'd25_2', 'd25_3', 'd25_4', 'd25_5', 'd25_6', 'i22_0', 'i22_1', 'i22_2', 'i22_3', 'i22_4', 'i22_5', 'i23_0', 'i23_1', 'i23_2', 'i23_3', 'i23_4', 'i23_5', 'i24_0', 'i24_1', 'i24_2', 'i24_3', 'i24_4', 'i24_5', 'i25_0', 'i25_1', 'i25_2', 'i25_3', 'i25_4', 'i25_5', 'mh2_1', 'mh3_1', 'mh3_2', 'mh4_1', 'mh4_2', 'mh4_3', 'mh5_1', 'mh5_2', 'mh5_3', 'mh5_4', 'mh5_5']
    output_file.write('\t'.join(head_list) + "\n")
    
    del_1bp = []
    ins_1bp = []
    del_2bp = []
    ins_2bp = []
    del_mh = []

    for i in range(0,12):
        del_1bp.append(0)
        ins_1bp.append(0)
    for i in range(0,24):    
        del_2bp.append(0)
        ins_2bp.append(0)
    for i in range(0,11):
        del_mh.append(0)

    input_line = input_file.readline().strip()
    input_line = input_file.readline().strip()

    while input_line:
        input_info = input_line.split('\t')[-1] ####
            
        info_split = input_info.split(';')

        if len(info_split[1]) == 1:
            if info_split[0] == 'deletion':
                homopolymer=min(int(info_split[2]),6)
                if info_split[1] == 'C' or info_split[1] == 'G':
                    del_1bp[0+homopolymer-1] += 1          
                else:
                    del_1bp[6+homopolymer-1] += 1  
            else:
                homopolymer=min(int(info_split[2]),5)
                if info_split[1] == 'C' or info_split[1] == 'G':
                    ins_1bp[0+homopolymer]+=1          
                else:
                    ins_1bp[6+homopolymer]+=1
        elif info_split[0] == 'insertion':
            ins_length = min(len(info_split[1]),5)
            homopolymer=min(int(info_split[2]),5)
            ins_2bp[6*(ins_length-2)+homopolymer]+=1

        else:
            if info_split[3] == 'noMH':
                del_length = min(len(info_split[1]),5)
                homopolymer=min(int(info_split[2]),6)        
                del_2bp[6*(del_length-2) + homopolymer-1] += 1
            else:
                del_length = min(len(info_split[1]),5)
                mh = min(int(info_split[3]),5)
                if del_length == 2:
                    del_mh[0] += 1
                elif del_length == 3:
                    del_mh[mh] += 1
                elif del_length == 4:
                    del_mh[3-1+mh] += 1
                else:
                    del_mh[6-1+mh] += 1
        input_line = input_file.readline().strip()

    temp_total = ""
    for i in del_1bp:
        temp_total = temp_total + '\t' + str(i)
    for i in ins_1bp:
        temp_total = temp_total + '\t' + str(i)
    for i in del_2bp:
        temp_total = temp_total + '\t' + str(i)
    for i in ins_2bp:
        temp_total = temp_total + '\t' + str(i)
    for i in del_mh:
        temp_total = temp_total + '\t' + str(i)

    output_file.write(input_fn.split('.indel')[0] +'\t'+ 'NA' + temp_total + '\n')        
    output_file.close()  
    
print 'THE END'

