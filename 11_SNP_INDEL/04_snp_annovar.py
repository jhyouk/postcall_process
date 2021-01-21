#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# SNP ANNOVAR

import sys, os
input_fn = sys.argv[1]
os.system ('python /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/annovar_script/annotate_using_annovar_v2.py -s -i "refGene,1000g2015aug_all,exac03,snp142,clinvar_20180603,dbnsfp35a,dbscsnv11,cosmic86_coding" %s' % input_fn)

input_file = open(input_fn + '.anv')
input_line = input_file.readline().strip() #input_line 은 input_file 의 첫 줄
output_file = open(input_fn.replace(".vcf",".anv.functional.vcf"),'w')
output_file.write(input_line+'\n') #첫 줄은 header 라서 그냥 써주기 (앞에서 .strip()으로 '\n'을 빼줬으니까 여기서 안넣어주면 계속 옆에 길게 붙는다.)
input_line = input_file.readline().strip() #input_line 은 input_file 의 그 다음줄

while input_line: #둘째줄부터는 아래 조건에 맞는 행만 추출
    input_split = input_line.split('\t')
    #print(input_split)
    #break
    if input_split[35+1] == 'exonic' and input_split[38+1] =='nonsynonymous SNV':
        output_file.write(input_line+'\n')
    else:
        if input_split[35+1] == 'exonic' and input_split[38+1] =='stopgain':
            output_file.write(input_line+'\n')
        else:
            if input_split[35+1] == 'splicing':
                output_file.write(input_line+'\n')
            else:
                'blank'
    input_line = input_file.readline().strip()

output_file.close()

