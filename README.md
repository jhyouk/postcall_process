sampleID="$1"

mttype="$2" # snv or indel

reference="$3" # mm10 or hg19

somaticbam="$4"

germlinebam="$5"

panelofnormal="$6"

sh /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/000_SNP_INDEL_annotation_filter_base.sh JHY-COL-HC09-06_blood-wgs-ILLUMINA.fmarked snp hg19 /home/users/team_projects/colon_LINE1/02_bam/JHY-COL-HC09-06-wgs-ILLUMINA.fmarked.bam /home/users/team_projects/colon_LINE1/02_bam/HC09_germline_30x.s.md.bam /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/korean36.36s.q0Q0.chr1.mpileup.snv.edit.gz

sh /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/000_SNP_INDEL_annotation_filter_base.sh JHY-COL-HC09-06_blood-wgs-ILLUMINA.fmarked indel hg19 /home/users/team_projects/colon_LINE1/02_bam/JHY-COL-HC09-06-wgs-ILLUMINA.fmarked.bam /home/users/team_projects/colon_LINE1/02_bam/HC09_germline_30x.s.md.bam /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/korean36.36s.q0Q0.chr1.mpileup.indel.edit.gz
