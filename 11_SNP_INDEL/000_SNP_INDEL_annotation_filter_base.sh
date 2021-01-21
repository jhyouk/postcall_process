#!/bin/bash
set -e

sampleID="$1"
mttype="$2" # snv or indel
reference="$3" # mm10 or hg19
somaticbam="$4"
germlinebam="$5"
panelofnormal="$6"

case $3 in
	mm10)	species="mouse";;
	hg19)	species="human";;
esac
case $2 in
	snp)	temp_mttype="snv";;
	indel)	temp_mttype="indel";;
esac


outDir=$(dirname $1)
log=$outDir/$1.$2.annot.log

echo $1 $2 $3 $4 $5 > $log
if false; then
echo varscan somatic filtering >> $log
(python /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/11_universe_filter/00_varscan_somaticfilter.py ../05_varscan/$1.varscan.$2.vcf) &>> $log || { c=$?;echo "Error";exit $c; }
echo done >> $log
echo "start: union of pass call in varscan2 somatic & strelka2" >>$log
(python /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/00_vcf_combination_by_Youk_$2.py $1 $species 2 ../02_Strelka2/$1/results/variants/somatic."$temp_mttype"s.vcf.gz ../05_varscan/$1.varscan.$2.somatic.vcf) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done" >>$log
echo "initial annotation" >> $log
(sh /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/PointMt_annot_filter/PointMt_annot.sh $1_$2_union_2.vcf $4 $5 /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/PointMt_annot_filter/src $3) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done" >>$log
echo "panel of normal annotation" >> $log
(python /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/02_AddNpanelToVCF_"$2".py $1_$2_union_2.readinfo.readc.rasmy.vcf $6 PanelofNormal $species) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done" >>$log
fi
echo "filter1 using sample and germline information"  >>$log
(python /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/03_$2_filter1.py $1_$2_union_2.readinfo.readc.rasmy_PanelofNormal.vcf) &>> $log || { c=$?;echo "Error";exit $c; }
echo "filter1 done" >> $log
if false; then
echo "Annovar and functional mutation selection"  >>$log
(python /home/users/jhyouk/82_post_calling_process_JYouk/11_SNP_INDEL/04_"$mttype"_annovar.py "$sampleID"_"$mttype"_filtered.vcf) &>> $log || { c=$?;echo "Error";exit $c; }
echo "annovar done" >> $log
echo "remove intermediate files" >> $log
#rm $1_$2_union_2.vcf
#rm $1_$2_union_2.readinfo.readc.rasmy.vcf
#rm $1_$2_union_2.readinfo.readc.rasmy_PanelofNormal.vcf
echo "Finish!!" >> $log
fi
