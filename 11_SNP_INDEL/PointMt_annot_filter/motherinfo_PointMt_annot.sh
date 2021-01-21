#!/bin/bash
input=$1
tbam=$2
nbam=$3
srcDir=$4
ref=$5

outDir=$(dirname $input)
log=$outDir/$input.annot.log

echo $input > $log
#echo "Starting: Readinfo annot"
#(python $srcDir/readinfo_anno_bwa_190314_Y.py $input $tbam) &>> $log || { c=$?;echo "Error";exit $c; }
#echo "done"
echo "Starting: N readcount annot"
(python $srcDir/readcount_only_anno_190314_clipinsdeladd_Y.py $input $nbam pairMom) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting: Local reassembly annot"
(python $srcDir/read_local_reassembly_RPsize_190311.py $input.readc $tbam $nbam $5) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
rm $input.readc
