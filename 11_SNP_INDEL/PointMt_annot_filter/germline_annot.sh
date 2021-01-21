#!/bin/bash
input=$1
tbam=$2
srcDir=$3

outDir=$(dirname $input)
log=$outDir/$input.annot.log

echo $input > $log
echo "Starting: Readinfo annot"
(python $srcDir/readinfo_anno_bwa_190314_Y.py $input $tbam) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
