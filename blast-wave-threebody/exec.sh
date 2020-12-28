#!/bin/bash

#source /storage/fdunphome/szhang/.conenv.sh

inDir=./

outDir=$inDir

target=H_3_qmax10_lmax4_Lmax0_angnum110_nmax5_nstart2_1700-scaled

export EXEH=$inDir/bin/analysis
export jobID=$1
export LogFile=$outDir/log-$jobID.log

$EXEH $outDir $jobID $target

exit
