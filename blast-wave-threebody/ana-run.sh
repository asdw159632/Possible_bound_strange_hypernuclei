#!/bin/bash

source /storage/fdunphome/szhang/.conenv.sh

inDir=/storage/fdunphome/szhang/work/cal_waveFunction/blast-wave-proton-deuteron

outDir=$inDir/results

export EXEH=$inDir/bin/analysis
export jobID=$1
export LogFile=$outDir/log-$jobID.log

$EXEH $inDir $jobID > $LogFile

exit
