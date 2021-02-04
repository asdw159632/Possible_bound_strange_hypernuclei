#!/bin/bash

#source /storage/fdunphome/szhang/.conenv.sh

inDir=./

outDir=$inDir

#target=H_3_qmax10_lmax4_Lmax0_angnum110_nmax5_nstart2_1700-scaled
#target=pnOmega_qmax10_lmax4_Lmax0_angnum33_nmax5_nstart2_2700-scaled
target_dibaryon=OmegaOmega-scaled
target_tri_body=OO-n-scaled

echo $target_dibaryon
echo $target_tri_body

export EXEH=$inDir/bin/analysis
export jobID=$1
export LogFile=$outDir/log-$jobID.log

$EXEH $outDir $jobID $target_dibaryon $target_tri_body

exit
