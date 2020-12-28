#!/bin/bash
#echo "ok?"

rootevn=$HOME"/.conenv.sh"

export in_dir=$PWD

rm -rf $in_dir/results
mkdir -p $in_dir/results

## for stOut
rm -rf $in_dir/stOut
mkdir -p $in_dir/stOut
stOut=$in_dir/stOut

##compile code
#cd $in_dir
#unpstar.sh src.des3
make clean;make
#rm -rf src

##calculate list number
#find ./ -name $1.con | xargs perl -pi -e 's|Queue .*|Queue '$i'|g'
find ./ -name $1.con | xargs perl -pi -e 's|Executable .*|Executable = '$1'.sh|g'
find ./ -name $1.con | xargs perl -pi -e 's|Initialdir .*|Initialdir = '$in_dir'|g'

##modify data dir in script
find ./ -name $1.sh | xargs perl -pi -e 's|source.*.sh|source '$rootevn'|g'
find ./ -name $1.sh | xargs perl -pi -e 's|HOME.*|HOME='$HOME'|g'
find ./ -name $1.sh | xargs perl -pi -e 's|inDir=.*|inDir='$in_dir'|g'

##submit the job
condor_submit $1.con

exit
