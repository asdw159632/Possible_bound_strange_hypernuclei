#!/bin/sh
Dir=$1
hadd result-all.root *.root
mkdir $Dir
mv *.root $Dir
