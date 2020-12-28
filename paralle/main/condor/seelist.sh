#!/bin/sh
for ((i=0;i<20;i++))
do
sed -n "138p" job_$i.out
done
