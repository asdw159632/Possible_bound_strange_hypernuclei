#!/bin/sh
list=*.txt
for i in $list
do
  echo $i
	sed -i 's/{//g' -- $i
	sed -i 's/}//g' -- $i
	sed -i 's/,/ /g' -- $i
	sed -i 's/\*\^/e/g' -- $i
done
