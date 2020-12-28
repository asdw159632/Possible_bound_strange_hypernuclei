#! /bin/bash
#
#cp jacobi_polynomial.hpp /$MYINC
#
g++ -c -Wall -I ./ jacobi_polynomial.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  #exit
fi
#
#mv jacobi_polynomial.o ~/mylib/jacobi_polynomial.o
#
echo "Normal end of execution."
