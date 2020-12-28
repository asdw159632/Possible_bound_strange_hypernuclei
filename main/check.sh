#!/bin/sh
cd bin
valgrind --leak-check=full ./analysis
cd ..
exit 0
