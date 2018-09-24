#!/bin/sh
NR="nr*"
FFT="fft"
HEAD=`echo $1 | sed -e "s/.f90//"`;
ifort -c ${NR}.f90 ${FFT}.f90 $1
ifort -o $2 ${NR}.o ${FFT}.o ${HEAD}.o

