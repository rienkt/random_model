#!/bin/sh
# This is a sub-program for exec2d.sh
#   Version 1
#     run $exe
#     output files of $exe should be fortran formatted double precision
#


exe=$1
fhead=$2
iseed=$3
nreal=$4

if [ ! -f $exe ] ; 
then
    echo "Main program file does not exist"
    exit
fi

#iseed=$2
#nreal=$3
#exe=$4
#exe=rnd1d_v0
#fhead=$1 #14june_1_
nhist=20
nx=`sed -n "1p" param2d.dat | awk '{print $2}'`; echo $nx
nz=`sed -n "2p" param2d.dat | awk '{print $2}'`; echo $nz
n=$nz
no2=$((nx/2+1)); echo $no2
no2m=$((nx/2-1)); echo $no2m
dt=1000;
dth=1000;
dts=1000;
#rm ${fhead}*.su

# Execute main program
echo $fhead
sed -e "s/SEED/$iseed/" param2d.dat | \
sed -e "s/REAL/$nreal/"  | $exe > ${fhead}.log 2>&1
cp param2d.dat ${fhead}param.dat

ftnstrip < tmp/s0.bin | recast in=double out=float | suaddhead ns=$n | \
    sushw key=dt a=$dts> ${fhead}s0.su


case $5 in
    0)
	niter=20
	i=1
	while [ $i -le $niter ]
	  do
	  ii=`perl -e "printf \"%2.2d\" ,${i}"`
	  #echo $ii
	  ftnstrip < tmp/sg${ii}.bin | recast in=double out=float |  suaddhead ns=$n | sushw key=dt a=$dts> ${fhead}${ii}sg.su 
	  cp tmp/sg${ii}.bin ${fhead}${ii}sg.bin
#	  ftnstrip < tmp/zge${ii}.bin | recast in=double out=float | suaddhead ns=$n | \
#	      sushw key=dt a=$dts> ${fhead}${ii}zge.su
#	  ftnstrip < tmp/zbe${ii}.bin | recast in=double out=float | suaddhead ns=$n | \
#	      sushw key=dt a=$dts> ${fhead}${ii}zbe.su
	  #sed -e "s/D/E/" tmp/sge${ii}.dat  | a2b n1=1 | suaddhead ns=$n | \
	  #    sushw key=dt a=$dts> ${fhead}${ii}sge.su
	  ftnstrip < tmp/sbe${ii}.bin  | recast in=double out=float | suaddhead ns=$n | \
	      sushw key=dt a=$dts> ${fhead}${ii}sbe.su
	  i=$((i+1))
#	  cp tmp/sg${ii}.dat ${fhead}sg${ii}.dat
	done
	;;
    1)
	i=1;
	while [ $i -le $nreal ]
	  do
	  ii=`perl -e "printf \"%5.5d\" ,${i}"`
	  ftnstrip < tmp/zz${ii}.bin | recast in=double out=float | suaddhead ns=$n | \
	      sushw key=dt a=$dt > ${fhead}zz${ii}.su
	  mv tmp/zz${ii}.bin ${fhead}zz${ii}.bin
#	  ftnstrip < tmp/zg${ii}.bin | recast in=double out=float | suaddhead ns=$n | \
#	      sushw key=dt a=$dt > ${fhead}zg${ii}.su
#	  mv tmp/zg${ii}.bin ${fhead}zg${ii}.bin
	  
#	  ftnstrip < tmp/zb${ii}.bin | recast in=double out=float | suaddhead ns=$n | \
#	      sushw key=dt a=$dt > ${fhead}zb${ii}.su



#	  ftnstrip < tmp/sb${ii}.bin | recast in=double out=float | suaddhead ns=$n | \
#	      sushw key=dt a=$dts> ${fhead}sb${ii}.su
	  cat ${fhead}sb${ii}.su >> ${fhead}.su

	  
#	  a2b < tmp/sg${ii}.dat n1=1 | suaddhead ns=$no2 | \
#	      sushw key=dt a=$dts> ${fhead}sg${ii}.su
#	  cat ${fhead}sg${ii}.su >> ${fhead}sg.su
	  
	  i=$((i+1))
	done

	cat ${fhead}s0.su >> ${fhead}sb.su

#	sed -e "s/D/E/"  tmp/zge.dat | a2b n1=1 | suaddhead ns=$n | \
#	    sushw key=dt a=$dts> ${fhead}zge.su
#	sed -e "s/D/E/"  tmp/zbe.dat  | a2b n1=1 | suaddhead ns=$n | \
#	    sushw key=dt a=$dts> ${fhead}zbe.su
#	sed -e "s/D/E/" tmp/sge.dat  | a2b n1=1 | suaddhead ns=$n | \
#	    sushw key=dt a=$dts> ${fhead}sge.su
#	sed -e "s/D/E/" tmp/sbe.dat  | a2b n1=1 | suaddhead ns=$n | \
#	    sushw key=dt a=$dts> ${fhead}sbe.su

	;;
esac
