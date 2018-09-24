#!/bin/sh

if [ ! -f param.dat ];
then
    echo " param.dat cannot be found.";
    exit;
fi

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
n=`sed -n "1p" param.dat | awk '{print $2}'`; echo $n
no2=$((n/2+1)); echo $no2
no2m=$((n/2-1)); echo $no2m
dt=1000;
dth=1000;
dts=1000;
#rm ${fhead}*.su

# Execute main program
sed -e "s/SEED/$iseed/" param.dat | \
sed -e "s/REAL/$nreal/"  | $exe > ${fhead}.log
cp param.dat ${fhead}param.dat

a2b < tmp/s0.dat  n1=1 | suaddhead ns=$no2 | \
    sushw key=dt a=$dts> ${fhead}s0.su


case $5 in
    0)
	niter=19
	i=1
	while [ $i -le $niter ]
	  do
	  ii=`perl -e "printf \"%2.2d\" ,${i}"`

	  sed -e "s/D/E/" tmp/sg${ii}.dat | a2b  n1=1 | suaddhead ns=$n | \
	      sushw key=dt a=$dts> ${fhead}${ii}sg.su
	  cp tmp/sg${ii}.dat ${fhead}${ii}sg.dat
	  sed -e "s/D/E/"  tmp/zge${ii}.dat | a2b n1=1 | suaddhead ns=$n | \
	      sushw key=dt a=$dts> ${fhead}${ii}zge.su
	  sed -e "s/D/E/"  tmp/zbe${ii}.dat  | a2b n1=1 | suaddhead ns=$n | \
	      sushw key=dt a=$dts> ${fhead}${ii}zbe.su
	  sed -e "s/D/E/" tmp/sge${ii}.dat  | a2b n1=1 | suaddhead ns=$no2 | \
	      sushw key=dt a=$dts> ${fhead}${ii}sge.su
	  sed -e "s/D/E/" tmp/sbe${ii}.dat  | a2b n1=1 | suaddhead ns=$no2 | \
	      sushw key=dt a=$dts> ${fhead}${ii}sbe.su
	  cp tmp/sg${ii}.dat ${fhead}psdsg${ii}.dat
	  cp tmp/zg${ii}.dat ${fhead}zg${ii}.dat
	  cp tmp/zb${ii}.dat ${fhead}zb${ii}.dat
	  cp tmp/sge${ii}.dat ${fhead}sge${ii}.dat
	  cp tmp/sbe${ii}.dat ${fhead}sbe${ii}.dat
	  i=$((i+1))
	done
	;;
    1)
	i=1;
	while [ $i -le $nreal ]
	  do
	  ii=`perl -e "printf \"%5.5d\" ,${i}"`
	  sed -e "s/D/E/"  tmp/zg${ii}.dat | a2b n1=1 | suaddhead ns=$n | \
	      sushw key=dt a=$dt > ${fhead}zg${ii}.su
	  sed -e "s/D/E/" tmp/zb${ii}.dat | a2b  n1=1 | suaddhead ns=$n | \
	      sushw key=dt a=$dt > ${fhead}zb${ii}.su
	  sed -e "s/D/E/" tmp/sb${ii}.dat | a2b  n1=1 | suaddhead ns=$no2 | \
	      sushw key=dt a=$dts> ${fhead}sb${ii}.su
	  cat ${fhead}sb${ii}.su >> ${fhead}.su
	  cp tmp/zg${ii}.dat ${fhead}zg${ii}.dat
	  cp tmp/zb${ii}.dat ${fhead}zb${ii}.dat

#	  a2b < tmp/sg${ii}.dat n1=1 | suaddhead ns=$no2 | \
#	      sushw key=dt a=$dts> ${fhead}sg${ii}.su
#	  cat ${fhead}sg${ii}.su >> ${fhead}sg.su
	  
	  i=$((i+1))
	done

	cat ${fhead}s0.su >> ${fhead}sb.su

	sed -e "s/D/E/"  tmp/zge.dat | a2b n1=1 | suaddhead ns=$n | \
	    sushw key=dt a=$dts> ${fhead}zge.su
	sed -e "s/D/E/"  tmp/zbe.dat  | a2b n1=1 | suaddhead ns=$n | \
	    sushw key=dt a=$dts> ${fhead}zbe.su
	sed -e "s/D/E/" tmp/sge.dat  | a2b n1=1 | suaddhead ns=$no2 | \
	    sushw key=dt a=$dts> ${fhead}sge.su
	sed -e "s/D/E/" tmp/sbe.dat  | a2b n1=1 | suaddhead ns=$no2 | \
	    sushw key=dt a=$dts> ${fhead}sbe.su

	;;
esac
