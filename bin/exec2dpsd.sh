#!/bin/sh
RHOME=$HOME/study/simulation/random
PBIN=$RHOME/bin
FDIR=$RHOME/f_src_2d
NDIR=`pwd`
compile="$FDIR/compile2d.sh"
est2df90='estpsd2d_v1.f90'
rnd2df90='rnd2d_v1.f90'


echo "Compile again? [y/n] (n) :"
read cyn
test -n "$cyn" && if [ $cyn = "y" -o $cyn = "Y"  ] ; then
    cd $FDIR
    $compile $est2df90 estpsd2d || exit
    $compile $rnd2df90 rnd2d || exit
    cp rnd2d $PBIN
    cp estpsd2d $PBIN
    cd $NDIR
fi

echo "Random seed : "
read iseed
echo "Number of Realization : "
read nreal
echo "Main Program to execute : (rnd2d,estpsd2d) "
read exe
test -n $exe && exe=$PBIN/estpsd2d
echo $exe
#----------------------------
# Estimate PSD
#----------------------------

for dir in ./*
do
  echo $dir
  if [ -d $dir ]; then
      cp $dir/param2d.dat ./param2d.dat
#    test -n $exe && exe=rnd2d
      $PBIN/exec2dsub.sh $exe $dir/psd $iseed $nreal 0
      mkdir $dir/psd ; mv $dir/psd*.dat $dir/psd*.su $dir/psd
  fi
done
