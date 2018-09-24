#!/bin/sh
RHOME=$HOME/study/simulation/random
PBIN=$RHOME/bin
FBIN=$RHOME/f_src_2d
PHOME=`pwd`
NDIR=`pwd`
compile='./compile2d.sh'
compile2='./compile2d_est.sh'
est2df90='estpsd2d_v1.f90'
rnd2df90='rnd2d_v1.f90'


echo "Compile again? [y/n] (n) :"
read cyn
test -n "$cyn" && if [ $cyn = "y" -o $cyn = "Y"  ] ; then
    cd $FBIN
    $compile $est2df90 estpsd2d || exit
    $compile $rnd2df90 rnd2d || exit
    cp rnd2d $FBIN
    cp estpsd2d $FBIN
    cd $NDIR
fi
#--------------------------
# Bi-modal or not
#--------------------------
echo "[1] Normal distribution [2] Bi-modal distribution (1)";
read dist_type


echo "Random seed : "
read iseed
echo "Number of Realization : "
read nreal
echo "Main Program to execute : (rnd2d,estpsd2d) "
read exe
test -n $exe && exe=$PBIN/rnd2d
echo $exe
#----------------------------
# Estimate PSD
#----------------------------

for dir in ./*
do
  echo $dir
  if [ -d $dir ]; then
      if [ "${dist_type}" == "2" ] ; then
	  echo "Determine which PSD model to use " ;
	  grep error  $dir/psd.log
	  echo "Now input PSD model number";
	  read mod_num; echo $mod_num
	  ii=`perl -e "printf \"%2.2d\" ,${mod_num}"` && echo $ii
	  cp $dir/psd${ii}sg.bin sg.bin
      fi

#	  cp $dir/sg.bin sg.bin

      cp $dir/param2d.dat ./param2d.dat
#    test -n $exe && exe=rnd2d
      $PBIN/exec2dsub.sh $exe $dir/rnd $iseed $nreal 1 
#      mkdir $dir/psd ; mv $dir/psd*.dat $dir/psd*.su $dir/psd
  fi
done
