#!/bin/sh
#
#       
#-------------------
# settings
#-------------------
RHOME=$HOME/study/simulation/random
PBIN=$RHOME/bin
FBIN=$RHOME/f_src_2d
NDIR=`pwd`
compile="$FBIN/compile2d.sh"
est2df90='estpsd2d_v1.f90'
rnd2df90='rnd2d_v1.f90'

#--------------------
# 
echo "Compile again? [y/n] (n) :"
read cyn
test -n "$cyn" && if [ $cyn = "y" -o $cyn = "Y"  ] ; then
    cd $FBIN
    $compile $est2df90 estpsd2d || exit
    $compile $rnd2df90 rnd2d || exit
    cp rnd2d estpsd2d $PBIN
    cd $NDIR
fi

#--------------------------
# Bi-modal or not
#--------------------------
echo "[1] Normal distribution [2] Bi-modal distribution (1)";
read dist_type


#--------------------------
# Set name of output files
#-------------------------

echo "Name of output directory"
read dir
if [ ! -d $dir ]; then
    mkdir $dir
fi

#----------------------------
# Estimate PSD
#----------------------------

if [ "${dist_type}" == "2" ] ; then

echo "Estimate Power Spectral Density ? [y/n] (y)"
read yn
if [ "$yn" != "n" -a "$yn" != "N"  ]; then
    echo "Paramter file name :  (param2d.dat) "
    read $paramdat
    cp $dir/param2d.dat param2d.dat
    test -n "${paramdat}" && cp $paramdat param2d.dat
    echo "Random seed : "
    read iseed
    echo "Number of Realization : "
    read nreal
    echo "Main Program to execute : (rnd2d,estpsd2d) "
    read exe
    test -n $exe && exe=$PBIN/estpsd2d
#    test -n $exe && exe=rnd2d
    $PBIN/exec2dsub.sh $exe $dir/psd $iseed $nreal 0
    mkdir $dir/psd ; mv $dir/psd*.dat $dir/psd*.su $dir/psd
fi

echo "PSD estimation was ended. Continue to generate random field?"
read yn
if [ $yn = "n" -o $yn = "N" ]; then
    exit
fi


echo "Determine which PSD model to use " ;
grep error  $dir/psd.log
echo "Now input PSD model number";
read mod_num; echo $mod_num
ii=`perl -e "printf \"%2.2d\" ,${mod_num}"` && echo $ii
cp $dir/psd${ii}sg.bin sg.bin

fi

echo " Paramter file name :  (param.dat) "
read paramdat
echo $paramdat
cp $dir/sg.bin ./
cp $dir/param2d.dat ./param2d.dat
test -n "${paramdat}" && cp ${paramdat} param2d.dat
echo "Random seed : "
read iseed
echo "Number of Realization : "
read nreal
echo "Main Program to execute : (rnd2d) "
read exe
test -n $exe && exe=$PBIN/rnd2d
echo "Do you want to add something to output flie header $dir/rnd [y/n] (n) ";
if [ "$yn" = "y" -o "$yn" = "Y" ]; then
    echo "input :"
read add 
fi
echo $exe
$PBIN/exec2dsub.sh $exe $dir/rnd$add $iseed $nreal 1
