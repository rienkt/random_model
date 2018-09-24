#!/bin/sh
#----------------
# settings
#----------------
PHOME=$HOME/study/simulation/random
PBIN=$PHOME/bin
FDIR=$PHOME/f_src_1d
# compiling script 
compile="$FDIR/compile1d.sh"
# fortran programs
est1df90='estpsd1d_v0.f90' # estimating input spectral
rnd1df90='rnd1d_v2.f90'    # creating random field

#---------------
# compiling
#---------------
echo "Compile again? [y/n] (n) :"
read cyn
test -n "$cyn" && if [ $cyn = "y" -o $cyn = "Y"  ] ; then
    cd $FDIR
    $compile $est1df90 estpsd1d || exit
    $compile $rnd1df90 rnd1d || exit
    cp estpsd1d rnd1d $PDIR
    cd $PHOME
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
    echo "Paramter file name :  (param.dat) "
    read $paramdat
    cp $dir/param.dat ./
    test -n "$paramdat" && cp $paramdat param.dat
    echo "Random seed : "
    read iseed
    echo "Number of Realization : "
    read nreal
#    echo "Main Program to execute : (estpsd1d) "
#    read exe
    test -n $exe && exe=$PBIN/estpsd1d
    $PDIR/exec1dsub.sh $exe $dir/psd $iseed $nreal 0
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
cp $dir/psd/psd${ii}sg.dat sg.dat

fi

echo " Paramter file name :  (param.dat) "
read paramdat
cp $dir/param.dat ./
test -n "$paramdat" && cp $paramdat param.dat
echo "Random seed : "
read iseed
echo "Number of Realization : "
read nreal
#echo "Main Program to execute : (estpsd1d) "
#read exe
test -n $exe && exe=$PBIN/rnd1d
echo "Do you want to add something to output flie header $dir/rnd [y/n] (n) ";
if [ "$yn" = "y" -o "$yn" = "Y" ]; then
    echo "input :"
    read add 
fi
$PDIR/exec1dsub.sh $exe $dir/rnd$add $iseed $nreal 1


