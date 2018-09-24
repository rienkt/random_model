#!/bin/bash
# Settings
F90=ifort
OPT_o=        # building options
OPT_c=        # compile options


# Mathmatical Libraries
#                    from Numerical Recipes, Tokyo-Univ, MT
nrlib='nrtype.f90 nr.f90 nrutil.f90' ; nrlibdir=..\\/f_lib_num_rec
sort='sort.f90'              ; sortdir=$nrlibdir
fft='modfft4g.f90 fft4g.f'     ; fftdir=..\\/f_lib_fft\\/fft1d
rnd='mt19937_mod.f90'        ; rnddir=..\\/f_lib_mt
stat='static.f90'            ; statdir=$nrlibdir
spfunc='sp_func.f90'         ; spfuncdir=$nrlibdir
fndroot='fnd_root.f90'       ; fndrootdir=$nrlibdir

# My modules
acf=acf.f90                  ; acfdir=../f_src_acf
head=header1d.f90
sub=subrnd1d.f90

# main program
main=$1
exe=$2


echo "compile libralies again? (n)"
read cyn
if [ $cyn = "y" ]; then
nrlibf=`eval echo $nrlibdir/$nrlib | sed -e "s/ / $nrlibdir\//g"`
fftf=`eval echo $fftdir/$fft | sed -e "s/ / $fftdir\//g"`
rndf=`eval echo $rnddir/$rnd | sed -e "s/ / $rnddir\//g"`
statf=`eval echo $statdir/$stat | sed -e "s/ / $statdir\//g"`
spfuncf=`eval echo $spfuncdir/$spfunc |  sed -e "s/ / $spfuncdir\//g"`
fndrootf=`eval echo $fndrootdir/$fndroot | sed -e "s/ / $fndrootdir\//g"`
sortf=`eval echo $sortdir/$sort | sed -e "s/ / $sortdir\//g"`

   $F90 $OPT_c -c $nrlibf
   $F90 $OPT_c -c $rndf
   $F90 $OPT_c -c $statf
   $F90 $OPT_c -c $spfuncf
   $F90 $OPT_c -c $fndrootf
   $F90 $OPT_c -c $fftf
   $F90 $OPT_c -c $sortf
fi


echo "compile subprograms again? (n)"
read cyn
if [ $cyn = "y" ]; then
   $F90 $OPT_c -c $acfdir/$acf
   $F90 $OPT_c -c ./$head 2>/dev/null || $F90 $OPT_c -c ./$head
   $F90 $OPT_c -c ./$sub
fi


$F90 $OPT_c -c $main


libs="$nrlib $sort $fft $rnd $stat $spfunc $fndroot $acf $head $sub "
libs=`echo $libs | sed -e "s/.f90 /.o /g" -e "s/.f /.o /g"` 
main=`echo $main | sed -e "s/.f90 /.o /"`

$F90 $OPT_o -o $exe $main $libs
