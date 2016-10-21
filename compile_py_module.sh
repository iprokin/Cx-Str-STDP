#!/bin/bash

# (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
# INRIA Rhone-Alpes
# STDP model : script to compile model with gfortran and f2py

LODEFOLDER='./odepack'
F2PY='python2 f2py2.py'
NAME="solve_py"

#compile odepack
cd $LODEFOLDER
for FILE in `ls -1 *.f`; do
    gfortran -c -Ofast -fPIC $FILE
    if [ $? -ne 0 ]; then
        echo "Errors compiling " $FILE
        exit
    fi
done
ar qc libodepack.a *.o
rm *.o
cd -

#compile python module linking it with lodepack
$F2PY -L$LODEFOLDER -lodepack -c --f90flags='-ffree-form -ffree-line-length-none' --fcompiler=gnu95 --opt='-Ofast' solve_py.pyf general_math.f95 statevars_mod.f95 pars_mod.f95 ghk_flux.f95 caL13.f95 TRPV1.f95 subcellular.f95 CaMKII_plast.f95 AMPA.f95 NMDA.f95 CB1R.f95 qsort_c_module.f95 stims.f95 comp_part.f95 solve_py.f95 > mp.txt 2>&1

#remove libodepack after it is linked
rm $LODEFOLDER/libodepack.a

echo "----------------------"
echo "Checking for errors..."
echo "----------------------"
ERRS=`grep -i "err" -C 1 mp.txt`
if [ -z "$ERRS" ]; then
    echo "no errors found"
else
    echo $ERRS
    echo '- - - - - - - - - - - - - - -'
    echo 'for details check file mp.txt'
fi
