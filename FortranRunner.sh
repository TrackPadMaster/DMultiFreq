#!/bin/bash
 

PHI="    0.0000000000000     0.3141592653590     0.6283185307180     0.9424777960769     1.2566370614359     1.5707963267949     1.8849555921539     2.1991148575129     2.5132741228718     2.8274333882308     3.1415926535898     3.4557519189488     3.7699111843078     4.0840704496667     4.3982297150257     4.7123889803847     5.0265482457437     5.3407075111026     5.6548667764616     5.9690260418206     6.2831853071795"

for phi in $PHI; do
    mydir='PHI_'$phi
    mkdir $mydir
    awk -v phi=$phi '{if($3 == "phi:") {print phi" ! phi: phase difference"} else {print $0}}' MultiIn0 > $mydir/input.dat 
    cd $mydir
    mpirun -n 8 ../a.out
    cd ..
done

