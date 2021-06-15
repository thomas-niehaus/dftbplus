#!/bin/sh

# Starting with a geometry file Monomer.geo.gen, creates a copy that is displaced and rotated
# with respect to the first monomer. Creates MO data files for both monomers and the 
# command file for the model Hamiltonian. 

# To be replaced by your executable (make sure the path in dftb_in.hsd
# is changed as well) 
exe=/data/niehaus/git_wrk/dftb+/master/_build/prog/dftb+/dftb+

# Hubbards for the LC-DFTB ob2 parameters (order the same as in gen file)
HUBBARDS='0.3494 0.3929 0.4223 0.4912'

# As defined by Liu et al. Used here to match gap of reference calculation
e0='0.0624'

# Displacement vector of second fragment with respect to the first
x1=3.4
x2=0.0 
x3=0.0

# Some care is needed to remove files not required for the task in question
rm -rf HamPar.dat 2by2.dat adia.p.dat adia.m.dat jointGeo.xyz
rm -rf *.pre.eigenvec.out sign.dat modelCT.cmd

# Monomer.geo.gen contains the monomer geometry properly oriented
cp Monomer.geo.gen geo.gen 
cp geo.gen A.geo.gen
nAtomFragA=`awk 'NR==1 {print $1}' geo.gen`
# Generate MO coefficient of fragment A
$exe > A.dftb.out

# This makes sure that phases are consistent 
cp eigenvec.out A.cur.eigenvec.out
cp eigenvec.out A.pre.eigenvec.out
cp eigenvec.out B.pre.eigenvec.out
echo "1 1 1 1" > sign.dat

# Determine fragment frontier orbitals and dimension
awk '$3 > 1.0 {print $0}' band.out | tail -1 > band.tmp
iHOMO_A=`awk '{print $1}' band.tmp`
awk '$3 < 1.0 {print $0}' band.out | head -1 > band.tmp
iLUMO_A=`awk '{print $1}' band.tmp`
nDimFragA=`tail -2 band.out | awk '{print $1}'`

# Loop over angle phi around x-axis 
for phi in `awk 'BEGIN { for( i=0; i<=90; i+=3 ) print i }'`; do
 echo "-> phi = " $phi

 # Important
 rm -rf modelCT.cmd

 # Create rotated and displaced fragment B, as well as geometry of complex
 cp Monomer.geo.gen geo.gen 
 awk 'NR<3 {print $0}' geo.gen > headerMonomer.gen
 awk 'BEGIN {pi = atan2(0,-1)} NF==5 {printf("%5d %5d %12.8f %12.8f %12.8f\n", $1, $2, $3 + x1, $4*cos(pi*phi/180.0) - $5*sin(pi*phi/180.0) + x2, $4*sin(pi*phi/180.0) + $5*cos(pi*phi/180.0) + x3)}' x1=$x1 x2=$x2 x3=$x3 phi=$phi geo.gen > tail.gen
 cat headerMonomer.gen tail.gen > B.geo.gen 
 nAtomDimer=`echo $nAtomFragA | awk '{print $1 * 2}'`
 cat A.geo.gen tail.gen > tmp.geo
 awk '{if(NR==1){print x, "C"} else {print}}' x=$nAtomDimer tmp.geo > Dimer.geo.gen
 rm -rf  headerMonomer.gen tail.gen tmp.geo
 gen2xyz Dimer.geo.gen 
 # Contencated dimer geometries can be viewed by jmol 
 cat Dimer.geo.xyz >> jointGeo.xyz

 # Create MO of fragment B. Strictly speaking this is not required in this example 
 # since the MOs do not change (they are invariant with respect to rigid rotation)
 cp B.geo.gen geo.gen
 $exe > B.dftb.out
 cp eigenvec.out B.cur.eigenvec.out

 # Setup of command file (in this example A==B)
 echo $nAtomFragA $nDimFragA $nAtomFragA $nDimFragA > modelCT.cmd
 echo $iHOMO_A $iLUMO_A $iHOMO_A $iLUMO_A >> modelCT.cmd
 echo $e0 >> modelCT.cmd
 echo $HUBBARDS >> modelCT.cmd
 cp Dimer.geo.gen geo.gen
 
 # The actual modelCT computations
 $exe > Dimer.dftb.out
 awk 'NR==2 {printf("%2d",phi);print $0}' phi=$phi modelCT.out >> HamPar.dat
 awk 'NR==4 {printf("%2d",phi);print $0}' phi=$phi modelCT.out >> 2by2.dat
 awk 'NR==6 {printf("%2d",phi);print $0}' phi=$phi modelCT.out >> adia.p.dat
 awk 'NR==8 {printf("%2d",phi);print $0}' phi=$phi modelCT.out >> adia.m.dat

 # Important for consistent MO signs (not relevant in this simple example)
 cp B.cur.eigenvec.out B.pre.eigenvec.out

done

# Clean up  
rm -rf band.out band.tmp detailed.out geo.gen *dftb.out eigenvec.bin charges.bin
rm -rf *eigenvec.out Dimer.geo.* ?.geo.gen sign.dat dftb_pin.hsd 
