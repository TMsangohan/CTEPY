#!/bin/bash

#setting base directory
base="/afs/cern.ch/work/t/tomerten/CTEPYlxplus"
mainprog=$1
# copying fortran files
cp $base/*.so .

# copying python files
cp $base/CTEsetparameters.py .
cp $base/CTEwrite.py .
cp $base/CTEraddamping.py .

# copying tfs fils
cp $base/*.tfs .

# copying main program
cp $base/$mainprog .

# run the program
/afs/cern.ch/user/t/tomerten/anaconda2/bin/python2.7 $mainprog >> "{$mainprog}test.out"

# copy back results
mv *.csv $base/results
mv *.out $base/results



