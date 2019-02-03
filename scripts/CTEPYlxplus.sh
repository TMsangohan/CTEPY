#!/bin/bash

#setting base directory
base="/afs/cern.ch/work/t/tomerten/CTEPYlxplus"

# copying fortran files
cp $base/*.so .

# copying python files
cp $base/CTEsetparameters.py .
cp $base/CTEwrite.py .
cp $base/CTEraddamping.py .

# copying tfs fils
cp $base/*.tfs .

# copying main program
cp $base/CTEPYv5.py .

# run the program
/afs/cern.ch/user/t/tomerten/anaconda2/bin/python2.7 CTEPYv5.py >> test.out

# copy back results
mv *.csv $base/results
mv *.out $base/results



