#!/bin/bash

outfn=$1
infn=$2

f2py -c --fcompiler=gnu95 -m $outfn $infn


