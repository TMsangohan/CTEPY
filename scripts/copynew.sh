#!/bin/bash

for file in *int_1p71.py
do 
  cp "$file" "${file/1p71/1p9}"
done
