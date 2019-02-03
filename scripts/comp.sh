for i in *.f95
do
  bash pythonfortrancompile.sh ${i%%.*} $i
done
