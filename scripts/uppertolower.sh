for i in *.tfs
do 
  mv $i `echo $i | tr [:upper:] [:lower:]`
done
