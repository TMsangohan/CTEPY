for i in *int_1p9.py
do
  bsub -q 1nd -J test  -o /dev/null -e /dev/null -G  u_SLAP CTEPYlxpluschamonix.sh $i
done
