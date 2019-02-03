for i in IP2_500 IP125_500_500_500 IP125_1000_500_1000
do
  for ex in 1p50 1p65 
   do
     for int in 1p6  1p85 2p1
       do
         cp "CTEPYv5.py" "CTEPY_${i}_${ex}_${int}.py" 
       done 
   done
done
