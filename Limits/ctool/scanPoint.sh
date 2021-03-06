#!/bin/sh

cmswd=$1
jobwd=$2
inp=$3
mass=$4
outp=$5

method=$6
toys=$7
comment=$8
special=$9
#seed=1789
seed=9130

source /local/cms/sw/cmsset_CMSSW6X.sh

cd ${cmswd}

cmsenv

cd ${jobwd}

rrange="--rMin 0.01 --rMax 10.0 --rAbsAcc 0.01"

if [ "${method}" == "MarkovChainMC" ]; then
   extra="-i 30000 --tries 50"
   extraObs="-i 30000 --tries 100"
fi

if [ "${method}" == "HybridNew" ]; then
   extra="--frequentist --testStat LHC --fork 8 "
   extraObs="--frequentist --testStat LHC --fork 4 "
fi


localOF=/tmp/sp_$$

rm -f ${outp}

echo "Limits for ${mass}" > ${outp}
echo "Method=${method}"  >> ${outp}
echo "Toys=${toys}"  >> ${outp}
echo "Special=${special}"  >> ${outp}
echo ${comment} >> ${outp}
echo " " >> ${outp}

echo "Observed Limit Pass" >> ${outp}
echo "==============================================" >> ${outp}

combine -v0  -n WRmuObs -m ${mass} -M ${method} -H ProfileLikelihood ${rrange} -s ${seed} ${extraObs} ${special} ${inp} > ${localOF} 2>&1 
 
cat ${localOF} >> ${outp}

echo " " >> ${outp}

mv ${localOF} check.txt
rm -f ${localOF}

echo "Expected Limit Pass" >> ${outp}
echo "==============================================" >> ${outp}

combine -v0 -t${toys} -n WRmu -m ${mass} -M ${method} -H ProfileLikelihood ${rrange} -s ${seed} ${extra} ${special} ${inp} > ${localOF} 2>&1 

cat ${localOF} >> ${outp}

rm ${localOF}




#echo combine -v0 -t${toys} -n WRmu -m ${mass} -M ${method} -H ProfileLikelihood -s ${seed} ${extra} ${inp} 

#-m 16001000 -s1002 -n WRmu -t100 -v0   
