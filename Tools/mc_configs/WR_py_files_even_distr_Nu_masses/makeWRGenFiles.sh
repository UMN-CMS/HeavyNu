#!/bin/bash

NuMasses=(100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900)

for i in {0..18}
do
	eval "sed 's/333/${NuMasses[$i]}/g' WR_2100_Nu_333_GEN.py > WR_2100_Nu_${NuMasses[$i]}_GEN.py"
	cmsRun WR_2100_Nu_${NuMasses[$i]}_GEN.py

done


