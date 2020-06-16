#!/bin/bash
allR="2.2 2.4"
allx="100 300 1000 3000"
allq="91 121 152 182"
alli="CI CI_HQ CI_HQ_SD CI_SD CI_HQ_SDOL70 MG PC PC_CI_HQ_SDOL70 PC_CI_SD PC_CI_HQ_SD MG_PC_CI_HQ_SDOL70"

for R in $allR
do
	rs=$(echo $R | awk '{print $1/2}')
	CovidSim /NR:10 /c:8 /PP:preGB_R0=2.0.txt /P:p_NoInt.txt /CLP1:100000 /CLP2:0 /O:MeanT8_NR10/NoInt_R0=${R} /D:../population/GB_pop2018.bin /L:../population/NetworkGB_8T.bin /R:${rs} 98798150 729101 17389101 4797132
	for x in $allx
	do
	   for q in $allq
	   do
		   qo=$(echo $q | awk '{print $1+30}')
		   #echo $qo
		   for i in $alli
		   do
			  CovidSim /NR:10 /c:8 /PP:preGB_R0=2.0.txt /P:p_${i}.txt /CLP1:${x} /CLP2:${q} /CLP3:${qo} /CLP4:${qo} /O:MeanT8_NR10/${i}_${x}_${q}_R0=${R} /D:../population/GB_pop2018.bin /L:../population/NetworkGB_8T.bin /R:${rs} 98798150 729101 17389101 4797132
		   done
	   done
	done
done

