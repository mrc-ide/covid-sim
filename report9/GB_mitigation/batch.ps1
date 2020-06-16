foreach ($R in "2.2","2.4")
{
$rs=$R/2
job submit ####flags#### .\CovidSim.exe /NR:10 /c:8 /PP:preGB_R0=2.0.txt /P:p_NoInt.txt /CLP1:100000 /CLP2:0 /O:MeanT8_NR10\NoInt_R0=${R} /D:..\population\GB_pop2018.bin /L:..\population\NetworkGB_8T.bin /R:${rs} 98798150 729101 17389101 4797132

foreach ($x in "100","300","1000","3000")
{
foreach ($q in 91,121,152,182)
{
$qo=$q+30
echo $qo
foreach ($i in "CI","CI_HQ","CI_HQ_SD","CI_SD","CI_HQ_SDOL70","MG","PC","PC_CI_HQ_SDOL70","PC_CI_SD","PC_CI_HQ_SD","MG_PC_CI_HQ_SDOL70")
{
job submit ####flags#### .\CovidSim.exe /NR:10 /c:8 /PP:preGB_R0=2.0.txt /P:p_${i}.txt /CLP1:${x} /CLP2:${q} /CLP3:${qo} /CLP4:${qo} /O:MeanT8_NR10\${i}_${x}_${q}_R0=${R} /D:..\population\GB_pop2018.bin /L:..\population\NetworkGB_8T.bin /R:${rs} 98798150 729101 17389101 4797132
}
}
}
}

