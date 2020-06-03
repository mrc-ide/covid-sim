foreach ($R in "2.0","2.2","2.4","2.6")
{
$rs=$R/2
job submit ###flags### .\CovidSim.exe /NR:10 /c:8 /PP:preGB_R0=2.0.txt /P:p_NoInt.txt /CLP1:100000 /CLP2:0 /O:meanT8_NR10\NoInt_R0=${R} /D:..\population\GB_pop2018.bin /L:..\population\NetworkGB_8T.bin /R:${rs} 98798150 729101 17389101 4797132

foreach ($x in 60,100,200,300,400)
{
foreach ($z in 0.75,0.5,0.25)
{
$q=1000
$y=$x * $z
echo $y
foreach ($i in "CI_HQ_SD","PC_CI_SD","PC_CI_HQ_SD")
{
job submit ###flags### .\CovidSim.exe /NR:10 /c:8 /PP:preGB_R0=2.0.txt /P:p_${i}.txt /CLP1:${x} /CLP2:${q} /CLP3:${q} /CLP4:${q} /CLP5:${y} /O:meanT8_NR10\${i}_${x}_${y}_R0=${R} /D:..\population\GB_pop2018.bin /L:..\population\NetworkGB_8T.bin /R:${rs} 98798150 729101 17389101 4797132
}
}
}
}

