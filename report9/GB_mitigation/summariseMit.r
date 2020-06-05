cur_path="S:/GB_mitigation/meanT8_NR10/"
f=TRUE
stats.names=c("R0","Trig","Dur","Scenario","variable","Total incidence","Maximum incidence","Mean time","Max time","Median time")
for(j in c("Mild","ILI","SARI","Critical","incMild","incILI","incSARI","incCritical","incDeath"))
{

for (R in c("2.4","2.2"))
{

	for(tr in c("100","300","1000","3000"))
	{
	cc=1
	for( v in c("91","121","152","182"))
		{ 
		  i="NoInt"
		  data1 <- read.delim(paste0(cur_path,i,"_R0=",R,".avNE.severity.xls"), header = TRUE)
		  if(cc==1) outd=data1[,c("t")]
		  cnm=paste0(i,":",v)
		  cc=cc+1
		  tmp=as.numeric(data1[[j]])
		  outd=cbind(outd,tmp)
		  colnames(outd)[cc]=cnm
		  si=sum(tmp)
		  mi=max(tmp)		  
		  t2=tmp*data1$t
		  ct=cumsum(tmp)
		  medi=si/2
		  tmed=which(ct>medi)[1]-1
		  mt=sum(t2)/si
		  tm=which(tmp==mi)[1]-1
		  if(f)
			{
		    stats=data.frame(R,tr,v,i,j,si,mi,mt,tm,tmed)
			names(stats)=stats.names
			f=FALSE
			}
		  else
			{
		    stats2=data.frame(R,tr,v,i,j,si,mi,mt,tm,tmed)
			names(stats2)=stats.names
			stats=rbind(stats,stats2)
			}
		   for (i in c("CI","CI_HQ","CI_HQ_SD","CI_SD","CI_HQ_SDOL70","MG","PC","PC_CI_HQ_SDOL70","PC_CI_SD","PC_CI_HQ_SD","MG_PC_CI_HQ_SDOL70"))
				{
				  data1 <- read.delim(paste0(cur_path,i,"_",tr,"_",v,"_R0=",R,".avNE.severity.xls"), header = TRUE)
				  cnm=paste0(i,":",v)
				  cc=cc+1
				  tmp=as.numeric(data1[[j]])
				  outd=cbind(outd,tmp)
				  colnames(outd)[cc]=cnm
				  si=sum(tmp)
				  mi=max(tmp)		  
				  t2=tmp*data1$t
				  ct=cumsum(tmp)
				  medi=si/2
				  tmed=which(ct>medi)[1]-1
				  mt=sum(t2)/si
				  tm=which(tmp==mi)[1]-1
				  stats2=data.frame(R,tr,v,i,j,si,mi,mt,tm,tmed)
				  names(stats2)=stats.names
				  stats=rbind(stats,stats2)
				}
		}	  
    write.csv(outd,paste0(cur_path,"R0=",R,"_",j,"_",tr,".csv"))
}
}
}
write.csv(stats,paste0(cur_path,"stats_mitigation.csv"))
