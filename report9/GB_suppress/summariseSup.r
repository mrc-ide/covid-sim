cur_path="S:/GB_suppress/mean8_10/"
f=TRUE
stats.names=c("R0","Trig","Dur","Scenario","variable","Total incidence","Maximum incidence","Mean time","Max time","Median time")
for(j in c("PropSocDist","Mild","ILI","SARI","Critical","incMild","incILI","incSARI","incCritical","incDeath"))
{

for (R in c("2.0","2.2","2.4","2.6"))
{

	for(tr in c(60,100,200,300,400))
	{
	cc=1
	for( v in c(0.75,0.5,0.25))
		{ 
		  y=tr*v
		  i="NoInt"
		  data1 <- read.delim(paste0(cur_path,i,"_R0=",R,".avNE.severity.xls"), header = TRUE)
		  data1=data1[data1$t<=730,]
		  if(cc==1) outd=data1[,c("t")]
		  cnm=paste0(i,":",v)
		  cc=cc+1
		  tmp=as.numeric(data1[[j]])
		  if(j=="PropSocDist")
			{
			outd=cbind(outd,tmp)
			colnames(outd)[cc]=cnm
			tmax=max(tmp)+1e-10
			tmp=tmp/tmax # unpopulated mcells never social distance
			tmed=mean(tmp)
			tmp=ifelse(tmp==0,0,1)
			tmp2=c(0,tmp[2:length(tmp)] - tmp[1:(length(tmp)-1)])
			mi=sum(ifelse(tmp2==1,1,0))
			ton=data1$t[which(tmp2==1)]
			toff=data1$t[which(tmp2==-1)]
			mt=toff[1]-ton[1]
			tm=ton[2]-toff[1]
			si=tmed*730/(730-80) # remove initial period before control starts
			}
		  else
			  {
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
			  }
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
		   for (i in c("CI_HQ_SD","PC_CI_SD","PC_CI_HQ_SD"))
				{
				  data1 <- read.delim(paste0(cur_path,i,"_",tr,"_",y,"_R0=",R,".avNE.severity.xls"), header = TRUE)
				  data1=data1[data1$t<=730,]
				  cnm=paste0(i,":",v)
				  cc=cc+1
				  tmp=as.numeric(data1[[j]])
				  if(j=="PropSocDist")
					{
					outd=cbind(outd,tmp)
					colnames(outd)[cc]=cnm
					tmax=max(tmp)+1e-10
					tmp=tmp/tmax # unpopulated mcells never social distance
					tmed=mean(tmp)
					tmp=ifelse(tmp==0,0,1)
					tmp2=c(0,tmp[2:length(tmp)] - tmp[1:(length(tmp)-1)])
					mi=sum(ifelse(tmp2==1,1,0))
					ton=data1$t[which(tmp2==1)]
					toff=data1$t[which(tmp2==-1)]
					mt=toff[1]-ton[1]
					tm=ton[2]-toff[1]
					si=tmed*730/(730-80) # remove initial period before controls start
					}
				  else
					{
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
					}
				stats2=data.frame(R,tr,v,i,j,si,mi,mt,tm,tmed)
				names(stats2)=stats.names
				stats=rbind(stats,stats2)

				}
		}	  
    write.csv(outd,paste0(cur_path,"R0=",R,"_",j,"_",tr,".csv"))
}
}
}
write.csv(stats,paste0(cur_path,"stats_contain.csv"))
