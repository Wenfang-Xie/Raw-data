library(survival)
library(survminer)
inputFile="km.txt"        
outFile="GZMA2.pdf"      
var="GZMA"                  

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=rt[,c("futime","fustat",var)]
colnames(rt)=c("futime","fustat","var")

group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 5,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf(file=outFile,onefile = FALSE,width = 6,height =5)
print(surPlot)
dev.off()


