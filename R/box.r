library('ggplot2')
library('RColorBrewer')
library(ggpubr)
data<-read.csv("box.csv",header=TRUE)
data$group<-factor(data$group,levels = c("num(T)=46","num(N)=16"))
a<-ggplot(data,aes(x=group,y=value,fill=group))+ #??fill=????????????ɫ
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ #?????Դ???????ͼû?к???ĩ??û?ж̺??ߣ?ʹ?????????ķ?ʽ????
  geom_boxplot(size=1.0,fill=c("#ffb6c1", "#87cefa"),outlier.fill="white",outlier.color="white")+ #size????????ͼ?ı߿??ߺͺ??????߿??ȣ?fill??????????ɫ??outlier.fill??outlier.color?????쳣????????
  geom_jitter(aes(fill=group),width =0.2,shape = 21,size=1.5)+ #????Ϊ??ˮƽ???򶶶???ɢ??ͼ??widthָ??????ˮƽ???򶶶??????ı???????ֵ
  scale_fill_manual(values = c("#ff4040", "#0000ff","#F0E442"))+  #????????????ɫ
  scale_color_manual(values=c("black","black","black"))+ #????ɢ??ͼ??ԲȦ????ɫΪ??ɫ
  ggtitle(" ")+ #?????ܵı???
  theme_bw()+ #??????Ϊ??ɫ
  theme(legend.position="none", #????Ҫͼ??
        axis.text.x=element_text(colour="black",family="Times",size=12), #????x???̶ȱ?ǩ??????????
        axis.text.y=element_text(family="Times",size=12,face="plain"), #????x???̶ȱ?ǩ??????????
        axis.title.y=element_text(family="Times",size = 12,face="plain"), #????y???ı?????????????
        axis.title.x=element_text(family="Times",size = 12,face="plain"), #????x???ı?????????????
        plot.title = element_text(family="Times",size=12,face="bold",hjust = 0.5), #?????ܱ?????????????
        panel.grid.major = element_blank(), #????ʾ??????
        panel.grid.minor = element_blank())+
  stat_compare_means()+
  ylab("GSDMB expression")+xlab("SKCM")
pdf('GSDMB_plot.pdf',width = 4,height = 8)
a
dev.off()
