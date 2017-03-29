
read.table("M700gene.MGS.117hostFree.relative.abun.pro")->Dat;
read.table("Sample.group.txt",header = F,sep = "\t")->label;
Dat[,match(label[,1],colnames(Dat))]->Dat
Result<-c();
Res_temp<-list(mgsID = character(),statistic = numeric(),p.value = numeric(),A=numeric(),B=numeric(),C=numeric(),
               #type = character(),
               A.occ=numeric(),B.occ=numeric(),C.occ=numeric())
t(Dat)->Data;
label<-as.factor(label[,2])
for(i in 1:dim(Data)[2]){
  colnames(Data)[i]->Res_temp$mgsID;
  kruskal.test(Data[1:dim(Data)[1],i], label)->Res;
  Res_temp$A<-mean(Data[which(label=="A"),i]);
  Res_temp$B<-mean(Data[which(label=="B"),i]);
  Res_temp$C<-mean(Data[which(label=="C"),i]);
  Res_temp$type<-"C";
  if(Res_temp$A==max(Res_temp$A,Res_temp$B,Res_temp$C)){Res_temp$type<-"A";}
  if(Res_temp$B==max(Res_temp$A,Res_temp$B,Res_temp$C)){Res_temp$type<-"B";}
  Res_temp$A.occ<-length(which(Data[which(label=="A"),i]!=0))/length(Data[which(label=="A"),i]);
  Res_temp$B.occ<-length(which(Data[which(label=="B"),i]!=0))/length(Data[which(label=="B"),i]);
  Res_temp$C.occ<-length(which(Data[which(label=="C"),i]!=0))/length(Data[which(label=="C"),i]);
  Res$statistic->Res_temp$statistic;
  Res$p.value->Res_temp$p.value;
  if(Res_temp$A.occ > 0.05 && (Res_temp$B.occ > 0.05 || Res_temp$C.occ> 0.05)){
    Result<-rbind(Result,Res_temp);
    }
}
unique(which(Result[,3]<0.05))->p05.index 
write.table(Result[p05.index,],"M700.Kruskal.p05.mgs.list.v2.txt",col.names = T,row.names = F,sep ="\t",quote = F)

Data<-Data[,match(as.character(Result[p05.index,1]),colnames(Data))]
read.table("Sample.group.EEN.txt",header = F,sep = "\t")->E.label
read.table("Sample.group.CDvsCT.txt",header = F,sep = "\t")->C.label
C.Data<-Data[match(as.character(C.label[,1]),rownames(Data)),]
C.Data<-cbind(as.character(C.label[match(rownames(C.Data),as.character(C.label[,1])),2]),as.character(C.label[match(rownames(C.Data),as.character(C.label[,1])),3]),C.Data)
E.Data<-Data[match(as.character(E.label[,1]),rownames(Data)),]
E.Data<-cbind(as.character(E.label[match(rownames(E.Data),as.character(E.label[,1])),2]),as.character(E.label[match(rownames(E.Data),as.character(E.label[,1])),3]),E.Data)
colnames(C.Data)[1:2]<-c("group","disease")
colnames(E.Data)[1:2]<-c("group","disease")

calculator<-function(Dat,index,g1,g2){
  res<-list(clinical=character(),pvalue=numeric(),enrichment=character(),mean1=numeric(),mean2=numeric(),occ1=numeric(),occ2=numeric())
  Result<-c()
  for(i in 3:ncol(Dat)){
    #if(sum(as.numeric(Dat[which(Dat[,index]==g1),i]))!=0 && sum(as.numeric(Dat[which(Dat[,index]==g2),i]))!=0){
    res$clinical<-colnames(Dat)[i]
    res$mean1<-mean(as.numeric(Dat[which(Dat[,index]==g1),i]))
    res$mean2<-mean(as.numeric(Dat[which(Dat[,index]==g2),i]))
    if(res$mean1>res$mean2){res$enrichment<-g1}
    if(res$mean1<res$mean2){res$enrichment<-g2}
    if(res$mean1==res$mean2){res$enrichment<-"NA"}
    res$occ1<-length(which(as.numeric(Dat[which(Dat[,index]==g1),i])!=0))/length(Dat[which(Dat[,index]==g1),i])
    res$occ2<-length(which(as.numeric(Dat[which(Dat[,index]==g2),i])!=0))/length(Dat[which(Dat[,index]==g2),i])
    res$pvalue<-wilcox.test(as.numeric(Dat[which(Dat[,index]==g1),i]),as.numeric(Dat[which(Dat[,index]==g2),i]))$p.value
    Result<-rbind(Result,res)
    #}
  }
  q.value<-p.adjust(Result[,2])
  Result<-cbind(Result[,c(1,2)],q.value,Result[,3:7])
  return(Result)
}
C.Data[which(C.Data[,"disease"]=="CT"),]->CT.Dat;
calculator(CT.Dat,1,"A","B")->CT.result;
C.Data[union(intersect(which(C.Data[,"group"]=="A"),which(C.Data[,"disease"]=="CT")),intersect(which(C.Data[,"group"]=="B"),which(C.Data[,"disease"]=="CD"))),]->A.CTvsB.CD.Dat;
calculator(A.CTvsB.CD.Dat,1,"A","B")->A.CTvsB.CD.result;
C.Data[union(intersect(which(C.Data[,"group"]=="A"),which(C.Data[,"disease"]=="CT")),intersect(which(C.Data[,"group"]=="C"),which(C.Data[,"disease"]=="CD"))),]->A.CTvsC.CD.Dat;
calculator(A.CTvsC.CD.Dat,1,"A","C")->A.CTvsC.CD.result;
C.Data[which(C.Data[,"group"]=="B"),]->B.Dat;
calculator(B.Dat,2,"CT","CD")->B.result;
C.Data[union(intersect(which(C.Data[,"group"]=="B"),which(C.Data[,"disease"]=="CT")),intersect(which(C.Data[,"group"]=="C"),which(C.Data[,"disease"]=="CD"))),]->B.CTvsC.CD.Dat;
calculator(B.CTvsC.CD.Dat,1,"B","C")->B.CTvsC.CD.result;
C.Data[which(C.Data[,"disease"]=="CD"),]->B.CDvsC.CD.Dat;
B.CDvsC.CD.Dat[-which(B.CDvsC.CD.Dat[,"group"]=="A"),]->B.CDvsC.CD.Dat
calculator(B.CDvsC.CD.Dat,1,"B","C")->B.CDvsC.CD.result;
calculator(E.Data,2,"1pre","2post")->PreEENvsPostEEN.result;

#read.table("/Users/simona/SGDev/BGI/project/IBD.v5/12.bacteria/abundance.list.txt",header = F,sep = "\t")->ref
#A.CTvsC.CD.result[which(A.CTvsC.CD.result[,3]<0.2),]->test
#cbind(as.character(ref[,2]),test[match(as.character(ref[,1]),as.character(test[,1])),1:3])

MGS<-as.character(unique(c(CT.result[which(CT.result[,2]<0.05),1],
                           A.CTvsB.CD.result[which(A.CTvsB.CD.result[,2]<0.05),1],
                           A.CTvsC.CD.result[which(A.CTvsC.CD.result[,2]<0.05),1],
                           B.result[which(B.result[,2]<0.05),1],
                           B.CTvsC.CD.result[which(B.CTvsC.CD.result[,2]<0.05),1],
                           B.CDvsC.CD.result[which(B.CDvsC.CD.result[,2]<0.05),1],
                           PreEENvsPostEEN.result[which(PreEENvsPostEEN.result[,2]<0.05),1])))
lda.data <- data.frame(Data[,match(MGS,colnames(Data))],
                        group = label) 
library("MASS")

means<-c()
for(i in 1:30){
  sample(1:nrow(lda.data),117,replace = T)->sl
  temp<-lda.data[sl,]
  z <- lda(group ~ ., temp, prior = c(1,1,1)/3)
  z$scaling[,1]->w
  w.unit<-w/sqrt(sum(w^2))
  ss<-temp[,-match("group",colnames(temp))]
  xy.matrix<-as.matrix(ss)
  LD<-xy.matrix%*%w.unit
  effect.size <- mean(abs(mean(LD[temp[,"group"]=="A"]) - mean(LD[temp[,"group"]=="B"])),
                      abs(mean(LD[temp[,"group"]=="A"]) - mean(LD[temp[,"group"]=="C"])),
                      abs(mean(LD[temp[,"group"]=="B"]) - mean(LD[temp[,"group"]=="C"])))
  wfinal <-  effect.size * w.unit
  mm <- (abs(z$means["A",]-z$means["B",])+abs(z$means["B",]-z$means["C",])+abs(z$means["A",]-z$means["C",]))/3
  coeff<-abs(wfinal)
  value<-log((((coeff+mm)*0.5)*(1e+6)+1),10)
  means<-cbind(means,value[match(colnames(test),names(value))])
}
means[which(rowMeans(means)>2),]->re
write.table(rowMeans(re),"lefse.list.p05.lda2.v2.txt",col.names = F,row.names = T,sep = "\t",quote = F)

#======table
read.table("lefse.list.p05.lda2.v2.txt",header = F,sep = "\t")->MGS.list
read.table("M700.Kruskal.p05.mgs.list.v2.txt",header = T,sep = "\t")->MGS.res
colnames(MGS.list)<-c("mgsID","lda.score")
merge(MGS.list,MGS.res,by="mgsID")->MGS.table
write.table(MGS.table,"85MGS.table.test.txt",col.names = T,row.names = F,sep = "\t",quote = F)

