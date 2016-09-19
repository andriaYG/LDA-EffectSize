####Kruskal.test
setwd("/Users/simona/SGDev/BGI/project/IBD.v2/fisher.lefse")
read.table("M100gene.MGS.123hostFree.relative.abun.pro")->Dat;
read.table("occ50.MGS100.list",header = F,colClasses = "character")->occ50;
Dat[match(occ50[,1],rownames(Dat)),]->Dat
read.table("Sample.group.txt",header = F,sep = "\t")->label;
Result<-c();
Res_temp<-list(mgsID = character(),statistic = numeric(),p.value = numeric(),A=numeric(),B=numeric(),C=numeric(),type = character(),A.occ=numeric(),B.occ=numeric(),C.occ=numeric())
t(Dat)->Data;
label<-as.factor(label[,2])
for(i in 1:dim(Data)[2]){
  colnames(Data)[i]->Res_temp$mgsID;
  kruskal.test(Data[1:dim(Data)[1],i], label)->Res;
  Res_temp$A<-mean(Data[which(label==2),i]);
  Res_temp$B<-mean(Data[which(label==1),i]);
  Res_temp$C<-mean(Data[which(label==3),i]);
  Res_temp$type<-"C";
  if(Res_temp$A==max(Res_temp$A,Res_temp$B,Res_temp$C)){Res_temp$type<-"A";}
  if(Res_temp$B==max(Res_temp$A,Res_temp$B,Res_temp$C)){Res_temp$type<-"B";}
  Res_temp$A.occ<-length(which(Data[which(label==2),i]!=0))/length(Data[which(label==2),i]);
  Res_temp$B.occ<-length(which(Data[which(label==1),i]!=0))/length(Data[which(label==1),i]);
  Res_temp$C.occ<-length(which(Data[which(label==3),i]!=0))/length(Data[which(label==3),i]);
  Res$statistic->Res_temp$statistic;
  Res$p.value->Res_temp$p.value;
  Result<-rbind(Result,Res_temp);
}
#unique(which(Result[,3]<0.01))->p01.index
unique(which(Result[,3]<0.05))->p05.index 
write.table(Result[p05.index,],"M100.Kruskal.p05.mgs.list",col.names = T,row.names = F,sep ="\t",quote = F)

##since our data cannot fit the wilconxon test,we skip this step
#An LDA model is finally built with the class as dependent variable and the remaining feature values, subclass, and subject values as independent variables. This model is used to estimate their effect sizes, which are obtained by averaging the differences between class means (using unmodified feature values) with the differences between class means along the first linear discriminant axis, which equally weights features' variability and discriminatory power. The LDA score for each biomarker is obtained computing the logarithm (base 10) of this value after being scaled in the [1,106] interval and, regardless of the absolute values of the LDA score, it induces the ranking of biomarker relevance. For robustness, LDA is additionally supported by bootstrapping (default 30-fold) and subsequent averaging.
##LDA score estimating the effect side
lda.data <- data.frame( Data[,match(Result[p01.index,1],colnames(Data))],
                        group = label)   
#colnames(lda.data)<-unlist(lda.data[1,])
#lda.data<-data.frame(lda.data[-1,])
#write.table(lda.data,"test",col.names = T,row.names = F,sep="\t",quote = F);
#read.table("test",header = T,sep = "\t")->lda.data
##train <- sample(1:150, 75)
##table(lda.data)
library("MASS")
#test<-lda.data[,-c(7,23,58,60,72,86,89,90,91,99,103,105,112)]
#test<-lda.data[,-c(3,9,12,32,37,53,64,69,70,72,93,94,97,104,107,109,114,124,127,133,146,149,152,153,156,172,173,174,175,182,188,192,193,195,196,198,199,202,203,204,223,224,232,234,237,242,245,258,273,300,324)]
#test<-lda.data[,-c(6,91,93,102,113,116,117,118,128,130,134)]
#test<-lda.data[,-c(46,48,57,70,73,74,75,87,88,90,96,194,202)]
z <- lda(group ~ ., lda.data, prior = c(1,1,1)/3)
##abs(z$scaling[,"LD1"])/sum(abs(z$scaling[,"LD1"]))->sca.LD1
##=============================================
#lda.re<-list(mgsID = character(),lda.score = numeric())
#test<-c();
#for(i in 1:dim(z$means)[2]){
#	dimnames(lda.data)[[2]][i]->lda.re$mgsID;
#	Dis1<-z$means[1,i]-z$means[2,i]
#	Dis2<-z$means[2,i]-z$means[3,i]
#	Dis3<-z$means[1,i]-z$means[3,i]
#	#log(sqrt((Dis1^2+Dis2^2+Dis3^2)*abs(z$scaling[i,1])/3),10)->lda.re$lda.score
#	sqrt((Dis1^2+Dis2^2+Dis3^2)*abs(z$scaling[i,1])/3)->lda.re$lda.score
#	test<-rbind(test,lda.re)
#}


lda.re<-list(mgsID = character(),lda.score = numeric())
lda.result<-c();
length(which(lda.data[,"group"]=="1"))->n1;
length(which(lda.data[,"group"]=="2"))->n2;
length(which(lda.data[,"group"]=="3"))->n3;
for(i in 1:dim(z$means)[2]){
  dimnames(lda.data)[[2]][i]->lda.re$mgsID;
  Dis1<-abs(z$means[1,i]-z$means[2,i]);
  Dis2<-abs(z$means[2,i]-z$means[3,i]);
  Dis3<-abs(z$means[1,i]-z$means[3,i]);
  sum((lda.data[which(lda.data[,"group"]=="1"),i]-z$means[1,i])^2)/(n1-1)->s1;
  if(s1==0){s1<-(1e-10)}
  sum((lda.data[which(lda.data[,"group"]=="2"),i]-z$means[2,i])^2)/(n2-1)->s2;
  if(s2==0){s2<-(1e-10)}
  sum((lda.data[which(lda.data[,"group"]=="3"),i]-z$means[3,i])^2)/(n3-1)->s3;
  if(s3==0){s3<-(1e-10)}
  s12<-sqrt(((n1-1)*s1+(n2-1)*s2)/(n1+n2-2))
  s23<-sqrt(((n2-1)*s2+(n3-1)*s3)/(n3+n2-2))
  s13<-sqrt(((n1-1)*s1+(n3-1)*s3)/(n1+n3-2))
  #log(sqrt((Dis1^2+Dis2^2+Dis3^2)*abs(z$scaling[i,1])/3),10)->lda.re$lda.score
  log(abs(z$scaling[i,1])*(Dis1/s12+Dis2/s23+Dis3/s13),10)->lda.re$lda.score
  lda.result<-rbind(lda.result,lda.re)
}


lda.result[which(lda.result[,2]>3),]->temp
write.table(temp,"temp>3",col.names = T,row.names = F,sep="\t",quote = F);
#read.table("temp1.5",header = T,sep = "\t")->temp
#temp[order(temp[,2],decreasing =T),] -> heatmap.data

