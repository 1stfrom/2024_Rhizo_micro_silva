library(data.table)
library(lme4)
options(scipen = 10)
d=read.table("ASV_GWAS_re839.txt",header = T,sep="\t",as.is=T)[,-c(2:6)]
inf=fread("sample_info.txt", header=T,data.table=F)
colnames(inf)[1]="ID"
d=merge(d,inf,by="ID")
d=d[d[,837]=="HN" & d[,840]=="Inbred",]
#d=d[d[,4]=="LN",]
he=NULL
bl=as.data.frame(unique(d[,836]))
colnames(bl)="id"
for(i in 2:835)
{
  na=colnames(d)[i]
data=d[,c(1,836:840,i)]
data$genotype=factor(data$genotype)
data$Block=factor(data$Block)
data$sublock=factor(data$sublock)
data$Type=factor(data$Type)
colnames(data)[7]="y"
fit=lmer(y~genotype+(1|Block)+ (1|Block:sublock),data=data)
isSingular(fit, tol = 1e-4)
x=summary(fit)
BLUE <- as.data.frame(fixef(fit))
rownames(BLUE)=gsub("genotype","",rownames(BLUE))
a=bl[,1][!bl[,1]%in%rownames(BLUE)]
da=na.omit(data)
da1=da[da[,2]%in%a,]
rownames(BLUE)[1]=as.character(unique(da1[,2]))
b=data.frame(id=rownames(BLUE),Phe=c(BLUE[1,1],(BLUE[-1,1]+BLUE[1,1])))
colnames(b)[2]=na
bl=merge(bl,b,by="id",all=T)
colnames(bl)[ncol(bl)]=na
###calculate H2
fit=lmer(y~(1|genotype)+(1|Block)+ (1|Block:sublock),data=data)
x=summary(fit)
v=as.data.frame(x$varcor)
vg=v[1,4]
ve=v[4,4]##Residual
H=vg/(vg+ve/2)
##H=Vg/(Vg+Vge/y+Ve/yb)  #y:year  b: block Ve:Residual
##reference: Hallauer, A.R. and Miranda, J.B. (1988) Quantitative Genetics in Maize Breeding. Ames, IA: Iowa State University, 283.
hee=c(na,H)
he=rbind(he,hee)
}
colnames(he)=c("Trait","Heritability")
he=data.frame(he,pop="Inbred",N="HN")
bl=data.frame(bl,pop="Inbred",N="HN")
write.table(bl,file="data_blue_HN_inb_result.txt",row.names = F,quote = F,sep="\t")
write.table(he,file="Heritability_HN_inb_result.txt",row.names = F,quote = F,sep="\t")

##############################################################
d=read.table("ASVre_count_GWAS.txt",header = T,sep="\t",as.is=T)[,-c(2:6)]
inf=fread("sample_info.txt", header=T,data.table=F)
colnames(inf)[1]="ID"
d=merge(d,inf,by="ID")
d=d[d[,2175]=="HN" & d[,2178]=="Hybrid",]
#d=d[d[,4]=="LN",]
he=NULL
bl=as.data.frame(unique(d[,13]))
colnames(bl)="id"
for(i in 2:12)
{
  na=colnames(d)[i]
  data=d[,c(1,13:17,i)]
  data$genotype=factor(data$genotype)
  data$Block=factor(data$Block)
  data$sublock=factor(data$sublock)
  data$Type=factor(data$Type)
  colnames(data)[7]="y"
  fit=lmer(y~genotype+(1 | Block)+ (1|Block:sublock),data=data)
  x=summary(fit)
  BLUE <- as.data.frame(fixef(fit))
  rownames(BLUE)=gsub("genotype","",rownames(BLUE))
  a=bl[,1][!bl[,1]%in%rownames(BLUE)]
  da=na.omit(data)
  da1=da[da[,2]%in%a,]
  rownames(BLUE)[1]=as.character(unique(da1[,2]))
  b=data.frame(id=rownames(BLUE),Phe=c(BLUE[1,1],(BLUE[-1,1]+BLUE[1,1])))
  colnames(b)[2]=na
  bl=merge(bl,b,by="id",all=T)
  colnames(bl)[ncol(bl)]=na

  ###calculate H2
  fit=lmer(y~(1|genotype)+(1 | Block)+ (1|Block:sublock),data=data)
  x=summary(fit)
  v=as.data.frame(x$varcor)
  vg=v[1,4]
  ve=v[4,4]##Residual
  H=vg/(vg+ve/2)
  ##H=Vg/(Vg+Vge/y+Ve/yb)  #y:year  b: block Ve:Residual
  ##reference: Hallauer, A.R. and Miranda, J.B. (1988) Quantitative Genetics in Maize Breeding. Ames, IA: Iowa State University, 283.
  hee=c(na,H)
  he=rbind(he,hee)
}
colnames(he)=c("Trait","Heritability")
he=data.frame(he,pop="Hybird",N="HN")
bl=data.frame(bl,pop="Hybird",N="HN")
write.table(bl,file="data_blue_HN_hyb_result.txt",row.names = F,quote = F,sep="\t")
write.table(he,file="Heritability_HN_hyb_result.txt",row.names = F,quote = F,sep="\t")

#############################################################
d=read.table("ASV_GWAS_re839.txt",header = T,sep="\t",as.is=T)[,-c(2:6)]
inf=fread("sample_info.txt", header=T,data.table=F)
colnames(inf)[1]="ID"
d=merge(d,inf,by="ID")
d=d[d[,837]=="LN" & d[,840]=="Inbred",]
#d=d[d[,4]=="LN",]
he=NULL
bl=as.data.frame(unique(d[,836]))
colnames(bl)="id"
for(i in 2:835)
{
  na=colnames(d)[i]
  data=d[,c(1,836:840,i)]
  data$genotype=factor(data$genotype)
  data$Block=factor(data$Block)
  data$sublock=factor(data$sublock)
  data$Type=factor(data$Type)
  colnames(data)[7]="y"
  fit=lmer(y~genotype+(1 | Block)+ (1|Block:sublock),data=data)
  x=summary(fit)
  BLUE <- as.data.frame(fixef(fit))
  rownames(BLUE)=gsub("genotype","",rownames(BLUE))
  a=bl[,1][!bl[,1]%in%rownames(BLUE)]
  da=na.omit(data)
  da1=da[da[,2]%in%a,]
  rownames(BLUE)[1]=as.character(unique(da1[,2]))
  b=data.frame(id=rownames(BLUE),Phe=c(BLUE[1,1],(BLUE[-1,1]+BLUE[1,1])))
  colnames(b)[2]=na
  bl=merge(bl,b,by="id",all=T)
  colnames(bl)[ncol(bl)]=na

  ###calculate H2
  fit=lmer(y~(1|genotype)+(1|Block)+ (1|Block:sublock),data=data)
  x=summary(fit)
  v=as.data.frame(x$varcor)
  vg=v[1,4]
  ve=v[4,4]##Residual
  H=vg/(vg+ve/2)
  ##H=Vg/(Vg+Vge/y+Ve/yb)  #y:year  b: block Ve:Residual
  ##reference: Hallauer, A.R. and Miranda, J.B. (1988) Quantitative Genetics in Maize Breeding. Ames, IA: Iowa State University, 283.
  hee=c(na,H)
  he=rbind(he,hee)
}
colnames(he)=c("Trait","Heritability")
he=data.frame(he,pop="Inbred",N="LN")
bl=data.frame(bl,pop="Inbred",N="LN")
write.table(bl,file="data_blue_LN_inb_836.txt",row.names = F,quote = F,sep="\t")
write.table(he,file="Heritability_LN_inb_836.txt",row.names = F,quote = F,sep="\t")

##############################################################
d=read.table("2022root_data.txt",header = T,sep="\t",as.is=T)[,-c(2:6)]
inf=fread("sample_info.txt", header=T,data.table=F)
colnames(inf)[1]="ID"
d=merge(d,inf,by="ID")
d=d[d[,14]=="LN" & d[,17]=="Hybrid",]
#d=d[d[,4]=="LN",]
he=NULL
bl=as.data.frame(unique(d[,13]))
colnames(bl)="id"
for(i in 2:12)
{
  na=colnames(d)[i]
  data=d[,c(1,13:17,i)]
  data$genotype=factor(data$genotype)
  data$Block=factor(data$Block)
  data$sublock=factor(data$sublock)
  data$Type=factor(data$Type)
  colnames(data)[7]="y"
  fit=lmer(y~genotype+(1 | Block)+ (1|Block:sublock),data=data)
  x=summary(fit)
  BLUE <- as.data.frame(fixef(fit))
  rownames(BLUE)=gsub("genotype","",rownames(BLUE))
  a=bl[,1][!bl[,1]%in%rownames(BLUE)]
  da=na.omit(data)
  da1=da[da[,2]%in%a,]
  rownames(BLUE)[1]=as.character(unique(da1[,2]))
  b=data.frame(id=rownames(BLUE),Phe=c(BLUE[1,1],(BLUE[-1,1]+BLUE[1,1])))
  colnames(b)[2]=na
  bl=merge(bl,b,by="id",all=T)
  colnames(bl)[ncol(bl)]=na

  ###calculate H2
  fit=lmer(y~(1|genotype)+(1 | Block)+ (1|Block:sublock),data=data)
  x=summary(fit)
  v=as.data.frame(x$varcor)
  vg=v[1,4]
  ve=v[4,4]##Residual
  H=vg/(vg+ve/2)
  ##H=Vg/(Vg+Vge/y+Ve/yb)  #y:year  b: block Ve:Residual
  ##reference: Hallauer, A.R. and Miranda, J.B. (1988) Quantitative Genetics in Maize Breeding. Ames, IA: Iowa State University, 283.
  hee=c(na,H)
  he=rbind(he,hee)
}
colnames(he)=c("Trait","Heritability")
he=data.frame(he,pop="Hybird",N="LN")
bl=data.frame(bl,pop="Hybird",N="LN")
write.table(bl,file="data_blue_LN_hyb_result.txt",row.names = F,quote = F,sep="\t")
write.table(he,file="Heritability_LN_hyb_result.txt",row.names = F,quote = F,sep="\t")
