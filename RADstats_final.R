setwd("D:/seafile/grayling/grayling/FrenchRAD")

parents=read.table("parents.merged.over9.counts.txt")
parents$length=parents$V3-parents$V2

#peaks 184 (and 94)
barplot(table(parents$length),main="parents combined",xlab="RAD length")   

#filter: must occur in both parents and length is 184 bp
parents=parents[parents$V4>1 & parents$length>=182 & parents$length<=186,]

#write.table(parents[,1:3],file="rad_locations.txt",sep="\t",quote=F,col.names=F,row.names=F)

##############################################
#
#   Get parental SNPs that have good coverage in both 
#    parents and have one het. and one hom. genotype or het. + het.


parent.snps.all=read.table("var.parent.flt.new.vcf",stringsAsFactors=F)

F1genos=unlist(lapply(strsplit(parent.snps.all$V10,":"),function(x){x[1]}))
M2genos=unlist(lapply(strsplit(parent.snps.all$V11,":"),function(x){x[1]}))

table(paste(F1genos,M2genos))

#ok genotypes
okgenos=c("0/0 0/1","0/1 0/0","0/1 1/1","1/1 0/1","1/1 1/2","1/2 1/1","0/1 0/1","1/2 1/2")

parent.snps.all=parent.snps.all[paste(F1genos,M2genos) %in% okgenos,]

F1cov=unlist(lapply(strsplit(parent.snps.all$V10,":"),function(x){x[3]}))
M2cov=unlist(lapply(strsplit(parent.snps.all$V11,":"),function(x){x[3]}))
F1cov=as.numeric(F1cov)
M2cov=as.numeric(M2cov)

par(mfrow=c(2,2))
barplot(tapply(F1cov,parent.snps.all$V1,mean),ylim=c(0,1000),main="F1 mean coverage")
barplot(tapply(M2cov,parent.snps.all$V1,mean),ylim=c(0,1000),main="M2 mean coverage")
barplot(tapply(F1cov,parent.snps.all$V1,var),main="F1 variance in coverage")
barplot(tapply(M2cov,parent.snps.all$V1,var),main="M2 variance in coverage")

#parents must have coverage between 9 and 300 bp
#write.table(file="parent_coverages.txt",data.frame(scaffold=parent.snps.all$V1,pos=parent.snps.all$V2,F1cov=F1cov,M2cov=M2cov),quote=F,row.names=F,sep="\t")
parent.snps.all=parent.snps.all[F1cov>9 & M2cov>9 & F1cov<300 & M2cov<300,]
F1cov=unlist(lapply(strsplit(parent.snps.all$V10,":"),function(x){x[3]}))
M2cov=unlist(lapply(strsplit(parent.snps.all$V11,":"),function(x){x[3]}))
F1cov=as.numeric(F1cov)
M2cov=as.numeric(M2cov)
summary(F1cov)
summary(M2cov)
boxplot(F1cov,M2cov)

#write.table(parent.snps.all[,1:2],file="snp_locations.txt",sep="\t",quote=F,col.names=F,row.names=F)

##############################################
#
#   Calculate how many scaffolds the SNPs cover

scaffoldswithSNPs=table(parent.snps.all$V1)
scaffoldswithSNPsMin2=scaffoldswithSNPs[table(parent.snps.all$V1)>=2]

#how many scaffolds have SNP(s)
length(scaffoldswithSNPs)
length(scaffoldswithSNPsMin2)

scafflens=read.table("Final_PilonAsm_ssalBroken_simplemask.fasta.fai",stringsAsFactors=F)

#how large fraction of scaffolds has SNPs

sum(scafflens$V2[scafflens$V1 %in% names(scaffoldswithSNPs)])/sum(scafflens$V2)
#55%

sum(scafflens$V2[scafflens$V1 %in% names(scaffoldswithSNPsMin2)])/sum(scafflens$V2)
#41%

##############################################
#
#   Get offspring genotypes of selected SNPs

offspring_geno=read.table("var.offspring.flt.new.vcf",stringsAsFactors=F)
##############################################
#
#   filter SNPs based on offspring coverage

#remove low coverage genotypes
for(i in 10:ncol(offspring_geno)) {
	offspring_geno[,i][as.numeric(unlist(lapply(strsplit(offspring_geno[,i],":"),function(x){x[3]})))<8]=NA
}

#count how many missing genotypes in each sample,
#remove offspring that have too low coverage

offspring_hascov=apply(offspring_geno[,10:ncol(offspring_geno)],2,function(x){sum(is.na(x))})
offspring_hascov
which(offspring_hascov>1000)

offspring_geno.flt=offspring_geno[,-c(9+which(offspring_hascov>1000))]

##############################################
#
#   filter SNPs based on odd genotypes and fix 1-2 encoded

#count how many missing genotypes in each sample after filtering
offspring_hascov.flt=apply(offspring_geno.flt[,10:ncol(offspring_geno.flt)],2,function(x){sum(is.na(x))})
offspring_hascov.flt

#parse genotypes
for(i in 10:ncol(offspring_geno.flt)) {
	offspring_geno.flt[,i]=unlist(lapply(strsplit(offspring_geno.flt[,i],":"),function(x){x[1]}))
}
table(unlist(offspring_geno.flt[,10:ncol(offspring_geno.flt)]))

#remove odd genotypes

#replace unknown genotypes with NA
offspring_geno.flt[,10:ncol(offspring_geno.flt)][offspring_geno.flt[,10:ncol(offspring_geno.flt)]=="./."]=NA

#change SNPs encoded with alleles 1 and 2 to be encoded with alleles 0 and 1
rowstochange=which(offspring_geno.flt$V10=="1/2" & offspring_geno.flt$V11=="1/1")
offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)][offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)]=="1/2"]="0/1"
offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)][offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)]=="2/2"]="0/0"
rowstochange=which(offspring_geno.flt$V10=="1/1" & offspring_geno.flt$V11=="1/2")
offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)][offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)]=="1/2"]="0/1"
offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)][offspring_geno.flt[rowstochange,10:ncol(offspring_geno.flt)]=="2/2"]="0/0"
table(unlist(offspring_geno.flt[,10:ncol(offspring_geno.flt)]))

#remove rows with remaining allele 1/2 (triallelics)
genostoremove=apply(offspring_geno.flt[,10:ncol(offspring_geno.flt)],1,function(x){sum(x=="1/2",na.rm=T)})
#from table command it can be seen that one full SNP and one odd individual should be removed
offspring_geno.flt=offspring_geno.flt[genostoremove<=1,]
table(unlist(offspring_geno.flt[,10:ncol(offspring_geno.flt)]))
offspring_geno.flt[offspring_geno.flt=="1/2"]=NA

#remove also 0/2 alleles (triallelics)
genostoremove=apply(offspring_geno.flt[,10:ncol(offspring_geno.flt)],1,function(x){sum(x=="0/2",na.rm=T)})
table(genostoremove)
offspring_geno.flt=offspring_geno.flt[genostoremove==0,]
table(unlist(offspring_geno.flt[,10:ncol(offspring_geno.flt)]))

#remove also one odd 2/2 geno
offspring_geno.flt[offspring_geno.flt=="2/2"]=NA

#double check parental combinations
table(paste(offspring_geno.flt$V10,offspring_geno.flt$V11))
offspring_geno.flt=offspring_geno.flt[paste(offspring_geno.flt$V10,offspring_geno.flt$V11)!="0/0 0/0",]
offspring_geno.flt=offspring_geno.flt[paste(offspring_geno.flt$V10,offspring_geno.flt$V11)!="1/1 1/1",]
offspring_geno.flt=offspring_geno.flt[!is.na(offspring_geno.flt$V11),]
table(paste(offspring_geno.flt$V10,offspring_geno.flt$V11))

##############################################
#
#   convert remaining genotypes to numeric
#

offspring_geno.num=offspring_geno.flt
for(i in 10:ncol(offspring_geno.num)) {
	offspring_geno.num[,i][offspring_geno.num[,i]=="0/0"]="0"
	offspring_geno.num[,i][offspring_geno.num[,i]=="0/1"]="1"
	offspring_geno.num[,i][offspring_geno.num[,i]=="1/1"]="2"
	offspring_geno.num[,i]=as.numeric(offspring_geno.num[,i])

}
table(is.na(offspring_geno.flt[,10:ncol(offspring_geno.flt)]))
summary(unlist(offspring_geno.num[,10:ncol(offspring_geno.num)]))

#sum of different genotypes
hom=apply(offspring_geno.num[,12:ncol(offspring_geno.num)],1,function(x){sum(x==0 |  x==2,na.rm=T)})
het=apply(offspring_geno.num[,12:ncol(offspring_geno.num)],1,function(x){sum(x==1,na.rm=T)})

boxplot(het/(hom+het),ylab="heterozygosity")
summary(het/(hom+het),ylab="heterozygosity")

##############################################
#
#   filter SNPs based on expected offspring genotypes

genofreqs=data.frame(type=rep("",nrow(offspring_geno.num)))

#Final filtering based on Mendelian segregation of the offspring genotypes

#in any case, 50% are heterozygotes. 
#"hom0het case" -> 50% 0-homozygote; "hom2het case" -> 50% 2-homozygote; "hethet case" -> 25% 0-homozygote + 25% 2-homozygote.

genofreqs$type=as.character(genofreqs$type)
genofreqs$type[offspring_geno.num$V10==0 & offspring_geno.num$V11==1]="hom0het"
genofreqs$type[offspring_geno.num$V10==1 & offspring_geno.num$V11==0]="hom0het"
genofreqs$type[offspring_geno.num$V10==2 & offspring_geno.num$V11==1]="hom2het"
genofreqs$type[offspring_geno.num$V10==1 & offspring_geno.num$V11==2]="hom2het"
genofreqs$type[offspring_geno.num$V10==1 & offspring_geno.num$V11==1]="hethet"
genofreqs$hom0=apply(offspring_geno.num[,12:ncol(offspring_geno.num)],1,function(x){sum(x==0,na.rm=T)})
genofreqs$het=apply(offspring_geno.num[,12:ncol(offspring_geno.num)],1,function(x){sum(x==1,na.rm=T)})
genofreqs$hom2=apply(offspring_geno.num[,12:ncol(offspring_geno.num)],1,function(x){sum(x==2,na.rm=T)})
hom0hetExp=c(0.4999,0.4999,0.0002)
hethetExp=c(0.25,0.5,0.25)
hom2hetExp=c(0.0002,0.4999,0.4999)

require(HardyWeinberg)
genofreqs$hom0hetChi=apply(genofreqs[,2:4],1,function(x){chisq.test(x, p=hom0hetExp)$p.value})
genofreqs$hethetChi=apply(genofreqs[,2:4],1,function(x){chisq.test(x, p=hethetExp)$p.value})
genofreqs$hom2hetChi=apply(genofreqs[,2:4],1,function(x){chisq.test(x, p=hom2hetExp)$p.value})

#genofreqs$ok=0
genofreqs$fdr[genofreqs$type=="hom0het"]=genofreqs$hom0hetChi[genofreqs$type=="hom0het"]
genofreqs$fdr[genofreqs$type=="hethet"]=genofreqs$hethetChi[genofreqs$type=="hethet"]
genofreqs$fdr[genofreqs$type=="hom2het"]=genofreqs$hom2hetChi[genofreqs$type=="hom2het"]
genofreqs$fdr=p.adjust(genofreqs$fdr)
(table(genofreqs$fdr>0.05))

#finally filter the genotypes not following mendelian segregation
offspring_geno.num=offspring_geno.num[genofreqs$fdr>0.05,]
dim(offspring_geno.num)

#Final heterozygosity for offspring
offspring_mean_Hobs=apply(offspring_geno.num[,12:120],2,function(x){sum(x==1,na.rm=T)/sum(!is.na(x))})
hist(offspring_mean_Hobs)
t.test(offspring_mean_Hobs,mu=0.5)  #P=0.1175
sd(offspring_mean_Hobs)
#write.table(offspring_geno.num[,-c(6,7,8)],"genos.new.txt",quote=F)

offspring_geno.lepmap=offspring_geno.num[,10:ncol(offspring_geno.num)]
for(i in 1:ncol(offspring_geno.lepmap)){offspring_geno.lepmap[,i]=as.character(offspring_geno.lepmap[,i]) }
offspring_geno.lepmap[offspring_geno.lepmap=="0"]="1 1"
offspring_geno.lepmap[offspring_geno.lepmap=="1"]="1 2"
offspring_geno.lepmap[offspring_geno.lepmap=="2"]="2 2"
offspring_geno.lepmap[is.na(offspring_geno.lepmap)]="0 0"
table(unlist(offspring_geno.lepmap))
indnames=c("F2","M1","F100","F105","F106","F111","F116","F12","F127","F135","F139","F147","F148","F150","F151","F152","F156","F162","F163","F164","F171","F172","F18","F186","F196","F198","F205","F212","F224","F227","F234","F235","F24","F245","F246","F247","F252","F253","F256","F26","F264","F272","F275","F277","F293","F301","F302","F305","F307","F309","F320","F323","F340","F344","F345","F350","F36","F39","F4","F48","F5","F54","F60","F65","F72","F76","F80","F83","F86","F91","F94","M10","M122","M131","M141","M15","M16","M17","M181","M19","M190","M191","M20","M201","M203","M213","M22","M222","M225","M231","M239","M240","M242","M258","M265","M276","M289","M3","M30","M304","M336","M34","M352","M37","M38","M47","M57","M62","M7","M73","M75","M8","M81","M95","M97")
indnames=indnames[-c(which(offspring_hascov>1000))]
colnames(offspring_geno.lepmap)=indnames
indsex=substr(indnames,1,1)  #1=MALE, 2=FEMALE
indsex[indsex=="M"]="1"
indsex[indsex=="F"]="2"
table(indsex)
offspring_geno.lepmap=t(offspring_geno.lepmap)
offspring_geno.lepmap=cbind(
	rep("FrenchRAD",nrow(offspring_geno.lepmap)),
	indnames,
	c(0,0,rep("M2",nrow(offspring_geno.lepmap)-2)),
	c(0,0,rep("F1",nrow(offspring_geno.lepmap)-2)),
	indsex,
	rep("0",nrow(offspring_geno.lepmap)),
	offspring_geno.lepmap
)
#write.table(offspring_geno.lepmap,file="lepmap.new.txt",quote=F,col.names=F,row.names=F,sep="\t")

##############################################
#
#   PCA

offspring_geno.num.full=offspring_geno.num[apply(offspring_geno.num[,10:ncol(offspring_geno.num)],1,function(x){sum(is.na(x))})==0,]
pca=prcomp(t(offspring_geno.num.full[,10:ncol(offspring_geno.num.full)]))
summary(pca)
plot(pca$x[,1],pca$x[,2],col=c("red","blue",rep("white",nrow(pca$x)-2)),cex=3,xlab="PC1 4.7%",ylab="PC2 4.2%")
points(pca$x[1:2,1],pca$x[1:2,2],,col=c("red","blue"),cex=3,lwd=2)
text(pca$x[,1],pca$x[,2],c("F","M",rep("o",nrow(pca$x)-2)))
