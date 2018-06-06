setwd("C:/Users/tmsavi/Desktop/NORWAY_tmp")

#blast columns (1-based):
#
#1. 	 qseqid 	 query (e.g., gene) sequence id
#2. 	 sseqid 	 subject (e.g., reference genome) sequence id
#3. 	 pident 	 percentage of identical matches
#4. 	 length 	 alignment length
#5. 	 mismatch 	 number of mismatches
#6. 	 gapopen 	 number of gap openings
#7. 	 qstart 	 start of alignment in query
#8. 	 qend 	 end of alignment in query
#9. 	 sstart 	 start of alignment in subject
#10. 	 send 	 end of alignment in subject
#11. 	 evalue 	 expect value
#12. 	 bitscore 	 bit score

#read TE-sequence lengths
elementlens=read.table("Tthy2_salmonFIltered_2.lib.fai",comment.char="%")
dim(elementlens)

#############################################################################################################
#
#			Non-TE host gene detection

#read known nucleotides per element
knownbases=read.csv2("Tthy2_salmonFIltered_2.acgt.txt",comment.char="%",sep="\t",header=T)
knownbases[is.na(knownbases)]=0
knownbases$tot=rowSums(knownbases[,-1])
knownbases$known=(knownbases$A+knownbases$C+knownbases$G+knownbases$T+knownbases$a+knownbases$c+knownbases$g+knownbases$t)
knownbases$unknown=knownbases$tot-knownbases$known
knownbases$elementlen=elementlens$V2[match(knownbases$sequence,elementlens$V1)]

#The following table contains the best hits of TE sequences to RepBase and/or Uniprot
#TE sequences with best hits to Uniprot are removed from further analysis
hostelements=read.table("best_multi_blast_score_grayling.txt",comment.char="%")
#Length of all hits
length(hostelements$V1[!duplicated(hostelements$V1)])
hostelements=hostelements[grep("swissprot_graylingrep.word7-blastx.txt",hostelements$V4),]
#Length of hits to exclude
length(hostelements$V1[!duplicated(hostelements$V1)])

#############################################################################################################
#
#			Repeat classification

#Read repet hits and select only ‘good’ BLASTN alignments (Wicker et al, trout repeat classification) that are not 
#detected as host genes in the previous step, have at least 80 bp long, have at least 80% sequence similarity 
#between query and subject and occupy at least 80% of the query repeat (after removing unknown nucleotides from query length)

repetnuc=read.table("repet_graylingrep.word7-blastn.txt",comment.char="%",stringsAsFactors=F)
repetnuc$elementlen=elementlens$V2[match(repetnuc$V1,elementlens$V1)]
repetnuc$unknown=knownbases$unknown[match(repetnuc$V1,knownbases$sequence)]
repetnuc$percelementlen=repetnuc$V4/(repetnuc$elementlen-repetnuc$unknown)
repetnuc$ishostelement=repetnuc$V1 %in% hostelements$V1
repetnuc=repetnuc[!duplicated(repetnuc$V1),] #best hit
dim(repetnuc) #1071 TE:s have hits
table(repetnuc$V4>=80)  #1012 are >=80 base pairs long alignment
table(repetnuc$V3>=80) #800 have >=80% identity
table(repetnuc$percelementlen>=0.8) #only 308 have >=80% of the sequence aligned, calculated by alignmentlength/(total seq. length-unknown bases)
repetnuc$isgood=repetnuc$V4>=80 & repetnuc$V3>=80 & repetnuc$percelementlen>=0.8 & repetnuc$ishostelement==FALSE
table(repetnuc$isgood) #in total, only 284 fulfill all the requirements
repetnuc$type="nuc"
repetnuc=repetnuc[repetnuc$isgood==TRUE & repetnuc$ishostelement==FALSE,]
table(repetnuc$isgood)
dim(repetnuc)

#If no BLASTN-based classification was found for a sequence, it was aligned RepBase using BLASTX and superfamily was assigned 
#if there was an existing hit with E-value less than 1e-10.

repetprot=read.table("repet_graylingrep.word7-blastx.txt",comment.char="%",stringsAsFactors=F) #all hits already with e-value < 1e-10
repetprot$elementlen=elementlens$V2[match(repetprot$V1,elementlens$V1)]
repetprot$unknown=knownbases$unknown[match(repetprot$V1,knownbases$sequence)]
repetprot$percelementlen=repetprot$V4/(repetprot$elementlen-repetprot$unknown)
repetprot$ishostelement=repetprot$V1 %in% hostelements$V1
repetprot$isgood=TRUE
repetprot=repetprot[!duplicated(repetprot$V1),]
repetprot$type="prot"
repetprot=repetprot[repetprot$ishostelement==FALSE,]
dim(repetprot)

#Combine annotations based on BLASTN and BLASTX searches (BLASTN annotations are prioritized)
TEanno=rbind(repetnuc,repetprot[!repetprot$V1 %in% repetnuc$V1,])
dim(TEanno) #total hits

#Parse TE class, order, superfamily and name from the RepBase name
TEanno$class=""
TEanno$order=""
TEanno$superfam=""
TEanno$name=""
for(i in 1:nrow(TEanno)) {
	splitname=strsplit(TEanno$V2[i],":")
	TEanno$name[i]=splitname[[1]][1]
	TEanno$class[i]=splitname[[1]][2]
	TEanno$order[i]=splitname[[1]][3]
	TEanno$superfam[i]=splitname[[1]][4]
}
TEanno$class=as.character(TEanno$class)
TEanno$order=as.character(TEanno$order)
TEanno$superfam=as.character(TEanno$superfam)
TEanno$name=as.character(TEanno$name)

TEanno$class[which(TEanno$class=="?")]="Unknown"
TEanno$order[which(TEanno$order=="?")]="Unknown"
TEanno$superfam[which(TEanno$superfam=="?")]="Unknown"
TEanno$name[which(TEanno$name=="?")]="Unknown"

write.table(TEanno,file="TEanno.txt",quote=F,sep="\t",row.names=F)
#############################################################################################################
#
#			Repeat abundancy in the European grayling genome



#Read RepeatMasker output and add repet classes
repeats=read.csv2("all.fasta.out.tab",skip=3, sep="\t",header=F,stringsAsFactors=F,comment.char="%")
#remove empty header related lines
repeats=repeats[!is.na(repeats$V15) &repeats$V15!="ID",]
#main and sub classes
#repeat coordinates starting from 0 are marked as (0) so convert to regular 0 
repeats$V6=sub("\\(","",repeats$V6)
repeats$V7=sub("\\(","",repeats$V7)
repeats$V12=sub("\\(","",repeats$V12)
repeats$V13=sub("\\(","",repeats$V13)
repeats$V14=sub("\\(","",repeats$V14)
repeats$V6=sub("\\)","",repeats$V6)
repeats$V7=sub("\\)","",repeats$V7)
repeats$V12=sub("\\)","",repeats$V12)
repeats$V13=sub("\\)","",repeats$V13)
repeats$V14=sub("\\)","",repeats$V14)
repeats$fullname=paste(repeats$V10,repeats$V11,sep="#")
#convert to numeric
for(i in c(1,2,3,4,6,7,12,13,14)) {repeats[,i]=as.numeric(repeats[,i])}
#calculate hit length
repeats$repeathitlength=unlist(apply(repeats[,c("V6","V7")],1,function(x){max(x)-min(x)}))

#Make a summary table of repeat coverages with added classification from RepBase
hitlentotals=aggregate(repeats$repeathitlength,list(repeats$fullname),sum)
hitlentotals$TEname=TEanno$name[match(hitlentotals$Group.1,TEanno$V1)]
hitlentotals$TEclass=TEanno$class[match(hitlentotals$Group.1,TEanno$V1)]
hitlentotals$TEorder=TEanno$order[match(hitlentotals$Group.1,TEanno$V1)]
hitlentotals$TEsuperfam=TEanno$superfam[match(hitlentotals$Group.1,TEanno$V1)]
hitlentotals$TEname[which(is.na(hitlentotals$TEname))]="Unknown"
hitlentotals$TEclass[which(is.na(hitlentotals$TEclass))]="Unknown"
hitlentotals$TEorder[which(is.na(hitlentotals$TEorder))]="Unknown"
hitlentotals$TEsuperfam[which(is.na(hitlentotals$TEsuperfam))]="Unknown"
write.table(hitlentotals,file="hitlentotals.txt",sep="\t",quote=F)

#############################################################################################################
#
#			Repeat abundancy in the Atlantic salmon genome

#Read RepeatMasker output and add repet classes
repeats.ssa=read.csv2("all.ssa.fasta.out.tab",skip=3, sep="\t",header=F,stringsAsFactors=F,comment.char="%")

#remove empty header related lines
repeats.ssa=repeats.ssa[!is.na(repeats.ssa$V15) &repeats.ssa$V15!="ID",]
#main and sub classes
#repeat coordinates starting from 0 are marked as (0) so convert to regular 0 
repeats.ssa$V6=sub("\\(","",repeats.ssa$V6)
repeats.ssa$V7=sub("\\(","",repeats.ssa$V7)
repeats.ssa$V12=sub("\\(","",repeats.ssa$V12)
repeats.ssa$V13=sub("\\(","",repeats.ssa$V13)
repeats.ssa$V14=sub("\\(","",repeats.ssa$V14)
repeats.ssa$V6=sub("\\)","",repeats.ssa$V6)
repeats.ssa$V7=sub("\\)","",repeats.ssa$V7)
repeats.ssa$V12=sub("\\)","",repeats.ssa$V12)
repeats.ssa$V13=sub("\\)","",repeats.ssa$V13)
repeats.ssa$V14=sub("\\)","",repeats.ssa$V14)
repeats.ssa$fullname=paste(repeats.ssa$V10,repeats.ssa$V11,sep="#")
#convert to numeric
for(i in c(1,2,3,4,6,7,12,13,14)) {repeats.ssa[,i]=as.numeric(repeats.ssa[,i])}

#Make a summary table of repeat coverages with added classification from RepBase
repeats.ssa$repeathitlength=unlist(apply(repeats.ssa[,c("V6","V7")],1,function(x){max(x)-min(x)}))
hitlentotals.ssa=aggregate(repeats.ssa$repeathitlength,list(repeats.ssa$fullname),sum)
hitlentotals.ssa$TEname=TEanno$name[match(hitlentotals.ssa$Group.1,TEanno$V1)]
hitlentotals.ssa$TEclass=TEanno$class[match(hitlentotals.ssa$Group.1,TEanno$V1)]
hitlentotals.ssa$TEorder=TEanno$order[match(hitlentotals.ssa$Group.1,TEanno$V1)]
hitlentotals.ssa$TEsuperfam=TEanno$superfam[match(hitlentotals.ssa$Group.1,TEanno$V1)]
hitlentotals.ssa$TEname[which(is.na(hitlentotals.ssa$TEname))]="Unknown"
hitlentotals.ssa$TEclass[which(is.na(hitlentotals.ssa$TEclass))]="Unknown"
hitlentotals.ssa$TEorder[which(is.na(hitlentotals.ssa$TEorder))]="Unknown"
hitlentotals.ssa$TEsuperfam[which(is.na(hitlentotals.ssa$TEsuperfam))]="Unknown"
#write.table(hitlentotals.ssa,file="hitlentotals.ssa.txt",sep="\t",quote=F)

#Total repeat content
sum(hitlentotals$x)/1485210005
sum(hitlentotals.ssa$x,na.rm=T)/2240221656

#############################################################################################################
#
#			Make a Repeat abundancy summary table based on RepBase classification (unknown elements are 
#			placed in "Unknown" category

hitlentotals.all=rbind(hitlentotals,hitlentotals.ssa)
hitlentotals.all=hitlentotals.all[!duplicated(paste(hitlentotals.all$TEclass,hitlentotals.all$TEorder,hitlentotals.all$TEsuperfam,hitlentotals.all$TEname)),]

graysums.tmp=aggregate(hitlentotals$x,list(hitlentotals$TEclass,hitlentotals$TEorder,hitlentotals$TEsuperfam,hitlentotals$TEname),sum)
salsums.tmp=aggregate(hitlentotals.ssa$x,list(hitlentotals.ssa$TEclass,hitlentotals.ssa$TEorder,hitlentotals.ssa$TEsuperfam,hitlentotals.ssa$TEname),sum,na.rm=T)
sum(graysums.tmp$x,na.rm=T)
sum(hitlentotals$x,na.rm=T)
sum(salsums.tmp$x,na.rm=T)
sum(hitlentotals.ssa$x,na.rm=T)

hitlentotals.all$gray.coverage.bp=graysums.tmp$x[match(paste(hitlentotals.all$TEclass,hitlentotals.all$TEorder,hitlentotals.all$TEsuperfam,hitlentotals.all$TEname),paste(graysums.tmp$Group.1,graysums.tmp$Group.2,graysums.tmp$Group.3,graysums.tmp$Group.4))]
hitlentotals.all$salmon.coverage.bp=salsums.tmp$x[match(paste(hitlentotals.all$TEclass,hitlentotals.all$TEorder,hitlentotals.all$TEsuperfam,hitlentotals.all$TEname),paste(salsums.tmp$Group.1,salsums.tmp$Group.2,salsums.tmp$Group.3,salsums.tmp$Group.4))]

hitlentotals.all$gray.coverage.bp[is.na(hitlentotals.all$gray.coverage.bp)]=0
hitlentotals.all$salmon.coverage.bp[is.na(hitlentotals.all$salmon.coverage.bp)]=0
hitlentotals.all$x=NULL
colnames(hitlentotals.all)[1]="namelong"

#check that total percentages are right
sum(hitlentotals.all$gray.coverage.bp)/1485210005
sum(hitlentotals.all$salmon.coverage.bp)/2240221656


hitlentotals.all$gray.coverage.perc=hitlentotals.all$gray.coverage.bp/1485210005
hitlentotals.all$salmon.coverage.perc=hitlentotals.all$salmon.coverage.bp/2240221656
hitlentotals.all=hitlentotals.all[order(-hitlentotals.all$gray.coverage.bp),] #order based on abundance in grayling

write.table(hitlentotals.all,file="wicker_classification_elements.txt",sep="\t",quote=F)

#############################################################################################################
#
#			Grayling-pike comparison figure
#

require(scales)
require(randomcoloR)

length(levels(as.factor(paste(hitlentotals.all$TEclass,hitlentotals.all$TEorder))))
orderpalette <- distinctColorPalette(11)

orders=levels(as.factor(paste(hitlentotals.all$TEclass,hitlentotals.all$TEorder)))
orders
orderpalette=c(
"darkgreen",
"darkolivegreen",
"darkolivegreen4",
"darkolivegreen2",
"gold",
"yellow",
"darkred",
"red",
"salmon3",
"pink",
"gray"
)
plot(1:11,1:11,col=orderpalette,pch=16,cex=3)
myxlims=range(log2(hitlentotals.all$gray.coverage.bp+1))
myxlims[1]=myxlims[1]-3
par(mfrow=c(1,1))
plot(log2(hitlentotals.all$gray.coverage.bp+1),log2(hitlentotals.all$salmon.coverage.bp+1),
	col=orderpalette[as.numeric(as.factor(paste(hitlentotals.all$TEclass,hitlentotals.all$TEorder)))],pch=16,
	xlab="log2(base pair abundance in European grayling +1)",ylab="log2(base pair abundance in Atlantic salmon +1)",xlim=myxlims,
	cex=1.5
)

#a lm with 0 intercept
actuallm=lm(log2(hitlentotals.all$salmon.coverage.bp+1)~-1+log2(hitlentotals.all$gray.coverage.bp+1))
namestoplot=abs(actuallm$resid)>(1.96*sd(actuallm$resid))
table(namestoplot)
abline(actuallm)

conf1y=actuallm$fit+1.96*sd(actuallm$resid)
conf2y=actuallm$fit-1.96*sd(actuallm$resid)
conf1x=actuallm$fit
conf1lm=lm(conf1y ~ conf1x)
conf2lm=lm(conf2y ~ conf1x)
#abline(conf1lm)
#abline(conf2lm)

polygon(
	c(-5,30,-5,-5),
	c(predict(conf1lm,newdata=data.frame(conf1x=-5)),
	predict(conf1lm,newdata=data.frame(conf1x=30)),
	predict(conf1lm,newdata=data.frame(conf1x=30)),
	predict(conf1lm,newdata=data.frame(conf1x=-5))),
	col=alpha("red",alpha=0.1),border=NA)
polygon(
	c(-5,30,30,-5),
	c(predict(conf2lm,newdata=data.frame(conf1x=-5)),
	predict(conf2lm,newdata=data.frame(conf1x=30)),
	predict(conf2lm,newdata=data.frame(conf1x=-5)),
	predict(conf2lm,newdata=data.frame(conf1x=-5))),
	col=alpha("blue",alpha=0.1),border=NA)


graysaldev=log2(hitlentotals.all$gray.coverage.bp+1)-log2(hitlentotals.all$salmon.coverage.bp+1)
text(log2(hitlentotals.all$gray.coverage.bp+1)[namestoplot],jitter(log2(hitlentotals.all$salmon.coverage.bp+1),amount=0.00001)[namestoplot],hitlentotals.all$TEname[namestoplot],pch=16,cex=0.8)
legend("bottomright",pch=16,col=orderpalette,legend=levels(as.factor(paste(hitlentotals.all$TEclass,hitlentotals.all$TEorder))),cex=1,bty="n")


#edít the output table (add residuals)
hitlentotals.all$residuals_coverage=actuallm$resid
hitlentotals.all$over1.96sd=abs(actuallm$resid)>(1.96*sd(actuallm$resid))

write.table(hitlentotals.all,file="wicker_classification_elements.txt",sep="\t",quote=F)

