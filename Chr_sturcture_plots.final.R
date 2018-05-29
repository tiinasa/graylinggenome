##########################################################################################################
##########################################################################################################
##########################################################################################################

			
			Chromosome structure plots for European grayling genome

			May 2018

			Tiina Sävilammi
			tmsavi@utu.fi


##########################################################################################################
##########################################################################################################
##########################################################################################################


#Set working directory

setwd("C:/Users/tmsavi/Desktop/NORWAY_tmp")

#Files used:
# - ncbi-acc.salmonchr.txt
# - Maker_filtered_m_filter2.gff
# - LGs_esox_converter.tabular
# - GraylingBED_Oct2017withGaps.interval
# - LepmapRemappedOct2017.txt
# - chromosome_conversion_table.txt
# - Tth_LGs_Oct2017.100KwindowACGT.txt
# - all.fasta.out.tab
# - grayling_recip_pep.sw.txt
# - grayling_pike_recip_pep.sw.txt
# - salmon_pike_recip_pep.sw.txt
# - Salmo_salar-annotation-cds.gff3
# - nucmerjelly_sal_mum.coords
# - nucmerjelly_pike_mum.pikesett.coords2
# - sdy_loci.txt
# - nucmer.[graylingchrbasename].[1|2|3].[salmonchrname].coords

#Parameters to set:
geneannowindow=1000000
gcwindow=1000000
repeatwindow=100000

#Linkage map parameters:
binmarkersby=6
labelstoshow=3
##########################################################################################################
#
#                 Read NCBI-to-ssa naming converter for salmon chromosomes (for mummer matching)

salmonconvert=read.table("ncbi-acc.salmonchr.txt",stringsAsFactors=F)

##########################################################################################################
#
#                 Read gene annotations

geneanno=read.csv2("Maker_filtered_m_filter2.gff",stringsAsFactors=F,sep="\t",header=F,skip=1)
dim(geneanno) #should be 1.106.688 annotation lines
geneanno=geneanno[which(geneanno$V3=="gene"),]
geneanno$ID=sub(".*ID=","",geneanno$V9)
geneanno$ID=sub(";.*","",geneanno$ID)
geneanno$note=sub(".*Note=","",geneanno$V9)
geneanno$note=sub(";.*","",geneanno$note)
geneanno=geneanno[which(geneanno$V1!="unmapped_scaffolds"),]
geneanno=geneanno[order(geneanno$V1,geneanno$V4),]
geneanno$win=cut(geneanno$V4,seq(1,max(geneanno$V4)+1,geneannowindow)) 
geneanno$winNum=as.numeric(geneanno$win)*geneannowindow-geneannowindow+1
annocounts=aggregate(geneanno$ID,list(geneanno$V1,geneanno$winNum),length)
colnames(annocounts)=c("chr","pos","genecount")
annocounts=annocounts[order(annocounts$chr,annocounts$pos),]
boxplot(annocounts$genecount~annocounts$chr,las=2,cex.axis=0.8)

##########################################################################################################
#
#                 Read homolog table for grayling, salmon and pike

homeologs=read.table("LGs_esox_converter.tabular",stringsAsFactors=F,sep="\t")
homeologs$ssa=sub("Ssa","ssa",substr(homeologs$V3,1,5))
homeologs$omy=sub("Omy","omy",substr(homeologs$V12,1,5))

homeologs$keysal=paste(homeologs$V1,homeologs$ssa) #keys with grayling-salmon chr
homeologs$keypike=paste(homeologs$V1,homeologs$V6) #keys with grayling-pike chr
homeologs$keytrout=paste(homeologs$V1,homeologs$omy) #keys with grayling-trout chr

##########################################################################################################
#
#                 Read bed file of scaffold ordering

bed=read.table("GraylingBED_Oct2017withGaps.interval",sep="\t")
dim(bed)
colnames(bed)[4]="chr"
bed$V3=bed$V3+100
bed=bed[bed$V1!="gap",]
bed$start=0		#V8=start coord, V9=end coord
bed$end=bed$V3[1]
#make incremental coordinates
for(i in 2:nrow(bed)) {
	#if chromosome continues then add to previous coordinates
	if(bed$chr[i]==bed$chr[i-1]) {
		bed$start[i]=bed$end[i-1]
		bed$end[i]=bed$V3[i]+bed$end[i-1]
	}
	#if chromosome changes then take values from scaffold coordinates
	else {
		bed$start[i]=0
		bed$end[i]=bed$V3[i]
	}
}
##########################################################################################################
#
#                  Linkage map (final)

map=read.table("LepmapRemappedOct2017.txt",header=T,sep="\t",stringsAsFactors=T)
colnames(map)[3]="Marker"

#generate marker-contig-chromosome name converter for new names
marker_contig_lg_pikelg=data.frame(marker=map$Marker,contig=map$Contig)
marker_contig_lg_pikelg$pikelg=bed$chr[match(map$Contig,bed$V1)]

head(marker_contig_lg_pikelg)

#Add representative marker for each bed entry (TODO needed?)
bed$Marker=map$Marker[match(bed$V1,map$Contig)]
table(!is.na(bed$Marker))

head(map)
table(map$Chr)
bed$femalepos=map$V5[match(bed$V1,map$Contig)]  #from final map 

##########################################################################################################
#
#                  get grayling chromosome names and make a cyto object for circlize

cytogray=bed[order(bed$chr,-bed$end),c("chr","V9")]
cytogray=cytogray[!duplicated(cytogray$chr),]
cytogray=data.frame(Chromosome=cytogray$chr,ChromStart=1,ChromEnd=cytogray$V9)
cytogray$Band=""
cytogray$Stain=""
cytogray=cytogray[order(as.character(cytogray$Chromosome)),]
head(cytogray)

##########################################################################################################
#
#                  get salmon chromosome names and make cyto object

cytosal=read.table("chromosome_conversion_table.txt",header=T,stringsAsFactors=F) #NC to ssa conversion
cytosal=data.frame(Chromosome=cytosal$Ssa,ChromStart=1,ChromEnd=cytosal$Length)
cytosal$Band=""
cytosal$Stain=""
cytosal=cytosal[order(as.character(cytosal$Chromosome)),]
head(cytosal)
##########################################################################################################
#
#                             Read in ACGT-content
#

nucleotidecontent=read.table("Tth_LGs_Oct2017.100KwindowACGT.txt",stringsAsFactors=F,header=T)
dim(nucleotidecontent)
nucleotidecontent$HEAD=sub(">","",nucleotidecontent$HEAD)
nucleotidecontent=nucleotidecontent[order(nucleotidecontent$HEAD,nucleotidecontent$WINDOWSTART),]

nucleotidecontent$Afreq=nucleotidecontent$A/rowSums(nucleotidecontent[,c("A","C","G","T")])
nucleotidecontent$Cfreq=nucleotidecontent$C/rowSums(nucleotidecontent[,c("A","C","G","T")])
nucleotidecontent$Gfreq=nucleotidecontent$G/rowSums(nucleotidecontent[,c("A","C","G","T")])
nucleotidecontent$Tfreq=nucleotidecontent$T/rowSums(nucleotidecontent[,c("A","C","G","T")])
nucleotidecontent$ATfreq=rowSums(nucleotidecontent[,c("A","T")])/rowSums(nucleotidecontent[,c("A","C","G","T")])
nucleotidecontent$CGfreq=rowSums(nucleotidecontent[,c("C","G")])/rowSums(nucleotidecontent[,c("A","C","G","T")])
nucleotidecontent$CG=nucleotidecontent$C+nucleotidecontent$G
nucleotidecontent$AT=nucleotidecontent$A+nucleotidecontent$T

nucleotidecontent$WINDOWEND=nucleotidecontent$WINDOWSTART+100000-1
nucleotidecontent$win=cut(nucleotidecontent$WINDOWSTART,seq(0,max(nucleotidecontent$WINDOWSTART)+1,gcwindow),include.lowest=T) #make adjustable windows
nucleotidecontent$winNum=as.numeric(nucleotidecontent$win)*gcwindow-gcwindow+1
gcfreqcounts=aggregate(list(nucleotidecontent$CG,nucleotidecontent$AT),list(nucleotidecontent$HEAD,nucleotidecontent$winNum),sum)
colnames(gcfreqcounts)[3:4]=c("CG","AT")
gcfreqcounts$cgfreq=gcfreqcounts$CG/(gcfreqcounts$CG+gcfreqcounts$AT)
head(gcfreqcounts)

boxplot(nucleotidecontent$CGfreq~nucleotidecontent$HEAD)
boxplot(gcfreqcounts$cgfreq~gcfreqcounts$Group.1)

##########################################################################################################################
#
#                             Read in Repeat Masker grayling-specific library element hits 
#

repeats=read.csv2("all.fasta.out.tab",skip=3, sep="\t",header=F,stringsAsFactors=F)

#remove empty header related lines
repeats=repeats[!is.na(repeats$V15) &repeats$V15!="ID",]
#main and sub classes
repeats$mainclass=sub("/.*","",repeats$V11)
repeats$subclass=sub(".*/","",repeats$V11)
repeats$subclass[repeats$subclass==repeats$mainclass]="Undefined"
repeats$subclass=paste(repeats$mainclass,repeats$subclass,sep="/")
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
#convert to numeric
for(i in c(1,2,3,4,6,7,12,13,14)) {repeats[,i]=as.numeric(repeats[,i])}
#calculate hit length
repeats$repeathitlength=unlist(apply(repeats[,c("V12","V13")],1,function(x){max(x)-min(x)}))
#count sub and mainclass occurences
repeatcounts=table(repeats$subclass)
repeatcounts=repeatcounts[order(-repeatcounts)]
repeatlength=tapply(repeats$repeathitlength,repeats$subclass,sum)
repeatlength=repeatlength[match(names(repeatcounts),names(repeatlength))]
repeatclasscounts=table(repeats$mainclass)
repeatclasscounts=repeatclasscounts[order(-repeatclasscounts)]
repeatclasslength=tapply(repeats$repeathitlength,repeats$mainclass,sum)
repeatclasslength=repeatclasslength[match(names(repeatclasscounts),names(repeatclasslength))]
table(paste(repeats$mainclass,repeats$subclass))
#merge LINE/RTE class to LINE/RTE-X
repeats$subclass[which(repeats$subclass=="LINE/RTE")]="LINE/RTE-X"
repeats$elementlen=(repeats$V7-repeats$V6)*(100-repeats$V2)  #element length corrected by divergence
#add windows and window starts to elements
repeats$win=cut(repeats$V6,seq(0,max(repeats$V6)+1,repeatwindow),includeLowest=T) #make windows
repeats$winNum=as.numeric(repeats$win)*repeatwindow-repeatwindow+1

##########################################################################################################################
#
#			Subset Tc1 mariner -elements
#

tc1mar=repeats[which(repeats$subclass=="DNA/TcMar-Tc1"),c("V5","V6","V7","elementlen","V2","win","winNum")]
colnames(tc1mar)=c("chr","start","end","elementlen","divergenceperc","win","winNum")
#count how many grayling chromosomes covered
#summed elementlength for each tc1 mariner bin
tc1marbins=aggregate(tc1mar$elementlen,list(tc1mar$chr,tc1mar$winNum),sum)
colnames(tc1marbins)=c("Chr","Pos","Len")
tc1marbins=tc1marbins[order(tc1marbins$Chr,tc1marbins$Pos),]
#make element length bins compatible for cytogray
tc1marbins$Chr=factor(tc1marbins$Chr,levels=cytogray$Chromosome)
tc1marbins=tc1marbins[order(tc1marbins$Chr,tc1marbins$Pos),]
rm(tc1mar)

##########################################################################################################################
#
#                             Read in RTE-X -elements
#

rtex=repeats[which(repeats$subclass=="LINE/RTE-X"),c("V5","V6","V7","elementlen","V2","win","winNum")]
colnames(rtex)=c("chr","start","end","elementlen","divergenceperc","win","winstart")
tail(rtex)
#summed elementlength for each tc1 mariner bin
rtexbins=aggregate(rtex$elementlen,list(rtex$chr,rtex$winstart),sum)
colnames(rtexbins)=c("Chr","Pos","Len")
head(rtexbins)
rtexbins=rtexbins[order(rtexbins$Chr,rtexbins$Pos),]
#make element length bins compatible for cytogray
rtexbins$Chr=factor(rtexbins$Chr,levels=cytogray$Chromosome)
rtexbins=rtexbins[order(rtexbins$Chr,rtexbins$Pos),]
rm(rtex)

##########################################################################################################################
#
#			Find peak positions
#

tc1marbinsMax=lapply(split(tc1marbins,tc1marbins$Chr),function(x){ 
	xlowess=lowess(x$Len~x$Pos,f=0.2); 
	xpositions=xlowess$x; 
	ypositions=xlowess$y; 
	#xpositions=x$Pos
	#ypositions=x$Len
	xpositions=xpositions[order(-ypositions)]; 
	names(xpositions[1])=x$Chr[1]
	return(xpositions[1])
})
tc1marbinsMax=data.frame(t(data.frame(tc1marbinsMax)))
tc1marbinsMax$Chrlen=cytogray$ChromEnd[match(rownames(tc1marbinsMax),cytogray$Chromosome)]
tc1marbinsMax$oneArm=tc1marbinsMax$Chrlen-tc1marbinsMax[,1]
tc1marbinsMax$longArmPerc=tc1marbinsMax$oneArm/tc1marbinsMax$Chrlen
tc1marbinsMax$longArmPerc[tc1marbinsMax$longArmPerc<0.5]=1-tc1marbinsMax$longArmPerc[tc1marbinsMax$longArmPerc<0.5]
tc1marbinsMax$longShortRatio=tc1marbinsMax$longArmPerc/(1-tc1marbinsMax$longArmPerc)
tc1marbinsMax$chrType="M"
tc1marbinsMax$chrType[tc1marbinsMax$longShortRatio>1.67]="SM"
tc1marbinsMax$chrType[tc1marbinsMax$longShortRatio>3]="T"
table(tc1marbinsMax$chrType)
tc1marbinsMax$longShortRatio[tc1marbinsMax$longShortRatio>10]=10
boxplot(tc1marbinsMax$longShortRatio~tc1marbinsMax$chrType)

rtexbinsMax=lapply(split(rtexbins,rtexbins$Chr),function(x){ 
	xlowess=lowess(x$Len~x$Pos,f=0.2); 
	xpositions=xlowess$x; 
	ypositions=xlowess$y;
	xpositions=xpositions[order(-ypositions)]; 
	names(xpositions[1])=x$Chr[1]
	return(xpositions[1:2])
})
rtexbinsMax=data.frame(t(data.frame(rtexbinsMax)))

##########################################################################################################################
#
#                             calculate dN/dS
#
ssa.grayling_genes=read.table("grayling_recip_pep.sw.txt",stringsAsFactors=F)
ssa.grayling_genes$V2=sub("grayling_","",ssa.grayling_genes$V2)
pike.grayling_genes=read.table("grayling_pike_recip_pep.sw.txt",stringsAsFactors=F)
pike.salmon_genes=read.table("salmon_pike_recip_pep.sw.txt",stringsAsFactors=F)
gff.sal=read.csv2("Salmo_salar-annotation-cds.gff3",sep="\t",stringsAsFactors=F,header=F)
gff.sal=gff.sal[which(gff.sal$V3=="CDS"),]
gff.sal$V9=sub(".*Genbank:","",gff.sal$V9)
gff.sal$V9=sub(";.*","",gff.sal$V9)
gff.sal=gff.sal[!duplicated(gff.sal$V9),]
gff.all=read.csv2("Maker_filtered_m_filter2.gff",sep="\t",stringsAsFactors=F,header=F)
gff.pike.gray=gff.all[which(gff.all$V3=="CDS"),]
gff.pike.gray$transcript=sub(".*Parent=","",gff.pike.gray$V9)
gff.pike.gray$transcript=sub(";.*","",gff.pike.gray$transcript)
gff.pike.gray$pike=pike.grayling_genes$V1[match(gff.pike.gray$transcript,pike.grayling_genes$V2)]
gff.pike.gray$score=pike.grayling_genes$V3[match(gff.pike.gray$transcript,pike.grayling_genes$V2)]
gff.pike.gray$hitlen=pike.grayling_genes$V4[match(gff.pike.gray$transcript,pike.grayling_genes$V2)]
gff.pike.gray$salmon=pike.salmon_genes$V2[match(gff.pike.gray$pike,pike.salmon_genes$V1)]
gff.pike.gray$salscore=pike.salmon_genes$V3[match(gff.pike.gray$pike,pike.salmon_genes$V1)]
gff.pike.gray$salhitlen=pike.salmon_genes$V4[match(gff.pike.gray$pike,pike.salmon_genes$V1)]
gff.pike.gray$salchr=gff.sal$V1[match(gff.pike.gray$salmon,gff.sal$V9)]
gff.pike.gray$salpos=gff.sal$V4[match(gff.pike.gray$salmon,gff.sal$V9)]

gff.pike.gray=gff.pike.gray[order(gff.pike.gray$transcript,-rowSums(gff.pike.gray[,c("salscore","score")],na.rm=T)),]
gff.pike.gray=gff.pike.gray[!duplicated(gff.pike.gray$transcript),]
gff.pike.gray=gff.pike.gray[order(gff.pike.gray$V1,gff.pike.gray$V4),]

gff.pike.gray$correcthomeolog=homeologs$ssa[match(gff.pike.gray$salchr,homeologs$ssa)]
gff.pike.gray=gff.pike.gray[which(gff.pike.gray$salchr==gff.pike.gray$correcthomeolog),]
gff.pike.gray=gff.pike.gray[which(!is.na(gff.pike.gray$pike) & !is.na(gff.pike.gray$salmon)),]

#make multiple sequence alignment & calculate divergence
###source("http://www.bioconductor.org/biocLite.R")
###biocLite("msa")
require("msa")
require(Biostrings)

#get grayling sequences

#cds
graycdsfile <- "grayling_cds_from_Tthy2_maker_filter2.fa"
salcdsfile <- "GCF_000233375.1_ICSASG_v2_cds_from_genomic.fna"
pikecdsfile <- "GCF_000721915.3_Eluc_V3_cds_from_genomic.fna"

#peps
graypepfile <- "Tthy2_maker_filter2_proteins.fasta"
salpepfile <- "GCF_000233375.1_ICSASG_v2_protein.faa"
pikepepfile <- "GCF_000721915.3_Eluc_V3_protein.faa"

graycds <- readAAStringSet(graycdsfile)
salcds <- readAAStringSet(salcdsfile)
pikecds <- readAAStringSet(pikecdsfile)


graypep <- readAAStringSet(graypepfile)
salpep <- readAAStringSet(salpepfile)
pikepep <- readAAStringSet(pikepepfile)

names(graycds)=sub("\\s.*","",names(graycds))
names(salcds)=sub(".*protein_id=","",names(salcds))
names(salcds)=sub("].*","",names(salcds))
names(pikecds)=sub(".*protein_id=","",names(pikecds))
names(pikecds)=sub("].*","",names(pikecds))
table(gff.pike.gray$pike %in% names(pikecds))
table(gff.pike.gray$salmon %in% names(salcds))

names(salpep)=sub("\\s.*","",names(salpep))
names(pikepep)=sub("\\s.*","",names(pikepep))
table(gff.pike.gray$pike %in% names(pikepep))
table(gff.pike.gray$salmon %in% names(salpep))
head(names(salpep))

pepgappedseq=thisseqs$seq[3]
cdsseq=as.character(allcds[[3]])
out=backtr(pepgappedseq,cdsseq)

paste(seqinr::translate(unlist(strsplit(as.character(tolower(salcdsseq)),split=""))),collapse="")
as.character(salseq)




require(seqinr)
require(ape)
require(msa)
dnds=list()

#for(i in 1:10) {

for(i in 1:nrow(gff.pike.gray)) {
	#old	
	#grayseqx=graypep[names(graypep)==gff.pike.gray$transcript[i]]
	#salseqx=salpep[names(salpep)==gff.pike.gray$salmon[i]]
	#pikeseqx=pikepep[names(pikepep)==gff.pike.gray$pike[i]]

	graycdsseq=graycds[names(graycds)==gff.pike.gray$transcript[i]]
	salcdsseq=salcds[names(salcds)==gff.pike.gray$salmon[i]]
	pikecdsseq=pikecds[names(pikecds)==gff.pike.gray$pike[i]]
	allcds=list(pikecdsseq,graycdsseq,salcdsseq)

	grayseq=seqinr::translate(unlist(strsplit(as.character(graycdsseq),split="")))
	salseq=seqinr::translate(unlist(strsplit(as.character(salcdsseq),split="")))
	pikeseq=seqinr::translate(unlist(strsplit(as.character(pikecdsseq),split="")))

	#old
	write.fasta(list(pikeseq,grayseq,salseq),names=c("pike","gray","sal"),file="temp.fa")
allpeps=readAAStringSet("temp.fa")
if(sum(grayseq=="*")>1) {print(paste("warning: many stops grayling",i))}
if(sum(salseq=="*")>1) {print(paste("warning: many stops salmon",i))}
if(sum(pikeseq=="*")>1) {print(paste("warning: many stops pike",i))}


	thismsa=msa(allpeps,method="Muscle",gapExtension=2,gapOpening=4,order="input")
	#thismsa=msa(c(pikeseq,grayseq,salseq),method="Muscle",gapExtension=2,gapOpening=4,order="input")
	thisseqs=msaConvert(thismsa,type="seqinr::alignment")

	gappedcds=list()
	ungappedcds=list()
	for(j in 1:length(thisseqs$seq)) {
		backtrseqs=backtr(thisseqs$seq[j],as.character(allcds[[j]]))
		gappedcds[[length(gappedcds)+1]]=backtrseqs
		
	}


	for(j in nchar(gappedcds[[1]]):1) {
		for(k in 1:length(gappedcds)) {
			if(substr(gappedcds[[k]],j,j)=="-" ) {
				for(kk in 1:length(gappedcds)) { substr(gappedcds[[kk]],j,j)="N"; }
 				j=j-1
			}
		}
		
	}
	for(k in 1:length(gappedcds)) { gappedcds[[k]]=gsub("N","",gappedcds[[k]])}
	

	thisseqsasvec=lapply(thisseqs$seq,function(x){unlist(strsplit(x, split=""))})
	dnds[[length(dnds)+1]]=kaks(seqinr::as.alignment(3, thisseqs$nam, unlist(gappedcds)))
	
}


(c(dnds[[1]]$ka))


alldn=unlist(lapply(dnds,function(x){if(is.na(x)) { return(c(-1,-1,-1))}; return(c(x$ka))}))
alldn=matrix(alldn,ncol=3,nrow=length(alldn)/3,byrow=T)
alldn[alldn==-1]=NA
rownames(alldn)=gff.pike.gray$transcript
colnames(alldn)=c("pikegray","pikesal","graysal")
head(alldn)

allds=unlist(lapply(dnds,function(x){if(is.na(x)) { return(c(-1,-1,-1))}; return(c(x$ks))}))
allds=matrix(allds,ncol=3,nrow=length(allds)/3,byrow=T)
allds[allds==-1]=NA
rownames(allds)=gff.pike.gray$transcript
colnames(allds)=c("pikegray","pikesal","graysal")
head(allds)

alldnds=alldn/allds
head(alldnds)
par(mfrow=c(1,3),mar=c(11,5,1,1))
boxplot(alldn,ylim=c(0,0.8),ylab="dN",cex.axis=1.5,cex.lab=1.5,names=c("pike - grayling","pike - salmon","grayling - salmon"),las=2)
boxplot(allds,ylim=c(0,0.8),ylab="dS",cex.axis=1.5,cex.lab=1.5,names=c("pike - grayling","pike - salmon","grayling - salmon"),las=2)
boxplot(alldnds,ylim=c(0,0.8),ylab="dN/dS",cex.axis=1.5,cex.lab=1.5,names=c("pike - grayling","pike - salmon","grayling - salmon"),las=2)
alldnds[alldnds==Inf]=NA

summary(aov(c(alldn[,1],alldn[,2],alldn[,3])~rep(c("P-G","P-S","G-S"),each=nrow(alldn))))
summary(aov(c(allds[,1],allds[,2],allds[,3])~rep(c("P-G","P-S","G-S"),each=nrow(alldn))))
summary(aov(c(alldnds[,1],alldnds[,2],alldnds[,3])~rep(c("P-G","P-S","G-S"),each=nrow(alldn))))

TukeyHSD(aov(c(alldn[,1],alldn[,2],alldn[,3])~rep(c("P-G","P-S","G-S"),each=nrow(alldn))))
TukeyHSD(aov(c(allds[,1],allds[,2],allds[,3])~rep(c("P-G","P-S","G-S"),each=nrow(alldn))))
TukeyHSD(aov(c(alldnds[,1],alldnds[,2],alldnds[,3])~rep(c("P-G","P-S","G-S"),each=nrow(alldn))))
plot(allds[,1],alldn[,1])
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(allds[which(alldn[,1]<2 & allds[,1]<4),1],alldn[which(alldn[,1]<2 & allds[,1]<4),1],col="white",xlab="dS",ylab="dN")
lines(lowess(alldn[which(alldn[,1]<2 & allds[,1]<4),1]~allds[which(alldn[,1]<2 & allds[,1]<4),1],f=1/60),col="gray",lwd=2)
lines(lowess(alldn[which(alldn[,2]<2 & allds[,2]<4),2]~allds[which(alldn[,2]<2 & allds[,2]<4),2],f=1/60),col="salmon",lwd=2)
lines(lowess(alldn[which(alldn[,3]<2 & allds[,3]<4),3]~allds[which(alldn[,3]<2 & allds[,3]<4),3],f=1/60),col="black",lwd=2)
legend("topleft",legend=c("pike-grayling","pike-salmon","grayling-salmon"),fill=c("gray","salmon","black"))


head(gff.pike.gray)
write.table(data.frame(alldn,chr=gff.pike.gray$V1[match(rownames(alldn),gff.pike.gray$transcript)],pos=gff.pike.gray$V4[match(rownames(alldn),gff.pike.gray$transcript)]),file="dn.txt",quote=F,sep="\t")
write.table(data.frame(allds,chr=gff.pike.gray$V1[match(rownames(alldn),gff.pike.gray$transcript)],pos=gff.pike.gray$V4[match(rownames(allds),gff.pike.gray$transcript)]),file="ds.txt",quote=F,sep="\t")
write.table(data.frame(alldnds,chr=gff.pike.gray$V1[match(rownames(alldn),gff.pike.gray$transcript)],pos=gff.pike.gray$V4[match(rownames(alldnds),gff.pike.gray$transcript)]),file="dnds.txt",quote=F,sep="\t")


allds=read.table("ds.txt",sep="\t",header=T)
alldnds=read.table("dnds.txt",sep="\t",header=T)
#cut the scale down to max 2
allds[,1:3][allds[,1:3]>2]=2
alldnds[,1:3][alldnds[,1:3]>2]=2


##########################################################################################################################
#
#                             Read in mummer for gray-sal and gray-pike
#

#read in mummer for salmon and pike (mummer was run using commands: 
#/homeappl/home/tmsavi/appl_taito/MUMmer3.23/nucmer --mum -c 100 -g 200 -b 250 -p ${name} ${genome} ${querygenome}
#and a summary of all the alignments was created using
#/homeappl/home/tmsavi/appl_taito/MUMmer3.23/show-coords -r -c -l ${name}.delta > ${name}.coords

mummer=read.table("nucmerjelly_sal_mum.coords",skip=5)
mummer$ssa=salmonconvert$V1[match(mummer$V18,salmonconvert$V2)]
mummer$key=paste(mummer$V19,mummer$ssa)  #contig+ssa
mummer=mummer[order(mummer$V19,-mummer$V10),]
#chr-chr hit matches betw. bed and mummer, order by contig and select best score
bed$keysal=paste(bed$V1,bed$chrinsalmon)
mummer=mummer[mummer$key %in% bed$keysal & !is.na(mummer$ssa),]
mummer=mummer[order(mummer$V19,-mummer$V10),]
mummer=mummer[!duplicated(mummer$V19),] 

mummer.pike=read.table("nucmerjelly_pike_mum.pikesett.coords2",skip=5)
mummer.pike$pike=homeologs$V6[match(mummer.pike$V18,homeologs$V17)]
mummer.pike$key=paste(mummer.pike$V19,mummer.pike$pike)  #contig+ssa
mummer.pike=mummer.pike[order(mummer.pike$V19,-mummer.pike$V10),]
#chr-chr hit matches betw. bed and mummer, order by contig and select best score
bed$keypike=paste(bed$V1,bed$chrinpike)
mummer.pike=mummer.pike[mummer.pike$key %in% bed$keypike & !is.na(mummer.pike$pike),]
mummer.pike=mummer.pike[order(mummer.pike$V19,-mummer.pike$V10),]
mummer.pike=mummer.pike[!duplicated(mummer.pike$key),]


mummer.trout=read.table("nucmerjelly_trout_mum.coords",skip=5)
dim(mummer.trout) #should be 629691
mummer.trout=mummer.trout[grep("omy",mummer.trout$V18),] #remove unmapped scaffolds
mummer.trout$trout=homeologs$omy[match(mummer.trout$V18,homeologs$omy)]
mummer.trout$key=paste(mummer.trout$V19,mummer.trout$trout)  #contig+ssa
mummer.trout=mummer.trout[order(mummer.trout$V19,-mummer.trout$V10),]
#chr-chr hit matches betw. bed and mummer, order by contig and select best score
bed$keytrout=paste(bed$V1,bed$chrintrout)
mummer.trout=mummer.trout[mummer.trout$key %in% bed$keytrout & !is.na(mummer.trout$trout),]
mummer.trout=mummer.trout[order(mummer.trout$V19,-mummer.trout$V10),]
mummer.trout=mummer.trout[!duplicated(mummer.trout$key),]

#keys in the bed file for matching contig+chr
bed$chrinsalmon=homeologs$ssa[match(bed$chr,homeologs$V1)]
bed$chrinpike=homeologs$V6[match(bed$chr,homeologs$V1)]
bed$chrintrout=homeologs$omy[match(bed$chr,homeologs$V1)]
#attach coordinates for contigs
bed$mummersalmonpos=mummer$V1[match(bed$V1,mummer$V19)]
bed$mummerpikepos=mummer.pike$V1[match(bed$V1,mummer.pike$V19)]
bed$mummertroutpos=mummer.trout$V1[match(bed$V1,mummer.trout$V19)]

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#
#                             make circlize chromosomes and colors

###install.packages("randomcoloR")
require(randomcoloR)
chr.colors.gray <- rep("gray",nrow(cytogray))
chr.colors.gray2 <- sample(distinctColorPalette(nrow(cytogray)),nrow(cytogray),replace=F)

##########################################################################################################################
#
#			read in sdY loci

sdy=read.table("sdy_loci.txt",header=T,sep="\t",stringsAsFactors=F)
head(sdy)
sdy$color="blue"
sdy$color[sdy$PercMale<0.5]="red"
sdy=sdy[grep("Tth",sdy$LG),]

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#
#                             Initialize circular plot
#
require(circlize)

pdf("circos5.pdf",width=12,height=12)
par(mfrow=c(1,1),mar=c(1,1,1,1))
circos.clear()
circos.par(canvas.xlim=c(-1, 1),canvas.ylim=c(-1.2, 1.2),cell.padding = c(0, 0, 0, 0),gap.after = c(rep(1,50),10), start.degree = 0)#,start.degree=-85)
circos.initializeWithIdeogram(cytogray,plotType =NULL)  #c("axis","labels") would print chr names also
circos.trackPlotRegion(ylim = c(0, 0.35), 
	track.height = 0.05,
	panel.fun = function(x, y) {
		chr = get.cell.meta.data("sector.index")
		xlim = get.cell.meta.data("xlim")
			ylim = get.cell.meta.data("ylim")
		#circos.rect(xlim[1], 0, xlim[2], 0.5)
		circos.text(mean(xlim), 1.3, chr, cex = 0.8,font=2, facing = "clockwise",col="black", niceFacing = TRUE)
 	},
	bg.lwd =2,
      bg.border = NA,
	bg.col = "gray"
)
#TE peaks
for(i in 1:nrow(tc1marbinsMax)) {
		thischr=rownames(tc1marbinsMax)[i]
		chrpos=tc1marbinsMax[i,1]
		chrmiddle=tc1marbinsMax$Chrlen[i]/2
		thislab=tc1marbinsMax$chrType[i]
		circos.points(sector.index=thischr,x=chrpos,y=0.2,pch=21,cex=0.7,bg="purple") #,col=chr.colors.gray[chr.colors.gray==thischr]
		circos.text(sector.index=thischr,x=chrmiddle,y=0.55,labels=thislab,cex=0.8,col="gray") #,col=chr.colors.gray[chr.colors.gray==thischr]
}
for(i in 1:nrow(rtexbinsMax)) {
		thischr=rownames(rtexbinsMax)[i]
		chrpos=rtexbinsMax[i,1]
		circos.points(sector.index=thischr,x=chrpos,y=0.2,pch=21,cex=0.7,bg="yellow") #,col=chr.colors.gray[chr.colors.gray==thischr]
}

##########################################################################################################################
#
#                             add new track for different data sets
#

#repeat tracks (100 kilobase pair bins)
addTrack(cytogray,tc1marbins,valuecol="Len",thiscol="purple",thislowessf=0.2,prob=T,thistrackheight=0.07,thisminy=0.4,thismaxy=1)
addTrack(cytogray,rtexbins,valuecol="Len",thiscol=adjustcolor("yellow", alpha.f = 0.7),thislowessf=0.2,add=T,prob=T,thisminy=0.4,thismaxy=1)

#gene density (1 megabase pair bins)
addTrack(cytogray,annocounts,valuecol="genecount",chrcol="chr",poscol="pos",thiscol=adjustcolor("red", alpha.f = 0.7),thisminy=20,thismaxy=60,thislowessf=0.2,add=F,prob=F)

#dN/dS
addTrack(cytogray,alldnds[!is.na(alldnds$graysal),],valuecol="graysal",thiscol="gray",poscol="pos",chrcol="chr",thisminy=0.15,thismaxy=0.4)

##########################################################################################################################
#
#			add new track for linkage map
#

circos.track(cytogray$Chromosome,ylim=c(0,max(bed$femalepos,na.rm=T)),track.height = 0.1)
#For each LG, plot bed

for(i in 1:nrow(cytogray)) {
		thischr=cytogray$Chromosome[i]
		chrpos=bed$end[bed$chr==thischr & !is.na(bed$femalepos)]
		salmonizedchrpos=bed$end[bed$chr==thischr & !is.na(bed$salmonizedfemalepos)]

		mappos=bed$femalepos[bed$chr==thischr & !is.na(bed$femalepos)]
		salmonizedmappos=bed$salmonizedfemalepos[bed$chr==thischr & !is.na(bed$salmonizedfemalepos)]

		circos.points(sector.index=thischr,x=chrpos,y=mappos,pch=16,cex=0.4) #,col=chr.colors.gray[chr.colors.gray==thischr]
		#circos.points(sector.index=thischr,x=salmonizedchrpos,y=salmonizedmappos,pch=16,cex=0.3,col="white") 
		
}

#sdY
circos.track(cytogray$Chromosome,ylim=c(0,4),track.height = 0.08)
circos.polygon(sector.index="Tth_03.1",x=c(36989262,36989262),y=c(0,1),cex=1,border="blue",lwd=3) #,col=chr.colors.gray[chr.colors.gray==thischr]
circos.polygon(sector.index="Tth_11.1",x=c(3485,1151185),y=c(0,3),cex=1,border="red",lwd=3) #,col=chr.colors.gray[chr.colors.gray==thischr]
circos.polygon(sector.index="Tth_11.1",x=c(3485,1151185),y=c(3,4),cex=1,border="blue",lwd=3) #,col=chr.colors.gray[chr.colors.gray==thischr]
circos.polygon(sector.index="Tth_11.2",x=c(36207966,36207966),y=c(0,1),cex=1,border="red",lwd=3) #,col=chr.colors.gray[chr.colors.gray==thischr]
circos.polygon(sector.index="Tth_11.2",x=c(36207966,36207966),y=c(1,4),cex=1,border="blue",lwd=3) #,col=chr.colors.gray[chr.colors.gray==thischr]
circos.polygon(sector.index="Tth_18.1",x=c(14366696,14366696),y=c(0,1),cex=1,border="blue",lwd=3) #,col=chr.colors.gray[chr.colors.gray==thischr]
	
dev.off()

#######################################################################################
#
#			Predict homology between salmon and grayling for comparative figure 

#generate file names based on grayling-salmon chr. pairs (chr to chr -mummer alignments)
files=list("Tth_01","Tth_02","Tth_03","Tth_04","Tth_05","Tth_06","Tth_07","Tth_08","Tth_09","Tth_10","Tth_11","Tth_12","Tth_13","Tth_14","Tth_15","Tth_16","Tth_17","Tth_18","Tth_19","Tth_20","Tth_21","Tth_22","Tth_23","Tth_24","Tth_25")
#salmon homeologs for the #1 tth  chromosomes
nc1=list("NC_027319.1","NC_027325.1","NC_027313.1","NC_027308.1","NC_027318.1","NC_027300.1","NC_027312.1","NC_027322.1","NC_027301.1","NC_027326.1","NC_027305.1","NC_027312.1","NC_027323.1","NC_027300.1","NC_027308.1","NC_027320.1","NC_027311.1","NC_027314.1","NC_027309.1","NC_027304.1","NC_027328.1","NC_027316.1","NC_027306.1","NC_027306.1","NC_027303.1")
#salmon homeologs for the #2 tth  chromosomes
nc2=list("NC_027308.1",
"NC_027310.1","NC_027302.1","NC_027304.1","NC_027327.1","NC_027317.1","NC_027303.1","NC_027309.1","NC_027311.1","NC_027313.1","NC_027302.1","NC_027314.1","NC_027319.1",
"NC_027310.1","NC_027300.1","NC_027324.1","NC_027321.1","NC_027305.1","NC_027315.1","NC_027301.1","NC_027318.1","NC_027315.1","NC_027316.1","NC_027317.1","NC_027307.1")
predictions.homolog=0
pdf("predictions_salmonvsgrayling.pdf")
for(i in 1:length(files)) {
	f=files[[i]]
	n1=nc1[[i]]
	n2=nc2[[i]]
	print(f)
	dat2=read.table(paste("nucmer.",f,".1.",n1,".coords",sep=""),skip=5)
	dat3=read.table(paste("nucmer.",f,".2.",n2,".coords",sep=""),skip=5)
	lowess2=data.frame(lowess(dat2$V10~dat2$V4,f=0.14))
	lowess3=data.frame(lowess(dat3$V10~dat3$V4,f=0.14))

	thischr1=paste(f,1,sep=".")
	thischr2=paste(f,2,sep=".")
	lowess2$chr=thischr1
	lowess3$chr=thischr2
	if(length(predictions.homolog)==1) {predictions.homolog=lowess2} else { predictions.homolog=rbind(predictions.homolog,lowess2) }
	if(length(predictions.homolog)==1) {predictions.homolog=lowess3} else { predictions.homolog=rbind(predictions.homolog,lowess3) }

	plot(lowess2$x,lowess2$y,col="red",lwd=2,type="l",main=f,xlab="grayling pos",ylab="predicted identity",ylim=range(c(lowess2$y,lowess3$y)))
	points(lowess3$x,lowess3$y,col="blue",type="l",lwd=2)
	legend("topleft",legend=c(paste(f,".1 vs. ",n1,sep=""),paste(f,".2 vs. ",n2,sep="")),fill=c("red","blue"))
}
dev.off()
#also for 13.3
dat2=read.table(paste("nucmer.Tth_13.3.NC_027323.1.coords",sep=""),skip=5)
lowess2=data.frame(lowess(dat2$V10~dat2$V4,f=0.14))
thischr2=paste("Tth_13.3")
lowess2$chr=thischr2
predictions.homolog=rbind(predictions.homolog,lowess2)
dim(predictions.homolog)
predictions.homolog=predictions.homolog[!duplicated(paste(predictions.homolog$x,predictions.homolog$y,predictions.homolog$chr)),]
predictions.homolog=predictions.homolog[order(predictions.homolog$chr,predictions.homolog$x),]
write.table(file="lowess_homolog_identities2.txt",predictions.homolog,sep="\t",quote=F,col.names=T,row.names=F)
homologidentities=read.table("lowess_homolog_identities2.txt",header=T,stringsAsFactors=F)

#######################################################################################
#
#			Make the comparative figure 

roundedidentities=table(round(homologidentities$y,1))
colorscalefunc=colorRampPalette(c("darkblue","lightblue","red","darkred","darkred"))
colorscale2=colorscalefunc(length(roundedidentities))
barplot(roundedidentities,col=colorscale2)
pdf("comparative9.pdf",width=5,height=10)
par(mfrow=c(6,5),mar=c(2,2,1,0.5))
for(i in 1:length(cytosal$Chromosome)) {
	thischr1=cytosal$Chromosome[i]
	plot(c(1,1.5),c(0,max(cytosal$ChromEnd)),col="white",xaxt="n",xlab="",ylab="bp",main=thischr1)
	grayx=1
	salx=1.5
	arrows(grayx,1,grayx,cytosal$ChromEnd[match(thischr1,cytosal$Chromosome)],code=3,angle=90,length=0.05)
	arrows(salx,1,salx,cytosal$ChromEnd[match(thischr1,cytosal$Chromosome)],code=3,angle=90,length=0.05)
	graychrs=bed$chr[which(bed$chrinsalmon==thischr1)]
	graypositions=bed$V8[which(bed$chrinsalmon==thischr1)]
	salmonchrs=bed$mummersalmonchr[which(bed$chrinsalmon==thischr1)]
	salmonpositions=bed$mummersalmonpos[which(bed$chrinsalmon==thischr1)]
	compositepositions=compositecoordinates(graychrs,graypositions,salmonpositions)
	compositepositions=adjustcoordinatesbyoverallorient(graychrs,compositepositions,salmonpositions)
	cols.opaque=chr.colors.gray2[match(graychrs,cytogray$Chromosome)]
	cols=adjustcolor(cols.opaque, alpha.f = 0.3)
	cols[match("Tth_13.2",cytogray$Chromosome)]=cols.opaque[match("Tth_13.2",cytogray$Chromosome)]
	
	#add grayling-salmon links
	for(j in 1:length(compositepositions)) {
		points(c(grayx,salx),c(compositepositions[j],salmonpositions[j]),col=cols[j],pch=16,type="l")
	}
	legend("topleft",legend=graychrs[!duplicated(graychrs)],fill=cols.opaque[!duplicated(graychrs)],cex=0.8,bty="n")
	
	#add homeolog identy for grayling
	tmphomologidentities=homologidentities[homologidentities$chr %in% graychrs[!duplicated(graychrs)],]
	tmphomologidentities=tmphomologidentities[order(-tmphomologidentities$y),]
	for(j in 1:nrow(tmphomologidentities)) {
		saladustedpos=closestfeat(tmphomologidentities$chr[j],tmphomologidentities$x[j],graychrs,graypositions,graychrs,salmonpositions,func="closest")
		
		roundedidentity=round(tmphomologidentities$y[j],1)
		tmpcolor=colorscale2[names(roundedidentities)==roundedidentity]
		tmpcolor=adjustcolor(tmpcolor, alpha.f = 0.3)
		#points(grayx,homeologadustedpos$pos[1],col=tmpcolor,pch="-",cex=1.4)
		points(salx,saladustedpos$pos[1],col=tmpcolor,pch="-",cex=1.4)
	}	
	#add centromere positions
	for(thisgraychr in graychrs[!duplicated(graychrs)]) {
		#grayling
		thispos=unlist(tc1marbinsMax[rownames(tc1marbinsMax)==thisgraychr,1])
		thispos2=unlist(rtexbinsMax[rownames(rtexbinsMax)==thisgraychr,1])
		if(length(thispos)>0) {
			ct1marapprox=closestfeat(thisgraychr,thispos,graychrs,graypositions,graychrs,compositepositions,func="closest")
			print(paste(thisgraychr,"ct1-mariner:",thispos,"->",ct1marapprox$pos))
			points(grayx,ct1marapprox$pos,pch=21,bg="purple",cex=1)
		}
		if(length(thispos)>0) {
			rtexapprox=closestfeat(thisgraychr,thispos2,graychrs,graypositions,graychrs,compositepositions,func="closest")
			print(paste(thisgraychr,"rte-x:",thispos2,"->",rtexapprox$pos))
			points(grayx,rtexapprox$pos,pch=21,bg="yellow",cex=1)
		}
	}
}
dev.off()
#################################################################################################################
#
#                  Inverted / non-inverted blocks (supplementary figure)
#
require(graphics) 
par( mfrow=c(2,1),mar=c(4,4,1,1))

#Manually add blocks of straight / inverted syntenic block. 
#click to add inversion start and end points
coords=list()
for(c in bed$chr[!duplicated(bed$chr)]) {
	thisxlim=c(0,max(as.numeric(bed$mummerpikepos[bed$chr==c]),na.rm=T)*1.25)
	plot(as.numeric(bed$mummerpikepos[bed$chr==c]),bed$V8[bed$chr==c],ylab="grayling pos",xlab="pike pos",main=c,xlim=thisxlim,col="green4")

	thisxlim=c(0,max(as.numeric(bed$mummersalmonpos[bed$chr==c]),na.rm=T)*1.25)
	plot(bed$V8[which(bed$chr==c)],as.numeric(bed$mummersalmonpos[which(bed$chr==c)]),ylab="grayling pos",xlab="salmon pos",col="salmon")
	p=identify(bed$V8[which(bed$chr==c)],as.numeric(bed$mummersalmonpos[bed$chr==c]),labels=bed$mummersalmonpos[bed$chr==c],tolerance=0.5,cex=0.7,col="black",pos=T)
	for(pos in bed$V8[which(bed$chr==c)][p$ind]) {coords[[length(coords)+1]]=pos; names(coords)[[length(coords)]]=c}
}
coords=data.frame(chr=names(coords),graylingpos=unlist(coords),c("start","end"))
#write.table(coords,file="salmon_final_blocks_using_interactive_plost.txt",sep="\t",quote=F,row.names=F)
rm(coords)
coords2=read.table("salmon_final_blocks_using_interactive_plost.txt",header=T,stringsAsFactors=F)
coords2$type=c("start","end")
coords2$blocknr=rep(1:sum(coords2$type=="start"),each=2)  #order in grayling
#write.table(coords2,file="salmon_final_blocks_using_interactive_plost_withblockid.txt",sep="\t",quote=F,row.names=F)
#add manually columns "inverted" (char) and "rearranged" (boolean)
coords2=read.table("salmon_final_blocks_using_interactive_plost_withblockid_edited.txt",header=T,stringsAsFactors=F,sep="\t")
#summary table
chrsummary=coords2[!duplicated(coords2$blocknr),c(1,4)]
chrsummary$chrlength=cytogray$ChromEnd[match(chrsummary$chr,cytogray$Chromosome)]
chrsummary$dsoutliercount=0
chrsummary$dndsoutliercount=0
chrsummary$totgenecount=0
chrsummary$meands=0
chrsummary$meandnds=0
chrsummary$meandsPS=0
chrsummary$meandndsPS=0
chrsummary$meandsGS=0
chrsummary$meandndsGS=0



blockcolors=rainbow(10)
#pdf("S:/2017_data_not_in_sync/NORWAY_MAY2017/Sep_chromosome_structure/salmon_final_blocks_using_interactive_plost_v5.pdf",width=11)
pdf("salmon_final_blocks_using_interactive_plost_v8.pdf",width=11)
for(c in cytogray$Chromosome) {
	par( mfrow=c(3,1),mar=c(4,4,1,4))
	thisxlim=c(0,max(as.numeric(bed$V8[bed$chr==c]),na.rm=T)*1.25)
	plot(bed$V8[bed$chr==c],bed$mummerpikepos[bed$chr==c],xlab="grayling pos",ylab="pike pos",main=c,xlim=thisxlim,col="green4")
	thisblocknr=1	
	for(p in seq(1,length(which(coords2$chr==c)),by=2)) {
		startpoint=coords2[coords2$chr==c,][p,"pos"] #pos=grayling pos.
		endpoint=coords2[coords2$chr==c,][p+1,"pos"] #pos=grayling pos.
		points(bed$V8[bed$chr==c& bed$V8>=startpoint & bed$V8<=endpoint],bed$mummerpikepos[bed$chr==c & bed$V8>=startpoint & bed$V8<=endpoint],col=blockcolors[thisblocknr],pch=16)
		thisblocknr=thisblocknr+1
	}
	legend("topright",legend=1:(thisblocknr-1),fill=blockcolors[1:(thisblocknr-1)] )
	thisxlim=c(0,max(as.numeric(bed$V8[bed$chr==c]),na.rm=T)*1.25)
	thisylim=range(bed$mummersalmonpos[which(bed$chr==c)],na.rm=T)
	plot(bed$V8[which(bed$chr==c)],bed$mummersalmonpos[which(bed$chr==c)],xlab="grayling pos",ylab="salmon pos",xlim=thisxlim,col="salmon",ylim=thisylim)
	thisblocknr=1
	for(p in seq(1,length(which(coords2$chr==c)),by=2)) {
		startpoint=coords2[coords2$chr==c,][p,"pos"] #pos=grayling pos.
		endpoint=coords2[coords2$chr==c,][p+1,"pos"] #pos=grayling pos.
		points(bed$V8[bed$chr==c& bed$V8>=startpoint & bed$V8<=endpoint],bed$mummersalmonpos[bed$chr==c & bed$V8>=startpoint & bed$V8<=endpoint],col=blockcolors[thisblocknr],pch=16)
		text(
			mean(bed$V8[bed$chr==c& bed$V8>=startpoint & bed$V8<=endpoint],na.rm=T),
			thisylim[2]*0.95,
			paste("blockid=",coords2$blocknr[coords2$chr==c][p],"\n",coords2$inverted[coords2$chr==c][p],"-",coords2$rearranged[coords2$chr==c][p],sep=""),
			cex=0.8)
		thisblocknr=thisblocknr+1
		numpoint=sum(bed$chr==c& bed$V8>=startpoint & bed$V8<=endpoint)
	}
	abline(v=tc1marbinsMax[which(rownames(tc1marbinsMax)==c),1],col="purple")
	#trout
	thisxlim=c(0,max(as.numeric(bed$V8[bed$chr==c]),na.rm=T)*1.25)
	thisylim=range(bed$mummertroutpos[which(bed$chr==c)],na.rm=T)
	if(sum(!is.na(bed$mummertroutpos[which(bed$chr==c)]))==0) { thisylim=c(0,0) }
		plot(bed$V8[which(bed$chr==c)],bed$mummertroutpos[which(bed$chr==c)],xlab="grayling pos",ylab="trout pos",xlim=thisxlim,col="red",ylim=thisylim)
		thisblocknr=1
		for(p in seq(1,length(which(coords2$chr==c)),by=2)) {
			startpoint=coords2[coords2$chr==c,][p,"pos"] #pos=grayling pos.
			endpoint=coords2[coords2$chr==c,][p+1,"pos"] #pos=grayling pos.
			points(bed$V8[bed$chr==c& bed$V8>=startpoint & bed$V8<=endpoint],bed$mummertroutpos[bed$chr==c & bed$V8>=startpoint & bed$V8<=endpoint],col=blockcolors[thisblocknr],pch=16)
			#text(mean(bed$V8[bed$chr==c& bed$V8>=startpoint & bed$V8<=endpoint],na.rm=T),thisylim[2],paste("blockid",coords2$blocknr[coords2$chr==c][p],sep="="),cex=0.8)
			thisblocknr=thisblocknr+1
			numpoint=sum(bed$chr==c& bed$V8>=startpoint & bed$V8<=endpoint)
		}
		abline(v=tc1marbinsMax[which(rownames(tc1marbinsMax)==c),1],col="purple")
	
}
dev.off()

#######################################################################################
#
#			Linear model inspecting if long:short arm ratio and chr length 
#			are related to linkage map size

maplens=data.frame(tapply(map$V5,map$Chr,max))
colnames(maplens)="femalepos"
maplens$chrlen=cytogray$ChromEnd[match(rownames(maplens),cytogray$Chromosome)]
maplens$longShortRatio=tc1marbinsMax$longShortRatio[match(rownames(maplens),rownames(tc1marbinsMax))]
maxposlm0=lm(maplens$femalepos~maplens$chrlen)
maxposlm1=lm(maplens$femalepos~maplens$longShortRatio+maplens$chrlen)
summary(maxposlm1)
anova(maxposlm0,maxposlm1,test="Chisq")

#######################################################################################
#
#			Draw linkage maps

par(mfrow=c(7,8),mar=c(0.1,2,1,0.1),bty="n")
pdf("Linkagemaps.pdf",width=14,height=32)
for(c in levels(map$Chr)) {
	plot(c(1,1),range(map$V5[map$Chr==c]),type="l",ylim=c(min(map$V5),max(map$V5)+10),col="gray",xaxt="n",cex.main=0.95,xlim=c(0.5,4),cex.axis=0.8)
	text(1,max(map$V5[map$Chr==c])+10,paste(c,"\nn=",sum(map$Chr==c),sep=""),cex=1)

	#bin points	
	pointnames=map$Marker[map$Chr==c]
	points=map$V5[map$Chr==c]
	pointbins=cut(points,seq(0,max(points)+binmarkersby,binmarkersby),include.lowest=T)
	pointbinstarts=seq(0,max(points)+binmarkersby,binmarkersby)[as.numeric(pointbins)]
	pointbinmedians=pointbinstarts+0.5*binmarkersby
	
	#combine text based on point bins
	pointlabels=rep("x",length(points))
	pointbinpositions=rep(0,length(points))

	for(i in 1:length(pointlabels)) { 
		theselabs=pointnames[pointbins==pointbins[i]]
		thesepoints=points[pointbins==pointbins[i]]
		if(length(theselabs)<=3) { pointlabels[i]=paste(theselabs,collapse=",")}
		if(length(theselabs)>3) { pointlabels[i]=paste(theselabs[1],"...",theselabs[length(theselabs)]," (",length(theselabs)-2," more)",sep="")}
		pointbinpositions[i]=mean(thesepoints)
	}
	#remove duplicate labs
	pointlabels_nondup=pointlabels[!duplicated(pointlabels)]
	pointbinpositions_nondup=pointbinpositions[!duplicated(pointlabels)]

	for(point in points) { points(c(0.85,1.15),c(point,point),col="blue",type="l")	}
	for(i in 1:length(pointlabels_nondup)) { text(1.1,pointbinpositions_nondup[i], pointlabels_nondup[i],pos=4,cex=0.7)	}

}
dev.off()

##########################################################################################################
##########################################################################################################	
##########################################################################################################
#
#			Functions
#

##########################################################################################################
#
#                 Function closestfeat takes two sets of aligned chromosomal coordinatas to predict
#			the corresponding position in set 2 based on query position in set 1  


#This function takes a query chromosome name and position as parameters, followed by lists of chr1 names, chr1 positions, 
#chr2 names and chr2 positions. The  The chr1 and chr2 coordinates are in the matching order. The query chromosomal
#coordinates belong to the chr1 set. The function predicts the query position location along the chr2 coordinates using
#function lm (linear regression), gam (local regression) or finding closest matching feature among chr1 coordinates. 
#The prediction method can be set using parameter func with values "lm", "gam" or "closest". 
#
#A prediction of the query position in chromosome 2 position list is returned as a list with chromosome name and 
#position as list items.  
closestfeat=function(chr1,pos1,chr1list,pos1list,chr2list,pos2list,func="lm") {
	chr2list_flt=chr2list[which(chr1list==chr1 & !is.na(chr2list))]
	pos2list_flt=pos2list[which(chr1list==chr1 & !is.na(chr2list))]
	pos1list_flt=pos1list[which(chr1list==chr1 & !is.na(chr2list))]
	chr1list_flt=chr1list[which(chr1list==chr1 & !is.na(chr2list))]
	
	pos1list_dev=abs(pos1list_flt-pos1)
	pos1list_minpos=which(pos1list_dev==min(pos1list_dev,na.rm=T))

	if(length(chr2list_flt)==0) { return(-1) }	
	if(length(chr2list_flt)==1) { return(list(chr=unlist(chr2list_flt),pos=unlist(pos2list_flt))) }
	if(length(chr2list_flt)>1) { 
		if(length(chr2list_flt)<4) { func="closest" }
		if(func=="lm") {posmodel=lm(pos2list_flt~pos1list_flt); predpos=predict(posmodel,data.frame(pos1list_flt=pos1)) }
		if(func=="gam") {posmodel=gam(pos2list_flt~s(pos1list_flt)); predpos=predict(posmodel,data.frame(pos1list_flt=pos1)) }
		if(func=="closest") { predpos=pos2list_flt[pos1list_minpos] }
		
		return(list(chr=chr2list_flt[1],pos=unlist(predpos)))
	}
}
x=closestfeat("chr1",200,c("chr0","chr1","chr1","chr1","chr2","chr3"),c(100,100,200,300,100,100),c("CHR0","CHR1","CHR1","CHR1","CHR2","CHR3"),1:6,func="lm")

##########################################################################################################################
#
#                             addTrack function to add new track
#

#valuedf contains:
# - factor Chr for chromosome with identical levels to Chromosome in cytodataframe
# - numerical Pos for positions in chromosomes
# - numerical Val to plot
#cytodataframe has to contain field
# - Chromosome, a factor with identical (and identically & alphabetically ordered) levels as valuedf
#valuecol is the column name with values to plot. Default: "Val"
addTrack=function(cytodataframe,valuedf,add=F,thistype="l",colorarea=T,prob=F,thiscol="gray",thisminy=NULL,thismaxy=NULL,valuecol="Val",thistrackheight=0.06,thislowessf=0.4,poscol="Pos",chrcol="Chr") {
	
	lowess_vals=list()
	lowess_vals_min=list()
	lowess_vals_max=list()
	lgs_with_data=list()

	for(i in 1:nrow(cytodataframe)) {
		thischr=as.character(cytodataframe$Chromosome[i])
		if(sum(valuedf[,chrcol]==thischr)>1) {
			lowess_vals[[length(lowess_vals)+1]]=lowess(c(valuedf[valuedf[,chrcol]==thischr,poscol]),c(valuedf[valuedf[,chrcol]==thischr,valuecol]),f=thislowessf)
			thisrange=range(lowess_vals[[length(lowess_vals)]]$y)
			lowess_vals_min[[length(lowess_vals_min)+1]]=thisrange[1]
			lowess_vals_max[[length(lowess_vals_max)+1]]=thisrange[2]
			lgs_with_data[[length(lgs_with_data)+1]]=thischr

		}
	}

	
	lowess_vals_min=min(unlist(lowess_vals_min))
	lowess_vals_max=max(unlist(lowess_vals_max))

	if(length(thisminy)>0) { lowess_vals_min=thisminy; }
	if(length(thismaxy)>0) { lowess_vals_max=thismaxy; }

	
	if(prob==T) {
		lowess_vals_min=0; lowess_vals_max=1;
		if(length(thisminy)>0) { lowess_vals_min=thisminy*lowess_vals_max; }
		if(length(thismaxy)>0) {  lowess_vals_max=thismaxy*lowess_vals_max;}
	}


	print(paste("min:",lowess_vals_min,"max:",lowess_vals_max))
	#Add track
	if(add==F) {circos.track(cytodataframe$Chromosome,ylim=c(lowess_vals_min,lowess_vals_max),track.height = thistrackheight) }


	for(i in 1:length(lgs_with_data)) {
		thischr=lgs_with_data[[i]]
		if(length(lowess_vals)>=i) {
			this_lowess=lowess_vals[[i]]
			thisx=this_lowess$x
			thisy=this_lowess$y
			if(prob==T) { thisy=thisy/max(thisy) }
			thisy[thisy<lowess_vals_min]=lowess_vals_min
			thisy[thisy>lowess_vals_max]=lowess_vals_max
		
			circos.lines(sector.index=thischr, type=thistype, area=colorarea, x=thisx, y=thisy, lwd=1, col=thiscol)
		}
	}
}#end function


#each parameter vector is equally long.
#names contains query chromosome names to be adjusted 
#originalpositions contains original positions along query chr
#referencepositions contains positions in a reference chr
#returns adjusted positions so query chromsoomes can be aligned
#against a single reference chromsome.
compositecoordinates=function(names, originalpositions, referencepositions) {
	medianpositions=c()   #median positions in the reference (for ordering)
	querynames=names[!duplicated(names)]
	for(c in querynames) {
		medianpositions=c(medianpositions,median(referencepositions[names==c],na.rm=T))
	}
	print(paste(querynames,"median:",medianpositions))
	offset=0
	compositepositions=list()
	for(c in querynames[order(medianpositions)]) {

		compositepositions[[length(compositepositions)+1]]=originalpositions[names==c]+offset
		names(compositepositions)[[length(compositepositions)]]=c
		print(paste(c,"median composite:",median(originalpositions[names==c]+offset,na.rm=T),"composite range:",paste(range(originalpositions[names==c]+offset),collapse="-")))
		
		offset=offset+max(originalpositions[names==c],na.rm=T)
		
	}
	resultpositions=originalpositions
	for(c in querynames[order(medianpositions)]) {
		resultpositions[names==c]=compositepositions[[c]]
	}
	return(resultpositions)
}
compositecoordinates(names=rep(paste("chr",1:3,sep=""),each=3),c(1:3,1:3,1:3),c(7,9,8,6:4,1:2,NA))

#if the whole originalpositions are in inverted order to the reference positions, flip all original positions
adjustcoordinatesbyoverallorient=function(names, originalpositions,referencepositions) {
	querynames=names[!duplicated(names)]
	overalldirection=list()
	for(c in querynames) {
		overalldirection[[length(overalldirection)+1]]=summary(lm(originalpositions[names==c]~referencepositions[names==c],na.action="na.omit"))$coeff[2,"Estimate"]
		names(overalldirection)[[length(overalldirection)]]=c
	}
	resultpositions=originalpositions
	for(c in querynames) {
		if(overalldirection[[c]]<0) {
			resultpositions[names==c]=rev(resultpositions[names==c])
		}
	}
	return(resultpositions)
}
adjustcoordinatesbyoverallorient(names=rep(paste("chr",1:3,sep=""),each=3),comp,c(7,9,8,6:4,1:2,NA))


#####################################################################
#
#  function for backtranslate

backtr=function(pepgappedseq,cdsseq) {
	pepgapped=unlist(strsplit(pepgappedseq, split=""))
	cds=unlist(strsplit(cdsseq, split=""))
	cds=c(cds,c("N","N","N"))
	cdsgapped=rep("-",(nchar(pepgappedseq)*3))
	cdsi=1
	cdsgappedi=1
	for(ii in 1:length(pepgapped)) {
		if(pepgapped[ii]!="-") { cdsgapped[cdsgappedi]=cds[cdsi]; cdsgapped[cdsgappedi+1]=cds[cdsi+1]; cdsgapped[cdsgappedi+2]=cds[cdsi+2]; cdsi=cdsi+3; }
		cdsgappedi=cdsgappedi+3;
	}
	return(paste(cdsgapped,collapse=""))
}
backtr("-A--BC--ABC-","aaabbbcccaaabbbccc")


