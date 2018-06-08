
############################################################
#
#       Configure
#

setwd("C:/Users/tmsavi/Desktop/NORWAY_tmp")
salmon_file="RADs_align_salmonizedTth_100.txt"
salmonizedorderfolder="salmonizedorders"

############################################################
#
#       Read maps with attached salmonized and troutized assembly data and check
#

#Assembly must have column names refstart, refend, querystart, queryend, femalepos, malepos, markeroffset, LG, refname
#####assembly_file=salmon_file
readAssembly=function(assembly_file){
	assembly=read.table(assembly_file,sep="\t",stringsAsFactors=F,header=T)
	assembly$refIsForward=unlist(apply(assembly,1,function(x) {if(x[colnames(assembly)=="refstart"]<x[colnames(assembly)=="refend"]) { return(T)} else return(F)}))
	assembly$queryIsForward=unlist(apply(assembly,1,function(x) {if(x[colnames(assembly)=="querystart"]<x[colnames(assembly)=="queryend"]) { return(T)} else return(F)}))
	
	assembly$femalepos=as.numeric(assembly$femalepos)
	assembly$malepos=as.numeric(assembly$malepos)
	assembly=assembly[!is.na(assembly$femalepos) & !is.na(assembly$malepos),]
	if(sum(assembly$refIsForward==F)>0) { return(NULL);}  #all ref coordinates must be forward
	
	#calculate snp position 

	assembly$snppos=-1
	#If query is in forw. strand, then snp position is ref start - query start + snp offset.
	#If query is in rev strand, then snp position is ref end + query start- snp offset.
	sel.refForQueFor=assembly$refIsForward==1&assembly$queryIsForward==1  #make selection vectors for forward/reverse ref. query combinations
	assembly$snppos[sel.refForQueFor]=assembly$refstart[sel.refForQueFor]-assembly$querystart[sel.refForQueFor]+assembly$markeroffset[sel.refForQueFor]
	assembly$snppos[!sel.refForQueFor]=assembly$refend[!sel.refForQueFor]+assembly$querystart[!sel.refForQueFor]-assembly$markeroffset[!sel.refForQueFor]

	#assembly is ordered first by LG, then female map, then male map.
	assembly=assembly[order(assembly$LG,assembly$femalepos,assembly$malepos),]
	return(assembly)
}


salmonassembly=readAssembly(salmon_file)
dim(salmonassembly)

############################################################
#
#       Salmonize the map: order first by LG, then by reference ID.
#       Only take the markers that match to the best hit chromosome of
#       that LG.
#

salmonizedAssembly=function(assembly) {
	assembly=assembly[order(assembly$LG,assembly$refstart),] #order by LG then salmon
	lgs=assembly$LG[!duplicated(assembly$LG)]

	#select only best hit chromosomes to salmonize
	bestrefsel=rep(0,nrow(assembly))
	for(i in lgs) { 
		bestref=names(sort(table(assembly$refname[assembly$LG==i]),decreasing=T))[1]
		bestrefsel[assembly$LG==i & assembly$refname==bestref]=1
	}
	assembly=assembly[bestrefsel==1,]
	return(assembly)
}


salmonizedassembly=salmonizedAssembly(salmonassembly)

plot(salmonizedassembly$femalepos[salmonizedassembly$LG=="1"],salmonizedassembly$snppos[salmonizedassembly$LG=="1"])

############################################################
#
#        Plot map positions against ref. positions
#



#Assembly is salmonassembly or troutassembly, lglist is list of LG names to plot
#assembly=salmonizedassembly
#lg="1"
#minId=100
plotAssemblyLGs=function(assembly,lg,malepos=F,minId=100){
	assembly=assembly[assembly$LG==lg,]   #subset to include only this LG
	lgpositions=assembly$femalepos
	if(malepos==T) { lgpositions=assembly$malepos }



	xlabel="Female map pos."
	if(malepos==T) { mapposlabel="Male map pos." }

	#Find out chromosome colors by frequency
	chromosomes=sort(table(assembly$refname),decreasing=T)
	chromosomecolors=1:length(chromosomes)
	assembly$color=chromosomecolors[match(assembly$refname,names(chromosomes))]

	#add space for legend
	ywidth=max(assembly$snppos)-min(assembly$snppos)

	#highlight low quality 
	assembly$pch=16
	assembly$pch[assembly$Id<=minId]=8 

	#plot
	plot(assembly$femalepos,assembly$snppos,xlab=xlabel,ylab="Marker pos. in reference",col=assembly$color,pch=assembly$pch,
		main=paste("LG",lg),ylim=c(min(assembly$snppos)-ywidth*0.05,max(assembly$snppos)+ywidth*0.10),cex=0.9,cex.axis=0.8)
	legend("topleft",horiz=TRUE, bty='n',legend=names(chromosomes),fill=chromosomecolors,cex=0.5)
	
	###loess
	#tra=data.frame(map=assembly$femalepos,pos=(assembly$snppos),chr=assembly$refname)   #if include only pops
	#tra=tra[which(!duplicated(tra$map) & tra$chr==names(sort(table(tra$chr),decreasing=T))[1]),]

	#loess50 <- loess(pos~ map, data=tra, span=0.40,normalize=F,family="symmetric")
	
	#smoothed50 <- predict(loess50)
	#lines(smoothed50, x=tra$map, col="blue")
	#slopecolor=c(0,as.numeric(diff(loess50$fitted)>0))+1
	#slopecolorx=seq(min(tra$map),max(tra$map),by=(max(tra$map)-min(tra$map))/(loess50$n-1))

	#points(slopecolorx,rep(min(tra$pos),length(slopecolor)),pch="|",col=slopecolor,cex=1.2)
	
}

plotAssemblyLGs(salmonassembly,"1")
plotAssemblyLGs(salmonizedassembly,"1")

############################################################
#
#       Loop through chromosomes and plot salmon assembly
#


#These are LGnames to plot
LGnames=as.character(1:51)
LGnums=1:length(LGnames)

#for each LG, plot salmonized LG assembly 

pdf("Salmonized_assembly.pdf",height=5,width=10)
for(i in 1:length(LGnames)) {
	par(mfrow=c(1,1))
	if(sum(salmonassembly$LG==LGnames[i])>1)  { 
		plotAssemblyLGs(salmonassembly,LGnames[i]) 
		mtext(side=1,cex=0.5,adj=1,"salmonized",line=3)
	}
}
dev.off()
plotAssemblyLGs(salmonassembly,"1")

############################################################
#
#       Loop through chromosomes and print 100% salmonized orders to file (for LepMAP)
#


for(i in 1:length(LGnames)) {
	bestref=names(sort(table(salmonizedassembly$refname[salmonizedassembly$LG==LGnames[i]]),decreasing=T))[1]
	x=data.frame(marker=salmonizedassembly[salmonizedassembly$LG==LGnames[i] & salmonizedassembly$refname==bestref,"markerid"])
	x$marker=sub("M","",x$marker)
	write.table(x,quote=F,row.names=F,col.names=F,file=paste(salmonizedorderfolder,"/ssa.",LGnames[i],".order.txt",sep=""))
}

############################################################
#
#       RUN LEPMAP HERE for 100% salmonized chromosomes and parsesep the files together
#

#Run Lep-MAP2 using shellscript norway_may2017_lepmap_salmonized.sh (give LG as param1)
#module load java/oracle/1.8
#lg=$1
#evaluate one LG (salmonized order)
#java -cp $USERAPPL/lepmap2/bin/ OrderMarkers improveOrder=0 numThreads=5 evaluateOrder=salmonizedorders/ssa.${lg}.order.txt removeDuplicates=0 polishWindow=30 filterWindow=10 useKosambi=1 minError=0.15 map=/wrk/tmsavi/DONOTREMOVE/RAD/lepmap/lepmap_lep/map_js.lep.lod9.txt data=/wrk/tmsavi/DONOTREMOVE/RAD/lepmap/lepmap_lep/data.linkage.lep>ssa.order.lod9.LG${lg}.txt

#Then combine the results using:
#perl $HOME/parse_lepmap_map_separately.pl ssa.order.lod9.LG1.txt >ssa.order.lod9.parsesep.txt

#(Repeat for troutized...)
#java -cp $USERAPPL/lepmap2/bin/ OrderMarkers improveOrder=0 numThreads=5 evaluateOrder=salmonizedorders/omy.${lg}.order.txt removeDuplicates=0 polishWindow=30 filterWindow=10 useKosambi=1 minError=0.15 map=/wrk/tmsavi/DONOTREMOVE/RAD/lepmap/lepmap_lep/map_js.lep.lod9.txt data=/wrk/tmsavi/DONOTREMOVE/RAD/lepmap/lepmap_lep/data.linkage.lep>omy.order.lod9.LG${lg}.txt


############################################################
#
#       Read Lep-MAP results to get 100% salmonized mapping distances
#
salmonized.remaps=read.table("ssa.order.lod9.parsesep.txt",sep="")
salmonized.remaps$V1=sub("ssa.order.lod9.LG","",salmonized.remaps$V1)
salmonized.remaps$V1=sub(".txt","",salmonized.remaps$V1)
salmonized.remaps$V3=paste("M",salmonized.remaps$V3,sep="")

#if there are LG splits in salmonized map, then add large gap 100cM
addgap=0
for(i in 2:nrow(salmonized.remaps)) {
	if(addgap>0 & salmonized.remaps$V1[i]==salmonized.remaps$V1[i-1]) {salmonized.remaps$V5[i]=salmonized.remaps$V5[i]+addgap}
	if(salmonized.remaps$V1[i]!=salmonized.remaps$V1[i-1]) { addgap=0}
	if(salmonized.remaps$V1[i]==salmonized.remaps$V1[i-1] & salmonized.remaps$V5[i]<salmonized.remaps$V5[i-1]) {
		print(paste("WARNING! SPLIT LG: ",salmonized.remaps$V1[i],"! 100 cM gap added!",sep=""))
		addgap=salmonized.remaps$V5[i-1]+100
		salmonized.remaps$V5[i]=addgap
	}
}

############################################################
#
#       Summarize new map distances of 100% salmonized data
#

salmon1=aggregate(femalepos~LG, salmonassembly[as.numeric(salmonassembly$LG)>=1&as.numeric(salmonassembly$LG)<=51,],max)
salmon2=aggregate(V5~V1, salmonized.remaps,max)
salmon_femalelen=salmon1[order(as.numeric(salmon1$LG)),]
colnames(salmon_femalelen)=c("LG","salmon1")
salmon_femalelen$salmon2=salmon2$V5[match(salmon_femalelen$LG,salmon2$V1)]

par(mfrow=c(1,1))
barplot(t(salmon_femalelen[,c(2:3)]),beside=T,cex.names=0.6,names=salmon_femalelen$LG)
legend("topright",title="Salmonized",legend=c("original","remapped"),fill=c("gray20","gray70"),bty="n")

############################################################
#
#       update the female positions of completely salmonized assemblies using LepMap 100% salmonized mapping positions
#

salmonizedassembly$femalepos=salmonized.remaps$V5[match(salmonizedassembly$markerid,salmonized.remaps$V3)]
salmonizedassembly$malepos=salmonized.remaps$V4[match(salmonizedassembly$markerid,salmonized.remaps$V3)]

plot(salmonizedassembly$femalepos[salmonizedassembly$LG=="1"],salmonizedassembly$snppos[salmonizedassembly$LG=="1"])

############################################################
#
#       Calculate distances to next marker in first (original) and second (salmonized) assemblies
#

for(i in 2:nrow(salmonassembly)) {
	if(salmonassembly$LG[i]!=salmonassembly$LG[i-1]) { salmonassembly$femalepos_dist_to_next[i]=0 }
	else { salmonassembly$femalepos_dist_to_next[i]=salmonassembly$femalepos[i]-salmonassembly$femalepos[i-1]}
	
}
for(i in 2:nrow(salmonizedassembly)) {
	if(salmonizedassembly$LG[i]!=salmonizedassembly$LG[i-1]) { salmonizedassembly$femalepos_dist_to_next[i]=0 }
	else { salmonizedassembly$femalepos_dist_to_next[i]=salmonizedassembly$femalepos[i]-salmonizedassembly$femalepos[i-1]}
}

salmonassembly$femalepos_dist_to_next[1]=0
salmonizedassembly$femalepos_dist_to_next[1]=0


############################################################
#
#       Plot assemblies and map distance differences between marker pairs
#
plotPairwiseFemaleDist=function(assembly,lg,distPlotRange=NULL,distSdLimToHighlight=8,distHardLimToHighlight=10) {
	require(Hmisc)

	this.distToNext=assembly$femalepos_dist_to_next[assembly$LG==lg]
	this.pos=assembly$femalepos[assembly$LG==lg]
	this.markers=assembly$markerid[assembly$LG==lg]

	this.xlab="Female map pos."
	if(length(distPlotRange)!=2) {ylim=c(range(this.distToNext))}
	plot(this.pos,this.distToNext,type="l",ylab="pairwise marker dist.",xlab=this.xlab,ylim=distPlotRange)
	mtext(paste("sd =",round(sd(this.distToNext),2),"mean =",round(mean(this.distToNext),2)),cex=0.5,adj=1)
	this.result=c()
	this.highlightpos.sd=c()
	this.highlightpos.hard=c()
	if(length(distSdLimToHighlight)>0) {
		this.highlightpos.sd=this.pos[this.distToNext>(distSdLimToHighlight*sd(this.distToNext))]
		this.result=this.markers[this.distToNext>(distSdLimToHighlight*sd(this.distToNext))]
		for(p in this.highlightpos.sd) {abline(v=p,col="red",lwd=2)}
	}
	if(length(distHardLimToHighlight)>0) {
		this.highlightpos.hard=this.pos[this.distToNext>(distHardLimToHighlight)]
		this.result=c(this.result,this.markers[this.distToNext>(distHardLimToHighlight)])
		for(p in this.highlightpos.hard) {abline(v=p,col="blue",lwd=1)}
	}
	#combine found markers,reorder and remove duplicates 

	this.result=this.result[order(c(this.highlightpos.sd,this.highlightpos.hard))]
	this.result=this.result[!duplicated(this.result)]


	return(c(this.result))
}
plotPairwiseFemaleDist(salmonizedassembly,lg="1",distSdLimToHighlight=8,distHardLimToHighlight=10)


pdf("Salmonized_assemblies_salmonizedremap.v2.pdf",height=5,width=10)
for(i in 1:length(LGnames)) {
	par(mfrow=c(2,2))
	if(sum(salmonassembly$LG==LGnames[i])>1)  { 
		plotAssemblyLGs(salmonassembly,LGnames[i]) 
		mtext(side=1,cex=0.5,adj=1,"salmonized",line=2)
	}
	if(sum(salmonizedassembly$LG==LGnames[i])>1)  {
		plotAssemblyLGs(salmonizedassembly,LGnames[i]) 
		mtext(side=1,cex=0.5,adj=1,"salmonized remapped",line=3)
	}
	ymax=max(c(salmonassembly$femalepos_dist_to_next[salmonassembly$LG==LGnames[i]],salmonizedassembly$femalepos_dist_to_next[salmonizedassembly$LG==LGnames[i]]))
	plotPairwiseFemaleDist(salmonassembly,lg=LGnames[i],distPlotRange=c(0,ymax))
	mtext(side=1,cex=0.5,adj=1,"salmonized",line=2)
	mtext(side=1,cex=0.5,adj=1,"salmonized remapped",line=3)
	plotPairwiseFemaleDist(salmonizedassembly,lg=LGnames[i],distPlotRange=c(0,ymax))
	tryCatch({subplot(barplot(unlist(salmon_femalelen[salmon_femalelen$LG==LGnames[i],c("salmon1","salmon2")]),names=c("ori","remap"),las=2),max(salmonizedassembly$femalepos,na.rm=T)*0.02,ymax,size=c(0.5,0.5))},
	error=print	)  #sometimes this says "plot region too large" and does not plot but that is ok

}
dev.off()

############################################################
#
#       Make list of candidate flipping intervals. The function plots intervals in both original and salmonized assembly,
#       and prints list of potential sites in the second assembly

findFlippingCandidates=function(assembly,lg,distSdLim=10,distHardLim=8){
	
	flipcandidateposdf=data.frame(LG=NULL,fromMarker=NULL,toMarker=NULL)
	#go through this LG and find out intervals that have map distance > distLim (by default 10)
	
		if(sum(assembly$LG==lg)>2)  { 
			par(mfrow=c(2,1))
			plotAssemblyLGs(assembly,lg) 
			#these markers are breakpoints that have high enough pairwise map distance
			#they start a new potential block
			candidates=plotPairwiseFemaleDist(assembly,lg=lg,distSdLimToHighlight=distSdLim, distHardLimToHighlight=distHardLim)

			mtext(side=1,line=3,adj=1,paste("sd limit",distSdLim,sep="="),cex=0.5,col="red")
			mtext(side=1,line=4,adj=1,paste("hard limit",distHardLim,sep="="),cex=0.5,col="blue")

			###RULES:
			#if only one candidate block position, then try to flip left half
			#one marker before the only candidate block marker
		#if(length(candidates)==1) { 
		#	flipcandidateposdf=data.frame(LG=lg,fromMarker=assembly$markerid[assembly$LG==lg][1],toMarker=(assembly$markerid[which(assembly$LG==lg & assembly$markerid==candidates)-1])) 
		#}
			#if two or more cb-positions, try to flip start-1st, 1st-2nd, ... 2nd-end 
		
		if(length(candidates)>0) {
															#Add first marker of the group to candidate flipping positions to make intervals
															#if first and last markers already in list then don't add them
			if(candidates[1]!=assembly$markerid[assembly$LG==lg][1]) { candidates=c(assembly$markerid[assembly$LG==lg][1],candidates)  } 
			if(candidates[length(candidates)]!=assembly$markerid[assembly$LG==lg][sum(assembly$LG==lg)]) { candidates=c(candidates,assembly$markerid[assembly$LG==lg][sum(assembly$LG==lg)])  } 	
			
			#Then, add every intermediate interval to the output list
			#flipcandidateposdf=data.frame(LG=rep(lg,length(candidates)-1),fromMarker=candidates[1:length(candidates)-1],toMarker=candidates[2:length(candidates)]) 
			flipcandidateposdf=data.frame(LG=rep(lg,length(candidates)-1),fromMarker=candidates[1:length(candidates)-1],toMarker=c(assembly$markerid[match(candidates[2:(length(candidates)-1)],assembly$markerid)-1],candidates[length(candidates)])) 


	
		}
	}
	return(flipcandidateposdf)
		
}
#Try it!
findFlippingCandidates(salmonizedassembly,"1",distSdLim=8,distHardLim=10)
sd(salmonizedassembly$femalepos_dist_to_next[salmonizedassembly$LG=="1"])*8
findFlippingCandidates(salmonizedassembly,"1",distSdLim=8,distHardLim=10)

#plots two assemblies and the flip positions
compareAssembliesAndFlipPos=function(assembly1,assembly2,lg,distSdLim=8,distHardLim=10) {
	par(mfrow=c(2,2))
	plotAssemblyLGs(assembly1,lg)
	plotAssemblyLGs(assembly2,lg)
	ymax=max(c(assembly1$femalepos[assembly2$LG==lg],assembly2$femalepos[assembly2$LG==lg]))
	plotPairwiseFemaleDist(assembly1,lg=lg,distPlotRange=c(0,ymax),distSdLimToHighlight=distSdLim, distHardLimToHighlight=distHardLim)
	plotPairwiseFemaleDist(assembly2,lg=lg,distPlotRange=c(0,ymax),distSdLimToHighlight=distSdLim, distHardLimToHighlight=distHardLim)

	totlen1=max(assembly1$femalepos[assembly1$LG==lg],na.rm=T)
	totlen2=max(assembly2$femalepos[assembly2$LG==lg],na.rm=T)
	tryCatch({subplot(barplot(c(totlen1,totlen2),names=c("ori","remap"),las=2),max(totlen1,totlen2)*0.02,ymax,size=c(0.5,0.5))},
		error=print	)
}

compareAssembliesAndFlipPos(salmonassembly,salmonizedassembly,"1")

#Salmonized assembly: flip positions

#For every chromosome: find candidate flip positions
pdf("Salmonized_assemblies_salmonizedremap_flipcoordinates.v2.pdf",height=5,width=10)
	candidateflips.salmonized=data.frame(LG=NULL,fromMarker=NULL,toMarker=NULL)
	for(lg in LGnames) {
		result.salmonized=findFlippingCandidates(salmonizedassembly,lg,distSdLim=8,distHardLim=10)
		candidateflips.salmonized=rbind(candidateflips.salmonized,result.salmonized)
		compareAssembliesAndFlipPos(salmonassembly,salmonizedassembly,lg,distSdLim=8,distHardLim=10)
	}
dev.off()
dim(candidateflips.salmonized)

#id blocks so they can be used later on in the analysis (when flipping salmon coordinates)
candidateflips.salmonized$blockid[1]=0
for(i in 2:nrow(candidateflips.salmonized)) {
	if(candidateflips.salmonized$LG[i-1]!=candidateflips.salmonized$LG[i]) {candidateflips.salmonized$blockid[i]=0}
	else { candidateflips.salmonized$blockid[i]=candidateflips.salmonized$blockid[i-1]+1 }
}

table(as.character(candidateflips.salmonized$fromMarker)==as.character(candidateflips.salmonized$toMarker))  #no single marker groups.
write.table(candidateflips.salmonized,"salmonized.candidateflips.v3.txt",col.names=T,row.names=F,quote=F,sep="\t")



############################################################
#
#       Run Lep-MAP2 with permutating candidate block order to find the best block order. 
#       Note: 3 first blocks of LG18 were merged because the original 8 blocks would have
#       resulted in 40K permutations


#FOR SALMON

#Run lepmap for each permutation
#S:\2017_work\NORWAY_MAY2017\salmonized>c:\Perl64\bin\perl printPotentiallyInvertedOrders2.0.pl -t o -p c:\lepmap2\bin -l salmonizedorders\ssa.18.order.txt -i lg.18.candidateflips.reduced >lg.18.candidateflips.reduced_out.tx
#S:\2017_work\NORWAY_MAY2017\salmonized>c:\Perl64\bin\perl printPotentiallyInvertedOrders2.0.pl -t o -p c:\lepmap2\bin -l salmonizedorders\ssa.18.order.txt -i troutized.candidateflips.v3.txt >troutized.candidateflips.v3_out.txt

#Combine best permutations to one file

#make reordered assembly

#read the reordered maps
salmonizedreorderedassembly=salmonizedassembly
reorders.salmon=read.table("salmonizedorders/ssa_best_ordered.txt")
reorders.salmon$order=1:nrow(reorders.salmon)
reorders.salmon$V1=paste("M",reorders.salmon$V1,sep="")
head(reorders.salmon)

#combine to the salmonized assembly by sorting by new order
salmonizedreorderedassembly$bestorder_ranks=reorders.salmon$order[match(salmonizedreorderedassembly$markerid,reorders.salmon$V1)]
salmonizedreorderedassembly$bestorder_ranks[is.na(salmonizedreorderedassembly$bestorder_ranks)]=0
salmonizedreorderedassembly$malepos[salmonizedreorderedassembly$markerid %in% reorders.salmon$V1]=reorders.salmon$V2[match(salmonizedreorderedassembly$markerid[salmonizedreorderedassembly$markerid %in% reorders.salmon$V1],reorders.salmon$V1)]
salmonizedreorderedassembly$femalepos[salmonizedreorderedassembly$markerid %in% reorders.salmon$V1]=reorders.salmon$V3[match(salmonizedreorderedassembly$markerid[salmonizedreorderedassembly$markerid %in% reorders.salmon$V1],reorders.salmon$V1)]
salmonizedreorderedassembly=salmonizedreorderedassembly[order(salmonizedreorderedassembly$LG,salmonizedreorderedassembly$bestorder_ranks),]
summary(salmonizedreorderedassembly)
salmonizedreorderedassembly$femalepos_dist_to_next[1]=0
for(i in 2:nrow(salmonizedreorderedassembly)) {
	salmonizedreorderedassembly$femalepos_dist_to_next[i]=salmonizedreorderedassembly$femalepos[i]-salmonizedreorderedassembly$femalepos[i-1]
}
pdf("Salmonized_assemblies_salmonizedremap_reordered.pdf",height=5,width=10)
par(mfrow=c(2,2))
for(lg in LGnames) {
	compareAssembliesAndFlipPos(salmonizedassembly,salmonizedreorderedassembly,lg)
}
dev.off()

compareAssembliesAndFlipPos(salmonassembly,salmonizedassembly,1)
compareAssembliesAndFlipPos(salmonizedassembly,salmonizedreorderedassembly,3,distSdLim=5,distHardLim=7)

tapply(salmonizedreorderedassembly$femalepos,salmonizedreorderedassembly$LG,max)

highlightPoints=function(assembly,pointsvector) {
	points(assembly$femalepos[assembly$markerid %in% pointsvector],assembly$snppos[assembly$markerid %in% pointsvector],col="red",pch=16)
	print(assembly$femalepos[assembly$markerid %in% pointsvector])
	print(assembly$snppos[assembly$markerid %in% pointsvector])
}

############################################################
#
#       Run Lep-MAP2 with best block orders and block inversions. Compare best order vs. each inversion and
#       accept inversions that shorten map lenght. Finally combine accepted inversions. Select salmonized, ordered
#       or inverted map for each LG. Read the new map in here.

#c:\Perl64\bin\perl printPotentiallyInvertedOrders3.0.pl -m map_js.lep.lod9.txt -d data.linkage.lep -t i -p c:\lepmap2\bin -l salmonizedorders\omy.1.order.txtpermutated_best_order.txt -i troutized.candidateflips.v3.txt >troutized.candidateflips_testflips.v3_out.txt

#make reordered assembly

#read the reordered maps
salmonizedreorderedassembly=salmonizedassembly
reorders=read.table("salmonizedorders/ssa_best_orders_checked/accepted_inversions/ssa_best_ordered_inverted.txt")
reorders$order=1:nrow(reorders)
reorders$V1=paste("M",reorders$V1,sep="")
head(reorders)

#combine to the salmonized assembly by sorting by new order
salmonizedreorderedassembly$bestorder_ranks=reorders$order[match(salmonizedreorderedassembly$markerid,reorders$V1)]
salmonizedreorderedassembly$bestorder_ranks[is.na(salmonizedreorderedassembly$bestorder_ranks)]=0
salmonizedreorderedassembly$malepos[salmonizedreorderedassembly$markerid %in% reorders$V1]=reorders$V2[match(salmonizedreorderedassembly$markerid[salmonizedreorderedassembly$markerid %in% reorders$V1],reorders$V1)]
salmonizedreorderedassembly$femalepos[salmonizedreorderedassembly$markerid %in% reorders$V1]=reorders$V3[match(salmonizedreorderedassembly$markerid[salmonizedreorderedassembly$markerid %in% reorders$V1],reorders$V1)]
salmonizedreorderedassembly=salmonizedreorderedassembly[order(salmonizedreorderedassembly$LG,salmonizedreorderedassembly$bestorder_ranks),]
summary(salmonizedreorderedassembly)
salmonizedreorderedassembly$femalepos_dist_to_next[1]=0
for(i in 2:nrow(salmonizedreorderedassembly)) {
	salmonizedreorderedassembly$femalepos_dist_to_next[i]=salmonizedreorderedassembly$femalepos[i]-salmonizedreorderedassembly$femalepos[i-1]
}
pdf("Salmonized_assemblies_salmonizedremap_reordered_flipped.pdf",height=5,width=10)
par(mfrow=c(2,2))
for(lg in LGnames) {
	
	compareAssembliesAndFlipPos(salmonizedassembly,salmonizedreorderedassembly,lg)
}
dev.off()
salmonizedreorderedassembly$femalepos[salmonizedreorderedassembly$LG=="2"]


############################################################
#
#       Flip salmon positions based on accepted reordering and flipping 

#command: perl ..\printPotentiallyInvertedOrders2.0_printPermutationFromOutfile.pl ssa_best_ordered.txt > ssa_best_ordered_order.txt

#read in ordering
acceptedorders=read.table("salmonizedorders/ssa_best_ordered_order.txt",header=T)
acceptedorders$FIRSTMARKER=paste("M",acceptedorders$FIRSTMARKER,sep="")
acceptedorders$LG=salmonassembly$LG[match(acceptedorders$FIRSTMARKER,salmonassembly$markerid)]
head(acceptedorders)
acceptedorders$FLIPPEDSTARTMARKER=candidateflips.salmonized$toMarker[match(paste(acceptedorders$LG,acceptedorders$BLOCK),paste(candidateflips.salmonized$LG,candidateflips.salmonized$blockid))]
acceptedorders$FLIPPEDENDMARKER=candidateflips.salmonized$fromMarker[match(paste(acceptedorders$LG,acceptedorders$BLOCK),paste(candidateflips.salmonized$LG,candidateflips.salmonized$blockid))]


#Go through each LG in the salmonizedreorderedassembly. Each time a start marker is found, mark down the corresponding salmon position. When end marker is met,
#mark down also this positon (in the acceptedorders data frame).
#After going through the whole chromosome, 

plotAssemblyLGs(salmonizedassembly,"2")
table(salmonizedreorderedassembly$refname)
compareAssembliesAndFlipPos(salmonizedassembly,salmonizedreorderedassembly,"2")

salmonizedreorderedmanualeditassembly=readAssembly("ssa_best_ordered_inverted_manualv1.txt")
compareAssembliesAndFlipPos(salmonizedreorderedassembly,salmonizedreorderedmanualeditassembly,"2")

#add dist to next 
salmonizedreorderedmanualeditassembly$femalepos_dist_to_next=0
for(i in 2:nrow(salmonizedreorderedmanualeditassembly)) {
	if(salmonizedreorderedmanualeditassembly$LG[i]!=salmonizedreorderedmanualeditassembly$LG[i-1]) { salmonizedreorderedmanualeditassembly$femalepos_dist_to_next[i]=0 }
	else { salmonizedreorderedmanualeditassembly$femalepos_dist_to_next[i]=salmonizedreorderedmanualeditassembly$femalepos[i]-salmonizedreorderedmanualeditassembly$femalepos[i-1]}
}
salmonizedreorderedmanualeditassembly$femalepos_dist_to_next[1]=0

##################################################################################