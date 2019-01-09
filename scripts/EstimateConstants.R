#EstimateConstants.R
#estimates genome-wide constants alpha1, alpha2, and beta (and prop. reads from signal)
#for ChIP-seq datasets with two replicates and one genomic input control
#requires bedtools (v2.26 or later), "parallel" package, and a chromosome size file (e.g. hg19.chromsizes.tbl, with chr name in column 1 and length in bp in column 2, tab delimited)
#parallelizes across chromosomes
#by Nicolas Altemose
#2015
####If you use this program, please cite Altemose et al. eLife 2017
####This is free software shared in the hope it may be of use; no warranty is given or implied
######

print(paste("Start!",format(Sys.time(), "%Y.%m.%d %H:%M:%S")))

##initialise inputs and outputs
library("parallel")
library(data.table)
source("functions.R")
options(scipen=20)
btpath = "bedtools" #full path to bedtools excutable file
genomesizefile = "hg38.sizes" #default path to chr size file
wide=100 #size of bins to do initial testing to estimate constants (100 is suitable)
slide=100 #distance between bin starting positions (100 is suitable)

args=commandArgs(TRUE)
datapath = args[1]
genomesizefile = args[2]
sample = args[3]
rep1suffix = args[4]
rep2suffix = args[5]
genomicsuffix = args[6]
autosomal_chrs = args[7]

#create vector of all chromosome names at which to call peaks (change if necessary)
chrs=c(1:autosomal_chrs ,"X")


#example hardwired input
# datapath="test/" #path to directory containing fragment position bed files
# sample = "223180_vs_220144" #sample name
# rep1suffix = "Fragment_Position_538916_223180.sorted.bed.PR1.sorted.bed" #filename for ChIP replicate 1 fragment position bed file
# rep2suffix = "Fragment_Position_538916_223180.sorted.bed.PR2.sorted.bed" #filename for ChIP replicate 2 fragment position bed file
# genomicsuffix = "Fragment_Position_538916_220144.sorted.bed" #filename for total chromatin input sample fragment position bed file
# autosomal_chrs = 22

#preprocessing: create bed file with window positions across genome
windowfile=paste0(datapath,"genome.windows.",wide,"wide.",slide,"slide.bed")

if(file.exists(windowfile)){
  print(paste("windowfile",windowfile,"allready exists, using this."))
}else{
  print("Generating window file")
  system(paste0(btpath," makewindows -g ",genomesizefile," -w ",wide," -s ",slide," >",windowfile))
}


#declare temporary intermediate filenames
#system(paste0("mkdir ",datapath,"EstimateConstants"),ignore.stdout = T, ignore.stderr = T)
infile1=paste0(datapath,"EstimateConstants_FragCount.",basename(rep1suffix),".",wide,"wide.",slide,"slide.bed")
infile2=paste0(datapath,"EstimateConstants_FragCount.",basename(rep2suffix),".",wide,"wide.",slide,"slide.bed")
infile3=paste0(datapath,"EstimateConstants_FragCount.",basename(genomicsuffix),".",wide,"wide.",slide,"slide.bed")

# Calculate Frag Count Overlaps
counts <- get_frag_overlap_counts(posfiles = c(rep1suffix, rep2suffix, genomicsuffix),
                                  infiles = c(infile1, infile2, infile3),
                                  windowfilepath = windowfile)

rep1=3
rep2=4
genomic=5

#declare function to call peaks in bins for one chromosome
getConstants=function(chr){

  chr0 = paste0("chr",chr)
  counts <- data.frame(counts[chr==chr0][,.(start, stop, countA, countB, countG)])

	#provide initial rough estimates for constants alpha1, alpha2, and beta
	alpha1.est = sum(counts[(counts[,rep2]==0),rep1])/sum(counts[(counts[,rep2]==0),genomic])
	alpha2.est = sum(counts[(counts[,rep1]==0),rep2])/sum(counts[(counts[,rep1]==0),genomic])
	beta.est0 = (mean(counts[,rep2])-alpha2.est*mean(counts[,genomic]))/(mean(counts[,rep1])-alpha1.est*mean(counts[,genomic]))

	#identify regions where one IP replicate is nonzero and the other is 0
	zeroregions1=sum((counts[,rep2]==0) & (counts[,rep1] + counts[,genomic])>0)
	zeroregions2=sum((counts[,rep1]==0) & (counts[,rep2] + counts[,genomic])>0)

	#remove regions with 0 coverage in all samples, and pseudocount regions with 0 genomic coverage but nonzero IP coverage
	counts[(counts[,genomic]==0 & (counts[,rep1]+counts[,rep2])>0),genomic]=0.5 #pseudocount regions with 0 genomic coverage and >0 IP coverage to have genomic coverage of 0.5
	counts=counts[counts[,genomic]>0,] #remove remaining regions with 0 genomic coverage (i.e. regions with 0 rep1, 0 rep2, and 0 genomic)

	#find initial set of p-values, identify confident set of peaks, re-do estimate of beta at these sites
	peaks = compute_likelihood(counts[,rep1], counts[,rep2], counts[,genomic], alpha1=alpha1.est, alpha2=alpha2.est, beta=beta.est0, return_slim=F)

	q=which(!is.na(peaks[,"p-value"]) & peaks[,"p-value"]<1e-10 & peaks[,"yhat_alt"]>0)
	if(length(q)<100){
		q=which(!is.na(peaks[,"p-value"]) & peaks[,"p-value"]<1e-5 & peaks[,"yhat_alt"]>0)
	}
	rm(peaks)
	beta.est = (mean(counts[q,rep2])-alpha2.est*mean(counts[q,genomic]))/(mean(counts[q,rep1])-alpha1.est*mean(counts[q,genomic]))

	#now redo p-value calls with new beta estimate
	peaks = compute_likelihood(counts[,rep1], counts[,rep2], counts[,genomic], alpha1=alpha1.est, alpha2=alpha2.est, beta=beta.est0, return_slim=F)
	gthresh = quantile(counts[,genomic],0.999) #peaks[,"cov_g"]


	#return constant values and estimates of signal/background in each replicate
	r1comb=mean(peaks[,"bhat_alt"]*(peaks[,"yhat_alt"]+alpha1.est))
	r1sig=mean(peaks[,"bhat_alt"]*peaks[,"yhat_alt"])
	r2comb=mean(peaks[,"bhat_alt"]*(beta.est*peaks[,"yhat_alt"]+alpha2.est))
	r2sig=mean(peaks[,"bhat_alt"]*peaks[,"yhat_alt"]*beta.est)
	signifbins = sum(peaks[,"p-value"]<1e-5 & peaks[,"yhat_alt"]>0,na.rm=TRUE)
	return(c(chr,alpha1.est,alpha2.est,beta.est,mean(counts[,genomic]),mean(counts[,rep1]),mean(counts[,rep2]),as.integer(dim(counts)[1]),zeroregions1,zeroregions2,length(q),r1sig/r1comb,r2sig/r2comb,gthresh,signifbins))
}


#declare final output column names
coln= c("chr","alpha1","alpha2","beta","meancovgenomic","meancovrep1","meancovrep2","totalnonzerobins","alpha1trainingregions","alpha2trainingregions","betatrainingregions","rep1signal","rep2signal","genomiccov999thpctile","significantbins1e-5")

#call functions on all chromosomes in parallel and combine results
# not actually nececcary to run in parallel (24 seconds vs 6 seconds)
print("Calculating Constants")
data=mclapply(chrs,getConstants,mc.preschedule=TRUE,mc.cores=length(chrs))
data2=t(simplify2array(data))
data2[,1]=unlist(lapply(data,function(x) as.character(x[[1]])))

#compute average values across autosomes
data2=rbind(data2,c("autosomal",rep("NA",14))) # 14 columns
for(m in c(2,3,4,5,6,7,12,13,14)){
	data2[nrow(data2),m]=weighted.mean(as.numeric(data2[1:autosomal_chrs,m]),as.numeric(data2[1:autosomal_chrs,8]))
}
for(m in c(8,9,10,11,15)){
	data2[nrow(data2),m]=sum(as.numeric(data2[1:autosomal_chrs,m]))
}

#write final output file with constant estimates
write.table(data2,file=paste0(datapath,"Constants.",sample,".tsv"),quote=FALSE,sep="\t",row.names=F,col.names=coln)

print(paste("printed results to",paste0(datapath,"Constants.",sample,".tsv")))
quit(save="no",runLast=FALSE)

