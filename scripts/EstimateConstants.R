#EstimateConstants.R
#estimates genome-wide constants alpha1, alpha2, and beta (and prop. reads from signal)
#for ChIP-seq datasets with two replicates and one genomic input control
#requires bedtools (v2.26 or later), "parallel" package, SplitChrBed.pl, and a chromosome size file (e.g. hg19.chromsizes.tbl, with chr name in column 1 and length in bp in column 2, tab delimited)
#parallelizes across chromosomes
#by Nicolas Altemose
#2015
####If you use this program, please cite Altemose et al. eLife 2017
####This is free software shared in the hope it may be of use; no warranty is given or implied
######


##initialise inputs and outputs
library("parallel")
library(data.table)
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
# datapath="" #path to directory containing fragment position bed files
# metadatapath="" #path to directory containing gneric files - chr size file & window files
# sample = "223180_vs_220144" #sample name
# rep1suffix = "Fragment_Position_538916_223180.sorted.bed.PR1.sorted.bed" #filename for ChIP replicate 1 fragment position bed file
# rep2suffix = "Fragment_Position_538916_223180.sorted.bed.PR2.sorted.bed" #filename for ChIP replicate 2 fragment position bed file
# genomicsuffix = "Fragment_Position_538916_220144.sorted.bed" #filename for total chromatin input sample fragment position bed file


#preprocessing: create bed file with window positions across genome
windowfile=paste0(datapath,"genome.windows.",wide,"wide.",slide,"slide.bed")

if(file.exists(windowfile)){
  print(paste("windowfile",windowfile,"allready exists, using this."))
}else{
  print("Generating window file")
  system(paste0(btpath," makewindows -g ",genomesizefile," -w ",wide," -s ",slide," >",windowfile))
}


# Fragment Position Files
posfileA = paste0(datapath,rep1suffix)
posfileB = paste0(datapath,rep2suffix)
posfileG = paste0(datapath,genomicsuffix)

#declare temporary intermediate filenames
#system(paste0("mkdir ",datapath,"EstimateConstants"),ignore.stdout = T, ignore.stderr = T)
infile1=paste0(datapath,"EstimateConstants_FragCount.",rep1suffix,".",wide,"wide.",slide,"slide.bed")
infile2=paste0(datapath,"EstimateConstants_FragCount.",rep2suffix,".",wide,"wide.",slide,"slide.bed")
infile3=paste0(datapath,"EstimateConstants_FragCount.",genomicsuffix,".",wide,"wide.",slide,"slide.bed")

# Calculate Frag Count Overlaps
recalculate_coverage=FALSE
if(all(file.exists(infile1, infile2, infile3)) & !recalculate_coverage){
  print("Skipping recalculation of coverage")
}else{

  print("Calculating coverage")

  posfiles = c(posfileA, posfileB, posfileG)
  infiles = c(infile1, infile2, infile3)

  calc_frag_overlap_counts <- function(i, posfilevec=posfiles, infilevec=infiles, windows=windowfile){
    system(paste0(btpath," coverage -a ",windows," -b ",posfiles[i]," -counts | cut -f4 > ",infiles[i]))
  }

  noreturn = mclapply(1:3, calc_frag_overlap_counts, mc.preschedule=TRUE, mc.cores=3)

  # non parallel version of above
  # system(paste0(btpath," coverage -a ",windowfile," -b ",posfileA," -counts | cut -f4 >",infile1))
  # system(paste0(btpath," coverage -a ",windowfile," -b ",posfileB," -counts | cut -f4 >",infile2))
  # system(paste0(btpath," coverage -a ",windowfile," -b ",posfileG," -counts | cut -f4 >",infile3))
}


#set output file name
outfile1 = paste0(datapath,"Constants.",sample,".tsv")
rep1=3
rep2=4
genomic=5

#read in and combine fragment coverage values into one dataframe "counts"
print("Reading in fragment count coverage")
counts = fread(windowfile, col.names=c('chr','start','stop'))
counts$countA <- fread(infile1)
counts$countB <- fread(infile2)
counts$countG <- fread(infile3)

setkey(counts, chr, start, stop)
counts <- counts[chr %in% paste0("chr",chrs)]

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
	peaks=makepeaks(counts,alpha1=alpha1.est,alpha2=alpha2.est,beta=beta.est0,r1=rep1,r2=rep2,g=genomic)
	q=which(!is.na(peaks[,"p-value"]) & peaks[,"p-value"]<1e-10 & peaks[,"yhat_alt"]>0)
	if(length(q)<100){
		q=which(!is.na(peaks[,"p-value"]) & peaks[,"p-value"]<1e-5 & peaks[,"yhat_alt"]>0)
	}
	rm(peaks)
	beta.est = (mean(counts[q,rep2])-alpha2.est*mean(counts[q,genomic]))/(mean(counts[q,rep1])-alpha1.est*mean(counts[q,genomic]))

	#now redo p-value calls with new beta estimate
	peaks=makepeaks(counts,alpha1=alpha1.est,alpha2=alpha2.est,beta=beta.est,r1=rep1,r2=rep2,g=genomic)
	gthresh = quantile(peaks[,"cov_g"],0.999)


	#return constant values and estimates of signal/background in each replicate
	r1comb=mean(peaks[,"bhat_alt"]*(peaks[,"yhat_alt"]+alpha1.est))
	r1sig=mean(peaks[,"bhat_alt"]*peaks[,"yhat_alt"])
	r2comb=mean(peaks[,"bhat_alt"]*(beta.est*peaks[,"yhat_alt"]+alpha2.est))
	r2sig=mean(peaks[,"bhat_alt"]*peaks[,"yhat_alt"]*beta.est)
	signifbins = sum(peaks[,"p-value"]<1e-5 & peaks[,"yhat_alt"]>0,na.rm=TRUE)
	return(c(chr,alpha1.est,alpha2.est,beta.est,mean(counts[,genomic]),mean(counts[,rep1]),mean(counts[,rep2]),as.integer(dim(counts)[1]),zeroregions1,zeroregions2,length(q),r1sig/r1comb,r2sig/r2comb,gthresh,signifbins))
}


#declare functions to find MLE values for each window

makepeaks=function(test=counts,alpha1,alpha2,beta,r1=rep1,r2=rep2,g=genomic){

	sumcov = test[,r1]+test[,r2]+test[,g]

	term1=(sumcov)*(beta+1)
	term2=1+alpha1+alpha2
	term3=beta+1
	term4=beta*alpha1*test[,r2]+alpha2*test[,r1]
	term5=beta*(test[,r1]+test[,r2])
	term6=alpha1*alpha2
	term7=alpha1*beta+alpha2

	aterm=term1*beta-term3*term5
	bterm=term1*term7-term2*term5-term3*term4
	cterm=term1*term6-term2*term4

	rm(term1,term2,term3,term4,term5,term6,term7)

	yvals=(-bterm+sqrt(bterm^2-4*aterm*cterm))/2/aterm
	rm(aterm,bterm,cterm)

	yvals[yvals<0]=0
	yvals[yvals>1e9]=1e9

	bvals=(sumcov)/(1+alpha1+alpha2+(beta+1)*yvals)

	bvalsnull=(sumcov)/(1+alpha1+alpha2)
	yvalsnull=rep(0,length(bvalsnull))

	lhooddiff=2*(lhood(yvals,bvals,test,alpha1,alpha2,beta,r1,r2,g)-lhood(yvalsnull,bvalsnull,test,alpha1,alpha2,beta,r1,r2,g))

	rm(bvalsnull,yvalsnull)

	signif=pchisq(lhooddiff,df=1,lower.tail=F)

	results=cbind(test[,1],test[,2],yvals,signif,lhooddiff,test[,r1],test[,r2],test[,g],bvals)
	colnames(results)=c("start","stop","yhat_alt","p-value","Lhood_diff","cov_r1","cov_r2","cov_g","bhat_alt")
	return(results)

}
lhood=function(yhat,bhat,test,alpha1,alpha2,beta,r1=rep1,r2=rep2,g=genomic){
	sumcov = test[,r1]+test[,r2]+test[,g]
	ourterm=sumcov*(log(bhat)-1)+test[,r1]*log(alpha1+yhat)+test[,r2]*log(alpha2+yhat*beta)
	return(ourterm)
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
write.table(data2,file=outfile1,quote=FALSE,sep="\t",row.names=F,col.names=coln)

print(paste("printed results to",outfile1))
quit(save="no",runLast=FALSE)
