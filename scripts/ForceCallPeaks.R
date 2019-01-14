#ForceCallPeaks.R
#implementation of an algorithm to perform LR testing of ChIP-seq data at a fixed set of sites
#given two ChIP replicates and one genomic input control
#call EstimateConstants.R before running
#by Nicolas Altemose
#2015
######
####If you use this program, please cite Altemose et al. eLife 2017
####This is free software shared in the hope it may be of use; no warranty is given or implied


##initialise inputs and outputs
library("parallel")
library(data.table)
source("functions.R")
options(scipen=20)
btpath = "bedtools" #path to bedtools executable
system(paste0("mkdir ",datapath,"forceCallCoverage"), ignore.stdout = T, ignore.stderr = T)

args=commandArgs(TRUE)
datapath = args[1] #path to folder to store intermediate files
posfileIP1base = args[1] #filename for ChIP replicate 1 fragment position bed file
posfileIP2base = args[2] #filename for ChIP replicate 2 fragment position bed file
posfileGbase = args[3] #filename for total chromatin input sample fragment position bed file
constfile = args[4] #full path to a file containing genome-wide estimates for constants alpha1/2 & beta (output of EstimateConstants.R)
bedfilebase = args[5] #full path to a 3-column bed file listing positions of windows in which to do force-calling
autosomal_chrs = args[6]
outfile = args[7] #path and filename of output file


##example hardwired input
#constfile = "Constants.YFP_HumanPRDM9.antiH3K4me3.100wide.100slide.txt"
#bedfilebase = "PromoterRegions.1kb.bed"
#posfileIP1base = "FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed.PR1.sorted.bed"
#posfileIP2base = "FragPos.YFP_HumanPRDM9.antiH3K4me3ProtocolN.bed.PR2.sorted.bed"
#posfileGbase = "FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed.sorted.bed"
#outfile = "ChIPseq_ForceCalling.HumanH3K4me3_on_YFP_HumanPRDM9_peaks.protocolN.p10e-5.sep1000.txt"

#create vector of all chromosome names at which to call peaks (change if necessary)
chrs=c(1:autosomal_chrs ,"X")

#read in constants
constdata=read.table(constfile,header=TRUE)
alpha1.est=constdata[which(constdata[,1]=="autosomal"),2]
alpha2.est=constdata[which(constdata[,1]=="autosomal"),3]
beta.est=constdata[which(constdata[,1]=="autosomal"),4]


tempfile1=paste0(datapath,"forceCallCoverage/FragCount.",basename(outfile),".temp1.bed")
tempfile2=paste0(datapath,"forceCallCoverage/FragCount.",basename(outfile),".temp2.bed")
tempfileG=paste0(datapath,"forceCallCoverage/FragCount.",basename(outfile),".tempG.bed")

counts <- get_frag_overlap_counts(posfiles = c(posfileIP1base, posfileIP2base, posfileGbase),
                                  infiles = c(tempfile1, tempfile2, tempfileG),
                                  windowfilepath = bedfilebase)


#declare function for each chromosome
getEnrichments=function(chr){

	chr0 = paste0("chr",chr)
	counts <- counts[chr==chr0][,.(start, stop, countA, countB, countG)]

	peaks = compute_likelihood(counts[,countA], counts[,countB], counts[,countG], alpha1=alpha1.est, alpha2=alpha2.est, beta=beta.est, return_slim=F)

	counts <- cbind("chr"=chr0, counts, peaks)

	setnames(counts, c("chr","center_start","center_stop","cov_r1","cov_r2","cov_input","enrichment","likelihood","pvalue","bhat"))

	return(counts)
}


#run all chromosomes in parallel and write output file
print(date())
data = mclapply(chrs,getEnrichments,mc.preschedule=TRUE,mc.cores=length(chrs))
data <- rbindlist(data)
data$bhat <- NULL
print(paste("done!:",date()))

write.table(data,file=outfile,quote=F,sep="\t",row.names=F,col.names=T)
