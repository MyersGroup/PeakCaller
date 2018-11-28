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
options(scipen=20)
btpath = "bedtools" #path to bedtools executable
system("mkdir tmp",ignore.stdout = T, ignore.stderr = T)

args=commandArgs(TRUE)
datapath = args[1] #path to folder containing fragment position bed files
posfileIP1base = args[2] #filename for ChIP replicate 1 fragment position bed file
posfileIP2base = args[3] #filename for ChIP replicate 2 fragment position bed file
posfileGbase = args[4] #filename for total chromatin input sample fragment position bed file
constfile = args[5] #full path to a file containing genome-wide estimates for constants alpha1/2 & beta (output of EstimateConstants.R)
bedfilebase = args[6] #full path to a 3-column bed file listing positions of windows in which to do force-calling
outfile = args[7] #path and filename of output file

##example hardwired input
#constfile = "Constants.YFP_HumanPRDM9.antiH3K4me3.100wide.100slide.txt"
#bedfilebase = "PromoterRegions.1kb.bed"
#posfileIP1base = "FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed.PR1.sorted.bed"
#posfileIP2base = "FragPos.YFP_HumanPRDM9.antiH3K4me3ProtocolN.bed.PR2.sorted.bed"
#posfileGbase = "FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed.sorted.bed"
#outfile = "ChIPseq_ForceCalling.HumanH3K4me3_on_YFP_HumanPRDM9_peaks.protocolN.p10e-5.sep1000.txt"

#create vector of all chromosome names at which to call peaks (change if necessary) 
chrs=seq(1,22,1) #change if chromosome number differs
chrs=c(chrs,"X")

#read in constants
constdata=read.table(constfile,header=TRUE)
alpha1.est=constdata[which(constdata[,1]=="autosomal"),2]
alpha2.est=constdata[which(constdata[,1]=="autosomal"),3]
beta.est=constdata[which(constdata[,1]=="autosomal"),4]
rep1=3
rep2=4
genomic=5

system(paste("perl SplitChrBed.pl ",bedfilebase," ","tmp/",bedfilebase,sep=""))


#declare function for each chromosome
getEnrichments=function(chr){
	posfileIP1 = paste(datapath,"/bychr/",posfileIP1base,".chr",chr,".bed",sep="")
	posfileIP2 = paste(datapath,"/bychr/",posfileIP2base,".chr",chr,".bed",sep="")
	posfileG= paste(datapath,"/bychr/",posfileGbase,".chr",chr,".bed",sep="")
	bedfile= paste("tmp/",bedfilebase,".chr",chr,".bed",sep="")
	tempfile1=paste("tmp/",outfile,".chr",chr,".temp1.bed",sep="")
	tempfile2=paste("tmp/",outfile,".chr",chr,".temp2.bed",sep="")
	tempfileG=paste("tmp/",outfile,".chr",chr,".tempG.bed",sep="")
	outfiletemp = paste("tmp/",outfile,".chr",chr,".bed",sep="")

	#now get coverage at windows specified in bedfile and combine into a dataframe "counts"
	system(paste(btpath," coverage -a ",posfileIP1," -b ",bedfile," -counts >",tempfile1,sep=""))
	system(paste(btpath," coverage -a ",posfileIP2," -b ",bedfile," -counts >",tempfile2,sep=""))
	system(paste(btpath," coverage -a ",posfileG," -b ",bedfile," -counts >",tempfileG,sep=""))
	
	counts = read.table(tempfile1,header=FALSE,colClasses=c('character','integer','integer','integer'))
	counts=counts[order(counts[,2]),]
	countstemp = read.table(tempfile2,header=FALSE,colClasses=c('NULL','integer','NULL','integer'))
	countstemp=countstemp[order(countstemp[,1]),]
	counts[,5]=countstemp[,2]
	countstemp = read.table(tempfileG,header=FALSE,colClasses=c('NULL','integer','NULL','integer'))
	countstemp=countstemp[order(countstemp[,1]),]
	counts[,6]=countstemp[,2]
	rm(countstemp)
	
	#remove regions with 0 coverage in all samples, and pseudocount regions with 0 genomic coverage but nonzero IP coverage
	countsfilt=counts[,2:6]
	countsfilt[(countsfilt[,genomic]==0 & (countsfilt[,rep1]+countsfilt[,rep2])>0),genomic]=0.5 #pseudocount regions with 0 genomic coverage and >0 IP coverage to have genomic coverage of 0.5
	countsfilt=countsfilt[countsfilt[,genomic]>0,] #remove remaining regions with 0 genomic coverage (i.e. regions with 0 rep1, 0 rep2, and 0 genomic)

	peaks=makepeaks(countsfilt,alpha1=alpha1.est,alpha2=alpha2.est,beta=beta.est,r1=rep1,r2=rep2,g=genomic)

	validcount = 0
	for (i in 1:dim(counts)[1]){
		startpos = counts[i,2]
		peakinfo = peaks[peaks[,1]==startpos,]
		if(length(peakinfo)>0){
			counts[i,7]=peakinfo[3]
			counts[i,8]=peakinfo[5]
			counts[i,9]=peakinfo[4]
			validcount=validcount+1
		}else{
			counts[i,7]="NA"
			counts[i,8]="NA"
			counts[i,9]="NA"
		}
	}

	colnames(counts)=c("chr","center_start","center_stop","cov_r1","cov_r2","cov_input","enrichment","likelihood","pvalue")

	#write final output file
	options(scipen=8)
	write.table(counts,file=outfiletemp,quote=FALSE,sep="\t",row.names=F,col.names=T)
	return(1)
}


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
	
	rm(term1,term4,term5,term6,term7)

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

	results=cbind(test[,1],test[,2],yvals,signif,lhooddiff,test[,r1],test[,r2],test[,g])
	colnames(results)=c("start","stop","yhat_alt","p-value","Lhood_diff","cov_r1","cov_r2","cov_g")
	return(results)
	
}
lhood=function(yhat,bhat,test,alpha1,alpha2,beta,r1=rep1,r2=rep2,g=genomic){
	sumcov = test[,r1]+test[,r2]+test[,g]
	ourterm=sumcov*(log(bhat)-1)+test[,r1]*log(alpha1+yhat)+test[,r2]*log(alpha2+yhat*beta)
	return(ourterm)
}


#run all chromosomes in parallel and write output file
print(date())
funfunc = mclapply(chrs,getEnrichments,mc.preschedule=TRUE,mc.cores=20)
print(paste("done!:",date()))

finaltable=read.table(paste("tmp/",outfile,".chr",chrs[1],".bed",sep=""), header=T)
for(chr in chrs[2:length(chrs)]){
	outfiletemp = paste("tmp/",outfile,".chr",chr,".bed",sep="")
	tempdata1=read.table(outfiletemp, header=T)
	finaltable=rbind(finaltable,tempdata1)
}
write.table(finaltable,file=outfile,col.names=T,row.names=F,quote=F)

quit(save="no",runLast=FALSE)
