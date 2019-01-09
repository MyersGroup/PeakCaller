#DeNovoPeakCalling.SingleBase.R
#A peak calling algorithm that computes single-base enrichment and likelihood values
#for ChIP-seq datasets with two replicates and one genomic input control
#call EstimateConstants.R before running
#by Nicolas Altemose
#2015
######
####If you use this program, please cite Altemose et al. eLife 2017
####This is free software shared in the hope it may be of use; no warranty is given or implied


##initialise inputs and outputs
library("parallel")
library("IRanges")
source("functions.R")
options("scipen"=8)
btpath = "bedtools"
batchsize=1000000 #how many bases to process in each batch of memory?
cithresh=0.99 #confidence interval threshold for output
searchwidth = 1000 #maximum CI width
genomesizefile = "hg38.sizes" #default path to chr size file


args=commandArgs(TRUE)
datapath = args[1] # intermediate files
genomesizefile = args[2]
outfileALL = args[3] #path to outfile
rep1suffix = args[4] #filename for ChIP replicate 1 fragment position bed file
rep2suffix = args[5] #filename for ChIP replicate 2 fragment position bed file
genomicsuffix = args[6] #filename for total chromatin input sample fragment position bed file
autosomal_chrs = as.integer(args[7]) #see above
pvalthresh = as.numeric(args[8]) #maximum p-value at peak centre
minsep = as.numeric(args[9]) #minimum separation between peak centres
constfile = args[10]


#example hardwired input
# datapath="test/"
# genomesizefile="hg38.sizes"
# rep1suffix = "Fragment_Position_538916_223180.sorted.bed.PR1.sorted.bed"
# rep2suffix = "Fragment_Position_538916_223180.sorted.bed.PR2.sorted.bed"
# genomicsuffix = "Fragment_Position_538916_220144.sorted.bed"
# pvalthresh = 0.000001
# minsep = 250
# autosomal_chrs=22


#create vector of all chromosome names at which to call peaks (change if necessary)
chrs=c(1:autosomal_chrs ,"X")

lthresh = qchisq(pvalthresh,df=1,lower.tail=F)
cilhood=qchisq(cithresh,df=2)

chrlengths = read.table(genomesizefile,header=F,colClasses=c('character','integer'))


#read in constants estimated previously by EstimateConstants.R, calculate constant terms used in likelihood calculations

constdata=read.table(constfile,header=TRUE)
alpha1.est=constdata[which(constdata[,1]=="autosomal"),2]
alpha2.est=constdata[which(constdata[,1]=="autosomal"),3]
beta.est=constdata[which(constdata[,1]=="autosomal"),4]

#declare output filenames
system(paste0("mkdir ",datapath,"bychr"), ignore.stdout = T, ignore.stderr = T)
system(paste0("mkdir ",datapath,"bychr/covzip"), ignore.stdout = T, ignore.stderr = T)
system(paste0("mkdir ",datapath,"bychr/covbin"), ignore.stdout = T, ignore.stderr = T)
system(paste0("mkdir ",datapath,"bychr/likelihoods"), ignore.stdout = T, ignore.stderr = T)
system(paste0("mkdir ",datapath,"bychr/enrichments"), ignore.stdout = T, ignore.stderr = T)
system(paste0("mkdir ",datapath,"bychr/peaks"), ignore.stdout = T, ignore.stderr = T)


singleBaseCoverageFP <- function(chr, FragPosFile, bedtools, covfile, genomesizefile){

  # 1) extract single chr from fragpos file (\\t to distinguish chr1 vs e.g chr11)
  # 2) compute bedgraph
  # 3) extract single chr
  # 4) epand to single base pair resolution

  command <- paste0("grep -P 'chr",chr,"\\t' ",FragPosFile," | \\
                            ",bedtools," genomecov -bga -i stdin -g ",genomesizefile," | \\
                            grep -P 'chr",chr,"\\t' | \\
                            perl -lane 'for ($F[1]+1..$F[2]) { print \"$F[0]\\t$_\\t$F[3]\" }' | \\
                            cut -f3 | \\
                            gzip > ",covfile)
  system(command)
}

#declare a function to perform calculations on one chromosome
getEnrichments=function(chr){

	#define input and output file names

  posfile = c("A"=paste0(rep1suffix),
              "B" = paste0(rep2suffix),
              "G" = paste0(genomicsuffix))

	covfile = c("A"= paste0(datapath,"bychr/covzip/",basename(rep1suffix),".FragDepth.chr",chr,".bed.gz"),
	            "B" = paste0(datapath,"bychr/covzip/",basename(rep2suffix),".FragDepth.chr",chr,".bed.gz"),
	            "G" = paste0(datapath,"bychr/covzip/",basename(genomicsuffix),".FragDepth.chr",chr,".bed.gz"))

	covfilebin = c("A"=paste0(datapath,"bychr/covbin/",basename(rep1suffix),".FragDepth.chr",chr,".binary.gz"),
	               "B"=paste0(datapath,"bychr/covbin/",basename(rep2suffix),".FragDepth.chr",chr,".binary.gz"),
	               "G"=paste0(datapath,"bychr/covbin/",basename(genomicsuffix),".FragDepth.chr",chr,".binary.gz"))

	outfileLhood = paste0(datapath,"bychr/likelihoods/SingleBaseLikelihood.",basename(outfileALL),".chr",chr,".binary.r")
	outfileEnrich = paste0(datapath,"bychr/enrichments/SingleBaseEnrichment.",basename(outfileALL),".chr",chr,".binary.r")
	outfilePeaks = paste0(datapath,"bychr/peaks/SingleBasePeaks.",basename(outfileALL),".p",pvalthresh,".sep",minsep,".chr",chr,".bed")

	chrlen = chrlengths[which(chrlengths[,1]==paste0("chr",chr)),2]


	#if not done already, compute single-base coverage values across chromosome and compress
	for (i in c("A","B","G")){
	  if(file.exists(covfilebin[i])){
	    print(paste(covfilebin[i], "already exists, skipping coverage calculation"))
	  }else{
	    singleBaseCoverageFP(chr, posfile[i], btpath, covfile[i], genomesizefile)

	    con=gzfile(covfile[i],open="r")
	    conb=gzfile(covfilebin[i],open="wb")
	    writeBin(scan(con,what=integer(1),nlines=chrlen,quiet=TRUE),conb)
	    close(con)
	    close(conb)

	    print(paste("done computing coverage for chr:",chr,format(Sys.time(), "%Y.%m.%d %H:%M:%S")))
	  }
	}

	#if not done already, compute single-base likelihood and enrichment values across chromosome, in batches
	if(all(file.exists(outfileLhood, outfileEnrich))){
	  print(paste("files",outfileLhood, outfileEnrich, "already exist, skipping recalculation of likelihood & enrichment"))
	}else{
		conAb=gzfile(covfilebin["A"],open="rb")
		conBb=gzfile(covfilebin["B"],open="rb")
		conGb=gzfile(covfilebin["G"],open="rb")
		conOUTL=gzfile(outfileLhood,open="wb")
		conOUTE=gzfile(outfileEnrich,open="wb")
		startpos=0
		while(startpos<chrlen){
			basenum=batchsize
			if((startpos+batchsize)>chrlen){
				basenum=chrlen-startpos
			}
			covA = readBin(conAb,what=integer(1),n=basenum)
			covB = readBin(conBb,what=integer(1),n=basenum)
			covG = readBin(conGb,what=integer(1),n=basenum)

			basepeaks = list(rep(0,basenum),rep(0,basenum))
			if(sum(covA+covB+covG)>0){
				basepeaks=compute_likelihood(r1=covA, r2=covB, g=covG, alpha1=alpha1.est, alpha2=alpha2.est, beta=beta.est)
			}

			writeBin(basepeaks[[1]],conOUTL)
			writeBin(basepeaks[[2]],conOUTE)

			startpos = startpos+batchsize
		}
		close(conAb)
		close(conBb)
		close(conGb)
		close(conOUTL)
		close(conOUTE)
		print(paste("done computing likelihoods for chr:",chr,format(Sys.time(), "%Y.%m.%d %H:%M:%S")))
	}


	#read in all single-base log likelihood values across chromosome
	con1=gzfile(outfileLhood,open="rb")
	vec=readBin(con1,what=numeric(1),n=chrlen)
	close(con1)

	#find bases above pvalue threshold that are local maxima (> minsep bases left and >= minsep bases to the right)
	q0=which(vec>=lthresh)
	q=q0[q0>searchwidth & q0<(chrlen-searchwidth)]
	rm(q0)
	qtest = vapply(q,function(x) (vec[x]>max(vec[(x-minsep):(x-1)]) && vec[x]>=max(vec[(x+1):(x+minsep)])),USE.NAMES=F,FUN.VALUE=logical(1))
	newq=q[qtest]
	rm(q)

	#around these bases find confidence intervals
	getci=function(i){
		left = i-searchwidth+max(0,which(vec[(i-searchwidth):(i-1)]<=(vec[i]-cilhood)))
		right = i + min(searchwidth,which(vec[(i+1):(i+searchwidth)]<=(vec[i]-cilhood)))
		return(c(i,vec[i],left,right))
	}

	confints = t(sapply(newq,getci,USE.NAMES=F,simplify="array"))
	rm(newq)

	#sort by CI width, find overlaps, indicate those not overlapping shorter interval
	confints=confints[order(confints[,4]-confints[,3]),]
	query = IRanges(confints[,3],confints[,4])
	overlaps = as.matrix(findOverlaps(query))

	includevec = rep(0,dim(confints)[1])
	includevec[1]=1
	for(i in 1:dim(confints)[1]){
		includevec[i]=1
		if(sum(includevec[overlaps[overlaps[,1]==i & overlaps[,2]<i,2]])>0){
			includevec[i]=0
		}
	}
	confints=cbind(confints,includevec)


	#refine centres within each CI (take mean position of all bases with llhood=max)
	refinecentres=function(i){
		floor(mean(range(confints[i,3]-1+which(vec[(confints[i,3]):(confints[i,4])]==confints[i,2]))))
	}
	newcentres=vapply(1:dim(confints)[1],refinecentres,USE.NAMES=F,FUN.VALUE=numeric(1))
	confints[,1]=newcentres
	confints=confints[order(confints[,1]),]
	rm(vec)

	#add in coverage and enrichment values and pvalues
	addtoConfints <- function(binfile, dtype=integer(1), condifence_intervals = confints){
	  conB = gzfile(binfile, open="rb")
	  cov = readBin(conB, what=dtype, n=chrlen)
	  close(conB)
	  condifence_intervals = cbind(condifence_intervals, cov[condifence_intervals[,1]])
	  rm(cov)
	  return(condifence_intervals)
	}

	confints = addtoConfints(covfilebin["A"])
	confints = addtoConfints(covfilebin["B"])
	confints = addtoConfints(covfilebin["G"])
	confints = addtoConfints(outfileEnrich, dtype = numeric(1))

	signif=pchisq(confints[,2],df=1,lower.tail=F)
	confints=cbind(confints,signif)

	#print final peaks to text file
	sub=confints[confints[,5]==1,]
	write.table(cbind(paste0("chr",chr),sub[,1]-1,sub[,1],sub[,3]-1,sub[,c(4,6,7,8,9,2,10)]),file=outfilePeaks,quote=F,row.names=F,col.names=F,sep="\t")
	print(paste("done getting peak intervals for chr:",chr,format(Sys.time(), "%Y.%m.%d %H:%M:%S")))
	return(1)
}


#call function on all chromosomes in parallel
print(format(Sys.time(), "%Y.%m.%d %H:%M:%S"))
funfunc = mclapply(chrs,getEnrichments,mc.preschedule=TRUE,mc.cores=length(chrs))

allpeaks = read.table(paste0(datapath,"bychr/peaks/SingleBasePeaks.",basename(outfileALL),".p",pvalthresh,".sep",minsep,".chr",chrs[1],".bed"),header=F)
for(i in 2:length(chrs)){
	chrpeaks = read.table(paste0(datapath,"bychr/peaks/SingleBasePeaks.",basename(outfileALL),".p",pvalthresh,".sep",minsep,".chr",chrs[i],".bed"),header=F)
	allpeaks=rbind(allpeaks,chrpeaks)
}
cnames=c("chr","center_start","center_stop","CI_start","CI_stop","cov_r1","cov_r2","cov_input","enrichment","likelihood","pvalue")
write.table(allpeaks,file=outfileALL,quote=F,row.names=F,col.names=cnames,sep="\t")


print(paste("done!:",format(Sys.time(), "%Y.%m.%d %H:%M:%S")))
quit(save="no",runLast=FALSE)
