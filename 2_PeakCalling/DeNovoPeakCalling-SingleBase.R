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
options("scipen"=8)
btpath = "bedtools"
batchsize=1000000 #how many bases to process in each batch of memory?
cithresh=0.99 #confidence interval threshold for output
searchwidth = 1000 #maximum CI width
genomesizefile = "hg19.chromsizes.tbl" #default path to chr size file


###set these all to 1 if preprocessing steps are needed, and to 0 if this has already been done for a given sample
#this allows for rapid re-processing if e.g. you would like to re-call peaks at a different p-value threshold later
#if these are set to 0 and corresponding temp files are not found, code will crash with an error
computecoverage="1,1,1" #preprocess read files to produce single-base coverage for r1,r2,genomic (0 or 1 for each)
createbin=1 #process coverage files into binary format for easy reading
computelikelihoods=1 #compute likelihood and enrichment values at each base

args=commandArgs(TRUE)
datapath = args[1] #path to folder containing fragment position bed files
sample = args[2] #a name for the sample
rep1suffix = args[3] #filename for ChIP replicate 1 fragment position bed file
rep2suffix = args[4] #filename for ChIP replicate 2 fragment position bed file
genomicsuffix = args[5] #filename for total chromatin input sample fragment position bed file
constfile = args[6] #path to file with constant estimates (output from EstimateConstants.R)
pvalthresh = as.numeric(args[7]) #maximum p-value at peak centre
minsep = as.numeric(args[8]) #minimum separation between peak centres
computecoverage = as.character(args[9]) #see above
createbin = as.integer(args[10]) #see above
computelikelihoods = as.integer(args[11]) #see above

#example hardwired input
#datapath="../1_Preprocessing"
#sample = "YFP_HumanPRDM9.antiGFP"
#rep1suffix = "FragPos.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.bed.sorted.bed"
#rep2suffix = "FragPos.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.bed.sorted.bed"
#genomicsuffix = "FragPos.YFP_HumanPRDM9.input.ProtocolN.bed.sorted.bed"
#constfile="Constants.YFP_HumanPRDM9.antiGFP.txt"
#pvalthresh = 0.000001
#minsep = 250
#computecoverage="1,1,1"
#createbin=1
#computelikelihoods=1


#create vector of all chromosome names at which to call peaks (change if necessary)
chrs=seq(1,22,1) #change if chromosome number differs
chrs=c(chrs,"X")

computecoverages=as.integer(strsplit(computecoverage,',')[[1]])
lthresh = qchisq(pvalthresh,df=1,lower.tail=F)
cilhood=qchisq(cithresh,df=2)

chrlengths = read.table(genomesizefile,header=F,colClasses=c('character','integer'))


#read in constants estimated previously by EstimateConstants.R, calculate constant terms used in likelihood calculations

constdata=read.table(constfile,header=TRUE)
alpha1=constdata[which(constdata[,1]=="autosomal"),2]
alpha2=constdata[which(constdata[,1]=="autosomal"),3]
beta=constdata[which(constdata[,1]=="autosomal"),4]

term2=1+alpha1+alpha2
term3=beta+1
term6=alpha1*alpha2
term7=alpha1*beta+alpha2

#declare output filenames
outfileALL = paste0("SingleBasePeaks.",sample,".p",pvalthresh,".sep",minsep,".ALL.bed")
system("mkdir tmp", ignore.stdout = T, ignore.stderr = T)
system("mkdir bychr", ignore.stdout = T, ignore.stderr = T)
system("mkdir bychr/covzip", ignore.stdout = T, ignore.stderr = T)
system("mkdir bychr/covbin", ignore.stdout = T, ignore.stderr = T)
system("mkdir bychr/likelihoods", ignore.stdout = T, ignore.stderr = T)
system("mkdir bychr/enrichments", ignore.stdout = T, ignore.stderr = T)
system("mkdir bychr/peaks", ignore.stdout = T, ignore.stderr = T)

singleBaseCoverageBG <- function(chr, bedgraph, path, covfile){

  # expand bedgraph to single base pair resolution for a single chr

  command <- paste0("grep -P 'chr",chr,"\\t' ",path,bedgraph," | \\
                            perl -lane 'for ($F[1]+1..$F[2]) { print\"$F[0]\\t$_\\t$F[3]\" }' | \\
                            cut -f3 | \\
                            gzip > ",path,"covzip/",covfile)
  system(command)
}

singleBaseCoverageFP <- function(chr, FragPosFile, path, covfile, genomesizefile){

  # 1) extract single chr from fragpos file
  # 2) compute bedgraph
  # 3) extract single chr
  # 4) epand to single base pair resolution

  command <- paste0("grep -P 'chr",chr,"\\t' ",path,FragPosFile," | \\
                            ",btpath," genomecov -bga -i stdin -g ",genomesizefile," | \\
                            grep -P 'chr",chr,"\\t' | \\
                            perl -lane 'for ($F[1]+1..$F[2]) { print \"$F[0]\\t$_\\t$F[3]\" }' | \\
                            cut -f3 | \\
                            gzip > ",path,"covzip/",covfile)
  system(command)
}

#declare a function to perform calculations on one chromosome
getEnrichments=function(chr){

	#define input and output file names

  posfileA = paste0(datapath,rep1suffix)
  posfileB = paste0(datapath,rep2suffix)
  posfileG = paste0(datapath,genomicsuffix)

	covfileA = paste0(datapath,"bychr/covzip/",rep1suffix,".FragDepth.chr",chr,".bed.gz")
	covfileB = paste0(datapath,"bychr/covzip/",rep2suffix,".FragDepth.chr",chr,".bed.gz")
	covfileG = paste0(datapath,"bychr/covzip/",genomicsuffix,".FragDepth.chr",chr,".bed.gz")

	covfileAb = paste0(datapath,"bychr/covbin/",rep1suffix,".FragDepth.chr",chr,".binary.gz")
	covfileBb = paste0(datapath,"bychr/covbin/",rep2suffix,".FragDepth.chr",chr,".binary.gz")
	covfileGb = paste0(datapath,"bychr/covbin/",genomicsuffix,".FragDepth.chr",chr,".binary.gz")

	outfileLhood= paste0(datapath,"bychr/likelihoods/SingleBaseLikelihood.",sample,".chr",chr,".binary.r")
	outfileEnrich= paste0(datapath,"bychr/enrichments/SingleBaseEnrichment.",sample,".chr",chr,".binary.r")
	outfilePeaks= paste0(datapath,"bychr/peaks/SingleBasePeaks.",sample,".p",pvalthresh,".sep",minsep,".chr",chr,".bed")

	chrlen = chrlengths[which(chrlengths[,1]==paste0("chr",chr)),2]


	#if not done already, compute single-base coverage values across chromosome and compress
	if(computecoverages[1]==1){
		#singleBaseCoverageBG() is ~2x faster, but requires bedgraph, not fragpos file
	  singleBaseCoverageFP(chr, posfileA, datapath, covfileA, genomesizefile)
	}
	if(computecoverages[2]==1){
	  singleBaseCoverageFP(chr, posfileB, datapath, covfileB, genomesizefile)
	}
	if(computecoverages[3]==1){
	  singleBaseCoverageFP(chr, posfileG, datapath, covfileG, genomesizefile)
	}

	if(createbin==1){

		conA=gzfile(covfileA,open="r")
		conAb=gzfile(covfileAb,open="wb")
		writeBin(scan(conA,what=integer(1),nlines=chrlen,quiet=TRUE),conAb)
		close(conA)
		close(conAb)

		conB=gzfile(covfileB,open="r")
		conBb=gzfile(covfileBb,open="wb")
		writeBin(scan(conB,what=integer(1),nlines=chrlen,quiet=TRUE),conBb)
		close(conB)
		close(conBb)

		conG=gzfile(covfileG,open="r")
		conGb=gzfile(covfileGb,open="wb")
		writeBin(scan(conG,what=integer(1),nlines=chrlen,quiet=TRUE),conGb)
		close(conG)
		close(conGb)

		print(paste("done converting coverage files to binary",chr,date()))

	}


	#if not done already, compute single-base likelihood and enrichment values across chromosome, in batches
	if(computelikelihoods==1){


		conAb=gzfile(covfileAb,open="rb")
		conBb=gzfile(covfileBb,open="rb")
		conGb=gzfile(covfileGb,open="rb")
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
				basepeaks=makepeaks.singlebase(r1=covA,r2=covB,g=covG)
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
		print(paste("done computing likelihoods",chr,date()))
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
	conAb=gzfile(covfileAb,open="rb")
	covA=readBin(conAb,what=integer(1),n=chrlen)
	confints=cbind(confints,covA[confints[,1]])
	rm(covA)
	close(conAb)

	conBb=gzfile(covfileBb,open="rb")
	covB=readBin(conBb,what=integer(1),n=chrlen)
	confints=cbind(confints,covB[confints[,1]])
	rm(covB)
	close(conBb)

	conGb=gzfile(covfileGb,open="rb")
	covG=readBin(conGb,what=integer(1),n=chrlen)
	confints=cbind(confints,covG[confints[,1]])
	rm(covG)
	close(conGb)

	con2=gzfile(outfileEnrich,open="rb")
	enrich=readBin(con2,what=numeric(1),n=chrlen)
	close(con2)
	confints=cbind(confints,enrich[confints[,1]])
	rm(enrich)

	signif=pchisq(confints[,2],df=1,lower.tail=F)
	confints=cbind(confints,signif)

	#print final peaks to text file
	sub=confints[confints[,5]==1,]
	write.table(cbind(paste("chr",chr,sep=""),sub[,1]-1,sub[,1],sub[,3]-1,sub[,c(4,6,7,8,9,2,10)]),file=outfilePeaks,quote=F,row.names=F,col.names=F,sep="\t")
	print(paste("done getting peak intervals",chr,date()))
	return(1)
}



#declare function used for calculating likelihoods and enrichments
makepeaks.singlebase=function(r1,r2,g){
	g[g==0]=0.5
	sumcov = r1+r2+g

	term1=sumcov*term3
	term4=beta*alpha1*r2+alpha2*r1
	term5=beta*(r1+r2)

	bterm=term1*term7-term2*term5-term3*term4
	cterm=term1*term6-term2*term4
	rm(term4)
	aterm=term1*beta-term3*term5
	rm(term5,term1)

	yvals=(-bterm+sqrt(bterm^2-4*aterm*cterm))/2/aterm
	rm(aterm,bterm,cterm)
	yvals[yvals<0]=0
	yvals[yvals>1e9]=1e9

	bvals=sumcov/(term2+(term3*yvals))

	bvalsnull=sumcov/term2
	yvalsnull=rep(0,length(bvalsnull))

	lhooddiff=2*( (sumcov*(log(bvals)-1)+r1*log(alpha1+yvals)+r2*log(alpha2+yvals*beta)) - (sumcov*(log(bvalsnull)-1)+r1*log(alpha1+yvalsnull)+r2*log(alpha2+yvalsnull*beta)) )
	rm(bvals, bvalsnull, yvalsnull)

	lhooddiff[is.na(lhooddiff)]=0
	lhooddiff[sumcov==0.5]=0
	yvals[sumcov==0.5]=0
	rm(sumcov)

	return(list(lhooddiff,yvals))
}


#call function on all chromosomes in parallel
print(date())
funfunc = mclapply(chrs,getEnrichments,mc.preschedule=TRUE,mc.cores=length(chrs))

allpeaks = read.table(paste("bychr/peaks/SingleBasePeaks.",sample,".p",pvalthresh,".sep",minsep,".chr",chrs[1],".bed",sep=""),header=F)
for(i in 2:length(chrs)){
	chrpeaks = read.table(paste("bychr/peaks/SingleBasePeaks.",sample,".p",pvalthresh,".sep",minsep,".chr",chrs[i],".bed",sep=""),header=F)
	allpeaks=rbind(allpeaks,chrpeaks)
}
cnames=c("chr","center_start","center_stop","CI_start","CI_stop","cov_r1","cov_r2","cov_input","enrichment","likelihood","pvalue")
write.table(allpeaks,file=outfileALL,quote=F,row.names=F,col.names=cnames,sep="\t")


print(paste("done!:",date()))
quit(save="no",runLast=FALSE)
