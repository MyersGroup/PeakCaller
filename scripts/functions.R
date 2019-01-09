
#' Calculate Enrichment and Likelihood at a single base pair
#'
#' @param r1 integer; read count in replicate 1
#' @param r2 integer; read count in replicate 2
#' @param g integer; read count in input/control
#' @param alpha1 numeric;
#' @param alpha2 numeric;
#' @param beta numeric;
#' @param return_slim logical; if TRUE (default) will only return likelihood difference and yvals, if FALSE will also return
#'
#' @example
#'
#' compute_likelihood(r1=rpois(100,3), r2=rpois(100,3), g=rpois(100,3), alpha1=0.5, alpha2=0.5, beta=0.5, return_slim = F)
#'

compute_likelihood = function(r1, r2, g, alpha1, alpha2, beta, return_slim=TRUE){
  g[g==0]=0.5
  sumcov = r1+r2+g

  term3=beta+1

  term1=sumcov*term3
  term2=1+alpha1+alpha2

  term4=beta*alpha1*r2+alpha2*r1
  term5=beta*(r1+r2)
  term6=alpha1*alpha2
  #term7=alpha1*beta+alpha2 only used once, bracketed term next line

  bterm=term1*(alpha1*beta+alpha2)-term2*term5-term3*term4
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

  lhooddiff=2*( (sumcov*(log(bvals)-1)    +  r1*log(alpha1+yvals)   +  r2*log(alpha2+yvals*beta)) -
                (sumcov*(log(bvalsnull)-1) + r1*log(alpha1+yvalsnull)+ r2*log(alpha2+yvalsnull*beta)) )
  rm(bvalsnull, yvalsnull) #bvals req. when return_slim=F

  lhooddiff[is.na(lhooddiff)]=0
  lhooddiff[sumcov==0.5]=0
  yvals[sumcov==0.5]=0
  rm(sumcov)

  if(return_slim){

    return(list(lhooddiff,yvals))

  }else{

    signif=pchisq(lhooddiff,df=1,lower.tail=F)

    return(cbind("yhat_alt"=yvals,
                  "Lhood_diff"=lhooddiff,
                  "p-value"=signif,
                  "bhat_alt"=bvals))
  }

}





#' Calculate Fragment Coverage Overlap
#'
#' @param recalculate_coverage logical; if infiles already exist should we use them or recalculate the coverage, (default FALSE)
#' @param posfiles string vector; vector of filenames of fragment position bed files
#' @param infiles string vector; vector of filenames of output filenames for the coverage files
#' @param windowfilepath string; filename/path of bedfile containing intervals over which to calculate coverage
#' @param chromosomes string vector; vector of chromosome names to return counts for
#' @param bedtools string; full path of bedtools executable
#'

get_frag_overlap_counts <- function(posfiles = c(posfileA, posfileB, posfileG), infiles = c(infile1, infile2, infile3), windowfilepath=windowfile, chromosomes=chrs, bedtools=btpath, recalculate_coverage=FALSE){

  if(all(file.exists(infiles)) & !recalculate_coverage){
    print("Skipping recalculation of coverage")
  }else{
    print("Calculating coverage")

    calc_frag_overlap_counts <- function(i, posfilevec=posfiles, infilevec=infiles, windows=windowfilepath, bedtoolspath=bedtools){
      system(paste0(bedtoolspath," coverage -a ",windows," -b ",posfiles[i]," -counts | cut -f4 > ",infiles[i]))
    }

    noreturn = mclapply(1:3, calc_frag_overlap_counts, mc.preschedule=TRUE, mc.cores=3)

    # non parallel version of above
    # system(paste0(btpath," coverage -a ",windowfile," -b ",posfileA," -counts | cut -f4 >",infile1))
    # system(paste0(btpath," coverage -a ",windowfile," -b ",posfileB," -counts | cut -f4 >",infile2))
    # system(paste0(btpath," coverage -a ",windowfile," -b ",posfileG," -counts | cut -f4 >",infile3))
  }

  #read in and combine fragment coverage values into one dataframe "counts"
  print("Reading in fragment count coverage")
  counts = fread(windowfilepath, col.names=c('chr','start','stop'))
  counts$countA <- fread(infiles[1])
  counts$countB <- fread(infiles[2])
  counts$countG <- fread(infiles[3])

  setkey(counts, chr, start, stop)
  counts <- counts[chr %in% paste0("chr",chromosomes)]
  return(counts)

}
