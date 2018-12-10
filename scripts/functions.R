
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
