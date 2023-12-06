#This is a tweaked version of the p.adjust function from the stats package, to which all the credit belongs. The change involves being able to use n which is not fixed to the length of the p vector to allow for setting this based on a PCA analysis for example.
p.adjust2 <- function (p, method = p.adjust.methods, n = length(p)) 
{
  method <- match.arg(method)
  if (method == "fdr") 
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  else p <- p[nna]
  lp <- length(p)
  #stopifnot(n >= lp) #+++JZ: Now we can run with n < lp !!!! Thank you Stack Overflow!
  if (n <= 1) 
    return(p0)
  if (n == 2 && method == "hommel") 
    method <- "hochberg"
  p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
    i <- seq_len(lp)
    o <- order(p)
    ro <- order(o)
    pmin(1, cummax((n + 1L - i) * p[o]))[ro]
  }, hommel = {
    if (n > lp) p <- c(p, rep.int(1, n - lp))
    i <- seq_len(n)
    o <- order(p)
    p <- p[o]
    ro <- order(o)
    q <- pa <- rep.int(min(n * p/i), n)
    for (j in (n - 1L):2L) {
      ij <- seq_len(n - j + 1L)
      i2 <- (n - j + 2L):n
      q1 <- min(j * p[i2]/(2L:j))
      q[ij] <- pmin(j * p[ij], q1)
      q[i2] <- q[n - j + 1L]
      pa <- pmax(pa, q)
    }
    pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
  }, hochberg = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin((n + 1L - i) * p[o]))[ro]
  }, BH = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin(n/i * p[o]))[ro]
  }, BY = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- sum(1/(1L:n))
    pmin(1, cummin(q * n/i * p[o]))[ro]
  }, none = p)
  
  p0 <- ifelse(p0 < p, p, p0) #+++JZ: Added this to avoid any p0 < p when using small n
  
  p0
}

getEffectiveNumberOfTests <- function(
    covarianceMatrix,
    symmetric = T,
    eigenSumLimit = 0.995
){
  pcaRes <- eigen(mvldsc$V_Stand,symmetric = symmetric)
  eigenSum<-sum(abs(pcaRes$values))
  for(iEV in 1:length(pcaRes$values)){
    if(sum(abs(pcaRes$values[1:iEV]))/eigenSum > eigenSumLimit) break
  }
  
  return(iEV)
}
# values1 = diag(p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S)[p$sumstats.sel.ssimp.code]
# standard_errors1 = diag(p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S.SE)[p$sumstats.sel.ssimp.code]
# values2 = p$sumstats$h2.LDSC.1kg.maf0_01[p$sumstats$code %in% p$sumstats.sel.ssimp.code]
# standard_errors2 = p$sumstats$h2.se.LDSC.1kg.maf0_01[p$sumstats$code %in% p$sumstats.sel.ssimp.code]
# n1 <- diag(p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$cov.blocks)[p$sumstats.sel.ssimp.code]
# n2 <- rep(200,length(p$sumstats.sel.ssimp.code))

#uses student's t test (rather than welch's); assumes ~equal variances (even though we clearly rely on them not being exactly equal)
#this test assumes per default that the variables are closely related/ NOT INDEPENDENT and have correlation ~ 1
test.meta.deltacovariance <- function(
    values1,
    standard_errors1,
    values2,
    standard_errors2,
    n1,
    n2,
    effectiveNumberOfTests=length(values1),
    fullyDependent=T #assumes practically fully dependent variables
    ){
  #test covariance difference
  mTest.values<-abs(values1-values2) #abs m because we want a one-sided test
  #old test was based on a normal
  if(fullyDependent){ # fully dependent samples
    sdeFullyDependentVariables<-sqrt(2*abs((n1-1)*standard_errors1^2-(n2-1)*standard_errors2^2)/(n1+n2-2)) #the intersect, from stochastic var. variance formulas
    pTest.values<-pt(
      q = ( mTest.values/(sdeFullyDependentVariables*sqrt((1/n1) + (1/n2)))),
      df = (n1+n2-2),
      lower.tail = F
    ) 
  } else { #independent samples
    pTest.values<-pt(
      q = ( mTest.values/(sqrt((n1-1)*standard_errors1^2+(n2-1)*standard_errors2^2/(n1+n2-2))*sqrt((1/n1) + (1/n2)))),
      df = (n1+n2-2),
      lower.tail = F
    )
  }
 
  pTest.values.adj<-p.adjust2(pTest.values, method = "fdr", n = effectiveNumberOfTests)

  #test standard error difference - assuming difference between two chi2
  t1<-n1*(n1-1)*(standard_errors1^2)/abs(values1) #we use the mean original variable in the denominator as in Pearsson's Chi2 test
  t2<-n2*(n2-1)*(standard_errors2^2)/abs(values2)
  
  #these are approximately normal with N(n-1,2(n-1)) with the same mean and variance as the Chi2
  
  #we need to do a double-sided/two-tailed test here, because of normal!!!
  mTest.standard_errors<-abs(t1 - t2)
  if(fullyDependent){
    sdeFullyDependentVariables<-sqrt(2*abs(2*(n1-1)-2*(n2-1))) #the intersect, again, from stochastic var. variance formulas
    pTest.standard_errors<-ifelse(mTest.standard_errors > abs(n1-n2),
                                  2*pnorm(q = mTest.standard_errors, mean =  abs(n1-n2), sd = sdeFullyDependentVariables,lower.tail = F),
                                  2*pnorm(q = mTest.standard_errors, mean =  abs(n1-n2), sd = sdeFullyDependentVariables,lower.tail = T))
  } else {
    sdeIndependentVariables <- sqrt(2*(n1-1)+2*(n2-1))
    pTest.standard_errors<-ifelse(mTest.standard_errors > abs(n1-n2),
                                  2*pnorm(q = mTest.standard_errors, mean =  abs(n1-n2), sd = sdeIndependentVariables, lower.tail = F),
                                  2*pnorm(q = mTest.standard_errors, mean =  abs(n1-n2), sd = sdeIndependentVariables, lower.tail = T))
  }
  
  pTest.standard_errors.adj<-p.adjust2(pTest.standard_errors, method = "fdr", n = effectiveNumberOfTests)

  return(
    list(
      pTest.values=pTest.values,
      pTest.values.adj=pTest.values.adj,
      pTest.standard_errors=pTest.standard_errors,
      pTest.standard_errors.adj=pTest.standard_errors.adj
      )
    )
  
}


# test.meta.deltacovariance for matrices
# mValues1 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.rblock.winfo.altcw$S
# mStandard_errors1 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.rblock.winfo.altcw$S.SE
# mValues2 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock$S
# mStandard_errors2 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock$S.SE
# mN1 <- p$mvLD$covstruct.mvLDSC.1kg.vbcs.rblock.winfo.altcw$cov.blocks
# mN2 <- p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock$cov.blocks
test.meta.deltacovariance.matrix <- function(
    mValues1,
    mStandard_errors1,
    mValues2,
    mStandard_errors2,
    mN1,
    mN2,
    symmetric=T,
    effectiveNumberOfTests=NULL,
    fullyDependent=T, #assumes practically fully dependent variables
    mValueCovariances=NULL,
    eigenSumLimit = 0.995
){
  
  if(!is.null(mValueCovariances)){
    effectiveNumberOfTests<-getEffectiveNumberOfTests(covarianceMatrix = mValueCovariances,symmetric = symmetric, eigenSumLimit = eigenSumLimit)
  }
  
  if(!symmetric){
    t<-test.meta.deltacovariance(
      values1 = mValues1,standard_errors1 = mStandard_errors1,values2 = mValues2,standard_errors2 = mStandard_errors2,n1 = mN1,n2 = mN2,effectiveNumberOfTests = effectiveNumberOfTests,fullyDependent = fullyDependent)
    mValues1
  }
  
  
}
