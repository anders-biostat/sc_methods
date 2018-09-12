# Negative Binomial Naive Bayesian Classifier


## Load packages
library(pbmcapply)




## Function: fitNB (R)
disp <- function(k) (var(k) - mean(k)) / mean(k) / mean(k) # var = mu + disp * mu^2
# remember to divide k by sf, and to prohibit negative values

fitNB <- function(k, sf, initialVals = c(mean(k), disp(k))) {   
  o = optim( 
    initialVals,  # mu, alpha
    function(x) {
      -sum( lgamma( k + 1/x[2] ) - lgamma( 1/x[2] ) - lgamma( k+1 ) - ( k + 1/x[2] ) * log( 1 + sf * x[2] * x[1] ) + k * log( sf * x[2] * x[1] ) )
    },
    
    function(x) c(
      -sum( ( k - sf * x[1] ) / ( x[1] + sf * x[2] * x[1]^2 ) ),
      -sum( ( x[2] * ( k - sf * x[1] ) / ( 1 + sf * x[2] * x[1] ) + log( 1 + sf * x[2] * x[1] ) - digamma( k + 1/x[2] ) + digamma( 1/x[2] ) ) / x[2]^2 ) ),
    hessian = TRUE,
    lower = c( 1e-10, 1e-10 ),
    method = "L-BFGS-B" )
  c( mean = o$par[1], disp = o$par[2], SE_m = 1 / sqrt(o$hessian[1,1]), SE_d = 1 / o$hessian[2,2] ) }


library(purrr) # for possibly
safe_fitNB <- possibly(fitNB, otherwise = c(mean=NA, disp=NA, SE_m = NA, SE_d = NA))



## Function: trainNB
trainNB <- function(countmatrix, isPositive, isNegative, sf = NULL) {
  ## FYI: faster for dense matrices and after excluding genes without a clear
  ##      expression signal to which we could fit (low-expression/dispersion).
  
  if(is.null(sf)) {
    sf <- colSums(countmatrix)
    names(sf) <- colnames(countmatrix)
  }
  
  dt <- pbmclapply(rownames(countmatrix), function(g) {
    if(sum(countmatrix[g, ]) == 0) return(
      data.frame(gene = NA, meanPos=NA, dispPos=NA, meanNeg=NA, dispNeg=NA))
    
    posCounts <- countmatrix[g, isPositive]
    negCounts <- countmatrix[g, isNegative]
    # if one group has 100% zeros this overemphasizes the genes, so remedy:
    if(sum(posCounts) == 0) posCounts <- c(1, rep(0, sum(isPositive)-1))
    if(sum(negCounts) == 0) negCounts <- c(1, rep(0, sum(isNegative)-1))
    
    pos =  safe_fitNB(posCounts, sf[isPositive])
    neg =  safe_fitNB(negCounts, sf[isNegative])
    data.frame(  gene = g, 
                 meanPos = pos["mean"], dispPos = pos["disp"],
                 meanNeg = neg["mean"], dispNeg = neg["disp"],
                 stringsAsFactors = F, row.names = NULL
    )
  })
  return(do.call(rbind, dt))
}
# Remember to remove NAs:
#    noNA <- colSums(apply(disptable, 1, is.na)) == 0






## Function: NBNB
NBNB <- function(countmatrix, dispersionTable, sf = NULL) {
  if(is.null(sf)) {
    sf <- colSums(countmatrix)
    names(sf) <- colnames(countmatrix)
  }
  sapply(colnames(countmatrix), function(cell) {
    
    sum(dnbinom(x  = countmatrix[dispersionTable$gene, cell],
                mu = dispersionTable$meanPos * sf[cell],
                size = 1 / dispersionTable$dispPos, log=T)) +
      -
      sum(dnbinom(x  = countmatrix[dispersionTable$gene, cell],
                  mu = dispersionTable$meanNeg * sf[cell],
                  size = 1 / dispersionTable$dispNeg, log=T))
  })
  
  
} # to do:
  #   - check dims of countmatrix and dispersiontable are the same -
  #     When excluding NAs, it can easily happen they're not
  #     (also genes should have same order). 
  #   - sf has to have names. Not sure if this is elegent, consider changing
  #     it but whatever you do this has to be documented!
  #
  #   - somehow make clear the size factors have to be the same as in trainNB.
  #     While this is obvious, some idiot WILL manage to mess it up (e.g.
  #     tired migrane felix on Friday, 10th August).



#  compute dispersion tables. This function is a wrapper for pbmclapply.
computeDT <- function(countmatrix) {
  sf <- colSums(countmatrix) / mean(colSums(countmatrix))
  dtList <- pbmclapply(rownames(countmatrix), function(g) {
  NBparams <- safe_fitNB(countmatrix[g, ], sf)
  data.frame(gene = g, t(NBparams))
})
  print("Converting list to matrix...")
  return(do.call(rbind,dtList))
}





# Example Usage -----------------------------------------------------------


# suppose we have raw scRNAseq counts (genes in rows, cells in cols) stored
# inside the object `rawC`. Then we can find super-poissonian genes like this
# (as of 8th August 18):

  #    disptable <- pbmclapply(rownames(rawC), function(g) {
  #      NBparams <- safe_fitNB(rawC[g, ], sf)
  #      data.frame(gene = g, t(NBparams))
  #    })  # typically runs about 20 min or so
  #    
  #    # even 1 NA is annoying and we kick it out:
  #    noNA <- colSums(apply(disptable, 1, is.na)) == 0
  #    # I pick superpoissonian genes as having dispersion clearly above 0:
  #    significant_dispersion <- (disptable$disp / disptable$SE_d) > 2
