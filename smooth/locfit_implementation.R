library(glmnet)
library( Matrix ) # recommended (necessary for sparse matrices)

# Locfit by hand ----------------------------------------------------------


dev <- function(counts, lambdas){
 # computes poisson deviance of counts from lambdas. Returns NA if count !=0 and lambda == 0.
 gives_inf <- lambdas == 0 & counts != 0
 deviances <- -2 * ( dpois(counts, lambdas, log=T) - dpois(counts, counts, log=T)) 
 if(sum(gives_inf)>0){
 deviances[gives_inf] <- NA
 cat(paste0("Replaced with NAs: ", sum(gives_inf), " expected infinite values (nonzero count but lambda==0).\n") ) 
 }
 return(deviances)
}


tricube <- function(x, halfwidth=1) {
  # tricube kernel function. All values outside c(-halfwidth, halfwidth) are 0,
  # integral is 1, variance is
  #     var(tricube) = 35/243 * halfwidth^2
  #
  # formula adapted from (wikipedia)[https://en.wikipedia.org/wiki/Kernel_(statistics)],
  # and expanded to contain `halfwidth`.
  tricube <- 70/81/halfwidth * (1 - abs(x/halfwidth)^3)^3
  # outside of kernel support (-halfwidth, +halfwidth), tricube is defined as 0:
  tricube[ abs(x) > halfwidth ] <- 0
  return(tricube)
}


manloc_smooth <- function(umis, totalUMI, featureMatrix, 
                          nn=NULL, h=NULL, leave_one_out=FALSE,
                          lambda_use = 10^-.5){ 
  # umis: raw umis for gene of interest. This should be an integer vector with 1
  #       element for each cell.
  # totalUMIs: totalUMIs of each cell (the colSums of gene-by-cell countmatrix)
  # featureMatrix: cells in rows, features in cols - typically first n PCs
  # nn: number of nearest neighbors (integer)
  
  
  # penalized regression is only meaningful on scaled dependent variable (very few exceptions).
  # glmnet does it by default but we opt to do it manually here, because 
  # coef(fit) will then return betas on the (interpretrable) scale where the
  # penalty was applied rather than reverse-standardizing them to fit the
  # dependent variable's original scale.
  # In any case, the distances 'ds' are computed on the _un_scaled features, as
  # high-variance features should count more for our application:
  ds <- as.matrix(dist( featureMatrix ))
  featureMatrix <- scale(featureMatrix)
  
  lambda_seq <- 10^(seq( log10(.001), log10(200), by = .1))

  fit_results <- do.call(rbind, lapply(1:length(umis), function(i){ 
    pos_distances <- ds[i, ]
   
   # compute weights for local regression - depending on how nn and/or h are set: 
   if(!is.null(nn)){ 
    hw_fit <- uniroot(
      f = function(root) { 
        # one_sd <- sqrt(35/243) * root
        sum( pos_distances < root) - nn},
      interval = c(1e-5, 20))}
    # weights for local regression:
    if(is.null(h)&!is.null(nn))  w <- tricube(pos_distances, halfwidth = hw_fit$root)
    if(!is.null(h)&is.null(nn))  w <- tricube(pos_distances, halfwidth = h)
    
    # leave-1-out crossvalidation:
    if( leave_one_out ) { w[i] <- 0}
    
    weighted_mean =  sum(umis*w/totalUMI) / sum(w) 
    
    
    # when all y are 0, we don't bother doing the fit and return 0 right away:
    # (in fact, glmnet throws warning otherwise)
    # cases (we reported it and by now it's fixed).
    if( max(umis[w>0]) == 0 ) {
      fit_results <- data.frame(
         smoothed  = 0,
         status    = "ok",
         intercept       = NA,
         coef_squaredSum = NA,
         coef_max        = NA,
         partial_linPred = NA,
         row.names=NULL)
    } else{
     fit_results <- tryCatch(
       {
        
        fit <- glmnet::glmnet(
                    x = featureMatrix[w>0, ],
                    y = umis[w>0],
               family = "poisson",
              weights = w[w>0],
               offset = log( totalUMI[w>0] ),
               lambda = lambda_seq,
                alpha = 0, # ridge regression
          standardize = FALSE  # featureMatrix was scaled already above
        )
        # glmnet returns linear predictor irrespective of whether type is "link"
        # or "response", so we have to exponentiate it to obtain the mean:
        fit_at_i <- exp( predict.glmnet(
          object    = fit,
          newx      = matrix(featureMatrix[i,], nrow=1),
          s         = lambda_use,
          newoffset = log(totalUMI[i]),
          type      = "link") )
        
        # convert from matrix to numeric:
        fit_at_i <- c(fit_at_i)
        
        # get coefficients:
        betas <- as.numeric(predict.glmnet(
          object    = fit,
          newx      = matrix(featureMatrix[i,], nrow=1),
          s         = lambda_use,
          newoffset = log(totalUMI[i]),
          type      = "coef") )
        
        # if fit runs _without_ error/warning, extract results:
        data.frame(
         smoothed        = fit_at_i,
         status          = "ok",
         intercept       = betas[1],
         coef_squaredSum = sum(betas[-1]^2),
         coef_max        = max(betas[-1][which.max(abs(betas[-1]))]),
         partial_linPred = sum(betas[-1] * featureMatrix[i,]),
         row.names=NULL)
        
       },
       error = function(e) {
         # print(paste0("error: ", i))
         return(data.frame(
         smoothed  = weighted_mean * totalUMI[i],
         status    = "error",
         intercept       = NA,
         coef_squaredSum = NA,        
         coef_max        = NA,
         partial_linPred = NA,
         row.names = NULL
         ))
          },
       warning= function(warn) {
         # print(paste0("warning: ", i))
         return(data.frame(
         smoothed  = weighted_mean * totalUMI[i],
         status    = "warning",
         intercept       = NA,
         coef_squaredSum = NA,        
         coef_max        = NA,
         partial_linPred = NA,
         row.names = NULL
         ))
         }
       
     )  # end of tryCatch
    }   # end of else
    
    # return value of cell i:
    data.frame(
         umis = umis[i],
         fit_results,
         weighted_mean = weighted_mean,
         neighbors = sum(w>0),
         nonzero_neighbors = sum(umis[w>0] != 0),
         stringsAsFactors = F )
    
  } ) )
  
}












# Other smoothing functions-----------------------------------------------------

topNN <- function(x, NN=20) { order(x)[1:NN] }

knn_smooth <- function(umis, totalUMI, featureMatrix) {
  # umis: raw umis for gene of interest. This should be a integer vector with 1
  #       element for each cell.
  # totalUMIs: the colSums of countmatrix
  # featureMatrix: cells in rows, features in cols - typically first n PCs
  
  ds <- as.matrix(dist( featureMatrix )) # cell-cell distances
  sapply(1:length(umis), function(i) {
    ids_knn <- topNN(ds[i, ], NN=20)
    mean(umis[ids_knn]/totalUMI[ids_knn]) })
}



# GAMs:
gam_form <- as.formula(
  paste0("gene ~ ",
         paste(paste0("s(PC", 1:10, ")"), collapse = " + "))
)



gam_smooth <- function(umis, totalUMI, featureMatrix) {
  # umis: raw umis for gene of interest. This should be a integer vector with 1
  #       element for each cell.
  # totalUMIs: the colSums of countmatrix
  # featureMatrix: cells in rows, features in cols - typically first n PCs
  fit <- gam(gam_form,
             data = data.frame(gene = umis, featureMatrix, totalUMI = totalUMI),
             offset = log( totalUMI ),
             family = "poisson")
  return( predict(fit, type = "response") )
}
# x <- gam_smooth(E13A$raw_umis["Hes5", ], E13A$totalUMI, E13A$PCA_embeddings)




