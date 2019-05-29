library(statmod)
library( locfit )
library( Matrix )

# Locfit by hand ----------------------------------------------------------


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

  fit_statmod <- do.call(rbind, lapply(1:length(umis), function(i){ 
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
    
    
    # when all y are 0, statmod::glmnb.fit explicitly returns 0 for all coefficients
    # which creates a bug in combination with an offset, so we handle it here explicitly:
    # [ specifically, fitted value=exp(log(s) + sum( coef(fit) * data)) would give s, not 0. ]
    if( max(umis[w>0]) == 0 ) {
      fit_results <- data.frame(
         smoothed  = 0,
         status    = "ok",
         row.names=NULL)
    } else{
     fit_results <- tryCatch(
       {
        # fit <- statmod::glmnb.fit(
        #   X = model.matrix(~featureMatrix[w>0,]),
        #   y = umis[w>0],
        #   dispersion = 0, # poisson regression
        #   weights = w[w>0],
        #   offset = log( totalUMI[w>0] ) )
        # 
        # if( leave_one_out ) {# for leave-1-out crossvalidation
        #   data_at_i <- c(Intercept = 1, featureMatrix[i,]) 
        #   fit_at_i <- exp( log(totalUMI[i]) + sum( coef(fit) * data_at_i ) )
        # }else{ # if we are not doing crossvalidation:
        #   fit_at_i <- fitted.values( fit )[ sum( (w>0)[1:i] ) ]
        # }
         
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
        
        # if fit runs _without_ error/warning, extract results:
        data.frame(
         smoothed  = fit_at_i,
         status    = "ok",
         row.names=NULL)
        
       },
       error = function(e) {
         # print(paste0("error: ", i))
         return(data.frame(
         smoothed  = weighted_mean * totalUMI[i],
         status    = "error",
         row.names = NULL
         ))
          },
       warning= function(warn) {
         # print(paste0("warning: ", i))
         return(data.frame(
         smoothed  = weighted_mean * totalUMI[i],
         status    = "warning",
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
















# Carter data from SDS ----------------------------------------------------

x <- readRDS("~/sds/sd17l002/p/smooth_carter/processed/E13A_Cerebella_exon_output.rds")





# Top2a -------------------------------------------------------------------



# Whenever locfit 'does not converge' (warning message "procv: parameters out of bounds"),
# it returns the totalUMIs (base aka offset) as fitted values:
gene <- "Top2a"
gene_locfit <- locfit.raw( x = do.call(lp, c(as.list(as.data.frame(x$PCA_embeddings[, 1:15])),
                                             nn = .1,
                                             deg=1)),
                           y = x$raw_umis[gene, ],
                           family = "poisson",
                           link = "log",
                           maxit = 50,
                           base = log(colSums(x$raw_umis)), # take totalUMIs into account
                           ev = dat() # no extrapolation between cells
)

# outlier values (thousands of UMIs, specifically: the cell's totalUMIs):
range(fitted.values(gene_locfit))
all.equal(unname(colSums(x$raw_umis)[fitted.values(gene_locfit) > 100]),
          fitted.values(gene_locfit)[fitted.values(gene_locfit) > 100])





# manloc_smooth without CV returns values between 0 and ~20, because unplausible
# values (> 1000 UMIs) are produced by non-converging glm fits, which
# manloc_smooth replaces with the tricube-weighted mean:
manfit <- manloc_smooth(x$raw_umis["Top2a", ],
                        colSums(x$raw_umis),
                        x$PCA_embeddings,
                        nn = 50)
range(manfit$smoothed)

# using leave-one-out crossvalidation, we observe a very large number of extreme outliers
# again:

# lambda is .3
l.3 <- manloc_smooth(
           umis = x$raw_umis["Top2a", ],
       totalUMI = colSums(x$raw_umis),
  featureMatrix = x$PCA_embeddings,
  leave_one_out = TRUE,
             nn = 50,
              h = NULL,
  lambda_use = 10^-.5)
range(l.3$smoothed)
table(l.3$smoothed > 1000)
plot(l.3$umis, l.3$smoothed, asp=1); abline(0,1)

l.200 <- manloc_smooth(
  umis = x$raw_umis["Top2a", ],
  totalUMI = colSums(x$raw_umis),
  featureMatrix = x$PCA_embeddings,
  leave_one_out = TRUE,
  nn = 50,
  h = NULL,
  lambda_use = 10^2.3)
range(l.200$smoothed)
plot(l.200$umis, l.200$smoothed, asp=1); abline(0,1)






dev <- function(counts, lambdas){
 gives_inf <- lambdas == 0 & counts != 0
 if(sum(gives_inf)>0) {
 counts <- counts[!gives_inf]; lambdas <- lambdas[!gives_inf]
 cat(paste0("Removed observations: ", sum(gives_inf), "\n") ) 
 cat( "(These would otherwise give -Inf as they have nonzero count but lambda==0)\n")
 }
 -2 * ( dpois(counts, lambdas, log=T) - dpois(counts, counts, log=T)) 
}

sum(dev(l.200$umis, l.200$smoothed))


