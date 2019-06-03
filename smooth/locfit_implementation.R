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

l.001 <- manloc_smooth(
           umis = x$raw_umis["Top2a", ],
       totalUMI = colSums(x$raw_umis),
  featureMatrix = x$PCA_embeddings,
  leave_one_out = TRUE,
             nn = 50,
              h = NULL,
  lambda_use = 1e-3)
range(l.001$smoothed)


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


plot(l.3$smoothed, l.3$weighted_mean * totalUMI, col = 1 + (l.3$status == "ok"), asp=1);abline(0,1)
plot(l.200$smoothed, l.200$weighted_mean * totalUMI, col = 1 + (l.200$status == "ok"), asp=1);abline(0,1)



# Top2a, 50 nn ------------------------------------------------------------



Top2a_50 <- lapply(10^(seq( log10(.001), log10(200), length.out = 9)),
       function(theLambda) {
         print(theLambda)
         res_df <- manloc_smooth(
           umis = x$raw_umis["Top2a", ],
           totalUMI = colSums(x$raw_umis),
           featureMatrix = x$PCA_embeddings,
           leave_one_out = TRUE,
           nn = 50,
           h = NULL,
           lambda_use = theLambda)
         res_df <- data.frame(lambda=theLambda,
                              res_df,
                              totalUMI = Matrix::colSums(x$raw_umis),
                              stringsAsFactors = F)
       })





# Top2a, 300 nn ------------------------------------------------------------



Top2a_300 <- lapply(10^(seq( log10(.001), log10(200), length.out = 9)),
       function(theLambda) {
         print(theLambda)
         res_df <- manloc_smooth(
           umis = x$raw_umis["Top2a", ],
           totalUMI = colSums(x$raw_umis),
           featureMatrix = x$PCA_embeddings,
           leave_one_out = TRUE,
           nn = 300,
           h = NULL,
           lambda_use = theLambda)
         res_df <- data.frame(lambda=theLambda,
                              res_df,
                              totalUMI = Matrix::colSums(x$raw_umis),
                              stringsAsFactors = F)
       })
tmp <- do.call(rbind, ll2)

tmp$lambda <- round(tmp$lambda, 3)
tmp$adjusted <- tmp$smoothed
tmp$adjusted[tmp$adjusted > 25] <- 25

# smoothing vs original
ggplot(tmp, aes(umis, adjusted))+geom_point()+geom_abline()+
  coord_fixed() + facet_wrap(~lambda)

# how much smoothing is different from weighted mean:
ggplot(tmp, aes(weighted_mean * totalUMI, adjusted))+geom_point()+geom_abline()+
  coord_fixed() + facet_wrap(~lambda)

# totalDeviance:
tmp %>% mutate( deviances = dev(umis, smoothed)) %>% group_by(lambda) %>% summarize(totalDeviance = sum(deviances, na.rm = T)) %>% ggplot(aes(lambda, totalDeviance))+geom_point() + scale_y_log10() + scale_x_log10()






# Hes5, 50 nn ------------------------------------------------------------



Hes5_50 <- lapply(10^(seq( log10(.001), log10(200), length.out = 9)),
             function(theLambda) {
               print(theLambda)
               res_df <- manloc_smooth(
                 umis = x$raw_umis["Hes5", ],
                 totalUMI = colSums(x$raw_umis),
                 featureMatrix = x$PCA_embeddings,
                 leave_one_out = TRUE,
                 nn = 50,
                 h = NULL,
                 lambda_use = theLambda)
               res_df <- data.frame(lambda=theLambda,
                                    res_df,
                                    totalUMI = Matrix::colSums(x$raw_umis),
                                    stringsAsFactors = F)
             })



# Hes5, 300 nn ------------------------------------------------------------



Hes5_300 <- lapply(10^(seq( log10(.001), log10(200), length.out = 9)),
             function(theLambda) {
               print(theLambda)
               res_df <- manloc_smooth(
                 umis = x$raw_umis["Hes5", ],
                 totalUMI = colSums(x$raw_umis),
                 featureMatrix = x$PCA_embeddings,
                 leave_one_out = TRUE,
                 nn = 300,
                 h = NULL,
                 lambda_use = theLambda)
               res_df <- data.frame(lambda=theLambda,
                                    res_df,
                                    totalUMI = Matrix::colSums(x$raw_umis),
                                    stringsAsFactors = F)
             })








# Plots -------------------------------------------------------------------
Hes5_50
Hes5_300
Top2a_50
Top2a_300



tmp <- do.call(rbind, Top2a_300)
tmp$lambda <- round(tmp$lambda, 3)
tmp$adjusted <- tmp$smoothed
# tmp$adjusted[tmp$adjusted > 60] <- 60

# smoothing vs original
ggplot(tmp, aes(umis, adjusted))+geom_point()+geom_abline()+
  coord_fixed() + facet_wrap(~lambda)

# how much smoothing is different from weighted mean:
ggplot(tmp, aes(weighted_mean * totalUMI, adjusted))+geom_point()+geom_abline()+
  coord_fixed() + facet_wrap(~lambda)

# totalDeviance:
tmp %>% mutate( deviances = dev(umis, smoothed)) %>% group_by(lambda) %>%
  summarize(totalDeviance = sum(deviances, na.rm = T)) %>%
  ggplot(aes(lambda, totalDeviance))+geom_point() + scale_y_log10() + scale_x_log10()



