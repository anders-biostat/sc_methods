library(pheatmap)



# Plausible means and sizefactors -----------------------------------------

# Find more detail in hvg.Rmd. Basically, below means/sfs are highly similar in
# magnitude and proportions to those observed in the CITEseq dataset:
plausible_means <- c(  
    runif(575, min = .5, max = 1),
    runif(267, min = 1,  max = 2),
    runif(190, min = 2,  max = 10),
    runif(70, min = 10, max = 40)  )
pm <- plausible_means   # for code

plausible_sf <- rf(1000,20,22, 0)


# simulate function -------------------------------------------------------

simulateFoldchange <- function(
                        mean_low  = 1,       # mean_low, mean_high: Poisson rate for
                        mean_high = 4,       #                      low/high expressors
                        n_low = 100,         # n_low, n_high: number of cells
                        n_high = 50,
                        sizefactors=NULL    # optional: provide arbitrary number of sizefactors from which to sample
                        ) {
   # Return value: vector with simulated counts. If sizefactors
   #               are provided, returns normalized counts.
   # function summary:
   # Simulate poisson counts for two populations of sizes n_low and n_high, with
   # means (Poisson rate lambda) mean_low and mean_high, respectively.
   # If sizefactors are provided, cells with a larger sf see a larger poisson rate
   # and are then later normalized again by this size factor. So small cells will
   # get more zeros overall, while large cels (large sf) will have small values instead
   # of zeros.
    if( is.null(sizefactors) ) {
      sizefactors <- sample(1, n_low+n_high, replace = TRUE)
    }else{ # also works if less sizefactors than cells due to `replace=T`:
      sizefactors <- sample(sizefactors, n_low+n_high, replace = TRUE)
    }
     
     c( rpois(n_low,  mean_low * sizefactors[1:n_low]),
        rpois(n_high, mean_high *sizefactors[ (n_low+1) : (n_low+n_high) ] )
       ) / sizefactors 
   
}






# Property of Variance (Wikipedia) ------------------------------
# If all values are scaled by a constant, the variance is scaled by 
# the square of that constant. In math:  var(c*X) <- c^2 * var(X)

p <- cbind(
p15 = rpois(1000, lambda = 15),
p5  = rpois(1000, lambda = 5),
p.1 = rpois(1000, lambda = .1)
)


apply(t(t(p) * c(1,2,3)), 2, var) # Poisson variance
apply(t(t(scale(p))), 2, var)     # scaled to unit variance
apply(t(t(scale(p)) * c(1,2,3)), 2, var)              # careful: var(c*X) <- c^2 * var(X)
apply(t(t(scale(p)) * c(1,sqrt(2),sqrt(3))), 2, var)  # inflate variance as I like!


# Understand PCA  -------------------------------------------------


# How PCA works:
counts <- cbind(fc2 = simulateFoldchange(2,4, n_low=500, n_high = 500),
                fc2high = simulateFoldchange(8, 16, n_low=500, n_high = 500),
                fc0 = simulateFoldchange(4,4, n_low=500, n_high = 500))
vmr <- function(vector) var(vector) / mean(vector)
apply(counts, 2, vmr)

# Observe how not scaling gives higher weight to the gene with the larger variance.
round(cbind(var =apply(counts, 2, var) ,prcomp(counts)$rotation), 2)
round(cbind(var = apply(scale(counts), 2, var), prcomp(scale(counts))$rotation), 2)
# We know variance is a stupid measure for a gene's importance in scRNAseq data,
# simply because the variance varies quite a lot over the range of expression magnitudes.
# Instead of stabilizing the variance, I will try to weight each gene by a better
# number; namely the VMR.

# Here is how one could weight one gene much stronger than all others:
rescaled <- t(  t(scale(counts)) * c(sqrt(15), 1, 1)  )
round(cbind(var = apply(rescaled, 2, var), prcomp(rescaled)$rotation), 2)

# We can see weighting by something other than the variance might be a good idea
# when we consider the following. Taking the log does
# somewhat stabilize variance, but in this case overshoots and gives too much
# weight to the lowly expressed gene:
round(cbind(var =apply(log1p(counts), 2, var) ,prcomp(log1p(counts))$rotation), 2)





# custom functions  ------------------------------------------------------
# useful to understand what's going on with top genes during PCA and weighted PCA.


ct_boxplot <- function(cellInfo = factor(sort(rep(1:2, 250))),
                       pca, pcs=1:6, plot_cols = 3, plot_rows=2) {
     par(mfrow = c(plot_rows, plot_cols))
    sapply(pcs, function(i) {
      scores <- pca$x[, i]
      plot(cellInfo, scores, main = paste0("PC", i))
    }); par(mfrow = c(1, 1))
}

markGeneRotations <- function(genesOfInterest.,
                              loadings. = pca$rotation[, 1],
                              n_pcTop. = 100,
                              interest_col. = "red",
                              ...) {
  topPCgenes <-head(loadings.[order(abs(loadings.), decreasing = T)], n_pcTop.)
  positionsOfInterest <- which( names(topPCgenes) %in% genesOfInterest.)
  set.seed(100)
  x <- rnorm(length(topPCgenes))
  plot(x, topPCgenes, pch = 20,
       ylab = paste0("Gene Loading  (Top ", n_pcTop., " genes", ")"),
       xlab = "jitter",
       ...)
  points(x[positionsOfInterest], topPCgenes[positionsOfInterest], col = interest_col., pch =20)
}

# Wrapper for markGeneRotations:
rotationPlot <- function(genesOfInterest,
                         pca,
                         pcs=1:6,
                         title = NULL,
                         n_pcTop=100,
                         interest_col = "red",
                         plot_cols=3,
                         plot_rows=2) {
    par(mfrow = c(plot_rows,plot_cols), oma = c(0, 0, 2, 0))
    sapply(pcs, function(i) {
      markGeneRotations(genesOfInterest. = genesOfInterest,
                        loadings. = pca$rotation[, i], n_pcTop. = n_pcTop ,
                        interest_col. =  interest_col,
                        main = paste0("PC", i))
    })
    mtext(title, outer = TRUE, cex = 1.5); par(mfrow = c(1,1), oma = c(0,0,0,0))
}






# Weight PCA by VMR -------------------------------------------------------

#
#    load CITEseq T cells (250 CD4, 250 CD8 T cells) using script
#    hvg.Rmd

# get data (after executing some sections in hvg.Rmd):
ct <- factor( c( rep("CD4", 250), rep("CD8", sum(isCD8Tcell)) ), levels = c("CD4", "CD8"))
vs <- VMRstats(rawC[, selecteven])
normeven <- log( t(t(rawC[, selecteven]) * median(colSums(rawC[, selecteven])) /
                  colSums(rawC[, selecteven])  ) + 1 ) 
normeven <- normeven[ rownames(vs), ]



# normal PCA (lognormalization, no scaling):
logPCA <- prcomp(t(normeven))

ct_boxplot(cellInfo = ct, pca = logPCA, 1:6)

rotationPlot(rownames(vs)[vs$aboveP > 5],
             logPCA,
             title = "Normal PCA on lognormalized counts\naboveP > 5")
rotationPlot(rownames(vs)[vs$aboveP < .5],
             logPCA,interest_col = "green",
             title = "Normal PCA on lognormalized counts\nNoisy genes (aboveP < .5")

logTSNE <- Rtsne::Rtsne(logPCA$x[, 1:6])
plot(logTSNE$Y, col = ct)


# weighted PCA (scale genes to unit variance, then inflate ~ `aboveP`)



unit_counts     <- t( scale(t(normeven[ rownames(vs)[vs$aboveP>0], ])))
weighted_counts <- unit_counts * sqrt(vs$aboveP[vs$aboveP>0])


vmrPCA <- prcomp(t(weighted_counts))

ct_boxplot(cellInfo = ct, pca = vmrPCA, 1:6)
#
#  We see that all PCs have one or two outlier cells that are very far away from
#  all others. This is because using 'aboveP' as weights overemphasizes a few
#  genes; the cell expressing such a gene the highest will be an extreme outlier.
#
# To do:   try empirical bayes to shrink 'aboveP' to something more moderate!





  
  

# Pearson correlation under Poisson noise ---------------------------------------
c1 <- cor(
     sapply(pm, function(mu) {
     simulateFoldchange(mean_low = mu, mean_high = 2 * mu, sizefactors = plausible_sf)
   })
)

noise1 <- cor(
  sapply(pm, function(mu) {
    simulateFoldchange(mean_low = mu, mean_high = mu, sizefactors = plausible_sf)
  })
)


rownames(c1) <- 1:length(pm)
pheatmap(c1, cluster_cols = F, cluster_rows = F, annotation_row = data.frame(GeneExpression = pm, row.names = 1:length(pm)))
# We can see the high correlation is observed best between highly expressed
# genes. The lower both genes are expressed, the less clear their (perfect) correlation
# is measured by Pearson coefficient.


logaverage_matrix <- outer(log(pm), log(pm), FUN = "+")
plot(
  logaverage_matrix,
  c1,
  pch = 20, cex= .2,
  main = "Poisson-simulated Gene-Gene correlations\nblack: with biological signal (fc=2, 100 low, 50 high cells)\norange: no biological signal",
  ylab = "Pearson Correlation",
  xlab = "log( GeneA_average ) + log( GeneB_average )"
)
points(
  logaverage_matrix,
  noise1,
  pch = 20, cex = .2, col = "orange"
)
# If at least one of the genes is expressed highly, correlations (black) lie above
# Poisson noise (orange).


















 