library( tidyverse )
library( Matrix )
library( locfit )
library( uwot )

colVars_spm <- function( spm ) {
  stopifnot( is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    mean <- sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}

rowVars_spm <- function( spm ) {
  colVars_spm( t(spm) )
}


# Load Carter data (pedestrian way) --------------------------------------------------------

dir <- "/home/felix/projects/Carter_Cerebellum/"


genes <- read_tsv(paste0(dir, "data/", "e13_A/genes.tsv"), col_names = c( "gene_id", "gene_name" ) )
# resolve duplicated symbols by concatenating ensemblID:
duplicate_gene_names <- genes$gene_name[ duplicated(genes$gene_name) ]
genes$gene_name <- ifelse(
  genes$gene_name %in% duplicate_gene_names,
  str_c( genes$gene_name, "_", genes$gene_id ),
  genes$gene_name )

barcodes <- readLines( paste0(dir, "data/", "e13_A/barcodes.tsv" ))
barcodes <- str_remove( barcodes, "-1" )

matrixTable <- read_delim( paste0(dir, "data/", "e13_A/matrix.mtx"), 
                           delim=" ", comment = "%", col_names = c( "gene", "barcode", "count" ) )

# first row contains number of rows|cols|totalEntries, 
# so it's removed from i/j/x and used as dims:
counts <-
  sparseMatrix( 
    i = matrixTable$gene[-1],
    j = matrixTable$barcode[-1], x = matrixTable$count[-1],
    dims = c( matrixTable$gene[1], matrixTable$barcode[1] ),
    dimnames = list( gene = genes$gene_name, barcode = barcodes ) )

rm( barcodes, genes, matrixTable )

# discard empty barcodes (kneeplot):
cell_barcodes <- names(which( colSums( counts>0 ) >= 1000 ))
counts <- counts[ , cell_barcodes ]

# discard empty genes (not detected in any cell):
counts <- counts[ rowSums(counts) > 0,  ]





# Normalization, PCA and UMAP -----------------------------------------------------------

norm_counts <- t(t(counts) / colSums(counts))


# informative genes:
gene_means <- rowMeans( norm_counts )
gene_vars <- rowVars_spm( norm_counts )
poisson_vmr <- mean( 1 / colSums( counts ) )

plot( gene_means, gene_vars / gene_means,
      log = "xy", cex = .3, col = adjustcolor("black", alpha=.3), 
      xlab = "mean", ylab = "variance / mean" )
abline(h = 2 * poisson_vmr, col = "red")

informative_genes <- names(which( 
  gene_vars / gene_means  >  2 * poisson_vmr ))

pcs_prcomp <- irlba::prcomp_irlba( t(sqrt(norm_counts[informative_genes,])),
                                   n = 20 )$x


umap_prcomp <- uwot::umap(pcs_prcomp,
                          n_neighbors = 30,
                          min_dist = .3,
                          metric = "cosine")





# select gene to fit:
gene <- "Cdkn1c"
# precomputations (outside loop, otherwise slow):
ds <- as.matrix(dist( pcs_prcomp[, 1:15] )) # cell-cell distances
y <- counts[gene, ]
cs <- colSums(counts)
x  <- pcs_prcomp[, 1:15]




# savepoint ---------------------------------------------------------------

# save.image(paste0("/home/felix/PhDother/scAnalysis/sc_methods/smooth/",
#                   "savepoint/",
#                   "locfit_by_hand_workspace.RData"))
load(paste0("/home/felix/PhDother/scAnalysis/sc_methods/smooth/",
                  "savepoint/",
                  "locfit_by_hand_workspace.RData"))
library( tidyverse )
library( Matrix )
library( locfit )
library( uwot )




# locfit package usage-----------------------------------------------------------



# using locfit package:
gene_locfit <- locfit.raw( x = do.call(lp, c(as.list(as.data.frame(pcs_prcomp[, 1:15])),
                                          nn = .1,
                                          deg=1)),
                        y = counts[gene, ],
                        family = "poisson",
                        link = "log",
                        maxit = 50,
                        base = log(colSums(counts)), # take totalUMIs into account
                        ev = dat() # no extrapolation between cells
)



gamma <-.5
data.frame(umap1 = umap_prcomp[,1],
           umap2 = umap_prcomp[,2],
           totalUMI = colSums(counts),
           umis  = counts[gene, ]) %>%
  mutate(locfit = fitted.values(gene_locfit)) %>%
  ggplot(aes(umap1, umap2, col = (locfit/totalUMI)^gamma)) +
  geom_point(size=1) + scale_color_gradientn(colours = rev(rje::cubeHelix(10))[2:10])










# locfit by hand ----------------------------------------------------------


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
# testing tricube function:
xseq <- seq(-10, 10, length.out = 500)
hw1 <- 2;plot(xseq, tricube(xseq, hw1), xlim = c(-10, 10)) # should scale with hw1
hw1 <- 2; integrate(tricube, halfwidth = hw1, lower = -20, upper = 20) # should be 1 for all hw1

rtricube <- function(n, halfwidth = 1) {
  # fun note: an improved implementation would use rejection sampling based on runif.
  possible_values <- seq(-halfwidth, halfwidth, length.out = 1e6)
  sample(possible_values,
         size = n,
         prob = tricube(possible_values, halfwidth = halfwidth))
}
# test rtricube:
hw1 <- 2; cat("runif var:  ", sd( rtricube(1000, hw1) )^2); cat("\ntheor. var:  ", 35/243 * hw1^2)















# h -----------------------------------------------------------------------

# loop:
cell_ids <- 1:ncol(counts)
manual.h2 <- sapply(cell_ids,function(i) { 
  pos_distances <- ds[i, ]
  fit <- glm(
    y ~ x,
    weights = tricube(pos_distances, halfwidth = .2),
    family = poisson("log"),
    offset = log( cs ) )
  return( fitted.values(fit)[i] )
  })
manual.h5 <- sapply(cell_ids,function(i) { 
  pos_distances <- ds[i, ]
  fit <- glm(
    y ~ x,
    weights = tricube(pos_distances, halfwidth = .5),
    family = poisson("log"),
    offset = log( cs ) )
  return( fitted.values(fit)[i] )
  })

# compute the same with locfit package:
package.h2 <- locfit.raw( x = do.call(lp, c(as.list(as.data.frame(pcs_prcomp[, 1:15])),
                              h = .2,
                              deg=1)),
            y = counts[gene, ],
            family = "poisson",
            link = "log",
            maxit = 50,
            base = log(colSums(counts)), # take totalUMIs into account
            ev = dat() # no extrapolation between cells
)
package.h5 <- locfit.raw( x = do.call(lp, c(as.list(as.data.frame(pcs_prcomp[, 1:15])),
                              h = .5,
                              deg=1)),
            y = counts[gene, ],
            family = "poisson",
            link = "log",
            maxit = 50,
            base = log(colSums(counts)), # take totalUMIs into account
            ev = dat() # no extrapolation between cells
)

# manual result is exactly the same as package output:
all.equal(fitted.values(package.h5), 
          unname(manual.h5))
all.equal(fitted.values(package.h2), 
          unname(manual.h2))
# sanity check: different bandwidths give different results:
plot(manual.h2, manual.h5, asp=1, main="Different bandwidths");abline(0,1)




# h: bug report -----------------------------------------------------------


# small h fail for some cells: 
pos_distances <- ds[109, ]
fit <- glm(
  y ~ x,
  weights = tricube(pos_distances, halfwidth = .1),
  family = poisson("log"),
  offset = log( cs ) )
# throws error "no valid set of coefficients has been found"



 




# nn ----------------------------------------------------------------------

# compute with package:
gene_loc_nn.1 <- locfit.raw( x = do.call(lp, c(as.list(as.data.frame(pcs_prcomp[, 1:15])),
                                            nn = .1,
                                            deg=1)),
                          y = counts[gene, ],
                          family = "poisson",
                          link = "log",
                          maxit = 50,
                          base = log(colSums(counts)), # take totalUMIs into account
                          ev = dat() # no extrapolation between cells
)

gene_loc_nn.7 <- locfit.raw( x = do.call(lp, c(as.list(as.data.frame(pcs_prcomp[, 1:15])),
                                            nn = .7,
                                            deg=1)),
                          y = counts[gene, ],
                          family = "poisson",
                          link = "log",
                          maxit = 50,
                          base = log(colSums(counts)), # take totalUMIs into account
                          ev = dat() # no extrapolation between cells
)

# compute manually:
man.1 <- sapply(1:ncol(counts), function(i) {
  pos_distances <- ds[i, ]
  hw_fit <- uniroot(
    f = function(root) { 
      # one_sd <- sqrt(35/243) * root
      sum( pos_distances < root) - ncol(counts) * .1},
    interval = c(1e-5, 20))
  fit <- glm(
    y ~ x,
    weights = tricube(pos_distances, halfwidth = hw_fit$root),
    family = poisson("log"),
    offset = log( cs ) )
  fitted.values(fit)[i]
})
man.7 <- sapply(1:ncol(counts), function(i) {
  pos_distances <- ds[i, ]
  hw_fit <- uniroot(
    f = function(root) { 
      # one_sd <- sqrt(35/243) * root
      sum( pos_distances < root) - ncol(counts) * .7},
    interval = c(1e-5, 20))
  fit <- glm(
    y ~ x,
    weights = tricube(pos_distances, halfwidth = hw_fit$root),
    family = poisson("log"),
    offset = log( cs ) )
  fitted.values(fit)[i]
})


# for the nn parameter, there are numerical differences between manual
# and package implementation:
all.equal(unname(man.1), fitted.values(gene_loc_nn.1))
all.equal(unname(man.7), fitted.values(gene_loc_nn.7))
# but these are not perceivable to the human eye:
plot(unname(man.1), fitted.values(gene_loc_nn.1), asp=1, pch=20);abline(0,1)

# This tiny but existent difference does __not__ go away when taking fewer
# cells, so it's not locfit's datapoint downsampling that causes it.
# I tried 100 and got similar relative differences (not shown here).








# number of PCs -----------------------------------------------------------
# Locfit package can NOT take more than 15 covariates and then throws a bug.
# Code provoking the bug:

#  15 is fine, 16 gives error:
locfit.raw( x = do.call(lp, c(as.list(as.data.frame(replicate(15, rnorm(100)))),
                              h = 5,
                              deg=1)),
            y = rnorm(100),
            ev=dat()
)
# Simon looked it up in locfit's source code, it probably is a bug:
# maxdim is a hard-coded constant set to 15, and throughout the code
# sometimes the constant is used, and sometimes a hard-coded 15.






# The manual implementation obviously compiles without errors:
x.16 <- pcs_prcomp[, 1:16]
pcs16 <- sapply(1:ncol(counts),function(i) { 
  pos_distances <- ds[i, ]
  fit <- glm(
    y ~ x.16,
    weights = tricube(pos_distances, halfwidth = .2),
    family = poisson("log"),
    offset = log( cs ) )
  return( fitted.values(fit)[i] )
})
# warnings:
#  32: glm.fit: algorithm stopped at boundary value
#  33: glm.fit: fitted rates numerically 0 occurred
#  34: step size truncated due to divergence










# speed benchmark ---------------------------------------------------------------------
ncells_halfwidth <- 50


glmfit_subsetting <- function() { 
  sapply(cell_ids,function(i) { 
    pos_distances <- ds[i, ]
    hw_fit <- uniroot(
      f = function(root) { 
        # one_sd <- sqrt(35/243) * root
        sum( pos_distances < root) - ncells_halfwidth},
      interval = c(1e-5, 20))
    weights <- tricube(pos_distances, halfwidth = hw_fit$root)
    nzw     <- weights != 0
    fit <- glm.fit(
      x[nzw,], y[nzw],
      weights = weights[nzw],
      family = poisson("log"),
      offset = log( cs[nzw] ) )
    return( fitted.values(fit)[i] )
    })
}
glmfit_no_subsetting <- function() { 
  sapply(cell_ids,function(i) { 
    pos_distances <- ds[i, ]
    hw_fit <- uniroot(
      f = function(root) { 
        # one_sd <- sqrt(35/243) * root
        sum( pos_distances < root) - ncells_halfwidth},
      interval = c(1e-5, 20))
    fit <- glm.fit(
      x, y,
      weights = tricube(pos_distances, halfwidth = hw_fit$root),
      family = poisson("log"),
      offset = log( cs ) )
    return( fitted.values(fit)[i] )
    })
}

glm_ <- function() {
  sapply(cell_ids,function(i) { 
  pos_distances <- ds[i, ]
  hw_fit <- uniroot(
    f = function(root) { 
      # one_sd <- sqrt(35/243) * root
      sum( pos_distances < root) - ncells_halfwidth},
    interval = c(1e-5, 20))
  fit <- glm(
    y ~ x,
    weights = tricube(pos_distances, halfwidth = hw_fit$root),
    family = poisson("log"),
    offset = log( cs ) )
  return( fitted.values(fit)[i] )
  })
}


cell_ids <- sample(1:ncol(counts), 10)
bm <- microbenchmark::microbenchmark(
glm_(),
glmfit_no_subsetting(),
glmfit_subsetting()
)

cell_ids <- sample(1:ncol(counts), 10)
bm2 <- microbenchmark::microbenchmark(
glm_(),
glmfit_no_subsetting(),
glmfit_subsetting()
)

cell_ids <- sample(1:ncol(counts), 10)
bm3 <- microbenchmark::microbenchmark(
glm_(),
glmfit_no_subsetting(),
glmfit_subsetting()
)






cell_ids <- 1:ncol(counts)
glmfit_subsetting()

locfit_package <- function() {
  
  locfit.raw( x = do.call(lp, c(as.list(as.data.frame(pcs_prcomp[, 1:15])),
                              nn = ncells_halfwidth / ncol(counts),
                              deg=1)),
            y = counts[gene, ],
            family = "poisson",
            link = "log",
            maxit = 50,
            base = log(colSums(counts)), # take totalUMIs into account
            ev = dat() # no extrapolation between cells
)
}


microbenchmark::microbenchmark(
  glmfit_subsetting(),
  locfit_package(),
  times = 10
)

