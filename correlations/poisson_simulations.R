



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
   sizefactos <- sample(sizefactors, n_low+nhigh, replace = TRUE)
 }
  
  c( rpois(n_low,  mean_low * sizefactors[1:n_low]),
     rpois(n_high, mean_high *sizefactors[ (n_low+1) : (n_low+n_high) ] )
    ) / sizefactors 
   
}







# Plausible gene means inspired by CITEseq dataset ----------------------------------
# When simulating scRNAseq data (e.g. as Poisson counts for homogeneous populations),
# it is useful to cover biologically relevant ranges of gene expression magnitude.
# Here are "plausible gene means", meaning they could actually be raw counts of ~ 1000 genes
# observed for a cell with average library size (nUMI aka colSums) in a 10X scRNAseq dataset.
# These can be used directly as Poisson rates to simulate data, or they can be
# adjusted with sizefactors centered around 1 (1 representing the cell with average sequencing depth).
  plausible <- c( 
    runif(575, min = .5, max = 1),
    runif(267, min = 1,  max = 2),
    runif(190, min = 2,  max = 10),
    runif(70, min = 10, max = 40)  )



  



# plausible_sizefactors ---------------------------------------------------

sf <- rf(1000,20,22, 0)
# turns out this is an acceptable simulation of sizefactors around 1!

 #  # load the citeseq dataset first, then you can get the sizefactors:
 #  sf <- colSums(rawC) / mean(colSums(rawC))
 #  hist(sf, breaks = 0:150 * .1)
 #  
 #  # let's fit an function to that. With guessing parameters, it'd look like this:
 #  x <- seq(0, 20, length.out = 1000)
 #  hist(sf, breaks = 0:150 * .1, freq  = F); lines(x, df(x, df1=3, df2=2, ncp = 0), col="blue")
 #  
 #  # lets find maximum likelihood estimates:
 #  optim(c(2, 5, 3), function(params) -sum(df(sf, df1 =  params[1],
 #                                             df2 =  params[2], ncp = 
 #                                               params[3], log = T) ) )
 #  
 #  hist(sf, breaks = 0:150 * .1, freq  = F); lines(x, df(x,20,22, 0), col="blue")


