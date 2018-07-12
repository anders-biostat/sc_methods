# Convenience functions to sample from lognormal and Poisson by specifiying distribution mean and variance

library( rootSolve )
rlognorm <- function( n, mean, var ) {
  a <- 
    multiroot(
      function(x) c( 
        exp( x[1] + x[2]/2 ) - mean ,
        ( exp(x[2]) - 1 ) * exp( 2*x[1] + x[2] ) - var  ),
      c( log(mean), log(mean+sqrt(var)) - log(mean) ) )
  exp( rnorm( n, a$root[1], sqrt(a$root[2]) ) )
}  

rgamma_mv <- function( n, mean, var ) 
  rgamma( n, scale = var/mean, shape = mean^2 / var )

# Get real cell size values

library( Seurat )
load( "~/Downloads/droplet_Kidney_seurat_tiss.Robj" )
kidney_cnts <- as.matrix(tiss@raw.data)
sf <- colSums(kidney_cnts) / mean(colSums(kidney_cnts))
rm( tiss )

plot( sf, log="y" )


# Simulate counts for a gene with given true mean and variance
# variance = mean^2 * dispersion
n <- length( sf )
truemean <- 2
truedisp <- .3

# Get rate parameters per cell by sampling for n cells from lognormal (or from gamma)
mu <- rlognorm( n, truemean, truemean^2 * truedisp  )
#mu <- rgamma_mv( n, 4, 3^2 )

# Check that we really get the requested means
c( mean(mu), sd(mu), var(mu) )

# Inspect them
hist(mu, 30 )

# Get counts by Poisson sampling with rate parameters as given by my, scaled up by size factors
k <- rpois( n, mu * sf )

# Tabulate and plot
plot( 0:29, tabulate( k, 30 ), type="h" )


# Find mean and dispersion by assuming a size-factor-ajusted negative binomial and find mean (mu) and dispersion (1/size) via maximum likelihood estimation.
optim(
  c( 1, 1 ),
  function(x) -sum( dnbinom( k, mu = x[1] * sf, size=1/x[2], log=TRUE  ) ) )

# The same optimization again, now with cost function (log likelihhod) written out explicitely, and gradient provided
optim( 
  c( 4, .3 ),  # mu, alpha
  function(x) 
    -sum( lgamma( k + 1/x[2] ) - lgamma( 1/x[2] ) - lgamma( k+1 ) - ( k + 1/x[2] ) * log( 1 + sf * x[2] * x[1] ) + k * log( sf * x[2] * x[1] ) ),
  function(x) c(
    -sum( ( k - sf * x[1] ) / ( x[1] + sf * x[2] * x[1]^2 ) ),
    -sum( ( x[2] * ( k - sf * x[1] ) / ( 1 + sf * x[2] * x[1] ) + log( 1 + sf * x[2] * x[1] ) - digamma( k + 1/x[2] ) + digamma( 1/x[2] ) ) / x[2]^2 ) ),
  lower = c( 1e-10, 1e-10 ),
  method = "L-BFGS-B" )

# Both optimizations find the same parameters

# Use this to run through all genes of the kidney data and get MLEs for mean and dispersion for them

kidney_gene_params <- 
  purrr::map_dfr( rownames(kidney_cnts)[1:100], function(g) {
    k <- kidney_cnts[g,]
    o = optim( 
      c( 4, .3 ),  # mu, alpha
      function(x) 
        -sum( lgamma( k + 1/x[2] ) - lgamma( 1/x[2] ) - lgamma( k+1 ) - ( k + 1/x[2] ) * log( 1 + sf * x[2] * x[1] ) + k * log( sf * x[2] * x[1] ) ),
      function(x) c(
        -sum( ( k - sf * x[1] ) / ( x[1] + sf * x[2] * x[1]^2 ) ),
        -sum( ( x[2] * ( k - sf * x[1] ) / ( 1 + sf * x[2] * x[1] ) + log( 1 + sf * x[2] * x[1] ) - digamma( k + 1/x[2] ) + digamma( 1/x[2] ) ) / x[2]^2 ) ),
      lower = c( 1e-10, 1e-10 ),
      method = "L-BFGS-B" )
    data.frame( gene = g, mean = o$par[1], disp = o$par[2] ) } )

plot( disp ~ mean, kidney_gene_params, log="xy")

# Next step: Get standard errors on the MLEs by calculating the second derivatives at the maximum
