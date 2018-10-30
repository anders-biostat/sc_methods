library( quantreg )

# Load test data
load("~/w/rlc_tutorial/citeseq.rda")

# size factors
sf <- colSums(counts)
sf <- sf / mean(sf)
xi <- mean( 1/sf )

# Calculate means and variances
means <- apply( counts, 1, function(x) mean( x / sf ) )
vars  <- apply( counts, 1, function(x) var( x / sf ) )

# Simulate Poisson counts with these
pois <- do.call( rbind, lapply( means, function(mu) {
  k <- rpois( length(sf), mu*sf )
  data.frame( mean=mean(k/sf), var=var(k/sf) ) } ) )

# fit 
fit <- rq( I(( var/mean - xi )^2) ~ I(1/mean), .9, pois[pois$mean > 1e-3,] )
rq_xi <- coef(fit)[1]
rq_b <- coef(fit)[2]

# Plot on the fit scales
mug <- 10^seq(-5,4,length.out=10000)
plot( log(pois$mean), log( (pois$var / pois$mean - xi )^2 ), pch=".")
lines( log(mug), log( rq_xi + rq_b / mug ), col="green" )

plot( pois$mean, (pois$var / pois$mean - xi ), pch=".", log="xy")
lines( mug, sqrt( rq_xi + rq_b / mug ), col="green" )



# Plot 
plot( means, vars/means, pch=".", log="xy" )
abline(h=xi)
points( pois$mean, pois$var / pois$mean, pch=".", col="orange" )
lines( mug, sqrt( rq_xi + rq_b / mug ) + xi, col="green", log="xy" )

# Check even quantile coverage
tapply( pois$var/pois$mean > sqrt( rq_xi + rq_b / pois$mean ) + xi, floor(log10(pois$mean)*2), mean )

qqnorm( ( ( pois$var/pois$mean - xi ) / sqrt( rq_xi + rq_b / pois$mean ) )[pois$mean>1e-3] )
qqline( ( ( pois$var/pois$mean - xi ) / sqrt( rq_xi + rq_b / pois$mean ) )[pois$mean>1e-3] )
plot( log(pois$mean), ( pois$var/pois$mean - xi ) / sqrt( rq_xi + rq_b / pois$mean ) )
