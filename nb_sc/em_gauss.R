# Code to demonstrate a simple Gaussian EM clustering

cl <- c( rep( "A", 200 ), rep( "B", 400 ) )
y <- c( rnorm( 200, 4, 1.5 ), rnorm( 400, 7, 1 ), rnorm( 500, 5, .5 ) )
hist(y,100)

plot( y, col = 1 + ( dnorm( y, 5, 2, log=TRUE ) - dnorm( y, 7, 1, log=TRUE ) + log(200/400) > 0 ) )

table( cl, dnorm( y, 5, 2, log=TRUE ) - dnorm( y, 7, 1, log=TRUE ) + log(200/400) > 0 )



icl <- 1 + ( runif( length(y) ) < .5 )
plot( y, col = icl )

mu1 <- mean( y[ icl==1 ] )
sd1 <- sd( y[ icl==1 ] )
mu2 <- mean( y[ icl==2 ] )
sd2 <- sd( y[ icl==2 ] )

c( mu1, sd1, mu2, sd2 )

icl <- 1 + ( dnorm( y, mu1, sd1, log=TRUE ) - dnorm( y, mu2, sd2, log=TRUE ) > 0 )

if( mean(y[icl==2]) < mean(y[icl==1]) )
  icl <- 3 - icl 

plot( y, col = icl )

table(icl)


icl <- sample( 1:3, length(y), replace=TRUE )

mus <- sapply( 1:3, function(i) mean( y[ icl==i ] ) )
sds <- sapply( 1:3, function(i) sd( y[ icl==i ] ) )

mus
sds

icl <- apply( sapply( 1:3, function(i) dnorm( y, mus[i], sds[i], log=TRUE ) ), 1, which.max )

plot( y, col = icl )

table(icl)
