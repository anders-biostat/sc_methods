rawCounts_RNA <- as.matrix( read.csv( "~/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", header=TRUE, row.names=1 ) )

# Remove mouse cells and mouse genes
gene_species <- sapply( strsplit( rownames(rawCounts_RNA), "_" ), `[`, 1 ) 
cell_counts_by_species <- sapply( unique( gene_species ), function(x) 
  colSums( rawCounts_RNA[ gene_species == x, , drop=FALSE ] ) )
rawCounts_RNA <- rawCounts_RNA[ 
  gene_species == "HUMAN",
  cell_counts_by_species[,"HUMAN"] / cell_counts_by_species[,"MOUSE"] > 10 ]
rownames(rawCounts_RNA) <- sub( "HUMAN_", "", rownames(rawCounts_RNA) )
counts <- rawCounts_RNA
rm( rawCounts_RNA )

# sum of UMIs per cell to nomalize
countsums <- colSums( counts )

means <- colMeans( t(counts) / countsums )
vars <- genefilter::rowVars(t( t(counts) / countsums ))

# let's use this as first plot
plot( means, vars, pch=".", log="xy" )
xg <- 10^seq( -9, -1, length.out=1000 )
lines( xg, xg * mean(1/countsums), col="red" )

# or better this
plot( means, vars / means / mean(1/countsums), pch=".", log="xy" )
abline( h=1, col="red")


# Now choose a gene

gene <- "CD3E"
k <- counts[gene,]

# Plot is against the countsum. Is this bimodal?
plot( countsums, counts[gene,] )

# Fit a mixure of two NBs on it, with these starting parameters
mu <- c( 1e-5, 1e-1 )
alpha <- c( 5, 5 )

# E step of EM algorithm
# What are the log odds of a cell coming from NB1 vs from NB2?
logodds <-
  dnbinom( k, mu = countsums * mu[1], size = 1/alpha[1], log=TRUE) - 
     dnbinom( k, mu = countsums * mu[2], size = 1/alpha[2], log=TRUE)

# Probabilities for class memberships, as calculated from log odds above
prob <- cbind( exp(logodds) / ( 1 + exp(logodds) ), NA )
prob[,2] <- 1 - prob[,1]

# plot it, colouring by probability
plot( countsums, k, col=rgb( prob[,1], prob[,2], 0 ) )

# M step
# now get MLEs for each class, by taking sums of log likelihoods, weighted
# by class probability
for( i in 1:2 ) {
  r <- optim( c( mu[i], alpha[i] ), function(x)
    -sum( dnbinom( k, mu = x[1] * countsums, size = 1/x[2], log=TRUE ) * prob[,i] ) )
  mu[i] <- r$par[1]
  alpha[i] <- r$par[2] }

# Show new parameters
mu
alpha

# Now iterate M and E steps
