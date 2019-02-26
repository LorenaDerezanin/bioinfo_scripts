###################################################################
###--- Script 1 of the procedure.
###---
###--- Here we first decide which reference markers are possible
###--- to use for estimating copy numbers.
###--- Step 1 - Sanity check on read coverage.
###--- Step 2 - Check on effect of BLAST-cutoff.
###--- Step 3 - Final check if copy number estimates are not 0.
###################################################################
options( stringsAsFactors = FALSE )
source( "copynumfun.R" )


### Loads data and print some summary
cat( "Loading data...\n" )
source( "hits_table_generator.R" )
uref <- gene_names
ref.hits <- hits
marker.length <- 80*3                            # marker length in bp
# G.size <- colMeans( genome.size )        
G.size <- c(2110508336, 2405352861, 2192934624)  # assembly size estimated from k=31-91, s=2 (kmergenie)
n.fish <- ncol( ref.hits )                       # number of species observed 
n.ref <- length( uref )                          # number of reference markers/genes
N.reads <- c(1121237882, 877225626, 891832240)   # N of concatenated raw PE reads
cuts <- c(-4,-6,-8,-10)
cat( "   we have data for", n.ref, "reference markers\n" )
cat( "   we have data for", n.fish, "species\n" )
cat( "   reference markers are", marker.length, "bases long\n" )
x <- sapply( 1:n.fish, function(i){
  cat(colnames(ref.hits)[i], "has", G.size[i],"basepairs and",N.reads[i],"reads\n")
} )


### First data sanity check - checking for too low coverage
### Some reference markers may, for some reason, have extremely low coverage
### We would like to exclude these since they may affect all other estimates
### of copy numbers later
low.cov <- matrix( 0, nrow=n.fish, ncol=n.ref )
rownames( low.cov ) <- colnames( ref.hits )
colnames( low.cov ) <- uref
coverage <- (N.reads*marker.length)/G.size    # This is the expected coverage per species
for( i in 1:n.ref ){
  cat( "Reference ", uref[i], "\n" )
  
  idd <- seq( 1, nrow(ref.hits), 4 )                          # hits under cutoff -4
  ohits_cutoff5 <- ref.hits[idd[i],]                          # observed number of hits for reference gene uref[i] at cutoff -4
  normalized_hits <- ohits_cutoff5/coverage                   # number of hits per unit of coverage (~percentage of hits compared to coverage)
  phi1 <- mean( normalized_hits )                             # very rough estimate of cutoff-factor Fs
  ehits1 <- coverage*phi1                                     # expected number of hits based on phi1
  # chi- square - for diff. between observed and estimated hits is high
  rr1 <- sign( ohits_cutoff5 - (ehits1-3*sqrt(ehits1)) )      # if rr1 is negative the coverage is very low
 
  
  idd <- seq( 2, nrow(ref.hits), 4 )                          # hits under cutoff -6
  phi2 <- mean( ref.hits[idd[i],]/coverage )                  # very rough estimate of cutoff-factor
  ehits2 <- coverage*phi2                                     # expected number of hits based on phi2
  rr2 <- sign( ref.hits[idd[i],] - (ehits2-3*sqrt(ehits2)) )  # if rr2 is negative the coverage is very low

  idd <- seq( 3, nrow(ref.hits), 4 )                          # hits under cutoff -8
  phi3 <- mean( ref.hits[idd[i],]/coverage )                  
  ehits3 <- coverage*phi3                                     # expected number of hits based on phi3 
  rr3 <- sign( ref.hits[idd[i],] - (ehits3-3*sqrt(ehits3)) )

  idd <- seq( 4, nrow(ref.hits), 4 )                          # hits under cutoff -10
  phi4 <- mean( ref.hits[idd[i],]/coverage )
  ehits4 <- coverage*phi4
  rr4 <- sign( ref.hits[idd[i],] - (ehits4-3*sqrt(ehits4)) )
  
  rmat <- matrix( c(rr1,rr2,rr3,rr4), ncol=4, byrow=F )       # The rr1,...,rr4 values
  low.cov[,i] <- rowSums( rmat )                              # summing the signs, i.e. need 3 out of 4 negative to
}                                                             # get a negative sum
ref.keep <- (low.cov >= 0 )
r.disc <- rowSums(!ref.keep)
f.disc <- colSums( !ref.keep )                                # number of species for which each ref.gene was discadred   
cat( "Number of discarded reference markers per species:\n" )
x <- sapply( 1:n.fish, function(i){
  cat(rownames(ref.keep)[i], "discards",r.disc[i],"reference markers\n")
} )
cat( "Number of discarded species per reference marker:\n" )
x <- sapply( 1:n.ref, function(i){
  cat(colnames(ref.keep)[i], "is discarded in",f.disc[i],"species\n")
} )


### Manual curation
### These two references were discarded in the end because they are
### probably too closely linked (physical position) to some other references
# ref.keep[,2] <- FALSE  # Ref 1b is OUT
# ref.keep[,7] <- FALSE  # Ref 2d is OUT



### Next, checking which reference genes have read-hits that
### cannot be explained by the cutoff-model.
### In short, we expect that as the blast-cutoff is made stricter, the
### number of hits should decrease in a smooth way.
### Here we look for cases where this is clearly violated
mse <- matrix( 0, nrow=n.fish, ncol=n.ref )       # mse represents residual variance (S2g) for a marker 
colnames( mse ) <- uref
rownames( mse ) <- colnames( ref.hits )
for( ss in 1:n.fish ){
  cat( rownames(mse)[ss], ":\n" )
  # same as C <- coverage[ss] 
  C <- (N.reads[ss]*marker.length)/G.size[ss]     # coverage for this species
  idx.in <- which( ref.keep[ss,] )                # indices of ref. genes used for this species
  n.in <- length( idx.in )                        # number of ref. genes used for this species
  # all observed hits for this species, columns = ref.genes, rows = 4 cutoffs
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F ) 
  y.mat <- y.mat[,idx.in]                         # hits for kept ref. genes for this species
  lst <- copynum.iter.iter( y.mat, C, rep( 1, n.in ), as.integer=FALSE )
  
  g.mat <- matrix( rep( lst$Gamma, 4 ), nrow=4, byrow=TRUE )
  p.mat <- matrix( rep( lst$Phi, n.in ), ncol=n.in, byrow=FALSE )
  # estimated number of ref. gene hits for this species for all 4 cutoffs
  y.hat <- C * g.mat * p.mat
  # difference between observed and estimated ref.gene hits 
  r.mat <- y.mat - y.hat
  # calculation of residual variance S^2g (formula 7. in SI) for each gene
  mse[ss,idx.in] <- apply( r.mat, 2, function(x){sum(x^2)/(4-1)} )
  # range returns min and max of given arguments
  rr <- range( c( as.vector( y.mat ), as.vector(y.hat) ) )
  plot( rr, rr, type="l", col="red", main=rownames(mse)[ss], 
        xlab="Observed number of reads", ylab="Predicted number of reads" )
  points( y.mat, y.hat, pch=16 )
  Sys.sleep( 2 )
}
# graph for observed and estimated ref. gene hits - closer the points are to the red line, smaller the diff.


sigma2 <- mean( mse, trim=0.01 )          # sigma2 is the variance of the error term, trimmed off 1% of extreme values
lambda <- (4-1)*mse/sigma2                # if lambda is large, it indicates that residual variance 
limit <- qchisq( 0.99, df=(4-1) )         # is much larger than we expected - data fits the model poorly
ref.keep <- (lambda <= limit)&ref.keep    # marker kept if lambda is lower than 99% quantile of the chi-square 
                                          # distribution with 3 degrees of freedom

cat( "Number of discarded reference markers per species:\n" )
x <- sapply( 1:n.fish, function(i){
  cat(rownames(ref.keep)[i], "discards",rowSums(!ref.keep)[i],"reference markers\n")
} )
cat( "Number of discarded species per reference marker:\n" )
x <- sapply( 1:n.ref, function(i){
  cat(colnames(ref.keep)[i], "is discarded in",colSums( !ref.keep )[i],"species\n")
} )


### Manual curation again
### Fish 21 suffers from having too few reference markers in the 
### analysis above. Inspect the plot for fish_21, it looks fine!
### We decide to keep reference 8 and 13 even if they are
### excluded in the procedure above.
# ref.keep[21,14] <- TRUE  # Ref 8 is IN
# ref.keep[21,20] <- TRUE  # Ref 13 is IN



### Third, eliminating reference genes that still estimate 
### to 0 copies, as these will blow the other copy number 
### estimates sky high later...
ref.hat <- matrix( 0, nrow= n.fish, ncol=n.ref )
for( ss in 1:n.fish ){
  cat( rownames(ref.keep)[ss], ":\n" )
  C <- (N.reads[ss]*marker.length)/G.size[ss]
  idx.in <- which( ref.keep[ss,] )
  n.in <- length( idx.in )
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F )
  y.mat <- y.mat[,idx.in]
  lst <- copynum.iter.iter( y.mat, C, rep( 1, n.in ), as.integer=FALSE )
  ref.hat[ss,idx.in] <- lst$Gamma
}
ref.keep <- (round(ref.hat) > 0)&ref.keep
cat( "Number of discarded reference markers per species:\n" )
x <- sapply( 1:n.fish, function(i){
  cat(rownames(ref.keep)[i], "discards",rowSums(!ref.keep)[i],"reference markers\n")
} )
cat( "Number of species per discarded reference marker:\n" )
x <- sapply( 1:n.ref, function(i){
  cat(colnames(ref.keep)[i], "is discarded in",colSums( !ref.keep )[i],"species\n")
} )


save( ref.keep, file="res/ref_keep.RData" )

