###################################################################
###--- Script 2 of the procedure.
###---
###--- You need to run script1 prior to this.
###--- Here we estimate the reference markers copy numbers
###--- Step 1 - Point estimates
###--- Step 2 - Bootstrapping 
###################################################################
options( stringsAsFactors = FALSE )
source( "copynumfun.R" )
source( "hits_table_generator.R" )

### Loads data and print some summary
cat( "Loading data...\n" )
#load( "data/Teleost_all_2015.RData" )
load( "res/ref_keep.RData" )
uref <- gene_names
ref.hits <- hits
marker.length <- 80*3                     # marker length in bp
#uref <- unique( gene )
#uref <- uref[3:length(uref)]             # two first are not reference markers
#r.rows <- 9:nrow( hits )                 # eight first row are not reference markers
#ref.hits <- hits[r.rows,]
#marker.length <- 90*3                    # marker length in bases
#G.size <- colMeans( genome.size )         # eight different estimates of genome size, using the mean
G.size <- c(2110508336, 2405352861, 2192934624)    # assembly size estimated from k=31-91, s=2 (kmergenie)
n.fish <- ncol( ref.hits )                         # N of species observed 
n.ref <- length( uref )                            # N of reference markers/genes
N.reads <- c(1121237882, 877225626, 891832240)     # N of concatenated raw PE reads
cuts <- c(-4,-6,-8,-10)
cat( "   we have data for", n.ref, "reference markers\n" )
cat( "   we have data for", n.fish, "species\n" )
cat( "   reference markers are", marker.length, "bases long\n" )
x <- sapply( 1:n.fish, function(i){
  cat(colnames(ref.hits)[i], "has", G.size[i],"basepairs and",N.reads[i],"reads\n")
} )


### Point estimates of reference marker copy number in each species
### based on read-counts and ALL reference markers for the specific species.
### The matrix ref.point contains one row for each species(__ rows)
### and one column for each marker (19), but since some
### markers have been excluded for the various species (see script1)
### there are some NA in the matrix.

# point estimation uses sample data to calculate a single value to serve as a best guess of copy numbers for sp.
ref.point <- matrix( NA, nrow=ncol(ref.hits), ncol=n.ref )
rownames( ref.point ) <- names( G.size )
colnames( ref.point ) <- uref
for( ss in 1:n.fish ){
  cat( rownames(ref.point)[ss], ":\n" )
  # same as C <- coverage[ss] 
  C <- (N.reads[ss]*marker.length)/G.size[ss]         # coverage for this species
  idx.in <- which( ref.keep[ss,] )                    # indices of ref. genes used for this species
  n.in <- length( idx.in )                            # number of ref. genes used for this sp.
  # all observed hits for this sp., columns = ref.genes, rows = 4 cutoffs
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F )
  y.mat <- y.mat[,idx.in]                             # hits for kept ref. genes for this sp.
  lst <- copynum.iter.iter( y.mat, C, rep( 1, n.in ), as.integer=FALSE )
  # estimated copy number for all kept ref. genes for this sp.
  ref.point[ss,idx.in] <- lst$Gamma
}
save( ref.point, file="res/ref_point.RData" )
boxplot(ref.point,las=2,pch=16,ylab="Reference marker copy number")
# ref_10 in fish_84 present in 3 copies, but kept in ref.keep?


### The bootstrap procedure

# The difference to the point estimate above is that
# estimates are based on a bootstrap-sample of the reference
# markers for current species (not all markers), and this is repeated 1000 times.
### When estimating the copy number for marker g, this marker must of course
# be included, but the OTHER markers are bootstrapped 1000 times, and this
# is repeated for every marker for each species.
### The 3-dimensional array ref.boot stores the results.

# phi.hat - normalized average number of observed ref. gene hits for each cutoff for this sp.
# aka estimated stringency factor Fs with a hat (^)
phi.hat <- matrix( -1, nrow=4, ncol=n.fish )
rownames( phi.hat ) <- paste( "log10(E)=", cuts, sep="" ) # gives names to stringency rows 
colnames( phi.hat ) <- colnames( ref.hits )
N.boot <- 1000
# 3-dim array
ref.boot <- array( NA, dim=c(n.fish, n.ref, N.boot), dimnames=list( Fish=names(G.size), REF=uref, Boot=1:N.boot ) )
for( ss in 1:n.fish ){
  cat( rownames(ref.point)[ss] )
  C <- (N.reads[ss]*marker.length)/G.size[ss]
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F )
  idx.in <- which( ref.keep[ss,] )
  y.mat <- y.mat[,idx.in]
  n.in <- length( idx.in )
  for( g in 1:n.in ){
    idb <- which( 1:n.in != g )    # randomly sampled g-1 genes from the ref.keep for this sp.
    for( b in 1:N.boot ){
      # indices of sampled g-1 genes with replacement
      idx.boot <- c( g, sample( idb, size=(n.in-1), replace=TRUE ) ) 
      y.boot <- y.mat[,idx.boot]    # randomized observed ref. hits regarding indices in idx.boot
      lst <- copynum.iter.iter( y.boot, C, rep( 1, n.in ), max.iter=10, as.integer=FALSE, verbose=FALSE )
      ref.boot[ss,idx.in[g],b] <- lst$Gamma[1]  # boot. estimate for this marker g for this sp.
    }
    cat( "." )
  }
  cat( "\n" )
}

save( ref.boot, file="res/ref_boot.RData" )
boxplot(t(ref.boot[3,,]), ylim=c(0,4), las=2,ylab="Reference marker copy number")
abline(h=1)

