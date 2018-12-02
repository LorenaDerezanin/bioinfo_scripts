###################################################################
###--- Script 3 of the procedure.
###---
###--- You need to run script1 and script2 prior to this.
###--- Here we estimate the copy numbers of the MHC genes.
###--- Only bootstrap-estimates are considered here. Point
###--- estimates are achieved by averaging over the bootstrap-results
###--- for each gene and species.
###################################################################
source( "copynumfun.R" )
options( stringsAsFactors=F )
source( "hits_table_generator.R" )

### Loads data and print som summary
cat( "Loading data...\n" )
#load( "data/Teleost_all_2015.RData" )
load( "res/ref_keep.RData" )
load( "res/ref_boot.RData" )
#uref <- unique( gene )
#uref <- uref[3:length(uref)]             # two first are not reference markers
#r.rows <- 9:nrow( hits )                 # eight first row are not reference markers
#ref.hits <- hits[r.rows,]
uref <- gene_names
ref.hits <- hits
marker.length <- 80*3                     # marker length in bases ## should be 210 for Hox copy number estimation
#G.size <- colMeans( genome.size )        # eight different estimates of genome size, using the mean
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


### Estimating copy number of the MHC-U gene in each species
### based on read-counts and reference marker copy numbers.
### The bootstrapping is as follows:
### For each bootstrap sample, the corresponding bootstrap-result
### for the reference markers is used. In addition, not all
### reference markers are used, but another bootstrap-sample
### selects which reference markers to use each time.
#gene.hits <- hits[1:4,]           # first 4 rows are MHC-U hits
N.boot <- 1000
MHC_U.boot <- matrix( 0, nrow=n.fish, ncol=N.boot )
rownames( MHC_U.boot ) <- colnames( gene.hits )
for( ss in 1:n.fish ){
  cat( "Bootstrapping", rownames(MHC_U.boot)[ss], "...\n" )
  yg.mat <- matrix( ref.hits[,ss], nrow=4, byrow=F )
  idx.in <- which( ref.keep[ss,] )
  yg.mat <- yg.mat[,idx.in]
  GGG <- length( idx.in )
  for( b in 1:N.boot ){
    ref.gamma <- ref.boot[ss,idx.in,b]               # using the reference bootstrap results
    M <- matrix( NA, nrow=4, ncol=GGG )
    idx <- sample( 1:GGG, size=GGG, replace=TRUE )   # sampling which reference markers to use
    for( i in 1:4 ){
      for( g in 1:GGG ){
        if( (yg.mat[i,idx[g]] != 0) & (ref.gamma[idx[g]] != 0) ){
          M[i,g] <- ref.gamma[idx[g]] * gene.hits[i,ss]/yg.mat[i,idx[g]] 
        }
      }
    }
    MHC_U.boot[ss,b] <- mean( M, na.rm=TRUE ) # averaging over all BLAST-cutoffs and selected references
  }
}
par( mar=c(5,6,1,1) )
boxplot( t(MHC_U.boot ),cex.axis=0.6,horizontal=T,las=2,xlab="MHC-U copy number")





### Estimating copy number of the MHC-Z gene in each species
### based on read-counts and reference marker copy numbers.
### The bootstrapping is as follows:
### For each bootstrap sample, the corresponding bootstrap-result
### for the reference markers is used. In addition, not all
### reference markers are used, but another bootstrap-sample
### selects which reference markers to use each time.
#gene.hits <- hits[5:8,]           # Rows 5-8 are MCH-Z hits
N.boot <- 1000
MHC_Z.boot <- matrix( 0, nrow=n.fish, ncol=N.boot )
rownames( MHC_Z.boot ) <- colnames( gene.hits )
for( ss in 1:n.fish ){
  cat( "Bootstrapping", rownames(MHC_Z.boot)[ss], "...\n" )
  yg.mat <- matrix( ref.hits[,ss], nrow=4, byrow=F )
  idx.in <- which( ref.keep[ss,] )
  yg.mat <- yg.mat[,idx.in]
  GGG <- length( idx.in )
  for( b in 1:N.boot ){
    ref.gamma <- ref.boot[ss,idx.in,b]               # using the reference bootstrap results
    M <- matrix( NA, nrow=4, ncol=GGG )
    idx <- sample( 1:GGG, size=GGG, replace=TRUE )   # sampling which reference markers to use
    for( i in 1:4 ){
      for( g in 1:GGG ){
        if( (yg.mat[i,idx[g]] != 0) & (ref.gamma[idx[g]] != 0) ){
          M[i,g] <- ref.gamma[idx[g]] * gene.hits[i,ss]/yg.mat[i,idx[g]] 
        }
      }
    }
    MHC_Z.boot[ss,b] <- mean( M, na.rm=TRUE ) # averaging over all BLAST-cutoffs and selected references
  }
}
par( mar=c(5,6,1,1) )
boxplot( t(MHC_Z.boot ),cex.axis=0.6,horizontal=T,las=2,xlab="MHC-Z copy number")




