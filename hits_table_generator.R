
# Create table with hits for each ref. marker for each e-value cutoff (-4, -6, -8, -10)

species = c(
  "brown_bear",
  "giant_panda",
  "american_black_bear"
)

get_species_files <- function(species) {
  species_path = paste0("/home/derezanin/species_comp/bears/", species, "/tblastn_hits/e2_hits_80aa")
  files <- list.files(path = species_path, full.names = TRUE)
  return(files)
}

get_gene_e_name <- function(gene_name, e) {
  return(paste0(gene_name, "_e_", e))
}

get_gene_names <- function() {
  first_species = species[1]
  files = get_species_files(first_species)
  gene_names <- c()
  for ( f in files ) {
    data <- read.csv(f, header=FALSE, sep="\t")
    gene_name <- strsplit(data[1,1], "_")[[1]][1]
    gene_names <- c(gene_names, gene_name)
  }
  
  return(gene_names)
}

initialize_hits_table <- function(species) {
  
  num_columns = length(species)
  num_genes = 19
  num_e_cutoffs = 4
  num_rows = num_genes * num_e_cutoffs

  gene_e_names = c()
  for (gene_name in get_gene_names()) {
    e_cutoffs <- c(4,6,8,10)
    for ( e in e_cutoffs ) {
      gene_e_name <- get_gene_e_name(gene_name, e)
      gene_e_names <- c(gene_e_names, gene_e_name)
    }
  }
  
  hits_table = matrix(0, nrow = num_rows, ncol = num_columns)
  colnames(hits_table) <- species
  rownames(hits_table) <- gene_e_names
  
  return(hits_table)
}

fill_hits <- function(species, hits) {
  for (s in species) {
    files <- get_species_files(s)
    
    for ( f in files ) {
      data <- read.csv(f, header=FALSE, sep="\t")
      data.e <- data[,11]
      gene_name <- strsplit(data[1,1], "_")[[1]][1]
      
      e_cutoffs <- c(4,6,8,10)
      for ( e in e_cutoffs ) {
        num_hits <- length(data.e[data.e<=10^-e])
        gene_e_name <- get_gene_e_name(gene_name, e)
        hits[gene_e_name, s] <- num_hits
      }
    } 
  }
  
  return(hits)
}

empty_hits <- initialize_hits_table(species)
gene_names <- get_gene_names()
hits <- fill_hits(species, empty_hits)
