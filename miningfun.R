filter.accession.string <- function(x) {
  splitted <- strsplit(x, ";|/| |:")[[1]]
  splitted <- splitted[grep(".*\\.[0-9]{1,2}", splitted)]
  underscore.search <- grep("[A-Z]{1,2}_", splitted)
  if (length(underscore.search)>0) {
    splitted <- splitted[-underscore.search]
  }
  return(c(splitted[1],
           paste(splitted[2:length(splitted)],
           collapse = ",")))
}
look.for.plasmid.hits <- function(hit) {
  plasmidhits <- grep("plasmid", hit)
  hits <- t(sapply(hit[plasmidhits],
                   filter.accession.string))
  colnames(hits) <- c("chromosomes", "plasmids")
  rownames(hits) <- NULL
  return(hits)
}
further.clean.up.accessions <- function(x) {
  a <- strsplit(x, ",")[[1]]
  return(paste(a[grep("^CP.*", a)], collapse = ","))
}
split.into.filenames <- function(idx, genome, hits) {
  accession <- strsplit(hits[idx,genome], ",")[[1]]
  filename <- paste("data/ncbi/", genome, "/",
                    accession, ".fasta", sep = "")
  names(filename) <- accession
  return(filename)
}
accessions.from.prokaryotes.csv <-
function(csvpath) {
  # From the provided csvpath, read the "Replicons" 
  # column from the csv, and list out the accessions in 
  # an "accession" column of a data.frame, where the 
  # "system" column indicate to which chromosome 
  # accession the accession belongs to, and the 
  # "folder" column indicate which class of either 
  # "chromosomes" or "plasmids" the accession belongs 
  # to.
  p <- read.csv(csvpath)
  accession <-
  apply(p[grep("plasmid",p[,"Replicons"]),
          c("Assembly","Replicons")],
        1, function(line) {
          accs <-
          gsub(".*(/|:)(.*)", "\\2",
               strsplit(line["Replicons"], ";")[[1]])
          return(c(line["Assembly"],
                   "chromosomes"=accs[1],
                   "plasmids"=accs[2:length(accs)]))
              })
  idx <- sapply(accession, function(x) {
                 grep("chromosomes|plasmids",names(x)) })
  system <- unlist(mapply(function(a, i) {
    rep(a[i][1], length(a[i])) }, accession, idx))
  assembly <- unlist(mapply(function(a, i) {
    rep(a["Assembly"], length(a[i])) }, accession, idx))
  folder <- unlist(mapply(function(a, i) {
    c("chromosomes", rep("plasmids", length(a[i])-1)) }, accession, idx))
  accession <- unlist(mapply(function(a, i) { a[i] }, accession, idx))
  a <- data.frame(assembly, system, accession, folder)
  rownames(a) <- NULL
  a
}
