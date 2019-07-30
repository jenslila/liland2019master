#!/usr/bin/env Rscript
options(width=190, stringsAsFactors = FALSE)

### Looking up what became simulated, and saving
#   the result to nanolocation.csv.
location <-
read.csv("data/ncbi/location.csv.xz")

paths <-
paste0("data/ncbi/nanopore/", location$folder, "/",
       location$accession, ".fastq")
idx <- file.exists(paths)
location <- location[idx,]
idx <- sapply(unique(location$system), function(system) {
  lines <- location[location$system==system,]
  length(table(lines$folder))})==2
location <- location[location$system %in%
                     unique(location$system)[idx],]
file <- "data/ncbi/nanolocation.csv"
write.csv(location, file = file, row.names = FALSE, quote = FALSE)
system(paste("xz -zf -9", file))

print("ok")
