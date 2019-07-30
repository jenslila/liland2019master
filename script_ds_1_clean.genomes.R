#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, width=170)

### Set a simple locale not knowing about funky chars
out <- Sys.setlocale("LC_ALL", "C")

source("knnfun.R")

### Read meta info defined in first mining script
artdir <- "data/ncbi/art"
meta <- read.csv(paste0(artdir, "/artmeta.csv"), sep="|")

location <- read.csv("data/ncbi/location.csv.xz")

# {{{ Taking in integer chunk number as the first 
#      argument from outside of the R script, i.e. the 
#      first command argument.  The argument is fed 
#      from SLURM into R in batch mode by the series of 
#      scripts created by script_8_slurm.R (the 
#      previous script to be run before this one, i.e. 
#      this script is not supposed to be run, only run 
#      through the shell scripts created by the 
#      previous SLURM script, i.e.  this script is not 
#      part of the main pipeline, it is a "side pipe"). 

first.arg <- commandArgs(trailingOnly=TRUE)[1]
if (!is.na(first.arg)) {
  chunk_num <- as.numeric(first.arg)
} else {
  stop("Specify number of count job index as first argument.")
}

second.arg <- commandArgs(trailingOnly=TRUE)[2]
if (!is.na(second.arg)) {
  n_countjobs <- as.numeric(second.arg)
} else {
  stop("Specify number of count jobs as second argument.")
}

if (n_countjobs>1 & !is.na(chunk_num)) {
  chunk_num <- as.numeric(chunk_num)
  chunks <- string.chunks(location$accession, n_countjobs)
  locidx <- chunks[chunk_num,"start"]:chunks[chunk_num,"stop"]
  location <- location[locidx,]
}
# }}}

ncbi.dir <- "data/ncbi/"
# {{{ Creating outdirs
dirs <-
paste0(paste0(ncbi.dir, "cleaned"), "/",
       c("plasmids", "chromosomes"))
if (sum(dir.exists(dirs))!=length(dirs)) {
  cmds <- paste("mkdir -p", dirs)
  out <- sapply(cmds, system)
}
# }}}

paths <- t(apply(location, 1, function(x) {
  c("orig"=
    paste0(ncbi.dir, x["folder"], "/",
           x["accession"], ".fasta.xz"),
    "cleaned"=
    paste0(ncbi.dir, "cleaned/", x["folder"], "/",
           x["accession"], ".fasta")
  )
}))

n_jobs <- meta$n_jobs

FUN <- function(x) {
  cmd <- paste("xz -cdT", n_jobs, x["orig"])
  Header <- system(paste(cmd, "| head -n1 | perl -pe 's/^>//g'"), intern=TRUE)
  cmd <- paste(cmd, "| awk 'NR > 1 { print }' | perl -pe 's/[^ACGT\n]/N/g'")
  Sequence <- paste(system(cmd, intern=TRUE), collapse="")
  microseq::writeFasta(out.file=x["cleaned"],
                       fdta=data.frame(Header, Sequence))
#   system(paste("xz -zf9T", n_jobs, x["cleaned"]))
}
apply(paths, 1, FUN)

print("ok")
