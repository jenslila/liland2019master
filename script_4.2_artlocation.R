#!/usr/bin/env Rscript
options(width = 175, stringsAsFactors = FALSE)

load(file = "data/ncbi/location.Rdata")
Location <- location
source("palindromicfun.R")
source("knnfun.R")

# {{{ Variables
artdir <- "data/ncbi/art"
tmpfile <- paste(artdir, "tmp", sep="/")
errs <- c("_errFree", "")

meta <- read.csv(paste0(artdir, "/artmeta.csv"), sep="|")

n_jobs <- meta$n_jobs
K.lengths.initial <- as.numeric(strsplit(meta$K.lengths.initial, ";")[[1]])
K.lengths.art <- as.numeric(strsplit(meta$K.lengths.art, ";")[[1]])
sequencing.technology <- meta$sequencing.technology
read.length <- meta$read.length
insert.size <- meta$insert.size
insert.size.std <- meta$insert.size.std
sequencing.depths <- as.numeric(strsplit(meta$sequencing.depths, ";")[[1]])
random.seed <- meta$random.seed
accession.suffix <- meta$accession.suffix
n_reads <- meta$n_reads
n_countjobs <- meta$n_countjobs
fqfilterbin <- meta$fqfilterbin
ftbin <- meta$ftbin
fsbin <- meta$fsbin
home <- meta$fsbin
fq2fa <- meta$fq2fa
spaces2zeroes <- meta$spaces2zeroes

# {{{ nucs
aa <- sapply(strsplit(meta$nucs, ",")[[1]], function(x) { strsplit(x, split=":")[[1]] })
nucs <- aa[2,]
names(nucs) <- aa[1,]
rm(aa)
# }}}
# }}}

# {{{ static vars for testing
K <- K.lengths.art[1]
accession <- location[1,"accession"]
err <- errs[1]
D <- sequencing.depths[1]
R <- 2
canonical <- T
# }}}

#  {{{ Taking in integer chunk number as the first 
#      argument from outside of the R script, i.e. the 
#      first command argument.  The argument is fed 
#      from SLURM into R in batch mode by the series of 
#      scripts created by script_8_slurm.R (the 
#      previous script to be run before this one, i.e. 
#      this script is not supposed to be run, only run 
#      through the shell scripts created by the 
#      previous SLURM script, i.e.  this script is not 
#      part of the main pipeline, it is a "side pipe"). 
chunk_num <- commandArgs(trailingOnly=TRUE)[1]

if (n_countjobs>1 & !is.na(chunk_num)) {
  chunk_num <- as.numeric(chunk_num)
  chunks <- string.chunks(location$accession, n_countjobs)
  locidx <- chunks[chunk_num,"start"]:chunks[chunk_num,"stop"]
  location <- location[locidx,]
}
# }}}

cl <- parallel::makeCluster(n_jobs)

# {{{ Determine which sequences are to be selected from 
#     ART fastq files.  First to a check to see if a 
#     set of sequence indecies are already determined, 
#     read them if they are, otherwise create new ones 
#     and save them.  This is important to make the 
#     forward and backward pipes work on the same data 
#     in the non-canonical case.  This is dog slow, and 
#     should be done in parallel once and for all 
#     before the for loops.

parallel::clusterExport(cl, "artdir")
parallel::clusterExport(cl, "sequencing.technology")
parallel::clusterExport(cl, "D")
X <- expand.grid(location[, "accession"], errs)

# {{{ Function: indecies.was.fetched 
indecies.was.fetched <-
function(x) {
  accession <- x[1]
  err <- x[2]
#   D <- as.numeric(x[3])
  rfiles <-
  paste0(artdir, "/", sequencing.technology, "/D=", D,
         "/", accession, err, ".R", 1:2, ".fq.xz.idx.xz")
  return(file.exists(rfiles[1]))
}
# }}}
out <- parallel::parApply(cl=cl, X=X, MARGIN=1,
                          FUN=indecies.was.fetched)
# {{{ New artlocation dataset consists only of genomes
#     with complete indecies
accessions.to.omit <- unique(as.character(X[!out,1]))

if (length(accessions.to.omit)>0) {
  systems.to.omit <-
  unique(location[location[,"accession"] %in%
                  accessions.to.omit,"system"])
  cat(paste0("Omitting ", length(systems.to.omit),
             " genomes from location dataset, because",
             " they could not be simulated by ART in a",
             " lean way.\n"))
  location <- location[!(location$system %in% systems.to.omit),]
  z <- xzfile("data/ncbi/artlocation.csv.xz", open = "w")
  write.csv(location, file = z, row.names = FALSE, quote = FALSE)
  close(z)
}
# }}}
# }}}

parallel::stopCluster(cl)

print("ok")
