#!/usr/bin/env Rscript
options(width=190, stringsAsFactors = FALSE)

source("palindromicfun.R")
source("knnfun.R")

# {{{ Variables
artdir <- "data/ncbi/art"
meta <- read.csv(paste0(artdir, "/artmeta.csv"), sep="|")

n_jobs <- meta$n_jobs
K.lengths.art <- as.numeric(strsplit(meta$K.lengths.art, ";")[[1]])

# {{{ nucs
aa <- sapply(strsplit(meta$nucs, ",")[[1]], function(x) { strsplit(x, split=":")[[1]] })
nucs <- aa[2,]
names(nucs) <- aa[1,]
rm(aa)
# }}}
# }}}

location <- read.csv("data/ncbi/nanolocation.csv.xz")

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

# {{{ Create outdir
outdir <- "data/ncbi/nanoporecounts/"
if (!dir.exists(outdir)) { dir.create(outdir) }
# }}}

loop.data <-
cbind(K=rep(K.lengths.art, each=nrow(location)),
  do.call(rbind, replicate(length(K.lengths.art),
    cbind(folder=location$folder, accession=location$accession),
    simplify=FALSE)))

FUN <- function(x) {
  K <- as.numeric(x["K"])
  folder <- x["folder"]
  accession <- x["accession"]

  path <- paste0("data/ncbi/nanopore/", folder,
                 "/", accession, ".fastq")
  outpath <- paste0(outdir, "a_K=", K, "_",
                    accession, ".csv")
  if (!file.exists(paste0(outpath, ".xz"))) {
    npk <- NamedPalindromicKmers(K=K, nucs=nucs)
    kmers <- sort(unique(unname(npk)))

    sequences <- microseq::readFastq(path)$Sequence
    if (length(sequences)>500) { sequences <- sequences[1:500] }
    m <- microclass::KmerCount(sequences=sequences,
                               K=K, col.names=TRUE)
    rm(sequences)
    data.table::fwrite(x=as.data.frame(m), file=outpath)
    system(paste("xz -zf -9 -T", n_jobs, outpath))

    m <- KmerCount.matrix.to.canonical.matrix(m=m, K=K, npk=npk, kmers=kmers)
    outpath <- paste0(outdir, "c_K=", K, "_",
                      accession, ".csv")
    data.table::fwrite(x=as.data.frame(m), file=outpath)
    system(paste("xz -zf -9 -T", n_jobs, outpath))
  }
}
cl <- parallel::makeCluster(n_jobs)
parallel::clusterExport(cl, "outdir")
parallel::clusterExport(cl, "n_jobs")
parallel::clusterExport(cl, "KmerCount.matrix.to.canonical.matrix")
parallel::clusterExport(cl, "NamedPalindromicKmers")
parallel::clusterExport(cl, "ExpectedKmers")
parallel::clusterExport(cl, "nucs")
out <- parallel::parApply(cl=cl, X=loop.data, MARGIN=1, FUN=FUN)
out
parallel::stopCluster(cl)

print("ok")
