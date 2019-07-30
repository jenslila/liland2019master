#!/usr/bin/env Rscript
options(width = 175, stringsAsFactors = FALSE)

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

D <- sequencing.depths[1]
location <- read.csv(z <- xzfile("data/ncbi/artlocation.csv.xz", open="r"))
close(z)

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
chunk_num <- commandArgs(trailingOnly=TRUE)[1]

if (n_countjobs>1 & !is.na(chunk_num)) {
  chunk_num <- as.numeric(chunk_num)
  chunks <- string.chunks(location$accession, n_countjobs)
  locidx <- chunks[chunk_num,"start"]:chunks[chunk_num,"stop"]
  location <- location[locidx,]
}
# }}}

cl <- parallel::makeCluster(n_jobs)

vars <- expand.grid(errs, K.lengths.art, location$accession)
colnames(vars) <- c("err", "K", "accession")

# {{{ get.seqs.from.fqxz
get.seqs.from.fqxz <-
function(rfile, fq2fa, n_jobs) {
  cmd <-
  paste("xz -cdT", n_jobs, rfile,
        " | ", fq2fa,
        " | sh script_multifasta_concat.sh",
        " | Rscript script_SelectSeqsFromIndecies.R", rfile)
  lines <- data.table::fread(cmd=cmd, header=FALSE)[[1]]
  return(lines[seq(1+1, length(lines)+1, 2)])
}
# }}}
# {{{ clusterExport lines
parallel::clusterExport(cl, "D")
parallel::clusterExport(cl, "NamedPalindromicKmers")
parallel::clusterExport(cl, "ExpectedKmers")
parallel::clusterExport(cl, "nucs")
parallel::clusterExport(cl, "get.seqs.from.fqxz")
parallel::clusterExport(cl, "artdir")
parallel::clusterExport(cl, "sequencing.technology")
parallel::clusterExport(cl, "fq2fa")
parallel::clusterExport(cl, "n_jobs")
parallel::clusterExport(cl, "KmerCount.matrix.to.canonical.matrix")
# }}}

# {{{ KmerCount loop
FUN <-
function(x) {
  err <- x["err"]
  K <- as.numeric(x["K"])
  accession <- x["accession"]
  # {{{ Required vectors to create the canonical set
  #     later on
  npk = NamedPalindromicKmers(K=K, nucs=nucs)
  kmers <- sort(unique(unname(npk)))
  # }}}
  # {{{ Construct the path to rfiles
  rfiles <-
  paste0(artdir, "/", sequencing.technology, "/D=", D,
          "/", accession, err, ".R", 1:2, ".fq.xz")
  # }}}
  # {{{ Read sequences
  seqs.R1 <-
  get.seqs.from.fqxz(rfile=rfiles[1], fq2fa=fq2fa, 
                     n_jobs=n_jobs)
  seqs.R2 <-
  get.seqs.from.fqxz(rfile=rfiles[2], fq2fa=fq2fa, 
                     n_jobs=n_jobs)
  # }}}
  # {{{ Perform the actual K-mer counting
  m1 <- microclass::KmerCount(sequences=seqs.R1, K=K, col.names=TRUE)
  m2 <- microclass::KmerCount(sequences=seqs.R2, K=K, col.names=TRUE)
  # }}}
  # {{{ Wrangle m1 and m2 into canonical and
  #     non-canonical sets
  for (canonical in c(T, F)) {
    outpath <- paste0(artdir, "/", sequencing.technology, "/D=", D,
                      "/xtest_", accession, err, "_",
                      if (canonical) { "c" } else { "a" },
                      "_K=", K, ".csv")
    counts <-
    if (canonical) {
      KmerCount.matrix.to.canonical.matrix(m=m1, K=K, npk=npk, kmers=kmers) +
      KmerCount.matrix.to.canonical.matrix(m=m2, K=K, npk=npk, kmers=kmers)
    } else {
      nuc.seqs.rev.R1 <- microseq::reverseComplement(nuc.sequences=seqs.R1)
      nuc.seqs.rev.R2 <- microseq::reverseComplement(nuc.sequences=seqs.R2)
      m1 + m2 +
      microclass::KmerCount(sequences=nuc.seqs.rev.R1, K=K, col.names=TRUE) +
      microclass::KmerCount(sequences=nuc.seqs.rev.R2, K=K, col.names=TRUE)
    }
    data.table::fwrite(x=as.data.frame(counts), file=outpath)
    system(paste("xz -zf -9 -T", n_jobs, outpath))
  }  # End of canonical loop
  # }}}
}
out <-
parallel::parApply(cl=cl, X=vars, MARGIN=1, FUN=FUN)
# }}}
parallel::stopCluster(cl)

print("ok")
