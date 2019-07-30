#!/usr/bin/env Rscript
options(width=125, stringsAsFactors = FALSE)

nucs <- c("A"="T", "C"="G", "G"="C", "T"="A")

### Read meta info defined in first mining script
artdir <- "data/ncbi/art"
meta <- read.csv(paste0(artdir, "/artmeta.csv"), sep="|")

### Use a few threads for handling count matrix (but 
#   not for counting, as unlike jellyfish(1), 
#   microclass::KmerCount is not paralellized.
n_jobs <- meta$n_jobs
K.lengths.initial <- as.numeric(strsplit(meta$K.lengths.initial, ";")[[1]])
source("palindromicfun.R")  # needs microclass lib

# {{{ Prepare directories
dirs <- paste("data/ncbi/", c("counts", "freq"), "/", sep="")
dirs <- paste(rep(dirs, each=2), c("npk", "all"), sep="")
for (dir in dirs) {  # Create dirs if nonexistent
  if (!dir.exists(dir)) { dir.create(dir, recursive=TRUE) }
}
# }}}

# {{{ Prepare for counting
read.sequence <- function(path) {
  zz <- xzfile(path, open = "r")
  fasta <- readLines(zz)
  close(zz)
  paste(fasta[2:length(fasta)], collapse = "")
}
GenericCountWrapper <-
function(path, K, kmers, nucs) {
  seq <- read.sequence(path=path)
  return(microclass::KmerCount(seq=seq, K=K,
                               col.names=T)[1,][kmers] +
         microclass::KmerCount(seq=ReverseComplement(seq=seq, nucs=nucs), K=K, 
                               col.names=T)[1,][kmers])
}
CountWrapper <- function(path, npk, cpk) {
  PalindromicCanonicalKmerCount(seq=read.sequence(path=path),
                                K=K, npk=npk)[cpk]
}

load(file = "data/ncbi/location.Rdata")
path <- paste("data/ncbi/", location$folder, "/",
              location$accession, ".fasta.xz", sep = "")

cl <- parallel::makeCluster(n_jobs)
parallel::clusterExport(cl, "read.sequence")
parallel::clusterExport(cl, "PalindromicCanonicalKmerCount")
parallel::clusterExport(cl, "ReverseComplement")
# }}}

for (K in K.lengths.initial) {
  parallel::clusterExport(cl, "K")
  cat(paste("Counting ", K, "-mers ... ", sep=""))
  ptm <- proc.time()

  # {{{ Complete design matrix
  counts <-
  t(parallel::parSapply(cl, path, GenericCountWrapper, K=K,
                        kmers=ExpectedKmers(K=K), nucs=nucs))
  dim.all <- dim(counts)
  file <- paste0("data/ncbi/counts/a_K=", K, ".csv")
  write.csv(counts, file = file,
            row.names = FALSE, quote = FALSE)
  system(paste("xz -zf -T", n_jobs, "-9", file))
  # }}}

  # {{{ Canonical design matrix
  counts <- 
  KmerCount.matrix.to.canonical.matrix(m=counts, K=K)
  write.csv(counts, file = file,
            row.names = FALSE, quote = FALSE)
  system(paste("xz -zf -T", n_jobs, "-9", file))
  dim.npk <- dim(counts)
  file <- paste0("data/ncbi/counts/c_K=", K, ".csv")
  # }}}

  # {{{ Timing and dimensions
  cat(paste("done after ",
            (proc.time() - ptm)["elapsed"],
            " seconds. ",
            paste('dim()=["a": (',
                  paste(dim.all,
                        collapse=", "),
                  '), "c": (',
                  paste(dim.npk,
                        collapse=", "),
                  ")].", sep=""),
            "\n", sep=""))
  # }}}
}

parallel::stopCluster(cl)
print("ok")
