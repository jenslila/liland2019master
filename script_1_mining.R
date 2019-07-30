#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, width=170)
source("miningfun.R")

location <-
accessions.from.prokaryotes.csv("data/ncbi/prokaryotes-2019-04-07.csv")

# {{{ Save meta variables for later scripts.  It is 
#     easier to define n_jobs and others here, and 
#     reference it from that csv in later scripts, than 
#     defining this locally in later scripts.
home <- Sys.getenv("HOME")
fsbin <- paste0(home, "/bin/fastasplitn")
nucs <- c("A"="T", "C"="G", "G"="C", "T"="A")

meta <- 
data.frame(
n_reads = 500,
n_countjobs = 4,
n_jobs = 5,
n_folds = 7,
pseudo_count = 0.5,
pca_lim = .9999,
K.lengths.initial = paste(1:8, collapse=";"),
K.lengths.art = paste(c(2, 5, 8), collapse=";"),
n_neighbors = paste(c(1, 3, 9), collapse=";"),
sequencing.technology = "HS25",  # HiSeq 2500
read.length = 150,
insert.size = 750,
insert.size.std = 100,
sequencing.depths = paste(c(10, 100), collapse=";"),  # X10, X100, and X1000 is equivalent to shallow, normal, and deep sequencing depth
random.seed = 777,
artdir = "data/ncbi/art",
accession.suffix = ".fasta.xz",
nucs = paste0(names(nucs), ":", nucs, collapse=","),
home = home,
ftbin = paste0(home, "/bin/fastatuples"),
fsbin = fsbin,
fq2fa = paste(fsbin, "-q4 -q2f -q -n 1 -p 1"),
spaces2zeroes = "awk '{gsub(/ /,\"0\");}1'",
fqfilterbin = "Rscript script_FilterFASTQRPairs.R",
art = "../art/art_src_MountRainier_Linux/art_illumina")
# }}}
write.table(x=meta, file=paste0(meta$artdir, "/artmeta.csv"),
            row.names = FALSE, quote = FALSE, na="", 
            sep="|")

# {{{ pdist
class.summary <-
t(sapply(unique(location$system),
function(system) {
  table(location[location$system==system, "folder"])
}))
cat("Statistics on plasmid counts:\n")
sapply(list(class.summary[,"plasmids"]), function(x) {
  c("X"=mean(x),
    "S"=sd(x),
    "n"=length(x),
    "quantile"=quantile(x))
})[,1]
p <- table(class.summary[,"plasmids"])
p <- cbind("n"=as.numeric(names(p)), "nn"=p)
file <- "data/ncbi/pdist.csv"
# }}}
write.csv(x=p, file=file,
          row.names=FALSE, quote=FALSE, na="")

# {{{ Echoing out number of plasmids and chromosomes from the
#     search filtering
cat(paste("From the NCBI search result a total of ",
          length(location[,"folder"]),
          " accession numbers, ",
          sum(location[,"folder"] == "plasmids"),
          " plasmids and ",
          sum(location[,"folder"] == "chromosomes"),
          " chromosomes, were filtered and structured.\n",
          sep=""))
# }}}

# {{{ Only looking at a toy dataset of the first few 
#     accessions for now, in the future comment this line to 
#     run the pipeline on all the accessions
location <- location[location$assembly %in% unique(location$assembly)[1:5],]
# }}}
save(location, file = "data/ncbi/location.Rdata")
z <- xzfile("data/ncbi/location.csv.xz", open = "w")
write.csv(location, file = z, row.names = FALSE, quote = FALSE)
close(z)

# {{{ A data.frame showing what files exists, to determine 
#     which potential new results to download and compress 
#     next, so as not to bother the NCBI server too much.  
#     Also create needed dirs, so nobody has to guess around 
#     as what paths to prepare before running this
start.path <- "data/ncbi/"
for (dir in paste0(start.path, unique(location$folder))) {  # Create dirs if nonexistent
  if (!dir.exists(dir)) { dir.create(dir) }
}
filepaths <-
paste0(start.path, location$folder, "/", 
       location$accession, ".fasta")
exists <- file.exists(paste0(filepaths, ".xz"))
overview <- data.frame(accession = location$accession,
                       filepath = filepaths,
                       exists = exists)
what.to.download <- overview[!overview[,"exists"],]
# }}}

# {{{ Feeding accession numbers and filepaths into 
#     micropan::entrezDownload, and compressing outfile 
#     using highest xz compression
if (nrow(what.to.download)>0) {
  out <-
  sapply(1:nrow(what.to.download), function(line) {
    accession <- what.to.download[line,"accession"]
    filepath <- what.to.download[line,"filepath"]
    micropan::entrezDownload(accession=accession,
                             out.file=filepath,
                             verbose=TRUE)
    system(paste("xz -zf -9 -T", meta$n_jobs, filepath))
  })
}
# }}}

print("ok")
