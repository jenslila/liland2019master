#!/usr/bin/env Rscript
options(width=220, stringsAsFactors = FALSE)

load(file = "data/ncbi/location.Rdata")
source("artfun.R")  # needs ART_Illumina and samtools binaries to be available in $PATH
source("knnfun.R")

# {{{ Variables
artdir <- "data/ncbi/art"
meta <- read.csv(paste0(artdir, "/artmeta.csv"), sep="|")

sequencing.depths = strsplit(meta$sequencing.depths, ";")[[1]]
n_countjobs <- meta$n_countjobs

tmpfile <- paste(meta$artdir, "tmp", sep="/")
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
  tmpfile <- paste(meta$artdir, paste0("tmp", chunk_num, "_"), sep="/")
}
# }}}

# {{{ Check for and create work directories
dirs <- c(meta$artdir,
          paste(meta$artdir, meta$sequencing.technology,
                paste0("D=", sequencing.depths), 
                sep="/"))
for (dir in dirs) {  # Create dirs if nonexistent
  if (!dir.exists(dir)) { dir.create(dir, recursive=TRUE) }
}
# }}}

# {{{ ART loop
loc <- 1
sequencing.depth <- 10

for (loc in 1:nrow(location)) {
# for (loc in 1:2) {
#   for (sequencing.depth in sequencing.depths) {
    folder <- location[loc, "folder"]
    accession <- location[loc, "accession"]
    xz.path <- paste("data/ncbi/", folder, "/",
                     accession, meta$accession.suffix,
                     sep = "")
    cmd <- paste("xz -dc ", xz.path, " > ", tmpfile, sep="")
    system(cmd)
#     print(cmd)
    read.ID.prefix <- fastq.prefix <- infile.path <- tmpfile

    run.ART_Illumina.and.perhaps.samtools(random.seed=meta$random.seed, 
      sequencing.technology=meta$sequencing.technology, 
      insert.size=meta$insert.size, 
      insert.size.std=meta$insert.size.std, 
      sequencing.depth=sequencing.depth, 
      read.length=meta$read.length, 
      infile.path=infile.path, 
      read.ID.prefix=read.ID.prefix, 
      fastq.prefix=fastq.prefix, 
      samout=FALSE,
      errfree=TRUE,
      art=meta$art)
    
    ### For any infile.path we are now always left with the 
    #   output in four FASTQ files tmp{1,2}.fq and 
    #   tmp_errFree{1,2}.fq.  We compress these at their 
    #   respective accession number.
    cmds <-
    paste0("xz -czf -T", meta$n_jobs, " ",
           tmpfile, c("", "_errFree"),
           rep(1:2, each=2), ".fq > ",
           meta$artdir, "/", meta$sequencing.technology,
           "/D=", sequencing.depth,
           "/",
           paste0(accession, c("", "_errFree")),
           ".R", rep(1:2, each=2), ".fq.xz")
    for (cmd in cmds) {
#       print(cmd)
      system(cmd)
    }
#   }
}

# system(paste("ls -lat", meta$artdir, "| head -n 15"))

### Remove the temporary files we already know exists
cmd <- paste0("rm ", meta$artdir, "/tmp", "{,*.fq,*.sam}")
system(cmd)
# print(cmd)

# }}}

print("ok")
