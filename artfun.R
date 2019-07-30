run.ART_Illumina.and.perhaps.samtools <-
function(# Wrapper function to run ART_Illumina version 
         # 2.5.8 and samtools version 1.9 on the system 
         # given a ton of knobs and switches listed 
         # here below.  These parameters are just a 
         # subset of all possible parameters of ART.
         # 
         # 
         # Integer.  Sets random seed, then randomness 
         # is toggled off, to turn randomness on again 
         # just paste in UNIX integer timestamp in 
         # here, or some other ever-changing integer 
         # you know of.
         random.seed,
         # String.  Which sequencing machine to use.  
         # E.g. HS25 for Hiseq 2500, MSv3 for MiSeq 
         # v3.
         sequencing.technology,
         # Integer.  Insert size, I, i.e.  expected 
         # value of the fragment length normal 
         # distribution.
         insert.size,
         # Integer.  Standard deviation in the fragment 
         # length normal distribution.
         insert.size.std,
         # Integer.  Sequencing depth is the D in the 
         # D=NL/G equation.
         sequencing.depth,
         # Integer.  Max read length on the HS25 
         # machine is 150bp.
         read.length,
         # String.  FASTA file with raw sequence.
         infile.path,
         # String.  Read ID prefix.
         read.ID.prefix,
         # String.  Prefix for FASTQ files.
         fastq.prefix,
         # Boolean.  Indicates whether a SAM file is 
         # output for the regular (error containing) 
         # reads, and after that running samtools fastq 
         # on the SAM file to yield fastq files like 
         # the regular (error containing) ones 
         # outputted by ART.
         samout,
         # Boolean.  Indicates whether a SAM file of 
         # error free reads is also going to be 
         # created. 
         errfree,
         # String.  Full or relative path to 
         # art_illumina binary, whatever suit your 
         # needs.
         art
         ) {
  if (errfree) { errfree.paste <- "--errfree" }
  else         { errfree.paste <- "" }

  if (samout)  { samout.paste <- "--samout" }
  else         { samout.paste <- "" }

  cmd <-
  paste(art,
        "--rndSeed", random.seed,
        "--paired",
        "--seqSys", sequencing.technology,
        "--noALN",
        samout.paste,
        errfree.paste,
        "--mflen", insert.size,
        "--sdev", insert.size.std,
        "--fcov", sequencing.depth,
        "--len", read.length,
        "--in",  infile.path,
        "--id", read.ID.prefix,
        "--out", fastq.prefix)
#   print(cmd)
  output <- system(cmd, intern=TRUE)
  cat(paste(output, collapse="\n"), "\n")

  if (errfree) {
    errname <- "_errFree"
    cmd <- paste("samtools fastq",
                 "-1", paste0(fastq.prefix, errname, "1.fq"),
                 "-2", paste0(fastq.prefix, errname, "2.fq"),
                 paste0(fastq.prefix, errname, ".sam"))
    ### samtools writes out a bunch of meaningsless 
    #   crap to stderr, so we set ignore.stderr here.
    output <- system(cmd, ignore.stderr=TRUE, intern=TRUE)
    cat(paste(output, collapse="\n"), "\n")
#     print(cmd)
  }
}
