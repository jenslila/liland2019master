# Programming appendix for Jens Rasmus Liland's 2019 master thesis *Recognizing plasmid-reads by machine learning and K-mer statistics*

**What's this?**  A pipeline of `R` and `python36` 
scripts to download, process, simulate reads from, and 
classify *E.coli* plasmid and chromosomal sequences for 
the 2019 master thesis of Rasmus Liland.

### Required `R` packages

R version 3.5.1.

* `micropan` 1.2
* `microclass` 1.1
* `data.table` 1.11.4

### Required `python3.6` packages

Python version 3.6.8.

* `numpy` 1.16.2
* `pandas` 0.24.2
* `scikit-learn` 0.20.2

Also listed out in `requirements.txt`.

*Note, August 2018:*  For the wide datasets, K>5, 
greater than 1024 columns, traditional 
[BLAS](http://www.netlib.org/lapack) seems to work with 
sklearn at least, so unfortunately 
[OpenBLAS](http://www.openblas.net/) can not be used as 
a substitute at the moment; BUT, R seems to be able to 
use both!

### Pipeline scripts

To be run one after the other.

1. **`script_1_mining.R`**  R code for mining the NCBI 
   Prokaryote database, which sources `miningfun.R` 
   library.
2. **s`cript_2_xtrain.R`**  R script for counting 
   K-mers of fully assembled sequences to create 
   training data sets.  It sources `palindromicfun.R`. 
3. **`script_3_art.R`**  R script for *in silico* 
   simulation of Illumina *Hiseq 2500* sequences using 
   ART art\_Illumina.  Is designed to be run in 
   parallell instances, reading from designated 
   portions of the indexed established in
   `mining.R`.  It sources `artfun.R` and `knnfun.R`.
4. **`script_ds_1_clean.genomes.R`**  R script for 
   cleaning genomes before running the 
   DeepSimulator.  It is essentially a shell script 
   calling `perl` and `awk`, which functions like they 
   should in a way that base R functions does not if 
   you do not branch out to using plyr functions.   It 
   sources `knnfun.R`.
5. **`script_ds_2_location.R`**  R script for indexing 
   genomes before running the DeepSimulator.
6. **`script_ds_3_xtest.R`**  R script for counting 
   K-mers on the Nanopore sequences.  Is designed to be 
   run as multiple scripts in paralell, reading 
   different chunks of the indexed established 
   previously.  It sources `palindromicfun.R` and 
   `knnfun.R`.
7. **`script_ds_4_skl.py`**  Python script for 
   classifying the K-mer profiles of the simulated 
   nanopore sequences using
   scikit-learn.  1:1 proportions of plasmid and 
   chromosomal reads are selected to make sure the 
   binary classification accuracy will never drop below 
   50%.
8. **`script_4.1_indecies.R`**  R script for sampling 
   indecies of simulated Illumina reads.  It sources 
   `palindromicfun.R` and `knnfun.R`.
9. **`script_4.2_artlocation.R`**  R script for 
   establishing an index of what has been simulated, 
   has sequences to be selected from, etc.
10. **`script_4.3_xtest.R`**  R script for counting 
    K-mer occurrence profiles of the sampled Illumina 
    sequences which exists based on checks done in 
    `script_4.1_indecies.R` and 
    `script_4.2_artlocation.R`.  Like previous scripts, 
    this one is also designed to be run in multiple 
    instances doing their designated tasks.
11. **`script_5_sklearn_train.py`**  Python script for 
    classifying the fully assembled sequences 
    (cross-validated training data) using 
    scikit-learn.
12. **`script_5_sklearn_art_partial.py`**  Python 
    script for classifying the K-mer occurrence 
    profiles of the *in silico* simulated Illumina 
    reads using scikit-learn.  1:1 proportions of 
    plasmid and chromosomal reads are selected to make 
    sure the binary classification accuracy will never 
    drop below 50%.
13. **`script_6_scores.py`**  Python script for 
    summarizing the three classification routines by 
    creating stacked confusion metrics for doing ANOVA 
    analyses later on.

Libraries `source`d by the R scripts:

* **`miningfun.R`**  R functions for mining.
* **`palindromicfun.R`**  R functions for K-mer counting.
* **`artfun.R`**  R functions for running `art_Illumina`.
* **`knnfun.R`**  R functions for running K-nearest neighbour algorithm.
