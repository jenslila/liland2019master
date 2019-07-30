#' @name ExpectedKmers
#' @aliases ExpectedKmers
#' @title Expected K-mers
#' 
#' @description Expected overlapping words of length K in
#' DNA/RNA sequences.
#' 
#' @param K Word length (integer).
#' @param nucs (optional) Named character vector
#' indicating how nucleotides are paired. (text)
#' 
#' @details Given a integer word length K and named
#' character vector indicating how nucleotides are paired,
#' a character vector of expected nucleotides are returned.
#' 
#' @return A named character vector of all possible K-mers.
#' 
#' @author Jens Rasmus Liland and Lars Snipen.
#' 
#' @seealso \code{\link{NamedPalindromicKmers}},
#' \code{\link{PalindromicCanonicalKmers}},
#' \code{\link{PalindromicCanonicalKmerCount}}.
#' 
#' @examples 
#' ExpectedKmers(K=2)
#' 
#' @export ExpectedKmers
#' 
ExpectedKmers <-
function(K,  # Word length (integer).
         # How nucleotides are paired.
         nucs=c("A"="T", "C"="G", "G"="C", "T"="A")) {
  sort(apply(expand.grid(rep(list(nucs), K)),
             1, paste, collapse = ""))
}

#' @name NamedPalindromicKmers
#' @aliases NamedPalindromicKmers
#' @title Named Palindromic K-mers
#' 
#' @description Expected palindromic overlapping words of
#' length K in DNA/RNA sequences.
#' 
#' @param K Word length (integer).
#' @param nucs (optional) Named character vector indicating
#' how nucleotides are paired. (text)
#' 
#' @details Given a integer word length K and named
#' character vector indicating how nucleotides are paired,
#' a named character vector of palindromic K-mers are
#' returned.
#' 
#' @return A named character vector of palindromic K-mers.
#' 
#' @author Jens Rasmus Liland and Lars Snipen.
#' 
#' @seealso \code{\link{ExpectedKmers}},
#' \code{\link{PalindromicCanonicalKmers}},
#' \code{\link{PalindromicCanonicalKmerCount}}.
#' 
#' @examples 
#' NamedPalindromicKmers(K=2)
#' 
#' @export NamedPalindromicKmers
#' 
NamedPalindromicKmers <-
function(K,  # Word length (integer).
         # How nucleotides are paired.
         nucs=c("A"="T", "C"="G", "G"="C", "T"="A")) {
  # Returns a named character vector of palindromic kmers, 
  # given word length and nucleotides

  k <- ExpectedKmers(K=K, nucs=nucs)

  # Takes in a vector of K-mers, kmer.vector, built up of 
  # nucleotides listed in nucs vector, and gives out 
  # reverse compliment sequences
  r <- sapply(lapply(lapply(strsplit(k, NULL), rev),
                     function(x) { return(nucs[x]) }),
              paste, collapse = "")

  # These are the named palindromic kmers
  mapply(function(y, z) { sort(unname(c(y, z)))[1] },
         r, k)
}

#' @name PalindromicCanonicalKmers
#' @aliases PalindromicCanonicalKmers
#' @title Palindromic Canonical K-mers
#' 
#' @description Expected palindromic canonical
#' overlapping words of length K in DNA/RNA sequences.
#' 
#' @param K Word length (integer).
#' 
#' @details Given a integer word length K and named
#' character vector indicating how nucleotides are paired,
#' a named character vector of palindromic K-mers are
#' returned; which might serve as the column names of
#' expected kmers to count in a count matrix
#' 
#' @return A named character vector of palindromic
#' canonical K-mers.
#' 
#' @author Jens Rasmus Liland and Lars Snipen.
#' 
#' @seealso \code{\link{ExpectedKmers}},
#' \code{\link{NamedPalindromicKmers}},
#' \code{\link{PalindromicCanonicalKmerCount}}.
#' 
#' @examples 
#' PalindromicCanonicalKmers(K=2)
#' 
#' @export PalindromicCanonicalKmers
#' 
PalindromicCanonicalKmers <-
function(K) {  # Word length (integer).
  sort(unique(unname(NamedPalindromicKmers(K))))
}

#' @name PalindromicCanonicalKmerCount
#' @aliases PalindromicCanonicalKmerCount
#' @title Palindromic Canonical K-mer Counter
#' 
#' @description Expected palindromic canonical overlapping
#' words of length K in DNA/RNA sequences.
#' 
#' @param seq Genomic sequence (text).
#' @param K Word length (integer).
#' @param npk (optional) Vector of named palindromic
#' K-mers (text).
#' 
#' @details Given a genomic sequence, integer word length,
#' and optional character vector of palindromic K-mers, a 
#' named numeric vector of the same K-mers jellyfish version 
#' 2.2.6 would have returned, e.g. for K=2, equivalent to 
#' the output in seq.txt of
#' 
#'   $ jellyfish count --mer-len 2 --canonical \
#'       --output jelly.cnt seq.fasta
#'   $ jellyfish dump --column --tab \
#'       --output seq.txt jelly.cnt
#' 
#' @return A named numeric vector of palindromic canonical
#' K-mer counts.
#' 
#' @author Jens Rasmus Liland and Lars Snipen.
#' 
#' @seealso \code{\link{ExpectedKmers}},
#' \code{\link{NamedPalindromicKmers}},
#' \code{\link{PalindromicCanonicalKmers}}.
#' 
#' @examples 
#' PalindromicCanonicalKmerCount(seq='ATTAG', K=2)
#' 
#' @export PalindromicCanonicalKmerCount
#' 
PalindromicCanonicalKmerCount <-
function(seq,  # A single sequence (text).
         K,    # Word length (integer).
         # Precomputed NamedPalindromicKmers vector
         npk = NamedPalindromicKmers(K=K)) {
  cnt <- microclass::KmerCount(seq=seq, K=K, col.names=T)[1,]
  idx <- npk != names(npk)
  canonical <- cnt[names(npk[!idx])]
  canonical[npk[idx]] <- cnt[npk[idx]] + cnt[names(npk[idx])]
  canonical
}
ReverseComplement <-
function(seq,
         nucs=c("A"="T", "C"="G", "G"="C", "T"="A")) {
  r <- nucs[rev(strsplit(seq, NULL)[[1]])]
  paste(r[!is.na(r)], collapse="")
}
KmerCount.matrix.to.canonical.matrix <-
function(m,
         K,
         npk = NamedPalindromicKmers(K=K),
         kmers=PalindromicCanonicalKmers(K=K)) {
  idx <- unname(npk) != names(npk)
  m[,npk[idx]] <- m[,npk[idx]] + m[,names(npk)[idx]]
  m[,kmers]
}
