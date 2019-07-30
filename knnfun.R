string.chunks <- function(txt,  # y, response vector 
                          n) {  # C, #chunks 
  #### Given response vector of unique strings, and 
  #    number of chunks to be sliced from this, this 
  #    function determines indecies for data groups for 
  #    doing unstratified K-fold cross validation. 
  df <- data.frame(c(), c())
  len <- length(txt)
  avg <- len/n
  last <- 1
  while (last < len) {
    df <- rbind(df, c(ceiling(last),
                      ceiling(last + avg - 1)))
    last <- last + avg
  }
  colnames(df) <- c("start", "stop")
  return(df)
}
genome.chunks <- function(location) {
  #### Given a dataframe with a character vector column 
  #    named system of sorted repeated character 
  #    elements, the function returns the start and 
  #    stop indecies of the character elements, so as 
  #    the number of returned index pairs is equal to 
  #    the number of unique strings within the column. 
  t(as.data.frame(sapply(unique(location$system), function(x, idx) {
    i <- idx[location$system == x]
    c("start"=i[1], "stop"=i[length(i)])
  }, idx=1:nrow(location))))
}
computeD <- function(q,   # row from train matrix Q
                     p) { # row from test matrix P
  a <- abs(log2(p/q))
  return(sum(a[!is.na(a)]))
}
min.idx <- function(Q,    # train matrix Q
                    FUN,  # distance routine
                    p) {  # row from test matrix P
  return(apply(X = Q, MARGIN = 1, FUN = FUN, p = p))
}
consensus <- function(p,  # row from test matrix, P
                      Q,  # train matrix, Q
                      k,  # K of KNN
                      y,  # response vector
                      dFUN) { # distance routine
  scores <- apply(Q, 1, dFUN, p = p)
  neighbours <- y[order(scores, decreasing = F)[1:k]]
  answer <- names(which.max(table(neighbours)))
  return(answer)
}
consensus.multi <- function(P,     # test matrix
                            Q,     # train matrix
                            K,     # K of KNN
                            dFUN,  # distance routine
                            y) {   # response vector
    return(apply(X = P, MARGIN = 1,
                 FUN = consensus,
                 Q = Q, K = K, dFUN = dFUN,
                 y = y))
}
