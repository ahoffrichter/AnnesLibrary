my_zscale <- function(x){
  M <- rowMeans(x, na.rm = TRUE)
  nsamples <- ncol(x)
  DF <- nsamples - 1L
  IsNA <- is.na(x)
  if (any(IsNA)) {
    mode(IsNA) <- "integer"
    DF <- DF - rowSums(IsNA)
    DF[DF == 0L] <- 1L
  }
  x <- x - M
  V <- rowSums(x^2L, na.rm = TRUE)/DF
  x <- x/sqrt(V + 0.01)
  out <- melt(x)
  out
}
