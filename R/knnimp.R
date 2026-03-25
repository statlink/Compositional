knnimp <- function(x, k = 2) {

  p <- dim(x)[2]  ## number of components
  a <- Rfast::rowsums(x)
  pou <- which( is.na(a) )  ## where are the missing values
  x1 <- x[-pou, , drop = FALSE]  ## rows with good data, i.e. no missing values
  ina <- 1:dim(x1)[1]
  x2 <- x[pou, , drop = FALSE]  ## rows with missing values (bad data)
  tot <- 1 - Rfast::rowsums(x2, na.rm = TRUE)  ## 1 - sum of the rows of the bad data

  if ( length(k) == 1 ) {
    for ( i in 1:length(pou) ) {  ## go to each row where missing values are
      id <- which( !is.na(x2[i, ]) )  ## find the good components that do not contain missing values
      w1 <- x1[, id, drop = FALSE] / Rfast::rowsums( x1[, id, drop = FALSE] )  ## choose the good components and normalise the good data
      w2 <- x2[i, id, drop = FALSE] / sum(x2[i, id, drop = FALSE])  ## choose the bad components and normalise the bad data
      #dis <- as.vector( Rfast::dista(w2, w1, type = "jensen_shannon") )  ## ESOV distance between good and bad bad data
      #dis <- cbind(dis, ina)
      #dis <- dis[order(dis[, 1]), ]  ## sort the ESOV distances
      #nn <- dis[1:k, 2]  ## choose the k smallest ones
      nn <- Rnanoflann::nn(w1, w2, k = k, method = "jensen_shannon", sorted = TRUE)$indices
      mu <- Rfast::colmeans( x1[nn, , drop = FALSE] )  ## take the average of the good data with the k smallest distances
      mu <- replace(mu, id, NA)  ## put NA in the place of the good components of the good data and keep the bad components only
      est <- mu / sum(mu, na.rm = TRUE) * tot[i]  ## normalize the above data and then multiply by the sum of the bad data
      x[pou[i], -id] <- est[-id]  ## replace the bad components in bad data with the imputed values
    }  ##  end for ( i in 1:length(pou) ) {
    y <- list()
    y[[ 1 ]] <- x
    res <- list(pou = pou, x = y)

  } else {

    names <- paste("kappa=", k, sep = "")
    y1 <- sapply(names, function(x) NULL)
    y <- y1
    z <- matrix(0, length(pou), p)
    for ( i in 1:length(k) )  {
      y1[[ i ]] <- z
      y[[ i ]] <- x
    }
    for ( i in 1:length(pou) ) {
      id <- which( !is.na(x2[i, ]) )
      w1 <- x1[, id, drop = FALSE] / Rfast::rowsums( x1[, id, drop = FALSE] )
      w2 <- x2[i, id, drop = FALSE] / sum(x2[i, id], drop = FALSE)
      NN <- Rnanoflann::nn(w1, w2, k = nrow(w1), method = "jensen_shannon", sorted = TRUE)$indices
      for ( j in 1:length(k) ) {
        #nn <- dis[1:k[j], 2]
        nn <- NN[ 1:k[j] ]
        mu <- Rfast::colmeans( x1[nn, , drop = FALSE] )
        mu <- replace(mu, id, NA)
        est <- mu / sum(mu, na.rm = TRUE) * tot[i]
        y[[ j ]][pou[i], -id] <- est[-id]
      }  ##  end  for ( j in 1:length(k) ) {
    }  ##  end for ( i in 1:length(pou) ) {
    res <- list(pou = pou, x = y)

  }  ## end if ( length(k) == 1 ) {
  res
}
