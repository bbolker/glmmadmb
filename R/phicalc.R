## recreate orthogonalization (phi) calculation from
##  glmmadmb.tpl, so we can retrieve real_beta from beta
##  in MCMC results ...
make_phi <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  phi <- matrix(0,p,p)
  diag(phi) <- 1
  TX <- t(X)
  for (i in seq(p)) {
    tmp <- sqrt(sum(TX[i,]^2))
    ## normalize 
    TX[i,] <- TX[i,]/tmp
    phi[i,] <- phi[i,]/tmp
    for (j in (i+1):p) {
      if (j>p) break
        a <- TX[j,] %*% TX[i,]
        TX[j,] <- TX[j,]-a*TX[i,]
        phi[j,] <- phi[j,]-a*phi[i,]
    }
  }
  phi
}

