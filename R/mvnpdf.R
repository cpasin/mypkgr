#' Fonction mvnpdf
#'
#' calcule la valeur de la densité d’une loi normale multivariée sur Rp en n points
#'
#' @param x : matrice d'observations en colonnes (n colonnes, p lignes)
#' @param mean : vecteur de moyenne
#' @param varcovM : matrice de variance covariance
#' @param Log : TRUE par défaut (calcul du log)
#'
#' @return matrice x et vecteur y de la valeur de la densité en chacun des vecteurs colonne de x
#' @export
#'
#' @examples
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))
  DetvarcoM<-det(varcovM)
  ## test changement
  ## commentaire
  ## test encore
  ## coucou Emilie
  ##hello darling, how are you today ?
  ## can I have a cup of tea please ?

  ## changement travis


  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  res<-(list(x=x,y=y))
  class(res)<-"mvnpdf"
  return(res)
}
