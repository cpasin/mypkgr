mvnpdftest<-function(x,mean,varcovM,Log=TRUE){
  p<-nrow(x)
  if (nrow(x)!=length(mean) | nrow(varcovM)!=ncol(varcovM) | nrow(x)!=nrow(varcovM) | length(mean)!=nrow(varcovM)){
    cat ("dimension issue")
  }else if (det(varcovM)==0){
    cat ("matrice non inversible")
  }else if (Log==FALSE){
    y<-1/((2*pi)^(p/2)*sqrt(det(varcovM)))*exp((-t(x-mean)%*%solve(varcovM)%*%(x-mean)))/2
  }else{
    y<-log(1/((2*pi)^(p/2)*sqrt(det(varcovM)))*exp((-t(x-mean)%*%solve(varcovM)%*%(x-mean)))/2)
  }
  out<-list(x,y)
  return (out)
}


mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))

  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  return(list(x=x,y=y))
}


### charger le package
devtools::load_all()


x<-matrix(rep(c(1,2),3),nrow=3,byrow = T)
mean<-rep(0,3)
varcovM<-diag(nrow(x))
