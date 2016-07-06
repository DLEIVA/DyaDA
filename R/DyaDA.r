
#######################################################################################
#                 ADDITIONAL R FUNCTIONS USED IN DYADA PACKAGE                        #
#######################################################################################

# Function to detect square matrices #

squaremat <- function (X){
  sqmat <- (nrow(X) == ncol(X));
  return(sqmat)
}

# Function to prepare matrices for rank orders without self-addressed cells #

prepmat <- function(X){
  dyad <- ncol(X)-1;
  X1 <- matrix(t(X)[row(X)!=col(X)],nrow=nrow(X),ncol=dyad,byrow=T);
  return(X1);
}


#######################################################################################
#                 R FUNCTIONS FOR ANALYZING LINEAR HIERARCHY IN GROUPS                #
#######################################################################################


# Function to transform original sociomatrix into a matrix of abilities #

freq2abil <- function (X,method=c('landau','matman','schmid')){
  method <- match.arg(method)
  if (method=='landau'){
    M1 <-  matrix(1*as.numeric(X>t(X)),nrow=nrow(X),byrow=F)
    M2 <-  matrix(-1*as.numeric(X<t(X)),nrow=nrow(X),byrow=F)
    M3 <-  matrix(0.*as.numeric(X==t(X)),nrow=nrow(X),byrow=T)
    V <- M1 + M2 + M3
  }

  if (method=='matman'){
    M1 <-  matrix(1*as.numeric(X>t(X)),nrow=nrow(X),byrow=F)
    M2 <-  matrix(1./2.*as.numeric(X==t(X)),nrow=nrow(X),byrow=T)
    V <- M1 +M2
  }
  
  if (method=='schmid'){
    M1 <-  matrix(1*as.numeric(X>t(X)),nrow=nrow(X),byrow=F)
    M2 <-  matrix(-1*as.numeric(X<t(X)),nrow=nrow(X),byrow=F)
    M3 <-  matrix(0.*as.numeric(X==t(X) & X+t(X)==0),nrow=nrow(X),byrow=T)
    V <- M1 + M2 + M3
  }
  diag(V) <- 0.
  dimnames(V) <- dimnames(X)
  return(V)
}

# Function to obtain the number of tied relationships #

gettied <- function(X,relfreq=TRUE){
  dyad <- nrow(X)*(nrow(X)-1)/2;
  tied <- sum(X==t(X) & X!=0.)/2
  if (relfreq == TRUE) list(tied.dyads=tied, perc.tied.dyads=(tied/dyad*100))
  else return(tied)  
}

# Function to obtain the number of unknown dyadic relationships #

getunknown <- function(X,relfreq=TRUE){
  dyad <- nrow(X)*(nrow(X)-1)/2;
  unknown <- (sum(X+t(X)==0.)-nrow(X))/2;
  if (relfreq == TRUE) list(unknown.dyads=unknown,
    perc.unknown.dyads=(unknown/dyad*100))
  else return(unknown)  
}

# Function to obtain the number of dyads with one-way relationships #

getonewayrel <- function (X, relfreq=TRUE){
  dyad <- nrow(X)*(nrow(X)-1)/2;
  c <- X + t(X);
  onewayrel <- sum(X-c==0. & c!=0)
  if (relfreq == TRUE) list(oneway.dyads=onewayrel, 
                            perc.oneway.dyads=(onewayrel/dyad*100))
  else return(onewayrel)
}

# Function to obtain the number of dyads with two-way relationships #

gettwowayrel <- function (X, relfreq=TRUE){
  dyad <- nrow(X)*(nrow(X)-1)/2;
  c <- X + t(X);
  twowayrel <- sum((X-c!= 0.) * (t(X) - c != 0.))/2
  if (relfreq == TRUE) list(twoway.dyads=twowayrel, perc.twoway.dyads=
                              (twowayrel/dyad*100))
  else return(twowayrel)
}

getsignificant <- function(X,alpha=0.05, relfreq=TRUE){
  dyad <- nrow(X)*(nrow(X)-1)/2;
  pvals <- mapply(x=X[row(X)!=col(X) & X+t(X)!=0 & !is.na(X+t(X))],
                  n=(X+t(X))[row(X)!=col(X) & X+t(X)!=0& !is.na(X+t(X))],
                  FUN=function(x,n) binom.test(x,n)$p.value)
  binommat <- matrix(NA,nrow=nrow(X),ncol(X))
  binommat[row(X)!=col(X) & X+t(X)!=0 & !is.na(X+t(X))]<-round(pvals,3)
  diag(binommat) <- 0
  binommat[upper.tri(binommat)] <- (X+t(X))[upper.tri(X)]
  ndyads <- sum(binommat[lower.tri(binommat)]<=alpha,na.rm=TRUE)
  if (!is.null(dimnames(X))) dimnames(binommat) <- dimnames(X)
  if (relfreq == TRUE) res<-list(original.mat=X,binomtest.matrix=binommat,significant.dyads=ndyads,
                            perc.significant.dyads=(ndyads/dyad*100))
  else res <- list(original.mat=X,binomtest.matrix=binommat,significant.dyads=ndyads)
  class(res) <- "sigdyads"
  res
}

plot.sigdyads <- function(x){
  #library(cowplot)
  #library(reshape2)
  #library(ggplot2)
  X <- x$original.mat
  pvals <- mapply(x=X[row(X)!=col(X) & X+t(X)!=0 & !is.na(X+t(X))],
                  n=(X+t(X))[row(X)!=col(X) & X+t(X)!=0 & !is.na(X+t(X))],
                  FUN=function(x,n) binom.test(x,n)$p.value)
  temp <- matrix(NA,nrow=nrow(X),ncol(X))
  temp[row(X)!=col(X) & X+t(X)!=0 & !is.na(X+t(X))]<-pvals
  myfun<-function(x){if (!is.na(x) & x >=0.01){res <-round(x,2)} else if(!is.na(x) & x <0.01)
    {res<-'<0.01'} else {res<-'NA'}; res}
  pvals2<- sapply(t(temp)[upper.tri(temp)],function(x) myfun(x))
  binommat <- matrix(NA,nrow=nrow(X),ncol(X))
  binommat[row(X)!=col(X) & X+t(X)!=0 & !is.na(X+t(X))]<-pvals
  binommat[row(binommat)>col(binommat)] <- (X+t(X))[lower.tri(X)]
  dyadn <- paste(t(X),'/',X+t(X),sep='')[lower.tri(X)]

  p <- ggplot(subset(melt(binommat)),
              aes(x=Var1,y=Var2))
  mi.ids <- subset(melt(binommat), Var1 == Var2)
  mi.low <- subset(melt(binommat), Var1 > Var2)
  mi.upp <- subset(melt(binommat), Var1 < Var2)
  p1 <- (p +  geom_tile(data=mi.ids, colour='black',alpha=0.2) + 
           geom_text(data=mi.ids,aes(label=if (is.null(rownames(X))){paste("Ind.",Var1)} else{
             rownames(X)}),colour='red',fontface=2)+
           geom_tile(data=mi.low,colour='black',alpha=0.1)+
           geom_text(data=mi.low,aes(label=dyadn))+
           geom_tile(data=mi.upp,aes(fill = value),colour = "black") +
           geom_text(data=mi.upp,aes(label = pvals2)) +
           scale_fill_gradient2(name="P-value",midpoint=0.5,low = "green",
           mid='orange', high = "red",na.value = 'grey') +
            xlab(NULL)+ylab(NULL)) + scale_y_reverse()
  p1 + theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank())
}


# Function to obtain several win and loss measures at individual level #

getwinloss <- function(X,names=NULL,method=c("Dij","Pij")){
  if (!squaremat(X))
    return("Error: Matrix X is not square and cannot be analyzed");
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric");
  method <- match.arg(method)
  dyadc <- X + t(X);
  if (method == "Dij"){
    Dij <- X/dyadc-(((X/dyadc)-0.5)/(dyadc+1));
    Dij[is.nan(Dij)] <- 0.;
    w1 <- rowSums(Dij);
    w2 <- Dij%*%w1;
    l1 <-colSums(Dij);
    l2 <- t(l1)%*%Dij;
  }
  if (method == "Pij"){
    Pij <- array(dim=c(nrow(X),ncol(X)),0.);
    Pij <- X/dyadc;
    Pij[is.nan(Pij)] <- 0.
    w1 <- rowSums(Pij);
    w2 <- Pij%*%w1;
    l1 <-colSums(Pij);
    l2 <- t(l1)%*%Pij;
  }
  wl<- data.frame(w=array(w1,dim=c(nrow(X),1)),weighted.w=w2,l=array(l1,dim=c(nrow(X),1)),
                  weighted.l=t(l2))
  if (is.null(names)) names <- paste('Ind.',1:nrow(X))
  rownames(wl) <- names
  return(wl)
}

# Function to obtain the David's scores (DS) based on dyadic dominance indices #

getDS <- function(X,names=NULL,method=c("Dij","Pij")){
  if (!squaremat(X))
    return("Error: Matrix X is not square and cannot be analyzed")
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric")
  method <- match.arg(method)
  dyadc <- X + t(X)
  if (method == "Dij"){
    Dij <- array(dim=c(nrow(X),ncol(X)),0.)
    Dij <- X/dyadc-(((X/dyadc)-0.5)/(dyadc+1))
    Dij[is.nan(Dij)] <- 0.
    w1 <- rowSums(Dij)
    w2 <- Dij%*%w1
    l1 <-colSums(Dij)
    l2 <- t(l1)%*%Dij
  }
  if (method == "Pij"){
    Pij <- array(dim=c(nrow(X),ncol(X)),0.)
    Pij <- X/dyadc
    Pij[is.nan(Pij)] <- 0.
    w1 <- rowSums(Pij)
    w2 <- Pij%*%w1
    l1 <-colSums(Pij)
    l2 <- t(l1)%*%Pij
  }
  DS <- w1 + w2 - l1 - t(l2);
  if (is.null(names)) names <- paste('Ind.',1:nrow(X))
  DS <- array(DS,dim=c(nrow(X),1),dimnames=c(list(names),"David's Scores"))
  return(DS)
}


# Function to obtain Landau's linear hierarchy index #

getlandau <- function(X){
  V <- freq2abil(X,'matman');
  Vvector <- rowSums(V)
  sumabil <- sum((Vvector - (length(Vvector) - 1)/2)^2)
  landauh <- (12/(nrow(X)^3-nrow(X)))*sumabil;
  return(landauh)
}

# Function to obtain improved Landau's linear hierarchy index #

getimplandau <- function(X){
  landauh <- getlandau(X);
  unknown <- getunknown(X)$unknown.dyads;
  improved.landauh <- landauh + (unknown*6/(nrow(X)^3-nrow(X)));
  return(improved.landauh)
}

# Function to obtain expected value of linear hierarchy index #

getexpech <- function(X){
  expectedh <- 3/(nrow(X)+1);
  return(expectedh)
}

# Function to obtain variance for linear hierarchy index #

getvarh <- function(X){
  varianceh <- 18/nrow(X)^3;
  return(varianceh)
}

# Function to compute the number of circular triads #

getcirc <- function (X){
  S <- freq2abil(X,'matman');
  Svector <- rowSums(S)
  sumabil <- sum(Svector^2)
  circular <- length(Svector)*(length(Svector)-1)*(2*length(Svector)-1)/
    12-(1/2*sumabil);
  return(circular)
}

# Function to compute the expected number of circular triads #

getexpeccirc <- function(X){
  expectedcirc <- nrow(X)*(nrow(X)-1)*(nrow(X)-2)/24;
  return(expectedcirc)
}

# Function to obtain the maximum number of circular triads #

getmaxcirc <- function(X){
  if (nrow(X) %% 2 == 0.) maxcirc <- 1/24*(nrow(X)^3-4*nrow(X))
    else maxcirc <- 1/24*(nrow(X)^3-nrow(X));
  return(maxcirc)
}

# Function to obtain Kendall's linearity index #

getkendall <- function(X){
  circular <- getcirc(X);
  maxcirc <- getmaxcirc(X);
  K <- 1 - circular/maxcirc;
  return(K)
}

# Function to obtain statistical significance for Kendall's linearity index #

getpvkendall <- function(X){
  circular <- getcirc(X);
  if (nrow(X) >= 10){
    degfr <-  nrow(X)*(nrow(X)-1)*(nrow(X)-2)/(nrow(X)-4)^2;
    chisq <- 8/(nrow(X)-4)*(nrow(X)*(nrow(X)-1)*(nrow(X)-2)/24-circular+1/2)+
             degfr;
    pval <- 1-pchisq(chisq,degfr)
    list(Chisq=chisq,df=degfr,pvalue=pval)
  }
  else cat("WARNING: ","\n","Matrix size < 10: ","\n",
         "Chi square test has not been carried out since it would 
          overestimate statistical significance","\n")
}


# Function to carry out the two-step randomization test for #
# testing linear hierarchy #

linear.hierarchy.test <- function (X,rep=9999){

# Is the matrix X square? #

  if (!squaremat(X))
    return("Error: Matrix X is not square and cannot be analyzed")
   
# Limit the number of replications #

  if ((rep < 1) | (rep > 1000000))
    return("Error: Number of replications must be between 1 and 1000000") 

# Matrix have to be transformed into vectors in order to pass to a .C routine #
 
vecX <- c(t(X));

# R wrapper of the .C routine #

out <- .C("linearh",
          as.double(vecX),
          as.integer(nrow(X)),
	        as.integer(rep),
          pvhright=double(1),
          pvhleft=double(1),
          PACKAGE="DyaDA")

h <- getlandau(X);
himp <- getimplandau(X);         
pvright <- out$pvhright;
pvleft <- out$pvhleft;
exph <- getexpech(X);
varh <- getvarh(X);
nind<-nrow(X);
dyads <- nind*(nind-1)/2
total<-sum(X);
circ <- getcirc(X);
expcirc <- getexpeccirc(X);
maxcirc <- getmaxcirc(X);
rho <- getkendall(X);
rho.test <- if (nrow(X)>=10){getpvkendall(X)} else {NULL};
desctable <- rbind(unlist(getunknown(X)),unlist(gettied(X)),
      unlist(getonewayrel(X)),unlist(gettwowayrel(X)),
      unlist(getsignificant(X)[3:4]))
colnames(desctable) <- c('Frequency','Percentage')
rownames(desctable) <- c('Unknown','Tied','One-way', 'Two-way', 'Significant')
res<-list(nind=nind,dyads=dyads,total=total,desctable=desctable,
          randomizations=rep,landauh=h,improved.landauh=himp,
          expected.h=exph,variance.h=varh,
          right.pvalue=pvright,left.pvalue=pvleft,circ=circ,
          expcirc=expcirc,maxcirc=maxcirc,rho=rho,rho.test=rho.test)
class(res) <- "hlineartest"
res
}

print.hlineartest <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Linearity test: Landau's h \n")
  cat("========================== \n\n")
  cat("\n Number of individuals:", x$nind," (",x$dyads,"dyads)\n")
  cat(" Total number of interactions:", x$total,"\n\n")
  cat(" Descriptive table of dyadic relationships: \n")
  print(x$desctable,digits)
  cat("\n")
  cat(" Landau's h index:", round(x$landauh,digits),"\n")
  cat(" E(h):", round(x$expected.h,digits),"\n")
  cat(" Var(h):", round(x$variance.h,digits),"\n")
  cat(" Improved Landau's h:", round(x$improved.landauh,digits),"\n")
  cat(" Number of randomizations:", x$randomizations,"\n")
  cat(" P-value (right):", round(x$right.pvalue,digits),"\n")
  cat(" P-value (left):", round(x$left.pvalue,digits),"\n")
  cat("\n\n")
  cat("Linearity test: Kendall's rho \n")
  cat("============================= \n\n")
  cat(" Number of circular dyads:", round(x$circ,digits),"\n")
  cat(" Expected number of circular dyads:", round(x$expcirc,digits),"\n")
  cat(" Maximum number of circular dyads:", round(x$maxcirc,digits),"\n")
  cat(" Kendall's rho index:", round(x$rho,digits),"\n")
  if (is.null(x$rho.test)){cat("Cannot computed since n<10.")} else{
  cat(" Chi-Squared:", round(x$rho.test[[1]],digits),"\n")
  cat(" df:", round(x$rho.test[[2]],digits),"\n")
  cat(" P-value:", round(x$rho.test[[3]],digits),"\n")
  }
  invisible(x)
}


#######################################################################################
#                                   I&SI METHOD IN R                                  #
#######################################################################################


DminS.comp <- function(X) {
  DminS <- rowSums(sign(X-t(X)))
  return(DminS)
}

ordering.matrix <- function(X,order) {
  N <- nrow(X);
  matord <- X[order,order];
  names1 <- dimnames(X)[[1]];
  dimnames(matord) <- list(c(names1[order]),c(names1[order]));
  return(matord);
}

ISI.comp <- function(X) {
  str.inc <- sum((col(X)*(X<t(X))-row(X)*(X<t(X)))[row(X)<col(X)]);
  numb.inc <- sum((X<t(X))[row(X)<col(X)]);
  list(number.inconsistencies=numb.inc,strength.inconsistencies=str.inc)
}

ISI.info <- function (X,names=NULL){
  N <- nrow(X);
  I <- sum((X<t(X))[row(X)<col(X)]);
  SI <- sum((col(X)*(X<t(X))-row(X)*(X<t(X)))[row(X)<col(X)]);
  if (is.null(names)) names <- paste('Ind.',c(1:N))
  else names <- dimnames(X)[[1]];
  info <- matrix(0,nrow=as.numeric(I),ncol=3,dimnames=list(c(1:as.numeric(I)),
          c("Individual i","Individual j","Strength Inconsistency")));
  inc.dyads <- matrix(0,nrow=as.numeric(I),ncol=2);
  inc.dyads[,1]<-(row(X)*(X<t(X)))[row(X)<col(X)][which((row(X)*(X<t(X)))[row(X)<col(X)]!=0)];
  inc.dyads[,2]<-(col(X)*(X<t(X)))[row(X)<col(X)][which((col(X)*(X<t(X)))[row(X)<col(X)]!=0)];
  info[,1]<-names[inc.dyads[,1]];
  info[,2]<-names[inc.dyads[,2]];
  info[,3]<-((col(X)*(X<t(X))-row(X)*(X<t(X)))[row(X)<col(X)])[which((col(X)*(X<t(X))-row(X)*
            (X<t(X)))[row(X)<col(X)]!=0)];
  info <- data.frame(info)[order(inc.dyads[,1]),];
  rownames(info) <- 1:nrow(inc.dyads);
  inc.dyads <- inc.dyads[order(inc.dyads[,1]),];
  list(inconsistent.dyads=inc.dyads,ISI.frame=info)
}

ISI.method <- function(X,names=NULL,tries=100){
  N <- nrow(X)  
  if (is.null(names)) names <- paste('Ind.',c(1:N))  # IF vector names is NULL rename it #
  # Matrices need to be transformed into vectors in order to be passed to a .C routine #
  DminS <- rowSums(sign(X-t(X)))
  ord <- rank(-DminS, ties.method='random')
  Xc <- X[ord,ord]
  namesord<- names[ord]
  vecX<-c(t(Xc))
  out <- .C("ISImethod",as.double(vecX),
            as.integer(N),
            as.character(namesord),
            as.integer(tries),
            matord=double(N*N),
            namord=character(N),
            PACKAGE="DyaDA")
  
  names.ordered <- out$namord
  BestXthusfar <- matrix(out$matord,nrow=nrow(X),ncol=ncol(X),byrow=TRUE,
                         dimnames=list(names.ordered,names.ordered))
  
  res <- list (call=match.call(),dataIni=X,dataFinal=BestXthusfar,initialI=ISI.comp(X)[[1]],
               initialSI=ISI.comp(X)[[2]], iniInfo=ISI.info(X)[[2]],finalI=ISI.comp(BestXthusfar)[[1]],
               finalSI=ISI.comp(BestXthusfar)[[2]],finalInfo=ISI.info(BestXthusfar)[[2]])
  class(res) <- "isi"
  res 
}

# A print method for isi object #

print.isi <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Ordering a dominance matrix: I&SI Method")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Initial order: ", rownames(x$dataIni),"\n\n")
  print.table(x$dataIni, zero.print = '.')
  cat("\n")
  cat(" Initial number of inconsistencies: ", x$initialI,"\n")
  cat(" Initial strength of incosistencies: ", x$initialI,"\n\n")
  cat(" Inconsistent dyads and strength of inconsistencies: ","\n")
  print(x$iniInfo,digits=digits)
  cat("\n")
  cat("**********************************************************","\n")
  cat(" Final order: ", rownames(x$dataFinal),"\n\n")
  print.table(x$dataFinal, zero.print = '.')
  cat("\n")
  cat(" Final number of inconsistencies: ", x$finalI,"\n")
  cat(" Final strength of inconsistencies: ", x$finalSI,"\n\n")
  cat(" Inconsistent dyads and strength of inconsistencies: ","\n")
  print(x$finalInfo, digits = digits)  
  cat("\n\n")
  invisible(x)
}


################################################################################
#                 R FUNCTIONS FOR  MATRIX CORRELATION METHODS                  #
################################################################################

# Function to obtain Mantel's Z statistic #

getZ <- function(X,Y){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

XY <- X*Y;
Z <- sum(XY);
return(Z);
}

# Function to compute expected value for Mantel's Z statistic #

getexpecZ <- function(X,Y){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

expecZ <- sum(X)*sum(Y)/(nrow(X)*(nrow(X)-1))
return(expecZ);
}

# Function to compute standar error for Mantel's Z statistic #

getSEZ <- function(X,Y){
  
  # Original sociomatrices should have the same size #
  
  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 
  
  if (squaremat(X) & squaremat(Y)){
    Bx <- sum(X*X);
    By <- sum(Y*Y);
    Cx <- sum(X*t(X));
    Cy <- sum(Y*t(Y));
    Dx <- sum(rowSums(X)^2);
    Dy <- sum(rowSums(Y)^2);
    Ex <- sum(colSums(X)^2);
    Ey <- sum(colSums(Y)^2);
    Fx <- c(rowSums(X)%*%colSums(X));
    Fy <- c(rowSums(Y)%*%colSums(Y));
    Gx <- sum(X)**2;
    Gy <- sum(Y)**2;
    Hx <- Dx - Bx;
    Hy <- Dy - By;
    Ix <- Ex - Bx;
    Iy <- Ey - By;
    Jx <- Fx - Cx;
    Jy <- Fy - Cy;
    Kx <- Gx - Bx - Cx - Hx - Ix - 2*Jx;
    Ky <- Gy - By - Cy - Hy - Iy - 2*Jy;
    VarZ <- ((Bx*By)+(Cx*Cy)+((Hx*Hy+Ix*Iy+2*Jx*Jy)/(nrow(X)-2))+(Kx*Ky/((nrow(X)-2)*(nrow(X)-3)))-
               (Gx*Gy/(nrow(X)*(nrow(X)-1))))/(nrow(X)*(nrow(X)-1)); 
    SEZ <- sqrt(VarZ);
    return(SEZ)
  }
  else
    return("Error: Matrices are not square and can not be analyzed")
}

# Function to compute Pearson's correlation for the two sociomatrices #

getpearson <- function(X,Y){
  
  # Original sociomatrices should have the same size #
  
  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 
  
  if (squaremat(X) & squaremat(Y)){
    X1 <- prepmat(X);
    Y1 <- prepmat(Y);
    vecX <- c(t(X1));
    vecY <- c(t(Y1));
    pearson <- cor(vecX,vecY);
  }
  else {
    vecX <- c(t(X));
    vecY <- c(t(Y));
    pearson <- cor(vecX,vecY);
  }
  return(pearson)
}

# Function to compute Dietz's R statistic #

getR <- function(X,Y){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (squaremat(X) & squaremat(Y)){
  X1 <- prepmat(X);
  Y1 <- prepmat(Y);
  Xrank <- rank(X1);
  Yrank <- rank(Y1);
  XYrank <- Xrank * Yrank;
  R <- sum(XYrank);
  return(R);}
  else {
      Xrank <- rank(X);
      Yrank <- rank(Y);
      XYrank <- Xrank * Yrank;
      R <- sum(XYrank);
      return(R);}
}

# Function to compute expected value for Dietz's R statistic #

getexpecR <- function(X,Y){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (squaremat(X) & squaremat(Y)){
  X1 <- prepmat(X);
  Y1 <- prepmat(Y);
  Xrank <- rank(X1);
  Yrank <- rank(Y1);
expecR <- sum(Xrank)*sum(Yrank)/(length(Xrank)*(length(Xrank)-1))
return(expecR);}
  else {
      Xrank <- rank(X);
      Yrank <- rank(Y);
      expecR <- sum(Xrank)*sum(Yrank)/(length(Xrank)*(length(Xrank)-1))
      return(expecR);}
}

# Function to compute standar error for Dietz's R statistic #

getSER <- function(X,Y){
  
  # Original sociomatrices should have the same size #
  
  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 
  
  if ((nrow(X) > 3) & (nrow(Y) > 3)){
    if (squaremat(X) & squaremat(Y)){
      X1 <- prepmat(X);
      Y1 <- prepmat(Y);
      Xrank <- rank(X1);
      Yrank <- rank(Y1);
      X[row(X)!=col(X)]<-Xrank;
      X<-t(X);
      Y[row(Y)!=col(Y)]<-Yrank;
      Y<-t(Y);
      Bx <- sum(X*X);
      By <- sum(Y*Y);
      Cx <- sum(X*t(X));
      Cy <- sum(Y*t(Y));
      Dx <- sum(rowSums(X)^2);
      Dy <- sum(rowSums(Y)^2);
      Ex <- sum(colSums(X)^2);
      Ey <- sum(colSums(Y)^2);
      Fx <- sum(rowSums(X)*colSums(X));
      Fy <- sum(rowSums(Y)*colSums(Y));
      Gx <- sum(X)**2;
      Gy <- sum(Y)**2;
      Hx <- Dx - Bx;
      Hy <- Dy - By;
      Ix <- Ex - Bx;
      Iy <- Ey - By;
      Jx <- Fx - Cx;
      Jy <- Fy - Cy;
      Kx <- Gx - Bx - Cx - Hx - Ix - 2*Jx;
      Ky <- Gy - By - Cy - Hy - Iy - 2*Jy;
      VarR <- ((Bx*By)+(Cx*Cy)+((Hx*Hy+Ix*Iy+2*Jx*Jy)/(nrow(X)-2))+(Kx*Ky/((nrow(X)-2)*(nrow(X)-3)))-
                 (Gx*Gy/(nrow(X)*(nrow(X)-1))))/(nrow(X)*(nrow(X)-1)); 
      SER <- sqrt(VarR);
      return(SER)}
    else stop("Matrices should be square");
  }
  else stop("Matrices should have size > 3x3");
}

# Function to compute Spearman's correlation for the two sociomatrices #

getspearman <- function(X,Y){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (squaremat(X) & squaremat(Y)){
  X1 <- prepmat(X);
  Y1 <- prepmat(Y);
  vecX <- array(dim=c(nrow(X1)*ncol(X1),1));
  m <- 0.;
  for (i in 1:nrow(X1))
     for (j in 1:ncol(X1)){
        m <- m+1;
        vecX[m] <- X1[i,j];
     }
  vecY <- array(dim=c(nrow(Y1)*ncol(Y1),1));
  m <- 0.;
  for (i in 1:nrow(Y1))
     for (j in 1:ncol(Y1)){
        m <- m+1;
        vecY[m] <- Y1[i,j];
     }
  spearman <- cor(vecX,vecY,method="spearman");
  return(spearman)
}
  else {
      vecX <- array(dim=c(nrow(X)*ncol(X),1));
      m <- 0.;
      for (i in 1:nrow(X))
         for (j in 1:ncol(X)){
            m <- m+1;
            vecX[m] <- X[i,j];
         }
      vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
      m <- 0.;
      for (i in 1:nrow(Y))
         for (j in 1:ncol(Y)){
            m <- m+1;
            vecY[m] <- Y[i,j];
         }
  spearman <- cor(vecX,vecY,method="spearman");
  return(spearman)      
  }
}

# Function to compute t statistic by means of Mantel's measures #

gett <- function(X,Y,method=c("mantel","dietz")){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

method <- match.arg(method)
if (squaremat(X) & squaremat(Y)){
  if (method == "mantel"){
  Z <- getZ(X,Y);
  expecZ <- getexpecZ(X,Y);
  SEZ <- getSEZ(X,Y);
  t <- (Z-expecZ)/SEZ;
  return(t)}
  if (method == "dietz"){
  R <- getR(X,Y);
  expecR <- getexpecR(X,Y);
  SER <- getSER(X,Y);
  t <- (R-expecR)/SER;
  return(t)}
  }
  else stop("Matrices should be square");
}

# Function to compute rowwise Zr statistic #

getZr <- function(X,Y,names=NULL){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (is.null(names)) {names <- paste('Ind.',LETTERS[c(1:nrow(X))])}
  else {names <- names};
wi <- array(dim=c(nrow(X),1),dimnames=list(c(names),'wi'));
Xterm <- array(dim=c(nrow(X),1));
Yterm <- array(dim=c(nrow(X),1));
if (squaremat(X) & squaremat(Y)){
  for (i in 1:nrow(X)){
     Xterm[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j > k){
             if ((i != j) & (i != k)) Xterm[i] <- Xterm[i] + (X[i,j]-X[i,k])**2
           }
        }
     Yterm[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j < k){
             if ((i != j) & (i != k)) Yterm[i] <- Yterm[i] + (Y[i,j]-Y[i,k])**2
           }
        }
     wi[i] <- sqrt(Xterm[i]*Yterm[i]);
  }
  ri <- array(dim=c(nrow(X),1),dimnames=list(c(names),'ri'));
  contrib <- array(dim=c(nrow(X),1),dimnames=list(c(names),'contrib'));
  Zr <- 0.;
  for (i in 1:nrow(X)){
     ri[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j > k){
             if ((i != j) & (i != k)) ri[i] <- ri[i] + ((X[i,j]-X[i,k])*(Y[i,j]-Y[i,k])/wi[i])
           }
        }
     Zr <- Zr + wi[i]*ri[i];
     contrib[i] <- wi[i]*ri[i];
  }
}
else {
    for (i in 1:nrow(X)){
       Xterm[i] <- 0.;
       for (j in 1:ncol(X))
          for (k in 1:ncol(X)){
             if (j > k){
               Xterm[i] <- Xterm[i] + (X[i,j]-X[i,k])**2}
          }
       Yterm[i] <- 0.;
         for (j in 1:ncol(X))
          for (k in 1:ncol(X)){
             if (j < k){
               Yterm[i] <- Yterm[i] + (Y[i,j]-Y[i,k])**2}
          }
       wi[i] <- sqrt(Xterm[i]*Yterm[i]);
    }
    ri <- array(dim=c(nrow(X),1),dimnames=list(c(names),'ri'));
    contrib <- array(dim=c(nrow(X),1),dimnames=list(c(names),'contrib'));
    Zr <- 0.;
    for (i in 1:nrow(X)){
       ri[i] <- 0.;
       for (j in 1:ncol(X))
          for (k in 1:ncol(X)){
             if (j > k){
               ri[i] <- ri[i] + ((X[i,j]-X[i,k])*(Y[i,j]-Y[i,k])/wi[i])}
          }
       Zr <- Zr + wi[i]*ri[i];
       contrib[i] <- wi[i]*ri[i];
    }
}
rrwav <- Zr/sum(wi);
rrw <- Zr/(sqrt(sum(Xterm)*sum(Yterm)));
list(Zr=Zr,rrwav=rrwav,rrw=rrw,rowpearson=ri,weighted_contributions=wi,contributions=contrib)
}
       
# Function to compute rowwise Rr statistic #

getRr <- function(X,Y,names=NULL){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (is.null(names)) {names <- paste('Ind.',LETTERS[c(1:nrow(X))])}
  else {names <- names};
if (squaremat(X) & squaremat(Y)){
  X1 <- prepmat(X);
  Y1 <- prepmat(Y);
  for (i in 1:nrow(X)){
     Xrank <- rank(X1[i,]);
     Yrank <- rank(Y1[i,])
     m <- 1.; 
     for (j in 1:ncol(X)){
        if (i == j) {X[i,j] <- 0.;
          Y[i,j] <- 0.}
          else {X[i,j] <- Xrank[m];
              Y[i,j] <- Yrank[m];
              m <- m+1}}
  }
  wi <- array(dim=c(nrow(X),1),dimnames=list(c(names),'wi'));
  Xterm <- array(dim=c(nrow(X),1));
  Yterm <- array(dim=c(nrow(X),1));
  for (i in 1:nrow(X)){
     Xterm[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j > k){
             if ((i != j) & (i != k)) Xterm[i] <- Xterm[i] + (X[i,j]-X[i,k])**2
           }
        }
     Yterm[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j < k){
             if ((i != j) & (i != k)) Yterm[i] <- Yterm[i] + (Y[i,j]-Y[i,k])**2
           }
        }
     wi[i] <- sqrt(Xterm[i]*Yterm[i]);
  }
  ri <- array(dim=c(nrow(X),1),dimnames=list(c(names),'ri'));
  contrib <- array(dim=c(nrow(X),1),dimnames=list(c(names),'contrib'));
  Rr <- 0.;
  for (i in 1:nrow(X)){
     ri[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j > k){
             if ((i != j) & (i != k)) ri[i] <- ri[i] + ((X[i,j]-X[i,k])*(Y[i,j]-Y[i,k])/wi[i])
           }
        }
     Rr <- Rr + wi[i]*ri[i];
     contrib[i] <- wi[i]*ri[i];
  }
}
  else {
      for (i in 1:nrow(X)){
         Xrank <- rank(X[i,]);
         Yrank <- rank(Y[i,])
         m <- 1.; 
         for (j in 1:ncol(X)){
            X[i,j] <- Xrank[m];
            Y[i,j] <- Yrank[m];
            m <- m+1}}
      wi <- array(dim=c(nrow(X),1),dimnames=list(c(names),'wi'));
      Xterm <- array(dim=c(nrow(X),1));
      Yterm <- array(dim=c(nrow(X),1));
      for (i in 1:nrow(X)){
         Xterm[i] <- 0.;
         for (j in 1:ncol(X))
            for (k in 1:ncol(X)){
               if (j > k){
                 Xterm[i] <- Xterm[i] + (X[i,j]-X[i,k])**2}
            }
         Yterm[i] <- 0.;
         for (j in 1:ncol(X))
            for (k in 1:ncol(X)){
               if (j < k){
                 Yterm[i] <- Yterm[i] + (Y[i,j]-Y[i,k])**2}
            }
         wi[i] <- sqrt(Xterm[i]*Yterm[i]);
      }
      ri <- array(dim=c(nrow(X),1),dimnames=list(c(names),'ri'));
      contrib <- array(dim=c(nrow(X),1),dimnames=list(c(names),'contrib'));
      Rr <- 0.;
      for (i in 1:nrow(X)){
         ri[i] <- 0.;
         for (j in 1:ncol(X))
            for (k in 1:ncol(X)){
               if (j > k){
                 ri[i] <- ri[i] + ((X[i,j]-X[i,k])*(Y[i,j]-Y[i,k])/wi[i])}
            }
         Rr <- Rr + wi[i]*ri[i];
         contrib[i] <- wi[i]*ri[i];
      }
  }
rhorwav <- Rr/sum(wi);
rhorw <- Rr/(sqrt(sum(Xterm)*sum(Yterm)));
list(Rr=Rr,rhorwav=rhorwav,rhorw=rhorw,rowrho=ri,weighted_contributions=wi,contributions=contrib)
}

# Function to compute rowwise Kr statistic #

getKr <- function(X,Y,names=NULL){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (is.null(names)) {names <- paste('Ind.',LETTERS[c(1:nrow(X))])}
  else {names <- names};
wi <- array(dim=c(nrow(X),1),dimnames=list(c(names),'wi'));
Xterm <- array(dim=c(nrow(X),1));
Yterm <- array(dim=c(nrow(X),1));
if (squaremat(X) & squaremat(Y)){
  for (i in 1:nrow(X)){
     Xterm[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j > k){
             if ((i != j) & (i != k)) Xterm[i] <- Xterm[i] + sign(X[i,j]-X[i,k])**2
           }
        }
     Yterm[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j < k){
             if ((i != j) & (i != k)) Yterm[i] <- Yterm[i] + sign(Y[i,j]-Y[i,k])**2
           }
        }
     wi[i] <- sqrt(Xterm[i]*Yterm[i]);
  }
  taui <- array(dim=c(nrow(X),1),dimnames=list(c(names),'taui'));
  contrib <- array(dim=c(nrow(X),1),dimnames=list(c(names),'contrib'));
  Kr <- 0.;
  for (i in 1:nrow(X)){
     taui[i] <- 0.;
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j > k){
             if ((i != j) & (i != k)) taui[i] <- taui[i] + (sign(X[i,j]-X[i,k])*sign(Y[i,j]-Y[i,k])/wi[i])
           }
        }
     Kr <- Kr + wi[i]*taui[i];
     contrib[i] <- wi[i]*taui[i];
  }
}
  else {
      for (i in 1:nrow(X)){
         Xterm[i] <- 0.;
         for (j in 1:ncol(X))
            for (k in 1:ncol(X)){
               if (j > k){
                 Xterm[i] <- Xterm[i] + sign(X[i,j]-X[i,k])**2}
            }
         Yterm[i] <- 0.;
         for (j in 1:ncol(X))
            for (k in 1:ncol(X)){
               if (j < k){
                 Yterm[i] <- Yterm[i] + sign(Y[i,j]-Y[i,k])**2}
            }
         wi[i] <- sqrt(Xterm[i]*Yterm[i]);
      }
      taui <- array(dim=c(nrow(X),1),dimnames=list(c(names),'taui'));
      contrib <- array(dim=c(nrow(X),1),dimnames=list(c(names),'contrib'));
      Kr <- 0.;
      for (i in 1:nrow(X)){
         taui[i] <- 0.;
         for (j in 1:ncol(X))
            for (k in 1:ncol(X)){
               if (j > k){
                 taui[i] <- taui[i] + (sign(X[i,j]-X[i,k])*sign(Y[i,j]-Y[i,k])/wi[i])}
            }
         Kr <- Kr + wi[i]*taui[i];
         contrib[i] <- wi[i]*taui[i];
      }
  }
taurwav <- Kr/sum(wi);
taurw <- Kr/(sqrt(sum(Xterm)*sum(Yterm)));
list(Kr=Kr,taurwav=taurwav,taurw=taurw,rowtau=taui,weighted_contributions=wi,contributions=contrib)
}

# Function to compute Kendall's statistics proposed by Dietz #
 
gettaukendall <- function(X,Y){

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (squaremat(X) & squaremat(Y)){
  Kr <- 0.;
  for (i in 1:nrow(X))
     for (j in 1:ncol(X))
        for (k in 1:ncol(X)){
           if (j > k){
     sig <- sign((X[i,j]-X[i,k])*(Y[i,j]-Y[i,k]));
     if ((i != j) & (i != k)) Kr <- Kr + sig
       else Kr <- Kr;}}
  Kc <- 0.;
  for (j in 1:ncol(X))
     for (i in 1:nrow(X))
        for (k in 1:nrow(X)){
           if (i > k ){
	     sig <- sign((X[i,j]-X[k,j])*(Y[i,j]-Y[k,j]));
	     if ((i !=j) & (j !=k)) Kc <- Kc + sig
	       else Kc <- Kc;}}
  Kc <- Kc + Kr;
  Ku <- 0.;
  for (i in 1:nrow(X))
     for (j in 1:ncol(X))
        for (k in 1:nrow(X))
           for (l in 1:ncol(X)){
	      sig <- sign((X[i,j]-X[k,l])*(Y[i,j]-Y[k,l]));
	      if ((i < j) & (i < k) & (k < l) & (j != k) & (j != l)) Ku <- Ku + sig
	        else Ku <- Ku;}
}
  else {
      Kr <- 0.;
      for (i in 1:nrow(X))
         for (j in 1:ncol(X))
            for (k in 1:ncol(X)){
               if (j > k){
                 sig <- sign((X[i,j]-X[i,k])*(Y[i,j]-Y[i,k]));
                 Kr <- Kr + sig}
            }
      Kc <- 0.;
      for (j in 1:ncol(X))
         for (i in 1:nrow(X))
            for (k in 1:nrow(X)){
               if (i > k ){
	         sig <- sign((X[i,j]-X[k,j])*(Y[i,j]-Y[k,j]));
	         Kc <- Kc + sig}
            }
      Kc <- Kc + Kr;
      Ku  <- 0.;
      for (i in 1:nrow(X))
         for (j in 1:ncol(X))
            for (k in 1:nrow(X))
               for (l in 1:ncol(X)){
	          sig <- sign((X[i,j]-X[k,l])*(Y[i,j]-Y[k,l]));
	          if ((i < j) & (i < k) & (k < l)) Ku <- Ku + sig}
  }
K <- Kc + Ku;
list(K=K,Kr=Kr,Kc=Kc,Ku=Ku);
}

# Function to compute partial rowwise Zr statistic #

getpartialZr <- function(X,Y,Z,names=NULL,perm=9999){
if (is.null(names)) {names <- paste('Ind.',LETTERS[c(1:nrow(X))])}
  else {names <- names}

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y)) | (nrow(X) != nrow(Z)) | (ncol(X) != ncol(Z))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

contrib <- array(dim=c(nrow(X),4),dimnames=list(c(names),c('rxy','rxz','ryz','Partial rxy.z')));
rrwxy <- getZr(X,Y)$rrw;
rrwxz <- getZr(X,Z)$rrw;
rrwyz <- getZr(Y,Z)$rrw;
rrwxy.z <- (rrwxy - rrwxz*rrwyz)/sqrt((1-rrwxz**2)*(1-rrwyz**2));
contrib[,1] <- getZr(X,Y)$rowpearson;
contrib[,2] <- getZr(X,Z)$rowpearson;
contrib[,3] <- getZr(Y,Z)$rowpearson;
term <- ((1-contrib[,2]**2)*(1-contrib[,3]**2))**0.5;
for (i in 1:nrow(X)){
if (term[i] == 0.) contrib[i,4] <- 0.
  else contrib[i,4] <- (contrib[i,1] - (contrib[i,2]*contrib[i,3]))/term[i];}

# Matrices have to be transformed into vectors in order to pass to a .C routine #
 
vecX <- array(dim=c(nrow(X)*ncol(X),1));
vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
vecZ <- array(dim=c(nrow(Z)*ncol(Z),1));
m <- 0.;
for (i in 1:nrow(X))
   for (j in 1:ncol(X)){
      m <- m+1;
      vecX[m] <- X[i,j];
      vecY[m] <- Y[i,j];
      vecZ[m] <- Z[i,j];
   }

if (squaremat(X) & squaremat(Y) & squaremat(Z)){

# R wrapper of the .C routine #
  if (nrow(X) < 9){
    out <- .C("partialexactsigZr",
           as.double(vecX),
           as.double(vecY),
           as.double(vecZ),
           as.integer(nrow(X)),
           pvpartialZrright=double(1),
           pvpartialZrleft=double(1),
           PACKAGE="DyaDA")}

    else {
        out <- .C("partialstatsigZr",
               as.double(vecX),
               as.double(vecY),
               as.double(vecZ),
               as.integer(nrow(X)),
               as.integer(perm),
               pvpartialZrright=double(1),
               pvpartialZrleft=double(1),
               PACKAGE="DyaDA")}
}
else {
    out <- .C("rectangularstatsigpartZr",
           as.double(vecX),
           as.double(vecY),
           as.double(vecZ),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           as.integer(perm),
           pvpartialZrright=double(1),
           pvpartialZrleft=double(1),
           PACKAGE="DyaDA")
}

pvright <- out$pvpartialZrright;
pvleft <- out$pvpartialZrleft;
list(rrwxy.z=rrwxy.z,p_value_right=pvright,p_value_left=pvleft,rrwxy=rrwxy,rrwxz=rrwxz,rrwyz=rrwyz,partial_rowwise_Z=contrib);
}

# Function to compute partial rowwise Rr statistic #

getpartialRr <- function(X,Y,Z,names=NULL,perm=9999){
if (is.null(names)) {names <- paste('Ind.',LETTERS[c(1:nrow(X))])}
  else {names <- names}

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y)) | (nrow(X) != nrow(Z)) | (ncol(X) != ncol(Z))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

contrib <- array(dim=c(nrow(X),4),dimnames=list(c(names),c('rhoxy','rhoxz','rhoyz','Partial rhoxy.z')));
rhorwxy <- getRr(X,Y)$rhorw;
rhorwxz <- getRr(X,Z)$rhorw;
rhorwyz <- getRr(Y,Z)$rhorw;
rhorwxy.z <- (rhorwxy - rhorwxz*rhorwyz)/sqrt((1-rhorwxz**2)*(1-rhorwyz**2));
contrib[,1] <- getRr(X,Y)$rowrho;
contrib[,2] <- getRr(X,Z)$rowrho;
contrib[,3] <- getRr(Y,Z)$rowrho;
term <- ((1-contrib[,2]**2)*(1-contrib[,3]**2))**0.5;
for (i in 1:nrow(X)){
if (term[i] == 0.) contrib[i,4] <- 0.
  else contrib[i,4] <- (contrib[i,1] - (contrib[i,2]*contrib[i,3]))/term[i];}

if (squaremat(X) & squaremat(Y) & squaremat(Z)){

# Matrices have to be transformed into vectors in order to pass to a .C routine #

  X1 <- prepmat(X);
  Y1 <- prepmat(Y);
  Z1 <- prepmat(Z);
  for (i in 1:nrow(X)){
     Xrank <- rank(X1[i,]);
     Yrank <- rank(Y1[i,]);
     Zrank <- rank(Z1[i,]);
     m <- 1.; 
     for (j in 1:ncol(X)){
        if (i == j) {X[i,j] <- 0.;
          Y[i,j] <- 0.;
          Z[i,j] <- 0.;}
          else {X[i,j] <- Xrank[m];
              Y[i,j] <- Yrank[m];
              Z[i,j] <- Zrank[m];
              m <- m+1}}
  }

  vecX <- array(dim=c(nrow(X)*ncol(X),1));
  vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
  vecZ <- array(dim=c(nrow(Z)*ncol(Z),1));
  m <- 0.;
  for (i in 1:nrow(X))
     for (j in 1:ncol(X)){
        m <- m+1;
        vecX[m] <- X[i,j];
        vecY[m] <- Y[i,j];
        vecZ[m] <- Z[i,j];
     }

# R wrapper of the .C routine #
  if (nrow(X) < 9){
    out <- .C("partialexactsigRr",
             as.double(vecX),
             as.double(vecY),
             as.double(vecZ),
             as.integer(nrow(X)),
             pvpartialRrright=double(1),
             pvpartialRrleft=double(1),
             PACKAGE=DYaDA)}

  else {
      out <- .C("partialstatsigRr",
             as.double(vecX),
             as.double(vecY),
             as.double(vecZ),
             as.integer(nrow(X)),
             as.integer(perm),
             pvpartialRrright=double(1),
             pvpartialRrleft=double(1),
             PACKAGE=DYaDA)}
}

else {
    for (i in 1:nrow(X)){
       Xrank <- rank(X[i,]);
       Yrank <- rank(Y[i,]);
       Zrank <- rank(Z[i,]);
       m <- 1.; 
       for (j in 1:ncol(X)){
          X[i,j] <- Xrank[m];
          Y[i,j] <- Yrank[m];
          Z[i,j] <- Zrank[m]; 
          m <- m+1}}

    vecX <- array(dim=c(nrow(X)*ncol(X),1));
    vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
    vecZ <- array(dim=c(nrow(Z)*ncol(Z),1));
    m <- 0.;
    for (i in 1:nrow(X))
       for (j in 1:ncol(X)){
          m <- m+1;
          vecX[m] <- X[i,j];
          vecY[m] <- Y[i,j];
          vecZ[m] <- Z[i,j];
       }

     out <- .C("rectangularstatsigpartRr",as.double(vecX),
            as.double(vecY),
            as.double(vecZ),
            as.integer(nrow(X)),
            as.integer(ncol(X)),
            as.integer(perm),
            pvpartialRrright=double(1),
            pvpartialRrleft=double(1),
            PACKAGE=DyaDA)
    }
pvright <- out$pvpartialRrright;
pvleft <- out$pvpartialRrleft;
list(rhorwxy.z=rhorwxy.z,p_value_right=pvright,p_value_left=pvleft,rhorwxy=rhorwxy,rhorwxz=rhorwxz,rhorwyz=rhorwyz,partial_rowwise_R=contrib);
}

# Function to compute partial rowwise Kr statistic #

getpartialKr <- function(X,Y,Z,names=NULL,perm=9999){
if (is.null(names)) {names <- paste('Ind.',LETTERS[c(1:nrow(X))])}
  else {names <- names}

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y)) | (nrow(X) != nrow(Z)) | (ncol(X) != ncol(Z))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

contrib <- array(dim=c(nrow(X),4),dimnames=list(c(names),c('tauxy','tauxz','tauyz','Partial tauxy.z')));
taurwxy <- getKr(X,Y)$taurw;
taurwxz <- getKr(X,Z)$taurw;
taurwyz <- getKr(Y,Z)$taurw;
taurwxy.z <- (taurwxy - taurwxz*taurwyz)/sqrt((1-taurwxz**2)*(1-taurwyz**2));
contrib[,1] <- getKr(X,Y)$rowtau;
contrib[,2] <- getKr(X,Z)$rowtau;
contrib[,3] <- getKr(Y,Z)$rowtau;
term <- ((1-contrib[,2]**2)*(1-contrib[,3]**2))**0.5;
for (i in 1:nrow(X)){
if (term[i] == 0.) contrib[i,4] <- 0.
  else contrib[i,4] <- (contrib[i,1] - (contrib[i,2]*contrib[i,3]))/term[i];}

# Matrices have to be transformed into vectors in order to pass to a .C routine #
 
vecX <- array(dim=c(nrow(X)*ncol(X),1));
vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
vecZ <- array(dim=c(nrow(Z)*ncol(Z),1));
m <- 0.;
for (i in 1:nrow(X))
   for (j in 1:ncol(X)){
      m <- m+1;
      vecX[m] <- X[i,j];
      vecY[m] <- Y[i,j];
      vecZ[m] <- Z[i,j];
   }

if (squaremat(X) & squaremat(Y) & squaremat(Z)){

# R wrapper of the .C routine #
  if (nrow(X) < 9){
    out <- .C("partialexactsiKr",
           as.double(vecX),
           as.double(vecY),
           as.double(vecZ),
           as.integer(nrow(X)),
           pvpartialKrright=double(1),
           pvpartialKrleft=double(1),
           PACKAGE="DyaDA")}

    else {
        out <- .C("partialstatsigKr",
               as.double(vecX),
               as.double(vecY),
               as.double(vecZ),
               as.integer(nrow(X)),
               as.integer(perm),
               pvpartialKrright=double(1),
               pvpartialKrleft=double(1),
               PACKAGE="DyaDA")}
}
else {
    out <- .C("rectangularstatsigpartKr",
           as.double(vecX),
           as.double(vecY),
           as.double(vecZ),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           as.integer(perm),
           pvpartialKrright=double(1),
           pvpartialKrleft=double(1),
           PACKAGE="DyaDA")
}

pvright <- out$pvpartialKrright;
pvleft <- out$pvpartialKrleft;
list(taurwxy.z=taurwxy.z,p_value_right=pvright,p_value_left=pvleft,taurwxy=taurwxy,taurwxz=taurwxz,taurwyz=taurwyz,partial_rowwise_Kr=contrib);
}

# Function to carry out the permutation test to estimate statistical significance for Mantel's Z statistic #

mantelZtest <- function (X,Y,perm=9999){
   
# Limit the number of replications #

  if ((perm < 1) | (perm > 1000000))
    return("Error: Number of permutations must be between 1 and 1000000");

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

# Matrices have to be transformed into vectors in order to pass to a .C routine #
 
vecX <- array(dim=c(nrow(X)*ncol(X),1));
m <- 0.;
for (i in 1:nrow(X))
for (j in 1:ncol(X))
{
m <- m+1;
vecX[m] <- X[i,j];
}
vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
m <- 0.;
for (i in 1:nrow(Y))
for (j in 1:ncol(Y))
{
m <- m+1;
vecY[m] <- Y[i,j];
}

if (squaremat(X) & squaremat(Y)){

# R wrapper of the .C routine #
  if (nrow(X) < 9){
    out <- .C("exactsigZ",as.double(vecX),
           as.double(vecY),
           as.integer(nrow(X)),
           pvZright=double(1),
           pvZleft=double(1),
           PACKAGE="DyaDA")}

    else {
        out <- .C("statsigZ",as.double(vecX),
               as.double(vecY),
               as.integer(nrow(X)),
               as.integer(perm),
               pvZright=double(1),
               pvZleft=double(1),
               PACKAGE="DyaDA")}
}
else {
    out <- .C("rectangularstatsigZ",as.double(vecX),
           as.double(vecY),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           as.integer(perm),
           pvZright=double(1),
           pvZleft=double(1),
           PACKAGE="DyaDA")
}
z <- getZ(X,Y);         
pvright <- out$pvZright;
pvleft <- out$pvZleft;
list(Mantel_Z=z,right_pvalue=pvright,left_pvalue=pvleft)      
}

# Function to carry out the permutation test to estimate statistical significance for Dietz's R statistic #

dietzRtest <- function (X,Y,perm=9999){
   
# Limit the number of replications #

  if ((perm < 1) | (perm > 1000000))
    return("Error: Number of permutations must be between 1 and 1000000");

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (squaremat(X) & squaremat(Y)){

# Matrices have to be transformed into vectors in order to pass to a .C routine #

  X1 <- prepmat(X);
  Y1 <- prepmat(Y);
  Xrank <- rank(X1);
  Yrank <- rank(Y1);

m <- 1.;
    for (i in 1:nrow(X))
       for (j in 1:ncol(X)){
          if (i == j) {X[i,j] <- 0.}
            else {X[i,j] <- Xrank[m];
                m <- m+1;}
       } 

m <- 1.;
    for (i in 1:nrow(Y))
       for (j in 1:ncol(Y)){
          if (i == j) {Y[i,j] <- 0.}
            else {Y[i,j] <- Yrank[m];
                m <- m+1;} 
       }

vecX <- array(dim=c(nrow(X)*ncol(X),1));
m <- 0.;
for (i in 1:nrow(X))
   for (j in 1:ncol(X)){
      m <- m+1;
      vecX[m] <- X[i,j];
   }
vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
m <- 0.;
for (i in 1:nrow(Y))
   for (j in 1:ncol(Y)){
      m <- m+1;
      vecY[m] <- Y[i,j];
   }

# R wrapper of the .C routine #
if (nrow(X) < 9){
  out <- .C("exactsigR",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             pvRright=double(1),
             pvRleft=double(1),
             PACKAGE="DyaDA")}

  else {
      out <- .C("statsigR",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             as.integer(perm),
             pvRright=double(1),
             pvRleft=double(1),
             PACKAGE="DyaDA")}
}

else {
    Xrank <- rank(X);
    Yrank <- rank(Y);

    m <- 1.;
    for (i in 1:nrow(X))
       for (j in 1:ncol(X)){
           X[i,j] <- Xrank[m];
           m <- m+1;}

    m <- 1.;
    for (i in 1:nrow(Y))
       for (j in 1:ncol(Y)){
          Y[i,j] <- Yrank[m];
                m <- m+1;}

    vecX <- array(dim=c(nrow(X)*ncol(X),1));
    m <- 0.;
    for (i in 1:nrow(X))
       for (j in 1:ncol(X)){
          m <- m+1;
          vecX[m] <- X[i,j];
       }
       vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
       m <- 0.;
       for (i in 1:nrow(Y))
          for (j in 1:ncol(Y)){
             m <- m+1;
             vecY[m] <- Y[i,j];
          }
     out <- .C("rectangularstatsigR",as.double(vecX),
            as.double(vecY),
            as.integer(nrow(X)),
            as.integer(ncol(X)),
            as.integer(perm),
            pvRright=double(1),
            pvRleft=double(1),
            PACKAGE="DyaDA")
    }
R <- getR(X,Y);         
pvright <- out$pvRright;
pvleft <- out$pvRleft;
list(Dietz_R=R,right_pvalue=pvright,left_pvalue=pvleft)
}

# Function to carry out the permutation test to estimate statistical significance for Mantel's Zr statistic #

rowwiseZrtest <- function (X,Y,perm=9999){
   
# Limit the number of replications #

  if ((perm < 1) | (perm > 1000000))
    return("Error: Number of permutations must be between 1 and 1000000");

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal");  

# Matrices have to be transformed into vectors in order to pass to a .C routine #

vecX <- array(dim=c(nrow(X)*ncol(X),1));
m <- 0.;
for (i in 1:nrow(X))
   for (j in 1:ncol(X)){
      m <- m+1;
      vecX[m] <- X[i,j];
   }
vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
m <- 0.;
for (i in 1:nrow(Y))
   for (j in 1:ncol(Y)){
      m <- m+1;
      vecY[m] <- Y[i,j];
   }

if (squaremat(X) & squaremat(Y)){

# R wrapper of the .C routine #

if (nrow(X) < 9){
  out <- .C("exactsigZr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             pvZrright=double(1),
             pvZrleft=double(1),
             PACKAGE="DyaDA")}

  else {
      out <- .C("statsigZr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             as.integer(perm),
             pvZrright=double(1),
             pvZrleft=double(1),
             PACKAGE="DyaDA")}
}

else {
    out <- .C("rectangularstatsigZr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             as.integer(ncol(X)),
             as.integer(perm),
             pvZrright=double(1),
             pvZrleft=double(1),
             PACKAGE="DyaDA")
}
  
Zr <- getZr(X,Y)$Zr;         
pvright <- out$pvZrright;
pvleft <- out$pvZrleft;
list(Zr=Zr,right_pvalue=pvright,left_pvalue=pvleft)
}

# Function to carry out the permutation test to estimate statistical significance for Dietz's Rr statistic #

rowwiseRrtest <- function (X,Y,perm=9999){

# Limit the number of replications #

  if ((perm < 1) | (perm > 1000000))
    return("Error: Number of permutations must be between 1 and 1000000");

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

if (squaremat(X) & squaremat(Y)){

# Matrices have to be transformed into vectors in order to pass to a .C routine #

X1 <- prepmat(X);
Y1 <- prepmat(Y);
for (i in 1:nrow(X)){
   Xrank <- rank(X1[i,]);
   Yrank <- rank(Y1[i,]);
   m <- 1.;
   for (j in 1:ncol(X)){
      if (i == j) {X[i,j] <- 0.;
        Y[i,j] <- 0.}
        else {X[i,j] <- Xrank[m];
            Y[i,j] <- Yrank[m];
            m <- m+1}}
}
vecX <- array(dim=c(nrow(X)*ncol(X),1));
m <- 0.;
for (i in 1:nrow(X))
for (j in 1:ncol(X))
{
m <- m+1;
vecX[m] <- X[i,j];
}
vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
m <- 0.;
for (i in 1:nrow(Y))
for (j in 1:ncol(Y))
{
m <- m+1;
vecY[m] <- Y[i,j];
}

# R wrapper of the .C routine #
if (nrow(X) < 9){
  out <- .C("exactsigRr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             pvRrright=double(1),
             pvRrleft=double(1),
             PACKAGE="DyaDA")}

  else {
      out <- .C("statsigRr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             as.integer(perm),
             pvRrright=double(1),
             pvRrleft=double(1),
             PACKAGE="DyaDA")}
}

else {
    for (i in 1:nrow(X)){
       Xrank <- rank(X[i,]);
       Yrank <- rank(Y[i,])
       m <- 1.; 
       for (j in 1:ncol(X)){
          X[i,j] <- Xrank[m];
          Y[i,j] <- Yrank[m];
          m <- m+1}}

    vecX <- array(dim=c(nrow(X)*ncol(X),1));
    m <- 0.;
    for (i in 1:nrow(X))
       for (j in 1:ncol(X)){
          m <- m+1;
          vecX[m] <- X[i,j];
       }

    vecY <- array(dim=c(nrow(Y)*ncol(Y),1));
    m <- 0.;
    for (i in 1:nrow(Y))
       for (j in 1:ncol(Y)){
          m <- m+1;
          vecY[m] <- Y[i,j];
       }
    out <- .C("rectangularstatsigRr",
             as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             as.integer(ncol(X)),
             as.integer(perm),
             pvRrright=double(1),
             pvRrleft=double(1),
             PACKAGE="DyaDA")
}

# R wrapper of the .C routine #
Rr <- getRr(X,Y)$Rr;         
pvright <- out$pvRrright;
pvleft <- out$pvRrleft;
list(Rr=Rr,right_pvalue=pvright,left_pvalue=pvleft)
}

# Function to carry out the permutation test to estimate statistical significance for Kendall's Kr statistic #

rowwiseKrtest <- function (X,Y,perm=9999){

# Limit the number of replications #

  if ((perm < 1) | (perm > 1000000))
    return("Error: Number of permutations must be between 1 and 1000000"); 

# Original sociomatrices should have the same size #

  if ((nrow(X) != nrow(Y)) | (ncol(X) != ncol(Y))) 
    return("Error: Sizes of original sociomatrices are expected to be equal"); 

# Matrices have to be transformed into vectors in order to pass to a .C routine #

vecX <- array(dim=c(nrow(X)*(ncol(X)-1)/2,1));
m <- 0.;
for (i in 1:nrow(X))
for (j in 1:ncol(X))
{
m <- m+1;
vecX[m] <- X[i,j];
}
vecY <- array(dim=c(nrow(Y)*(ncol(Y)-1)/2,1));
m <- 0.;
for (i in 1:nrow(Y))
for (j in 1:ncol(Y))
{
m <- m+1;
vecY[m] <- Y[i,j];
}

if (squaremat(X) & squaremat(Y)){

# R wrapper of the .C routine #
if (nrow(X) < 9){
  out <- .C("exactsigKr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             pvKrright=double(1),
             pvKrleft=double(1),
             PACKAGE="DyaDA")}

  else {
      out <- .C("statsigKr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             as.integer(perm),
             pvKrright=double(1),
             pvKrleft=double(1),
             PACKAGE="DyaDA")}
}
  else {  
      out <- .C("rectangularstatsigKr",as.double(vecX),
             as.double(vecY),
             as.integer(nrow(X)),
             as.integer(ncol(X)),
             as.integer(perm),
             pvKrright=double(1),
             pvKrleft=double(1),
             PACKAGE="DyaDA")
  }

Kr <- getKr(X,Y)$Kr;         
pvright <- out$pvKrright;
pvleft <- out$pvKrleft;
list(kendall_Kr=Kr,right_pvalue=pvright,left_pvalue=pvleft)
}

################################################################################
#                 R FUNCTIONS FOR  DYADIC INTERDEPENDENCE                      #
################################################################################

# A function for estimating dyadic interdependence in standard dyadic designs with distinguishable dyad members and interval responses #

distinguishable.dyad <- function(dataset,conf.int=0.95){
  N<-length(na.omit(dataset)[,2])
  na1<- sum(is.na(dataset[,2]))
  na2<- sum(is.na(dataset[,3]))
  dataset<-na.omit(dataset)
  summary.statistics <- array(c(N,na1,min(dataset[,2]),max(dataset[,2]),mean(dataset[,2]),
                                sd(dataset[,2]),quantile(dataset[,2],.25,names=F),quantile(dataset[,2],.50,names=F),
                                quantile(dataset[,2],.75,names=F),N,na2,min(dataset[,3]),max(dataset[,3]),
                                mean(dataset[,3]),sd(dataset[,3]),quantile(dataset[,3],.25,names=F),quantile(dataset[,3],.50,names=F),
                                quantile(dataset[,3],.75,names=F)),dim=c(9,2),dimnames=list(c("N:", "NAs:","Min:","Max:","Mean:","Sd:",
                                                                                              "25th Pctl:","50th Pctl:", "75th Pctl:"),c(colnames(dataset[,2:3]))))
  
  r <- cor(dataset[,2],dataset[,3])
  df <- N-2
  alpha <- 1 - conf.int
  t.statistic <- r*sqrt(N-2)/(sqrt(1-r**2))
  p.value <- 1 - pt(t.statistic,df)
  zr.value <- 0.5*log((1+r)/(1-r))
  zr.low <- zr.value - qnorm((1-alpha/2),0,1)/sqrt(N-3)
  zr.upper <- zr.value + qnorm((1-alpha/2),0,1)/sqrt(N-3)
  r.low <- (exp(1)**(2*zr.low)-1)/(exp(1)**(2*zr.low)+1)
  r.upper <- (exp(1)**(2*zr.upper)-1)/(exp(1)**(2*zr.upper)+1)
  res <- list(call=match.call(),data=dataset,stats=summary.statistics,r=r,rCI=c(r.low,r.upper),
              t=t.statistic,df=df,pval=2*p.value)
  class(res) <- 'dyadr'
  res
}

print.dyadr <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Distinguishable Members: Pearson's Correlation")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Data",  if (length(x$data[,1]) <= 5) ": " else " (5 first rows shown): ", "\n")
  print( if (length(x$data[,1]) >= 5) x$data[1:5,] else x$data[1:length(x$data),])
  cat("\n\n")
  cat("Descriptive Statistics","\n")
  print.table(x$stats, digits=digits)
  cat("\n\n")
  cat("t Test: ","\n")
  results <- cbind(x$t,as.integer(x$df),x$pval)
  colnames(results) <- c("t Value", "Df", "Pr(>|t|)")
  rownames(results) <- ""
  print.table(results,digits=digits)
  cat("\n\n")
  cat(paste((1-x$alpha)*100,"% Confidence Interval: ",sep=''),"\n")
  results <- t(x$rCI)
  colnames(results) <- c("Lower","Upper")
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  cat("Sample estimate","\n")
  results <- cbind(x$r)
  colnames(results) <- "cor"
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  invisible(x)
}

# A function for estimating dyadic interdependence in standard dyadic designs with indistinguishable dyad members and interval responses #

indistinguishable.dyad <- function(dataset,conf.int=0.95){
  N<-length(na.omit(dataset)[,2])
  na1<- sum(is.na(dataset[,2]))
  na2<- sum(is.na(dataset[,3]))
  dataset<-na.omit(dataset)
  summary.statistics <- array(c(N,na1,min(dataset[,2]),max(dataset[,2]),mean(dataset[,2]),
                                sd(dataset[,2]),quantile(dataset[,2],.25,names=F),quantile(dataset[,2],.50,names=F),
                                quantile(dataset[,2],.75,names=F),N,na2,min(dataset[,3]),max(dataset[,3]),
                                mean(dataset[,3]),sd(dataset[,3]),quantile(dataset[,3],.25,names=F),quantile(dataset[,3],.50,names=F),
                                quantile(dataset[,3],.75,names=F)),dim=c(9,2),dimnames=list(c("N:", "NAs:","Min:","Max:","Mean:","Sd:",
                                                                                              "25th Pctl:","50th Pctl:", "75th Pctl:"),c(colnames(dataset[,2:3]))))
  MSb <- 2*var(apply(dataset[,2:3],1,mean))
  MSw <- sum((dataset[,2]-dataset[,3])**2)/(2*N)
  ICC <- (MSb-MSw)/(MSb+MSw)
  alpha <- 1 - conf.int
  if ( MSb > MSw ) {
    F.statistic <- MSb/MSw
    df1 <- N-1
    df2 <- N}
  if (MSw > MSb ){
    F.statistic <- MSw/MSb
    df1 <- N
    df2 <- N-1}
  p.value <- pf(F.statistic,df1,df2,lower.tail=FALSE)
  fl <- MSb/MSw/qf(1-alpha/2,df1,df2)
  fu <- MSb/MSw*qf(1-alpha/2,df1,df2)
  icc.low <- (fl-1)/(fl+1)
  icc.upper <- (fu-1)/(fu+1)
  res=list(call=match.call(),data=dataset,stats=summary.statistics,
           intracor=ICC,MSb=MSb,MSw=MSw,Fstat=F.statistic,df1=df1,df2=df2,pval=2*p.value,
           alpha=alpha,iccCI=c(icc.low,icc.upper))
  class(res)="dyadicc"
  res
}

print.dyadicc <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Indistinguishable Members: Intraclass Correlation")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Data",  if (length(x$data[,1]) <= 5) ": " else " (5 first rows shown): ", "\n")
  print( if (length(x$data[,1]) >= 5) x$data[1:5,] else x$data[1:length(x$data),])
  cat("\n\n")
  cat("Descriptive Statistics","\n")
  print.table(x$stats, digits=digits)
  cat("\n\n")
  cat("F Test: ","\n")
  results <- cbind(rbind(x$MSb,x$MSw),rbind(x$df1,x$df2),rbind(x$Fstat,NA),rbind(x$pval,NA))
  colnames(results) <- c("Mean Sq", "Df", "F value", "Pr(>|F|)")
  rownames(results) <- c("Between-Dyads","Within-Dyads")
  print.table(results,digits=digits)
  cat("\n\n")
  cat(paste("Intraclass Correlation and ",(1-x$alpha)*100,"% Confidence Interval: ",sep=''),"\n")
  results <- cbind(x$intracor,t(x$iccCI))
  colnames(results) <- c("ICC", "Lower","Upper")
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  invisible(x)
}

# A function for estimating dyadic interdependence in standard dyadic designs with indistinguishable dyad members and interval responses by means of the pairwise correlation #

pairwise.correlation <-
  function (dataset, conf.int = 0.95){
    alpha <- 1 - conf.int
    double.entry.dataset <- array(c(c(1:(2 * length(dataset[,1]))), c(dataset[, 2], dataset[, 3]),
                                  c(dataset[, 3],dataset[, 2])), dim = c((2 * length(dataset[, 1])), length(dataset)))
    colnames(double.entry.dataset) <- colnames(dataset)
    dataset <- double.entry.dataset
    N <- length(na.omit(dataset)[, 2])
    na1 <- sum(is.na(dataset[, 2]))
    na2 <- sum(is.na(dataset[, 3]))
    dataset <- na.omit(dataset)
    summary.statistics <- array(c(N, na1, min(dataset[, 2]),
                                  max(dataset[, 2]), mean(dataset[, 2]), sd(dataset[, 2]),
                                  quantile(dataset[, 2], 0.25, names = F), quantile(dataset[,2],0.5, names = F),
                                  quantile(dataset[, 2], 0.75,names = F), N, na2, min(dataset[, 3]), max(dataset[,
                                  3]), mean(dataset[, 3]), sd(dataset[, 3]), quantile(dataset[,
                                          3], 0.25, names = F), quantile(dataset[, 3], 0.5,
                                  names = F), quantile(dataset[, 3], 0.75, names = F)),                                                                               
                                  dim = c(9, 2), dimnames = list(c("N:", "NAs:", "Min:",
                                                                 "Max:", "Mean:", "Sd:", "25th Pctl:", "50th Pctl:",
                                                                 "75th Pctl:"), c(colnames(dataset[, 2:3]))))
    rp <- cor(dataset[, 2], dataset[, 3])
    standard.error <- 1/sqrt(N/2)
    z <- rp/standard.error
    p.value <- 1 - pnorm(abs(z), 0, 1)
    rp.low <- rp - (qnorm((1 - alpha/2), 0, 1)/sqrt(N/2))
    rp.upper <- rp + (qnorm((1 - alpha/2), 0, 1)/sqrt(N/2))
    res <- list(call = match.call(), data = dataset, stats = summary.statistics,
                rp = rp, alpha = alpha, rpCI = c(rp.low, rp.upper), zVal = z,
                stdErr = standard.error, pval = 2 * p.value)
    class(res) <- "rp"
    res
  }

print.rp <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Indistinguishable Members: Pairwise Correlation")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Data",  if (length(x$data[,1]) <= 5) ": " else " (5 first rows shown): ", "\n")
  print( if (length(x$data[,1]) >= 5) x$data[1:5,] else x$data[1:length(x$data),])
  cat("\n\n")
  cat("Descriptive Statistics","\n")
  print.table(x$stats, digits=digits)
  cat("\n\n")
  cat("z Test: ","\n")
  results <- cbind(x$zVal,x$stdErr,x$pval)
  colnames(results) <- c("Z Value", "Std Error", "Pr(>|Z|)")
  rownames(results) <- ""
  print.table(results,digits=digits)
  cat("\n\n")
  cat(paste((1-x$alpha)*100,"% Confidence Interval: ",sep=''),"\n")
  results <- t(x$rpCI)
  colnames(results) <- c("Lower","Upper")
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  cat("Sample estimate","\n")
  results <- cbind(x$rp)
  colnames(results) <- "rp"
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  invisible(x)
}

# A function for estimating dyadic interdependence in standard dyadic designs with distinguishable dyad members and categorical responses #

categorical.distinguishable.dyad <- function(dataset,conf.int=0.95){
  tabulated.data <- table(dataset[,2],dataset[,3])
  alpha <- 1 - conf.int
  no <- sum(diag(tabulated.data))
  ne <- sum(margin.table(tabulated.data,1)*margin.table(tabulated.data,2)/sum(tabulated.data))
  kappa <- (no - ne)/(sum(tabulated.data)-ne)
  table.props <- prop.table(tabulated.data)
  prop.margins <-array(c(margin.table(table.props,1),margin.table(prop.table(tabulated.data),2)),dim=c(nrow(tabulated.data),2))
  sumterm1 <- 0.
  sumterm2 <- 0.
  sumterm3 <- 0.
  pe <- ne/sum(tabulated.data)
  for (i in 1:nrow(table.props)){
    sumterm1 <- sumterm1+(table.props[i,i]*(1-sum(prop.margins[i,])*(1-kappa))**2+(1-kappa)**2)
    sumterm3 <- sumterm3+(prop.margins[i,1]*prop.margins[i,2]*sum(prop.margins[i,]))} 
  for (i in 1:nrow(table.props)){
    for (j in 1:ncol(table.props)){
      sumterm2 <- sumterm2+(table.props[i,j]*(prop.margins[i,2]+prop.margins[j,1])**2-(kappa-pe*(1-kappa))**2)}} 
  se.kappa2 <- sqrt((sumterm1*sumterm2)/(sum(tabulated.data)*(1-pe)**2))
  se.kappa <- sqrt((pe+pe**2-sumterm3)/(sum(tabulated.data)*(1-pe)**2))
  tabulated.data <- addmargins(tabulated.data)
  z.value <- kappa/se.kappa
  p.value <- 1-pnorm(abs(z.value),0,1)
  kappa.low <- kappa-qnorm((1-alpha/2),0,1)*se.kappa
  kappa.upper <- kappa+qnorm((1-alpha/2),0,1)*se.kappa
  res <-list(call=match.call(),contingencyTab=tabulated.data,kappa=kappa,alpha=alpha,kappaCI=c(kappa.low,kappa.upper),
             zval=z.value,stdError=se.kappa,pval=2*p.value)
  class(res) <- "dyadkappa"
  res
}


print.dyadkappa <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Distinguishable Members with categorical data: Kappa Index")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Contingency table: \n")
  print.table(x$contingencyTab)
  cat("\n\n")
  cat("Agreement Test: ","\n")
  results <- cbind(x$kappa,x$stdError,x$zval,x$pval)
  colnames(results) <- c("Kappa Value", "Asymp. Std Error", "z Value", "Approx. Sig.")
  rownames(results) <- ""
  print.table(results,digits=digits,scientific=FALSE)
  cat("\n\n")
  cat(paste((1-x$alpha)*100,"% Confidence Interval: ",sep=''),"\n")
  results <- t(x$kappaCI)
  colnames(results) <- c("Lower","Upper")
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  invisible(x)
}

# A function for estimating dyadic interdependence in standard dyadic designs with indistinguishable dyad members and categorical responses #

categorical.indistinguishable.dyad <- function(dataset,conf.int=0.95){
  tabulated.data <- table(dataset[,2],dataset[,3])
  tabulated.data <- (tabulated.data+t(tabulated.data))/2
  alpha <- 1 - conf.int
  no <- sum(diag(tabulated.data))
  ne <- sum(margin.table(tabulated.data,1)*margin.table(tabulated.data,2)/sum(tabulated.data))
  kappa <- (no - ne)/(sum(tabulated.data)-ne)
  table.props <- prop.table(tabulated.data)
  prop.margins <-array(c(margin.table(table.props,1),margin.table(prop.table(tabulated.data),2)),dim=c(nrow(tabulated.data),2))
  sumterm1 <- 0.
  sumterm2 <- 0.
  sumterm3 <- 0.
  pe <- ne/sum(tabulated.data)
  for (i in 1:nrow(table.props)){
    sumterm1 <- sumterm1+(table.props[i,i]*(1-sum(prop.margins[i,])*(1-kappa))**2+(1-kappa)**2)
    sumterm3 <- sumterm3+(prop.margins[i,1]*prop.margins[i,2]*sum(prop.margins[i,]))} 
  for (i in 1:nrow(table.props)){
    for (j in 1:ncol(table.props)){
      sumterm2 <- sumterm2+(table.props[i,j]*(prop.margins[i,2]+prop.margins[j,1])**2-(kappa-pe*(1-kappa))**2)}} 
  se.kappa2 <- sqrt((sumterm1*sumterm2)/(sum(tabulated.data)*(1-pe)**2))
  se.kappa <- sqrt((pe+pe**2-sumterm3)/(sum(tabulated.data)*(1-pe)**2))
  tabulated.data <- addmargins(tabulated.data)
  z.value <- kappa/se.kappa
  p.value <- 1-pnorm(abs(z.value),0,1)
  kappa.low <- kappa-qnorm((1-alpha/2),0,1)*se.kappa
  kappa.upper <- kappa+qnorm((1-alpha/2),0,1)*se.kappa
  res <-list(call=match.call(),contingencyTab=tabulated.data,kappa=kappa,alpha=alpha,kappaCI=c(kappa.low,kappa.upper),
             zval=z.value,stdError=se.kappa,pval=2*p.value)
  class(res) <- "dyadkappa2"
  res
}

print.dyadkappa2 <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Indistinguishable Members with categorical data: Kappa Index")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Contingency table: \n")
  print.table(x$contingencyTab)
  cat("\n\n")
  cat("Agreement Test: ","\n")
  results <- cbind(x$kappa,x$stdError,x$zval,x$pval)
  colnames(results) <- c("Kappa Value", "Asymp. Std Error", "z Value", "Approx. Sig.")
  rownames(results) <- ""
  print.table(results,digits=digits)
  cat("\n\n")
  cat(paste((1-x$alpha)*100,"% Confidence Interval: ",sep=''),"\n")
  results <- t(x$kappaCI)
  colnames(results) <- c("Lower","Upper")
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  invisible(x)
}

################################################################################
#                 R FUNCTIONS FOR  SRM ROUND ROBIN DESIGNS                     #
################################################################################

# Function for formatting data according a round robin design

prepare.SRM.matrix <- function(X,numb.individuals,numb.times,numb.groups,names=NULL,times=NULL,groups=NULL){
  if (is.null(names)) names <- paste('Ind.',c(1:numb.individuals))
  if (is.null(times)) times <- paste('Time',c(1:numb.times))
  if (is.null(groups)) groups <- paste('Group',c(1:numb.groups))
  X <- array(X,c(numb.individuals,numb.individuals,numb.times,numb.groups))
  
  temp <- array(0,c(numb.individuals,numb.individuals,numb.times,numb.groups))
  for (i in 1:(dim(X)[1])){
    for (j in 1:(dim(X)[1])){
      for (k in 1:(dim(X)[3])){
        for (l in 1:(dim(X)[4])){
          temp[i,j,k,l] <- X[j,i,k,l]}}}}
  count <- 1;
  count2 <- 1;
  for (k in 1:(dim(X)[3])){
    for (i in 1:(dim(X)[1])){
      if (count <= (dim(X)[3])){ 
        X[count2,,count,] <- temp[i,,k,]
        count <- count + 1}
      else {
        count <- 1
        count2 <- count2 + 1
        X[count2,,count,] <- temp[i,,k,]
        count <- count + 1}}}
  res<-list(data=X,individuals=numb.individuals,ntimes=numb.times,ngroups=numb.groups,
            names=names,times=times,groups=groups)
  class(res) <- 'srmRRMatrix'
  res
}


# Function for obtaining SRM effects #

SRM.effects <- function (X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  # Compute row means #
  row.means <- array(dim=c(N,G),0.)
  for (l in 1:G){
    row.means[,l] <- rowSums(X$data[,,,l])/(G*(N-1));}
  time.rowmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      time.rowmeans[,k,l] <- rowSums(X$data[,,k,l])/(N-1)}}
  # Compute column means #
  column.means <- array(dim=c(N,G),0.)
  if (T <= 1) {
    for (l in 1:G){
      column.means[,l] <- colSums(X$data[,,,l])/(T*(N-1))}}
  if (T > 1) {
    for (l in 1:G){
      column.means[,l] <- rowSums(colSums(X$data[,,,l]))/(T*(N-1))}}
  time.colmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      time.colmeans[,k,l] <- colSums(X$data[,,k,l])/(N-1)}}
  # Compute grand mean #
  grand.mean <- array(dim=c(G,1),0.)
  for (l in 1:G){
    grand.mean[l] <- sum(X$data[,,,l])/(T*(N*(N-1)))}
  time.grandmean <- array(dim=c(T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      time.grandmean[k,l] <- sum(X$data[,,k,l])/((N*(N-1)))}}
  # Estimate actor effects for all individuals #
  actor <- array(dim=c(N,T,G),dimnames=list(X$names,X$times,X$groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      for (i in 1:N){
        actor[i,k,l] <- (N-1)*((N-1)*time.rowmeans[i,k,l]+time.colmeans[i,k,l]-
                                 N*time.grandmean[k,l])/(N*(N-2))}}}
  # Estimate partner effects for all individuals #
  partner <- array(dim=c(N,T,G),dimnames=list(X$names,X$times,X$groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      for (i in 1:N){
        partner[i,k,l] <- (N-1)*((N-1)*time.colmeans[i,k,l]+time.rowmeans[i,k,l]-
                                   N*time.grandmean[k,l])/(N*(N-2))}}}
  # Estimate relationships effects for all individuals #
  relationship <- array(dim=c(N,N,T,G),
                        dimnames=list(X$names,X$names,X$times,X$groups),0.)
  for (i in 1:N){
    for (j in 1:N){
      for (k in 1:T){
        for (l in 1:G){
          if (i != j)
          {relationship[i,j,k,l] <- X$data[i,j,k,l]-actor[i,k,l]-partner[j,k,l]-time.grandmean[k,l]}}}}}
  res<-list(actor.effects=actor,partner.effects=partner,relationship.effects=relationship)
  class(res)="srmRREffects"
  res
}

print.srmRREffects <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("SRM Effects: Round Robin Design")
  cat("\n\n")
  N <- dim(x$actor.effects)[1]
  T <- dim(x$actor.effects)[2]
  G <- dim(x$actor.effects)[3]
  names <- unlist(dimnames(x$actor.effects)[1])
  times <- unlist(dimnames(x$actor.effects)[2])
  groups <- unlist(dimnames(x$actor.effects)[3])
  group <- rep(groups,2*N*T)
  time <- rep(times,2*N*G)
  individual <- rep(names,2)
  type <- c(rep("actor",N*T*G),rep("partner",N*T*G))
  effect <- c (x$actor.effects,x$partner.effects)
  results<-data.frame(group,time,individual,type,effect)
  cat("Actor and Partner Effects: \n")
  print(results)
  cat("\n\n")
  results <- vector()
  
  for (j in 1:G)
  {  
    for (i in 1:T)
    {  
      mat <- x$relationship.effects[,,i,j]
      actornames <- c(rownames(mat)[(row(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)]],
                      rownames(mat)[(row(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]])
      partnernames <- c(colnames(mat)[(col(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)]],
                        colnames(mat)[(col(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]])
      effects <- mat[cbind(c((row(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)],
                             (row(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]),
                           c((col(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)],
                             (col(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]))]
      group <- rep(groups[j],N*(N-1)/2)
      time <- rep(times[i],N*(N-1)/2)
      res <- data.frame(group,time,actor=actornames,partner=partnernames,relationship=effects)
      results<- rbind(results,res)
    }
  }
  results<-results[order(results$actor,results$partner),]
  rownames(results) <- NULL
  cat("Relationship Effects: \n")
  print(results)
  cat("\n\n")
  invisible(x)
}

# Function to estimate SRM variances #

SRM.variances <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups 
  effects <- SRM.effects(X)
  actor <- effects$actor.effects
  partner <- effects$partner.effects
  relationship <- effects$relationship.effects
  # Estimate mean squares for actor, partner and cross products actor-partner #
  MSactor <- array(dim=c(T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      MSactor[k,l] <- sum(actor[,k,l]**2)/(N-1)}}
  MSpartner <- array(dim=c(T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      MSpartner[k,l] <- sum(partner[,k,l]**2)/(N-1)}}
  MCP <- array(dim=c(T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      MCP[k,l] <- sum(actor[,k,l]*partner[,k,l])/(N-1)}}
  # Estimate mean squares between and within dyads #
  sum.relationship <- array(dim=c(N,N,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      sum.relationship[,,k,l] <- relationship[,,k,l] + t(relationship[,,k,l])}}
  MSbetween <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      MSbetween[k,l] <- (sum(sum.relationship[,,k,l]**2)/2)/((N-1)*(N-2)-2)}}
  subs.relationship <- array(dim=c(N,N,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      subs.relationship[,,k,l] <- relationship[,,k,l] - t(relationship[,,k,l]);}}
  MSwithin <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      MSwithin[k,l] <- (sum(subs.relationship[,,k,l]**2)/2)/((N-1)*(N-2))}}
  # Estimate variance for actor effect #
  variance.actor <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      variance.actor[k,l] <- MSactor[k,l] - 0.5*((MSbetween[k,l]/(N-2))+(MSwithin[k,l]/N))
      if (variance.actor[k,l]<0) variance.actor[k,l]<- NA}}
  # Estimate variance for partner effect #
  variance.partner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      variance.partner[k,l] <- MSpartner[k,l] - 0.5*((MSbetween[k,l]/(N-2))+(MSwithin[k,l]/N))
      if (variance.partner[k,l]<0) variance.partner[k,l]<- NA}}
  # Estimate covariance for actor-partner effect #
  covariance.actorpartner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      covariance.actorpartner[k,l] <- MCP[k,l] - 0.5*((MSbetween[k,l]/(N-2))-(MSwithin[k,l]/N))}}
  # Estimate variance for relationship effect #
  variance.relationship <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      variance.relationship[k,l] <- 0.5*(MSbetween[k,l]+MSwithin[k,l])
      if (variance.relationship[k,l]<0) variance.relationship[k,l]<- NA}}
  # Estimate covariance for relationship effect #
  covariance.relationship <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      covariance.relationship[k,l] <- 0.5*(MSbetween[k,l]-MSwithin[k,l])}}
  res <- list(MSbetween=MSbetween,MSwithin=MSwithin,actor.variance=variance.actor,partner.variance=variance.partner,
              relationship.variance=variance.relationship,actorpartner.covariance=covariance.actorpartner,
              dyadic.covariance=covariance.relationship)
  
  class(res) <- "srmRRVars"
  res
}

print.srmRRVars <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("SRM Variances: Round Robin Design")
  cat("\n\n")
  T <- dim(x$actor.variance)[1]
  G <- dim(x$actor.variance)[2]
  times <- unlist(dimnames(x$actor.variance)[1])
  groups <- unlist(dimnames(x$actor.variance)[2])
  group <- rep(groups,T)
  time <- rep(times,G)
  actor <- array(x$actor.variance)
  partner <- array(x$partner.variance)
  relationship <- array(x$relationship.variance)
  actor.partner <- array(x$actorpartner.covariance)
  dyadic <- array(x$dyadic.covariance)
  results<-data.frame(group,time,actor,partner,relationship,actor.partner,dyadic)
  print(results)
  cat("\n\n")
  invisible(x)
}

# Function to compute relative variance due to actor, partner and relationship effects #

SRM.relative.variances <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  variance.relationship <- variances[[5]]
  prop.variance.actor <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  prop.variance.partner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  prop.variance.relationship <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  total.variance <- rowSums(data.frame(variance.actor,variance.partner,variance.relationship),na.rm=TRUE)
  prop.variance.actor <- variance.actor/total.variance
  prop.variance.partner <- variance.partner/total.variance
  prop.variance.relationship <- variance.relationship/total.variance
  res<-list(actor.variance=prop.variance.actor,partner.variance=prop.variance.partner,
            relationship.variance=prop.variance.relationship)
  class(res) <- "srmRRrelVars"
  res
}

print.srmRRrelVars <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("SRM Relative Variances: Round Robin Design")
  cat("\n\n")
  T <- dim(x$actor.variance)[1]
  G <- dim(x$actor.variance)[2]
  times <- unlist(dimnames(x$actor.variance)[1])
  groups <- unlist(dimnames(x$actor.variance)[2])
  group <- rep(groups,T)
  time <- rep(times,G)
  actor <- array(x$actor.variance)
  partner <- array(x$partner.variance)
  relationship <- array(x$relationship.variance)
  results<-data.frame(group,time,actor,partner,relationship)
  rownames(results) <- ''
  print(results)
  cat("\n\n")
  invisible(x)
}

# Function to compute reliability for the estimates of the actor and partner variances #

SRM.reliability <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  variance.relationship <- variances[[5]]
  covariance.relationship <- variances[[7]]
  actor.reliability <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  actor.reliability <- variance.actor/(variance.actor + variance.relationship/(N-1)-
                                         covariance.relationship/(N-1)**2)
  partner.reliability <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  partner.reliability <- variance.partner/(variance.partner + variance.relationship/(N-1)-
                                             covariance.relationship/(N-1)**2)
  list(actor.reliability=actor.reliability,partner.reliability=partner.reliability)
}

# Testing SRM effects by means of Jackknife tests #

SRM.jackknife <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  # Is group(s) size greater than or equal to 5 individuals? #
  if (N < 5)
    return("Error: Group(s) size must be greater than or equal to 5 individuals for carrying out the jackknife.")
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  jackknife <- array(dim=c(N,5,T,G),dimnames=list(names,source,times,groups),0.)
  jackknife.t.statistic <- array(dim=c(T,dim(jackknife)[2],G),dimnames=list(times,source,groups))
  jackknife.p.value <- array(dim=c(T,dim(jackknife)[2],G),dimnames=list(times,source,groups))
  jackknife.mean <- array(dim=c(T,dim(jackknife)[2],G),0.)
  jackknife.variance <- array(dim=c(T,dim(jackknife)[2],G),0.)
  Xtemp <- array(dim=c(N-1,N-1,T,G),0.)
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  variance.relationship <- variances[[5]]
  actorpartner.covariance <- variances[[6]]
  dyadic.covariance <- variances[[7]]
  for (i in 1:N){
    Xtemp <- as.numeric(X$data[-i,-i,,])
    Xtemp <- prepare.SRM.matrix(Xtemp,(N-1),T,G)
    variancestemp <- SRM.variances(Xtemp)
    for (k in 1:T){
      for (l in 1:G){
        jackknife[i,1,k,l] <- N*variance.actor[k,l]-variancestemp$actor.variance[k,l]*(N-1)
        jackknife[i,2,k,l] <- N*variance.partner[k,l]-variancestemp$partner.variance[k,l]*(N-1)
        jackknife[i,3,k,l] <- N*variance.relationship[k,l]-variancestemp$relationship.variance[k,l]*(N-1)
        jackknife[i,4,k,l] <- N*actorpartner.covariance[k,l]-variancestemp$actorpartner.covariance[k,l]*(N-1)
        jackknife[i,5,k,l] <- N*dyadic.covariance[k,l]-variancestemp$dyadic.covariance[k,l]*(N-1)}}}
  for (k in 1:T){
    for (l in 1:G){
      jackknife.mean[k,,l] <- apply(jackknife[,,k,l],2,mean)
      jackknife.variance[k,,l] <- apply(jackknife[,,k,l],2,var)/N
      jackknife.t.statistic[k,,l] <- jackknife.mean[k,,l]/sqrt(jackknife.variance[k,,l])}}
  df <- N-1
  for (k in 1:T){
    for (m in 1:(dim(jackknife)[2])){
      for (l in 1:G){
        if ( !is.na(jackknife.t.statistic[k,m,l])){
          if (sign(jackknife.t.statistic[k,m,l])>0)
            jackknife.p.value[k,m,l] <- pt(jackknife.t.statistic[k,m,l],df,lower.tail=FALSE)
          else jackknife.p.value[k,m,l] <- pt(jackknife.t.statistic[k,m,l],df)}}}}
  list(t.statistic=jackknife.t.statistic,df=df,two.tailed.p.value=2*jackknife.p.value)
}

# Function to estimate generalized reciprocity #

SRM.generalized.reciprocity <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  effects <- SRM.effects(X)
  actor <- effects$actor.effects
  partner <- effects$partner.effects
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  MSbetween <- variances[[1]]
  MSwithin <- variances[[2]]
  generalized.reciprocity <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      if (is.na(variance.actor[k,l]) | is.na(variance.partner[k,l])) generalized.reciprocity[k,l] <- NA
      else{
        
        generalized.reciprocity[k,l] <- (sum(actor[,k,l]*partner[,k,l])/(N-1) - MSbetween[k,l]/(2*(N-2)) +
                                           MSwithin[k,l]/(2*N))/sqrt(variance.actor[k,l]*variance.partner[k,l])
        if ((generalized.reciprocity[k,l] > 1.0) | (generalized.reciprocity[k,l] < -1.0))
          generalized.reciprocity[k,l] <- sign(generalized.reciprocity[k,l])*1.0}}}
  return(generalized.reciprocity)
}

# Compute dyadic.reciprocity #

SRM.dyadic.reciprocity <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  variances <- SRM.variances(X)
  MSbetween <- variances[[1]]
  MSwithin <- variances[[2]]
  dyadic.reciprocity <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      dyadic.reciprocity[k,l] <- (MSbetween[k,l] - MSwithin[k,l])/(MSbetween[k,l] + MSwithin[k,l])}}
  return(dyadic.reciprocity)
}

# Testing dyadic covariance by means of a F ratio #

SRM.dyadic.covariance.Ftest <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  variances <- SRM.variances(X)
  MSbetween <- variances[[1]]
  MSwithin <- variances[[2]]
  F.value <- array(dim=c(T,G),dimnames=list(times,groups))
  for (k in 1:T){
    for (l in 1:G){
      F.value[k,l] <- MSbetween[k,l]/MSwithin[k,l]}}
  df1 <- (N-1)*(N-2)/2-1
  df2 <- (N-1)*(N-2)/2
  F.p.value <- array(dim=c(T,G),dimnames=list(times,groups))
  for (k in 1:T){
    for (l in 1:G){
      if (sign(F.value[k,l])>0)
        F.p.value[k,l] <- pf(F.value[k,l],df1,df2,lower.tail=FALSE)
      else F.p.value[k,l] <- pf(F.value[k,l],df1,df2)}}
  res=list(MSb=MSbetween,MSw=MSwithin,F.value=F.value,df1=df1,df2=df2,p.value=F.p.value)
  class(res) <- "srmCovFTest"
  res
}

print.srmCovFTest <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("SRM Dyadic Covariance F Test: Round Robin Design")
  cat("\n\n")
  T <- dim(x$MSb)[1]
  G <- dim(x$MSb)[2]
  times <- unlist(dimnames(x$MSb)[1])
  groups <- unlist(dimnames(x$MSb)[2])
  group <- rep(groups,T)
  time <- rep(times,G)
  MSb <- x[[1]]
  MSw <- x[[2]]
  df1 <- x[[4]]
  df2 <- x[[5]]
  F <- x[[3]]
  pval <- x[[6]]
  results<-data.frame(group,time,MSb,MSw,df1,df2,F, pval)
  colnames(results)<-c("Group","Time","MSb","MSw","df1","df2","F value",if (sign(F)>0) "Pr(> F)" else "Pr(< F)")
  rownames(results) <- ''
  print(results)
  cat("\n\n")
  invisible(x)
}

# Testing SRM effects by means of a between-group t Test #

SRM.between.groups.tTest <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  # Are there more than 1 group? #
  if (G < 2)
    return("Error: Number of groups analyzed should be greater than 1.")
  
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  variance.relationship <- variances[[5]]
  actorpartner.covariance <- variances[[6]]
  dyadic.covariance <- variances[[7]]
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  between.group.t.statistic <- array(dim=c(T,5),dimnames=list(times,source))
  between.group.p.value <- array(dim=c(T,5),dimnames=list(times,source))
  actor.variance.mean <- apply(variance.actor,1,mean)
  partner.variance.mean <- apply(variance.partner,1,mean)
  relationship.variance.mean <- apply(variance.relationship,1,mean)
  actorpartner.covariance.mean <- apply(actorpartner.covariance,1,mean)
  dyadic.covariance.mean <- apply(dyadic.covariance,1,mean)
  actor.variance.var <- apply(variance.actor,1,var)
  partner.variance.var <- apply(variance.partner,1,var)
  relationship.variance.var <- apply(variance.relationship,1,var)
  actorpartner.covariance.var <- apply(actorpartner.covariance,1,var)
  dyadic.covariance.var <- apply(dyadic.covariance,1,var)
  for (k in 1:T){
    between.group.t.statistic[k,1] <- actor.variance.mean[k]/sqrt(actor.variance.var[k]/G)
    between.group.t.statistic[k,2] <- partner.variance.mean[k]/sqrt(partner.variance.var[k]/G)
    between.group.t.statistic[k,3] <- relationship.variance.mean[k]/sqrt(relationship.variance.var[k]/G)
    between.group.t.statistic[k,4] <- dyadic.covariance.mean[k]/sqrt(dyadic.covariance.var[k]/G)
    between.group.t.statistic[k,5] <- actorpartner.covariance.mean[k]/sqrt(actorpartner.covariance.var[k]/G)}
  df <- G-1
  for (k in 1:T){ 
    for (m in 1:(dim(between.group.t.statistic)[2])){
      if ( !is.na(between.group.t.statistic[k,m])){
        if (sign(between.group.t.statistic[k,m]) > 0)
          between.group.p.value[k,m] <- pt(between.group.t.statistic[k,m],df,lower.tail=FALSE)
        else between.group.p.value[k,m] <- pt(between.group.t.statistic[k,m],df)}}}
  list(t.statistic=between.group.t.statistic,df=df,two.tailed.p.value=2*between.group.p.value)
}

# Testing SRM effects by means of a within-group t Test #  

SRM.within.groups.tTest <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  # Computing some terms that will be used below #
  h1 <- 2*(G**2*N**6-7*G**2*N**5+18*G**2*N**4+10*G*N**4-20*G**2*N**3-
             46*G*N**3+8*G**2*N**2+70*G*N**2+24*N**2-36*G*N-48*N+32)
  h2 <- 2*(G**2*N**5-5*G**2*N**4+2*G*N**4+8*G**2*N**3-6*G*N**3-
             4*G**2*N**2+14*G*N**2+8*N**2-20*G*N-16*N+32)
  h3 <- 4*(G**2*N**5-5*G**2*N**4+8*G**2*N**3+12*G*N**3-
             4*G**2*N**2-38*G*N**2+28*G*N+32*N-32)
  h4 <- (N**2)*((N-2)**2)*(G*N-G+2)*(G*N**2-3*G*N+4)*(G*N**2-3*G*N+
                                                        2*G+4)
  h5 <- G**3*N**7-7*G**3*N**6+19*G**3*N**5+10*G**2*N**5-25*G**3*N**4-
    54*G**2*N**4-4*G*N**4+16*G**3*N**3+120*G**2*N**3+44*G*N**3-
    4*G**3*N**2-132*G**2*N**2-124*G*N**2-16*N**2+56*G**2*N+
    168*G*N+32*N-64*G-64
  h6 <- G**3*N**7-7*G**3*N**6-2*G**2*N**6+19*G**3*N**5+26*G**2*N**5-
    25*G**3*N**4-100*G**2*N**4-20*G*N**4+16*G**3*N**3+176*G**2*N**3+
    124*G*N**3-4*G**3*N**2-156*G**2*N**2-236*G*N**2-48*N**2+56*G**2*N+
    200*G*N+96*N-64*G-64
  # Estimating mean and variance of round robin parameters #
  variances <- SRM.variances(X)
  variance.actor <- unlist(variances[[3]])
  variance.partner <- unlist(variances[[4]])
  variance.relationship <- unlist(variances[[5]])
  actorpartner.covariance <- unlist(variances[[6]])
  apcor <- unlist(SRM.generalized.reciprocity(X))
  dyadcor <- unlist(SRM.dyadic.reciprocity(X))
  dyadic.covariance <- unlist(variances[[7]])
  actor.variance.mean <- apply(variance.actor,1,function(x) mean(x,na.rm=TRUE))
  partner.variance.mean <- apply(variance.partner,1,function(x) mean(x,na.rm=TRUE))
  relationship.variance.mean <- apply(variance.relationship,1,function(x) mean(x,na.rm=TRUE))
  actorpartner.covariance.mean <- apply(actorpartner.covariance,1,function(x) mean(x,na.rm=TRUE))
  dyadic.covariance.mean <- apply(dyadic.covariance,1,function(x) mean(x,na.rm=TRUE))
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  SRM.variance.parameters <- array(dim=c(T,5),dimnames=list(times,source),0.)
  within.group.t.statistic <- array(dim=c(T,5),dimnames=list(times,source),0.)
  within.group.p.value <- array(dim=c(T,5),dimnames=list(times,source),0.)
  for (k in 1:T){
    SRM.variance.parameters[k,1] <- ((2*(actor.variance.mean[k])**2)/(G*N-G+2))+(((4*actor.variance.mean[k])*
                                                                                    ((N-1)*relationship.variance.mean[k]+dyadic.covariance.mean[k]))/(N*(N-2)*(G*N-G+2)))+
      ((h1*(relationship.variance.mean[k])**2+h2*(dyadic.covariance.mean[k])**2+h3*
          (relationship.variance.mean[k]*dyadic.covariance.mean[k]))/h4)
    SRM.variance.parameters[k,2] <- ((2*(partner.variance.mean[k])**2)/(G*N-G+2))+(((4*partner.variance.mean[k])*
                                                                                      ((N-1)*relationship.variance.mean[k]+dyadic.covariance.mean[k]))/(N*(N-2)*(G*N-G+2)))+
      ((h1*(relationship.variance.mean[k])**2+h2*(dyadic.covariance.mean[k])**2+h3*
          (relationship.variance.mean[k]*dyadic.covariance.mean[k]))/h4)
    SRM.variance.parameters[k,3] <- SRM.variance.parameters[k,5] <- ((2*((G*N**2)-(3*G*N)+G+4)*
                                                                        (relationship.variance.mean[k]**2+dyadic.covariance.mean[k]**2))/(((G*N**2)-(3*G*N)+4)*
                                                                                                                                            ((G*N**2)-(3*G*N)+2*G+4)))+((4*G*relationship.variance.mean[k]*dyadic.covariance.mean[k])/(((G*N**2)-(3*G*N)+4)*
                                                                                                                                                                                                                                         ((G*N**2)-(3*G*N)+2*G+4)))
    SRM.variance.parameters[k,4] <- ((G*N-G-2)*actorpartner.covariance.mean[k]**2+G*(N-1)*
                                       actor.variance.mean[k]*partner.variance.mean[k])/((G*N-G+2)*(G*N-G-1))+(h5*relationship.variance.mean[k]**2+
                                                                                                                 h6*dyadic.covariance.mean[k]**2+h3*relationship.variance.mean[k]*dyadic.covariance.mean[k])/
      ((G*N-G-1)*h4)+(G*(N-1)*(actor.variance.mean[k]+partner.variance.mean[k])*((N-1)*
                                                                                   relationship.variance.mean[k]+dyadic.covariance.mean[k]))/(N*(N-2)*(G*N-G+2)*(G*N-G-1))+
      (2*(G*N-G-2)*actorpartner.covariance.mean[k]*(relationship.variance.mean[k]+(N-1)*dyadic.covariance.mean[k]))/
      (N*(N-2)*(G*N-G+2)*(G*N-G-1))}
  SRM.relative.parameters <- array(dim=c(T,5),dimnames=list(times,source),0.)
  SRM.mean.parameters <- array(dim=c(T,5),dimnames=list(times,source),c(actor.variance.mean,
                                                                        partner.variance.mean,relationship.variance.mean,actorpartner.covariance.mean,dyadic.covariance.mean))
  SRM.relative.parameters[1:3] <- SRM.mean.parameters[1:3]/sum(unlist(SRM.mean.parameters[1:3]),na.rm=TRUE)
  SRM.relative.parameters[4] <- apply(apcor,1,function(x) mean(x,na.rm=TRUE))
  SRM.relative.parameters[5] <- apply(dyadcor,1,function(x) mean(x,na.rm=TRUE))  
  SRM.stderror.parameters <- sqrt(SRM.variance.parameters)
  for (k in 1:T){
    within.group.t.statistic[k,1] <- actor.variance.mean[k]/sqrt(SRM.variance.parameters[k,1])
    within.group.t.statistic[k,2] <- partner.variance.mean[k]/sqrt(SRM.variance.parameters[k,2])
    within.group.t.statistic[k,3] <- relationship.variance.mean[k]/sqrt(SRM.variance.parameters[k,3])
    within.group.t.statistic[k,4] <- actorpartner.covariance.mean[k]/sqrt(SRM.variance.parameters[k,4])
    within.group.t.statistic[k,5] <- dyadic.covariance.mean[k]/sqrt(SRM.variance.parameters[k,5])}
  df <- G*(N-1)
  for (k in 1:T){ 
    for (m in 1:(dim(within.group.t.statistic)[2])){
      if ( !is.na(within.group.t.statistic[k,m])){
        if (sign(within.group.t.statistic[k,m]) > 0)
          within.group.p.value[k,m] <- pt(within.group.t.statistic[k,m],df,lower.tail=FALSE)
        else within.group.p.value[k,m] <- pt(within.group.t.statistic[k,m],df)}}}
  res <- list(estimates=SRM.mean.parameters,standard.values=SRM.relative.parameters,standard.error=SRM.stderror.parameters,
              t.value=within.group.t.statistic,df=df,p.value=within.group.p.value)
  class(res) <-"srmWithintTest"
  res
}

print.srmWithintTest <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Within group t-Test: Round Robin Design")
  cat("\n\n")
  estimates <- as.numeric(x$estimates)
  std.value <- as.numeric(x$standard.values)
  std.error <- as.numeric(x$standard.error)
  t <- as.numeric(x$t.value)
  df <- as.numeric(x$df)
  p.value<- as.numeric(x$p.value)
  results<-data.frame(estimates,std.value,std.error,t,df,p.value)
  results[is.na(results)]<-NA
  results[is.na(results[,2]),3:6]<-NA
  rownames(results) <- c('actor var','partner var','relationship var','actor-partner cov','dyadic cov')
  print(results)
  cat("\n\n")
  invisible(x)
}

################################################################################
#                      R FUNCTIONS FOR  SRM BLOCK DESIGNS                     #
################################################################################

# Function for preparing data in order to be analized as a SRM block design #

prepare.block.matrix <- function (X,subgroup1,subgroup2,numb.times,numb.groups,
                                  names1=NULL,names2=NULL,times=NULL,groups=NULL,
                                  subgroups=NULL){
  numb.individuals <- subgroup1 + subgroup2
  if (is.null(names1)) names1 <- paste('Ind.',c(1:(subgroup1)))
  else names1 <- names1
  if (is.null(names2)) names2 <- paste('Ind.',c((numb.individuals-subgroup2+1):numb.individuals))
  else names2 <- names2
  if (is.null(times)) times <- paste('Time',c(1:numb.times))
  else times <- times
  if (is.null(groups)) groups <- paste('Group',c(1:numb.groups))
  else groups <- groups
  if (is.null(subgroups)) subgroups <- paste('Subgroup',c(1,2))
  else subgroups <- subgroups
  X <- array(X,c(numb.individuals,numb.individuals,numb.times,numb.groups),
             dimnames=list(c(names1,names2),c(names1,names2),times,groups))
  temp <- array(0,c(numb.individuals,numb.individuals,numb.times,numb.groups))
  for (i in 1:(dim(X)[1])){
    for (j in 1:(dim(X)[2])){
      for (k in 1:(dim(X)[3])){
        for (l in 1:(dim(X)[4])){
          temp[i,j,k,l] <- X[j,i,k,l]}}}}
  count <- 1;
  count2 <- 1;
  for (k in 1:(dim(X)[3])){
    for (i in 1:(dim(X)[1])){
      if (count <= (dim(X)[3])){ 
        X[count2,,count,] <- temp[i,,k,]
        count <- count + 1}
      else {
        count <- 1
        count2 <- count2 + 1
        X[count2,,count,] <- temp[i,,k,]
        count <- count + 1}}} 
  res<-list(data=X,individuals1=subgroup1,individuals2=subgroup2,
            individuals=numb.individuals,ntimes=numb.times,ngroups=numb.groups,
            names1=names1,names2=names2,times=times,groups=groups,
            subgroups=subgroups)
  class(res) <- 'srmBlockMatrix'
  res
}

# Function for obtaining SRM effects in Block designs #

block.SRM.effects <- function (X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  groups<- X$groups
  times <- X$times
  # Compute row means #
  time.rowmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      for (i in 1:N){
        if (i <= G1) time.rowmeans[i,k,l] <- sum(X$data[i,,k,l])/G2
        else time.rowmeans[i,k,l] <- sum(X$data[i,,k,l])/G1}}}
  # Compute column means #
  time.colmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      for (j in 1:N){
        if (j <= G1) time.colmeans[j,k,l] <- sum(X$data[,j,k,l])/G2
        else time.colmeans[j,k,l] <- sum(X$data[,j,k,l])/G1}}}
  # Compute grand mean #
  time.grandmean <- array(dim=c(2,T,G),0.)
  for (i in 1:N){
    for (j in 1:N){
      for (k in 1:T){
        for (l in 1:G){
          if ((i <= G1) && (j > G1)){
            time.grandmean[1,k,l] <- time.grandmean[1,k,l] + X$data[i,j,k,l]}
          if ( (i > G1) && (j <= G1)){
            time.grandmean[2,k,l] <- time.grandmean[2,k,l] + X$data[i,j,k,l]}}}}}
  time.grandmean <- time.grandmean/(G1*G2)
  # Estimate actor effects for all individuals #
  actor <- array(dim=c(N,T,G),dimnames=list(c(names1,names2),times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      for (i in 1:N){
        if (i <= G1){
          actor[i,k,l] <- time.rowmeans[i,k,l]-time.grandmean[1,k,l]}
        else actor[i,k,l] <- time.rowmeans[i,k,l]-time.grandmean[2,k,l]}}}
  # Estimate partner effects for all individuals #
  partner <- array(dim=c(N,T,G),dimnames=list(c(names1,names2),times,groups),0.)
  for (j in 1:N){
    for (k in 1:T){
      for (l in 1:G){
        if (j > G1){
          partner[j,k,l] <- time.colmeans[j,k,l]-time.grandmean[1,k,l]}
        if (j <= G1){
          partner[j,k,l] <- time.colmeans[j,k,l]-time.grandmean[2,k,l]}}}}
  # Estimate relationships effects for all individuals #
  relationship <- array(dim=c(N,N,T,G),
                        dimnames=list(c(names1,names2),c(names1,names2),times,groups),0.)
  for (i in 1:N){
    for (j in 1:N){
      for (k in 1:T){
        for (l in 1:G){
          if ((i <= G1) && (j > G1))
          {relationship[i,j,k,l] <- X$data[i,j,k,l]-actor[i,k,l]-partner[j,k,l]-time.grandmean[1,k,l]}
          if ( (i > G1) && (j <= G1))
          {relationship[i,j,k,l] <- X$data[i,j,k,l]-actor[i,k,l]-partner[j,k,l]-time.grandmean[2,k,l]}}}}}
  list(actor.effects=actor,partner.effects=partner,relationship.effects=relationship)
}

# Function to estimate SRM variances in block designs #

block.SRM.variances <- function(X){ 
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  groups <- X$groups
  subgroups <- X$subgroups
  times <- X$times
  effects <- block.SRM.effects(X)
  actor <- effects$actor.effects
  partner <- effects$partner.effects
  relationship <- effects$relationship.effects 
  # Estimate variance for relationship effect #
  variance.relationship <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (i in 1:N){
    for (j in 1:N){
      for (k in 1:T){
        for (l in 1:G){
          if ((i <= G1) && (j > G1))
          {variance.relationship[1,k,l] <- variance.relationship[1,k,l]+relationship[i,j,k,l]**2}
          if ( (i > G1) && (j <= G1)) 
          {variance.relationship[2,k,l] <- variance.relationship[2,k,l]+relationship[i,j,k,l]**2}}}}}
  variance.relationship <- variance.relationship/((G1-1)*(G2-1))
  # Estimate variance for actor effect #
  variance.actor <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (i in 1:N){
    for (k in 1:T){
      for (l in 1:G){
        if ( i <= G1 )
        {variance.actor[1,k,l] <- variance.actor[1,k,l]+actor[i,k,l]**2}
        if ( i > G1 ) 
        {variance.actor[2,k,l] <- variance.actor[2,k,l]+actor[i,k,l]**2}}}}
  variance.actor <- (variance.actor/(G1-1))-(variance.relationship/G2)
  for (k in 1:T){
    for (l in 1:G){
      for (m in 1:length(subgroups)){
        if ( sign(variance.actor[m,k,l]) < 0. ) variance.actor[m,k,l] <- 0.}}} 
  # Estimate variance for partner effect #
  variance.partner <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (j in 1:N){
    for (k in 1:T){
      for (l in 1:G){
        if ( j <= G2 )
        {variance.partner[1,k,l] <- variance.partner[1,k,l]+partner[j,k,l]**2}
        if ( j > G1 ) 
        {variance.partner[2,k,l] <- variance.partner[2,k,l]+partner[j,k,l]**2}}}}
  variance.partner[1,,] <- (variance.partner[1,,]/(G2-1))-(variance.relationship[2]/G2)
  variance.partner[2,,] <- (variance.partner[2,,]/(G2-1))-(variance.relationship[1]/G2)
  for (k in 1:T){
    for (l in 1:G){
      for (m in 1:length(subgroups)){
        if ( sign(variance.partner[m,k,l]) < 0. ) variance.partner[m,k,l] <- 0.}}} 
  # Estimate covariance for actor-partner effect #
  dyads1 <- array(dim=c(G1*G2,T,G),0.)
  dyads2 <- array(dim=c(G1*G2,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      count1 <- 1
      count2 <- 1
      for (i in 1:N){
        for (j in 1:N){
          if ((i <= G1) && (j > G1))
          {dyads1[count1,k,l] <- relationship[i,j,k,l]
          count1 <- count1 +1}
          if ( (i > G1) && (j <= G1))
          {dyads2[count2,k,l] <- relationship[i,j,k,l]
          count2 <- count2 +1}}}}}
  covariance.actorpartner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      covariance.actorpartner[k,l] <- cov(dyads1[,k,l],dyads2[,k,l])}}
  # Estimate covariance for relationship effect #
  covariance.relationship <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  crosspod.relationship <- array(dim=c(N,N,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      crosspod.relationship[,,k,l] <- relationship[,,k,l] * t(relationship[,,k,l])}}
  for (i in 1:N){
    for (k in 1:T){
      for (l in 1:G){
        if ( i <= G1 )
        {covariance.relationship[1,k,l] <- covariance.relationship[1,k,l]+actor[i,k,l]*partner[i,k,l]}
        if ( i > (N-G2) ) 
        {covariance.relationship[2,k,l] <- covariance.relationship[2,k,l]+actor[i,k,l]*partner[i,k,l]}}}}
  for (k in 1:T){
    for (l in 1:G){
      covariance.relationship[1,k,l] <- (covariance.relationship[1,k,l]/(G2-1))-(sum(crosspod.relationship[,,k,l])/
                                                                                   2/((G2-1)*(G1-1)))/G2
      covariance.relationship[2,k,l] <- (covariance.relationship[2,k,l]/(G1-1))-(sum(crosspod.relationship[,,k,l])/
                                                                                   2/((G2-1)*(G1-1)))/G2}}
  list(actor.variance=variance.actor,partner.variance=variance.partner,relationship.variance=variance.relationship,
       actorpartner.covariance=covariance.actorpartner,dyadic.covariance=covariance.relationship)
}

# Function to compute relative variance due to actor, partner and relationship effects i block designs #

block.SRM.relative.variances <- function(X){ 
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  groups <- X$groups
  subgroups <- X$subgroups
  times <- X$times
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  prop.variance.actor <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  prop.variance.partner <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  prop.variance.relationship <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  total.variance <- variance.actor + variance.partner + variance.relationship
  prop.variance.actor <- variance.actor/total.variance
  prop.variance.partner <- variance.partner/total.variance
  prop.variance.relationship <- variance.relationship/total.variance
  list(relative.actor.variance=prop.variance.actor,relative.partner.variance=prop.variance.partner,
       relative.relationship.variance=prop.variance.relationship)
}

# Function to compute reliability for the estimates of the actor and partner variances in block designs #

block.SRM.reliability <- function(X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups <- X$subgroups
  groups <- X$sgroups
  times <- X$times  
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  actor.reliability <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0)
  partner.reliability <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      if ( (variance.actor[1,k,l] != 0) | (variance.relationship[1,k,l] != 0) )
        actor.reliability[1,k,l] <- variance.actor[1,k,l]/(variance.actor[1,k,l] + 
                                                             variance.relationship[1,k,l]/G1)
      if ( (variance.actor[2,k,l] != 0) | (variance.relationship[2,k,l] != 0) )
        actor.reliability[2,k,l] <- variance.actor[2,k,l]/(variance.actor[2,k,l] + 
                                                             variance.relationship[2,k,l]/G2)
      if ( (variance.partner[1,k,l] != 0) | (variance.relationship[1,k,l] != 0) )
        partner.reliability[1,k,l] <- variance.partner[1,k,l]/(variance.partner[1,k,l] + 
                                                                 variance.relationship[1,k,l]/G1)
      if ( (variance.partner[2,k,l] != 0) | (variance.relationship[2,k,l] != 0) )
        partner.reliability[2,k,l] <- variance.partner[2,k,l]/(variance.partner[2,k,l] + 
                                                                 variance.relationship[1,k,l]/G2)}}
  list(actor.reliability=actor.reliability,partner.reliability=partner.reliability)
}

# Function to estimate generalized reciprocity in block designs #

block.SRM.generalized.reciprocity <- function(X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups<- X$subgroups
  groups <- X$groups  
  times <- X$times  
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  covariance.relationship <- variances$dyadic.covariance
  generalized.reciprocity <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (m in 1:length(subgroups)){ 
    for (k in 1:T){
      for (l in 1:G){
        if ((variance.actor[m,k,l] == 0.) | (variance.partner[m,k,l] == 0.)) generalized.reciprocity[m,k,l] <- 0.
        else 
          generalized.reciprocity[m,k,l] <- covariance.relationship[m,k,l]/sqrt(variance.actor[m,k,l]*variance.partner[m,k,l])}}
    if ((generalized.reciprocity[m,k,l] > 1.0) | (generalized.reciprocity[m,k,l] < -1.0))
      generalized.reciprocity[m,k,l] <- sign(generalized.reciprocity[m,k,l])*1.0}
  return(generalized.reciprocity)
}

# Compute dyadic.reciprocity in block designs #

block.SRM.dyadic.reciprocity <- function(X,symmetric=FALSE){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups<- X$subgroups
  groups <- X$groups  
  times <- X$times
  effects <- block.SRM.effects(X)
  relationship <- effects$relationship.effects 
  variances <- block.SRM.variances(X)
  variance.relationship <- variances$relationship.variance
  # Estimate covariance for actor-partner effect #
  dyads1 <- array(dim=c(G1*G2,T,G),0.)
  dyads2 <- array(dim=c(G1*G2,T,G),0.)
  for (k in 1:T){
    for (l in 1:G){
      count1 <- 1
      count2 <- 1
      for (i in 1:N){
        for (j in 1:N){
          if ((i <= G1) && (j > G1))
          {dyads1[count1,k,l] <- relationship[i,j,k,l]
          count1 <- count1 +1}
          if ( (i > G1) && (j <= G1))
          {dyads2[count2,k,l] <- relationship[i,j,k,l]
          count2 <- count2 +1}}}}}
  dyadic.reciprocity <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  if (symmetric == FALSE) {
    for (k in 1:T){
      for (l in 1:G){
        dyadic.reciprocity[k,l] <- cor(dyads1[,k,l],dyads2[,k,l])}}}
  else {
    for (k in 1:T){
      for (l in 1:G){
        dyadic.reciprocity[k,l] <- cor(dyads1[,k,l],dyads2[,k,l])*(sqrt(variance.relationship[1,k,l]*
                                                                          variance.relationship[2,k,l])/((variance.relationship[1,k,l]+variance.relationship[2,k,l])/2))}}}
  return(dyadic.reciprocity)
}

# Testing SRM effects in block designs by means of Jackknife tests #

block.SRM.jackknife <-  function(X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups<- X$subgroups
  groups <- X$groups  
  times <- X$times
  # Is group(s) size greater than or equal to 6 individuals? #
  if (N < 6)
    return("Error: Group(s) size must be greater than or equal to 6 individuals for carrying out the jackknife for a block design.")
  if (G1 < 3 | G2 < 3)
    return("Error: Subgroups must be greater than or equal to 3 individuals for carrying out the jackknife for a block design.")
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  jackknife <- array(dim=c(length(subgroups),N,5,T,G),dimnames=list(subgroups,c(names1,names2),
                                                                    source,times,groups),0.)
  jackknife.t.statistic <- array(dim=c(length(subgroups),T,dim(jackknife)[3],
                                       G),dimnames=list(subgroups,times,source,groups),0.)
  jackknife.p.value <- array(dim=c(length(subgroups),T,dim(jackknife)[3],G),
                             dimnames=list(subgroups,times,source,groups),NA)
  jackknife.mean <- array(dim=c(length(subgroups),T,dim(jackknife)[3],G),0.)
  jackknife.variance <- array(dim=c(length(subgroups),T,dim(jackknife)[3],G),0.)
  Xtemp <- array(dim=c(N-1,N-1,T,G),0.)
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  actorpartner.covariance <- variances$actorpartner.covariance
  dyadic.covariance <- variances$dyadic.covariance
  for (m in 1:length(subgroups)){
    count <- 1
    for (i in 1:N){
      if (count <= G1) {
        ind.temp.1 <- G1-1
        ind.temp.2 <- G2}
      else {
        ind.temp.1 <- G1
        ind.temp.2 <- G2-1}
      if (count <= G1){ 
        term1 <- G1 
        term2 <- ind.temp.1}
      else {term1 <- G2 
      term2 <- ind.temp.2}
      count <- count + 1
      Xtemp <- prepare.block.matrix(as.numeric(X$data[-i,-i,,]),N-1,N-1,T,G)
      variancestemp <- block.SRM.variances(Xtemp)
      for (k in 1:T){
        for (l in 1:G){
          jackknife[m,i,1,k,l] <- term1*variance.actor[m,k,l]-variancestemp$actor.variance[m,k,l]*term2
          jackknife[m,i,2,k,l] <- term1*variance.partner[m,k,l]-variancestemp$partner.variance[m,k,l]*term2
          jackknife[m,i,3,k,l] <- term1*variance.relationship[m,k,l]-variancestemp$relationship.variance[m,k,l]*term2
          jackknife[m,i,4,k,l] <- term1*dyadic.covariance[m,k,l]-variancestemp$dyadic.covariance[m,k,l]*term2
          jackknife[m,i,5,k,l] <- term1*actorpartner.covariance[k,l]-variancestemp$actorpartner.covariance[k,l]*term2}}}}
  for (m in 1:length(subgroups)){
    for (k in 1:T){
      for (l in 1:G){
        jackknife.mean[m,k,,l] <- apply(jackknife[m,,,k,l],2,mean)
        jackknife.variance[m,k,,l] <- apply(jackknife[m,,,k,l],2,var)/(N)
        jackknife.t.statistic[m,k,,l] <- jackknife.mean[m,k,,l]/sqrt(jackknife.variance[m,k,,l])}}}
  df1 <- G1-1
  df2 <- G2-2
  df <- array(dim=c(2,1),dimnames=list(subgroups),c(df1,df2))
  for (m in 1:length(subgroups)){    
    for (k in 1:T){
      for (n in 1:(dim(jackknife)[3])){
        for (l in 1:G){
          if (sign(jackknife.t.statistic[m,k,n,l])>0){
            if (m == 1)
              jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df1,lower.tail=FALSE)
            if (m == 2)
              jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df2,lower.tail=FALSE)}
          else{
            if (m == 1)
              jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df1)
            if (m == 2)
              jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df2)}}}}}
  list(t.statistic=jackknife.t.statistic,df=df,two.tailed.p.value=2*jackknife.p.value)
}

# Testing SRM effects in block designs by means of a between-group t Test #

block.SRM.between.groups.tTest <-function(X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups<- X$subgroups
  groups <- X$groups  
  times <- X$times
  # Are there more than 1 group? #
  if (G < 2)
    return("Error: Number of groups analyzed should be greater than 1.")
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  actorpartner.covariance <- variances$actorpartner.covariance
  dyadic.covariance <- variances$dyadic.covariance
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  between.group.t.statistic <- array(dim=c(2,T,5),dimnames=list(subgroups,times,source))
  between.group.p.value <- array(dim=c(2,T,5),dimnames=list(subgroups,times,source))
  actor.variance.mean <- array(dim=c(2,T))
  partner.variance.mean <- array(dim=c(2,T))
  relationship.variance.mean <- array(dim=c(2,T)) 
  actorpartner.covariance.mean <- array(dim=c(T))
  dyadic.covariance.mean <- array(dim=c(2,T))
  actor.variance.var <- array(dim=c(2,T))
  partner.variance.var <- array(dim=c(2,T))
  relationship.variance.var <- array(dim=c(2,T)) 
  actorpartner.covariance.var <- array(dim=c(T))
  dyadic.covariance.var <- array(dim=c(2,T))
  for (k in 1:T){
    actor.variance.mean[,k] <- apply(variance.actor[,k,],1,mean)
    partner.variance.mean[,k] <-  apply(variance.partner[,k,],1,mean)
    relationship.variance.mean[,k] <-  apply(variance.relationship[,k,],1,mean)
    actorpartner.covariance.mean[k] <-  mean(actorpartner.covariance[k,])
    dyadic.covariance.mean[,k] <-  apply(dyadic.covariance[,k,],1,mean)
    actor.variance.var[,k] <- apply(variance.actor[,k,],1,var)
    partner.variance.var[,k] <-  apply(variance.partner[,k,],1,var)
    relationship.variance.var[,k] <-  apply(variance.relationship[,k,],1,var)
    actorpartner.covariance.var[k] <-  var(actorpartner.covariance[k,])
    dyadic.covariance.var[,k] <-  apply(dyadic.covariance[,k,],1,var)}
  for (m in 1:length(subgroups)){
    for (k in 1:T){
      between.group.t.statistic[m,k,1] <- actor.variance.mean[m,k]/sqrt(actor.variance.var[m,k]/G)
      between.group.t.statistic[m,k,2] <- partner.variance.mean[m,k]/sqrt(partner.variance.var[m,k]/G)
      between.group.t.statistic[m,k,3] <- relationship.variance.mean[m,k]/sqrt(relationship.variance.var[m,k]/G)
      between.group.t.statistic[m,k,4] <- dyadic.covariance.mean[m,k]/sqrt(dyadic.covariance.var[m,k]/G)
      between.group.t.statistic[m,k,5] <- actorpartner.covariance.mean[k]/sqrt(actorpartner.covariance.var[k]/G)}}
  df <- G-1
  for (m in 1:(dim(between.group.t.statistic)[1])){ 
    for (k in 1:T){ 
      for (n in 1:(dim(between.group.t.statistic)[3])){
        if ( !is.na(between.group.t.statistic[m,k,n]) ){
          if ( sign(between.group.t.statistic[m,k,n]) > 0 )
            between.group.p.value[m,k,n] <- pt(between.group.t.statistic[m,k,n],df,lower.tail=FALSE)
          else between.group.p.value[m,k,n] <- pt(between.group.t.statistic[m,k,n],df)}}}}
  list(t.statistic=between.group.t.statistic,df=df,two.tailed.p.value=2*between.group.p.value)
}

# Testing dyadic covariance in block designs by means of r Pearson correlation test #

block.SRM.dyadic.reciprocity.tTest <- function(X,symmetric=FALSE){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups<- X$subgroups
  groups <- X$groups  
  times <- X$times
  dyadic.reciprocity <- block.SRM.dyadic.reciprocity(X,symmetric)
  t.value <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
    for (l in 1:G){
      t.value[k,l] <- dyadic.reciprocity[k,l]*sqrt(G1*(G2-1)-1)/
        sqrt(1-dyadic.reciprocity[k,l]**2)}}
  df <- (G1-1)*(G2-1)-1
  t.p.value <- array(dim=c(T,G),dimnames=list(times,groups),NA)
  for (k in 1:T){
    for (l in 1:G){
      if (sign(t.value[k,l])>0)
        t.p.value[k,l] <- pt(t.value[k,l],df,lower.tail=FALSE)
      else t.p.value[k,l] <- pt(t.value[k,l],df)}}
  list(t.value=t.value,two.tailed.p.value=2*t.p.value)
}

################################################################################
#     R FUNCTIONS FOR  MEASURING RECIPROCITY WITH DYADIC DISCREPANCIES         #
################################################################################

# Function to obtain DC index by means of the observed sociomatrix #

getdc <- function(X){
  N = sum(X);
  dif_n = abs(X-t(X));
  sumdif <- sum(dif_n)
  dc = sumdif/(2*N)
  return(dc)
}

# Function to obtain delta index #

getdelta <- function(X){
  psi <- getpsi(X);
  phi <- getphi(X);
  delta = phi/psi;
  return(delta)
}

# Function to compute generalized reciprocity index -epsilon- #

getepsilon <- function(X){
  nu <- getnu(X);
  nuj <- getnuj(X);
  nui <- colSums(nu);
  epsilon = 0.;
  for (i in 1:nrow(X)){
    epsilon = epsilon + abs(nuj[i] - nui[i]);}
  if (sum(nuj) != 0)
  {epsilon = 1 - (epsilon/(sum(nuj)))}
  else
  {epsilon=1.}
  return(epsilon)
}

# Function to compute dyadic reciprocity index -kappa- #

getkappa <- function(X){
  rationu <- getrationu(X);
  dif_rationu = abs(rationu - t(rationu));
  kappa =1 - ((sum(dif_rationu)/2)/nrow(X));
  return(kappa)
}

# Function to compute matrix of unweighted contributions to symmetry -lambda- #

getlambda <- function(X){
  S = (X + t(X))/2;
  trX = sum(diag(t(X)%*%X));
  Srow = S*S;
  vectss = rowSums(Srow);
  K = (X - t(X))/2;
  Krow = K*K;
  vectkk = rowSums(Krow);
  lambda = array(dim=c(nrow(X),ncol(X)));
  sumvect = vectss+vectkk;
  for(i in 1:nrow(X))
    for(j in 1:ncol(X))
    {if (sumvect[i] == 0)
    {lambda[i,j] = 0.}
      else
      {lambda[i,j] = (S[i,j]*S[i,j])/sumvect[i]}
    }
  lambda;
}

# Function to compute individual contributions to symmetry -lambdaj- #

getlambdaj <- function(X){
  lambda <- getlambda(X);
  lambdaj = rowSums(lambda);
  lambdaj;
}

# Function to compute matrix of unweighted contributions to skew-symmetry -nu- #

getnu <- function(X){
  S = (X + t(X))/2;
  trX = sum(diag(t(X)%*%X));
  Srow = S*S;
  vectss = rowSums(Srow);
  K = (X - t(X))/2;
  Krow = K*K;
  vectkk = rowSums(Krow);
  nu = array(dim=c(nrow(X),ncol(X)));
  sumvect = vectss+vectkk;
  for(i in 1:nrow(X))
    for(j in 1:ncol(X))
    {if (sumvect[i] == 0)
    {nu[i,j] = 0.}
      else
      {nu[i,j] = (K[i,j]*K[i,j])/sumvect[i]}
    }
  nu;
}

# Function to compute individual contributions to skew-symmetry -nuj- #

getnuj <- function(X){
  nu <- getnu(X);
  nuj = rowSums(nu);
  nuj;
}

# Function to compute matrix of ratios skew-symmetry/symmetry -omega- #

getomega <- function(X){
  nu <- getnu(X);
  lambda <- getlambda(X);
  omega = array(dim=c(nrow(X),ncol(X)));
  for(i in 1:nrow(X))
    for(j in 1:ncol(X))
    {if (lambda[i,j] == 0.)
    {omega[i,j]=0.}
      else
      {omega[i,j] = nu[i,j]/lambda[i,j]}
    }
  omega;
}

# Function to obtain phi index by means of the observed sociomatrix #

getphi <- function(X){
  K = (X - t(X))/2;
  trX = sum(diag(t(X)%*%X));
  trK = sum(diag(t(K)%*%K));
  phi = trK / trX
  return(phi)
}

# Function to obtain psi index by means of the observed sociomatrix #

getpsi <- function(X){
  S = (X + t(X))/2;
  trX = sum(diag(t(X)%*%X));
  trS = sum(diag(t(S)%*%S));
  psi = trS / trX
  return(psi)
}

# Function to compute dyadic decomposition of symmetry -ratiolambda- #

getratiolambda <- function(X){
  lambda <- getlambda(X);
  lambdaj <- getlambdaj(X);
  ratiolambda = array(dim=c(nrow(X),ncol(X)));
  for(i in 1:nrow(X))
    for(j in 1:ncol(X))
    {if (lambdaj[i] == 0)
    {ratiolambda[i,j] = 0.}
      else
      {ratiolambda[j,i] = lambda[i,j]/lambdaj[i]}
    }
  ratiolambda;
}

# Function to compute dyadic decomposition of skew-symmetry -rationu- #

getrationu <- function(X){
  nu <- getnu(X);
  nuj <- getnuj(X);
  rationu = array(dim=c(nrow(X),ncol(X)));
  for(i in 1:nrow(X))
    for(j in 1:ncol(X))
    {if (nuj[i] == 0)
    {rationu[i,j] = 0.}
      else
      {rationu[j,i] = nu[i,j]/nuj[i]}
    }
  rationu;
}

# Function to estimate p values for social reciprocity statistics at different levels of analysis under the specified null hypothesis #

reciptest1<-function(X,pi,rep=9999,names=NULL,label=FALSE){
  
  if (!squaremat(X))
    return("Error: Matrix X is not square and cannot be analyzed")
  
  if (!squaremat(pi))
    return("Error: Matrix Pi is not square and cannot be analyzed")
  
  # Limit the number of replications #
  
  if ((rep < 1) | (rep > 1000000))
    return("Error: Number of replications must be between 1 and 1000000");
  
  # Compute the social reciprocity statistics at different levels for the original matrix #
  
  phi <- getphi(X);
  psi <- getpsi(X);
  dc <- getdc(X);
  delta <- getdelta(X);
  epsilon <- getepsilon(X);
  kappa <- getkappa(X);
  nuj <- getnuj(X);
  lambdaj <- getlambdaj(X);
  rationu <- getrationu(X);
  ratiolambda <- getratiolambda(X);
  omega <- getomega(X);
  
  # Matrices have to be transformed into vectors in order to pass to a .C routine #
  
  vecX=array(dim=c(nrow(X)*(ncol(X)-1)/2,1)); 
  m=0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      m =m+1;
      vecX[m]=X[i,j];}
  vecpi=array(dim=c(nrow(X)*(ncol(X)-1)/2,1));
  m=0.;
  for (i in 1:nrow(pi))
    for (j in 1:ncol(pi))
    {
      m =m+1;
      vecpi[m]=pi[i,j];}
  
  # R wrapper of the .C routine #
  
  out <- .C("recip0",
            as.double(vecX),
            as.double(vecpi),
            as.integer(nrow(X)),
            as.integer(rep),
            pvphi=double(1),
            pvpsi=double(1),
            pvdc=double(1),
            pvdelta=double(1),
            pvepsilon=double(1),
            pvkappa=double(1),
            pvnuj=double(nrow(X)),
            pvlambdaj=double(nrow(X)),
            pvrationu=double(length(vecX)),
            pvratiolambda=double(length(vecX)),
            pvomega=double(length(vecX)),
            PACKAGE="DyaDA")
  
  pvphi = out$pvphi;
  pvpsi = out$pvpsi;
  pvdelta = out$pvdelta;
  pvdc = out$pvdc;
  pvepsilon = out$pvepsilon;
  pvkappa = out$pvkappa;
  pvrationu = array(dim=c(nrow(X),ncol(X)),0.);
  m = 1;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      pvrationu[i,j] = out$pvrationu[m];
      m = m+1;
    }
  pvratiolambda = array(dim=c(nrow(X),ncol(X)),0.);
  m = 1;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      pvratiolambda[i,j] = out$pvratiolambda[m];
      m = m+1;
    }
  pvomega = array(dim=c(nrow(X),ncol(X)),0.);
  m = 1;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      pvomega[i,j] = out$pvomega[m];
      m = m+1;
    }
  
  pvnuj = out$pvnuj;
  pvlambdaj = out$pvlambdaj;
  
  # Computation of some summary statistics and print results #
  
  if (label == TRUE){dimnames(rationu)<- list(c(names),c(names))}
  if (label == TRUE){dimnames(pvrationu)<- list(c(names),c(names))}
  if (label == TRUE){dimnames(ratiolambda)<- list(c(names),c(names))}
  if (label == TRUE){dimnames(pvratiolambda)<- list(c(names),c(names))}
  if (label == TRUE){dimnames(omega)<- list(c(names),c(names))}
  if (label == TRUE){dimnames(pvomega)<- list(c(names),c(names))}
  
  res <- list(label=label,names=names,phi=phi,pvphi=pvphi,psi=psi,pvpsi=pvpsi,delta=delta,pvdelta=pvdelta,
              dc=dc,pvdc=pvdc,epsilon=epsilon,pvepsilon=pvepsilon,kappa=kappa,pvkappa=pvkappa,nuj=nuj,
              pvnuj=pvnuj,lambdaj=lambdaj,pvlambdaj=pvlambdaj,ratiolambda=ratiolambda,
              pvratiolambda=pvratiolambda,rationu=rationu,pvrationu=pvrationu,omega=omega,
              pvomega=pvomega)
  class(res) <- "reciptest1"
  res
}

print.reciptest1 <- function(x,digits=max(4,getOption("digits")-4),...)
{
    cat("SOCIAL RECIPROCITY ANALYSIS AT DIFFERENT LEVELS OF ANALYSIS","\n")
    cat("===========================================================","\n")
    cat("    ","\n")    
    cat("OVERALL MEASURES OF SOCIAL RECIPROCITY","\n")
    cat("======================================","\n")
    cat("Overall skew-symmetry index -phi-","\n")
    cat("=================================","\n")
    cat(round(x$phi,digits),"\n");
    cat("p value","\n")
    cat("=======","\n")
    cat(round(x$pvphi,digits),"\n")
    cat("    ","\n")
    cat("Overall symmetry index -psi-","\n")
    cat("============================","\n")
    cat(round(x$psi,digits),"\n");
    cat("p value","\n")
    cat("=======","\n")
    cat(round(x$pvpsi,digits),"\n")
    cat("    ","\n")
    cat("Overall skew-symmetry/symmetry ratio -delta-","\n")
    cat("============================================","\n")
    cat(round(x$delta,digits),"\n");
    cat("p value","\n")
    cat("=======","\n")
    cat(round(x$pvdelta,digits),"\n")
    cat("    ","\n") 
    cat("Directional consistency index -dc-","\n")
    cat("==================================","\n")
    cat(round(x$dc,digits),"\n");
    cat("p value","\n")
    cat("=======","\n")
    cat(round(x$pvdc,digits),"\n")
    cat("    ","\n")
    cat("Generalized reciprocity index -epsilon-","\n")
    cat("=======================================","\n")
    cat(round(x$epsilon,digits),"\n");
    cat("p value","\n")
    cat("=======","\n")
    cat(round(x$pvepsilon,digits),"\n")
    cat("    ","\n")
    cat("Dyadic reciprocity index -kappa-","\n")
    cat("================================","\n")
    cat(round(x$kappa,digits),"\n");
    cat("p value","\n")
    cat("=======","\n")
    cat(round(x$pvkappa,digits),"\n")
    cat("    ","\n")
    cat("    ","\n")
    cat("INDIVIDUAL MEASURES OF SOCIAL RECIPROCITY","\n")
    cat("=========================================","\n")
    cat("    ","\n")
    nujmat = array(dim=c(nrow(X),2));
    if (x$label == TRUE){dimnames(nujmat)<- list(c(x$ames),c("nuj","p value"))}
    if (x$label == FALSE){colnames(nujmat)=c("nuj","p value")}
    nujmat[,1] = round(x$nuj,digits);
    nujmat[,2] = round(x$pvnuj,digits);
    cat("Individual contribution to skew-symmetry","\n")
    cat("========================================","\n")
    print(nujmat);
    cat("    ","\n")
    lambdajmat = array(dim=c(nrow(X),2));
    if (x$label == TRUE){dimnames(lambdajmat)<- list(c(x$names),c("lambdaj","p value"))};
    if (x$label == FALSE){colnames(lambdajmat)=c("lambdaj","p value")};
    lambdajmat[,1] = round(x$lambdaj,digits);
    lambdajmat[,2] = round(x$pvlambdaj,digits);
    cat("Individual contribution to symmetry","\n")
    cat("===================================","\n")
    print(lambdajmat);
    cat("    ","\n")
    cat("    ","\n")
    cat("DYADIC MEASURES OF SOCIAL RECIPROCITY","\n")
    cat("=====================================","\n")
    cat("    ","\n")
    cat("Dyadic contributions to skew-symmetry","\n")
    cat("=====================================","\n")
    print.table(x$rationu,digits);
    cat("p value","\n")
    cat("=======","\n")    
    print.table(x$pvrationu,digits);
    cat("    ","\n")
    cat("Dyadic contributions to symmetry","\n")
    cat("================================","\n")
    print.table(x$ratiolambda,digits);
    cat("p value","\n")
    cat("=======","\n")
    print.table(x$pvratiolambda,digits);
    cat("    ","\n")
    cat("Matrix of dyadic balanced reciprocity","\n")
    cat("=====================================","\n")
    print.table(x$omega,digits);
    cat("p value","\n")
    cat("=======","\n")
    cat("    ","\n")
    print.table(x$pvomega,digits);
    cat("    ","\n")
    cat("    ","\n")
    invisible(x)
  }
  
# Function to compute the dyadic contribution to asymmetry - PHIij - #

getPHIij <- function (X)
{
  phirmat <- getphirmat(X);
  PHIij = phirmat + t(phirmat);
  return(PHIij)
}
# Function to compute the overall asymmetry index - PHIr - #

getPHIr <- function (X)
{
  phirmat <- getphirmat(X);
  PHIr = sum(phirmat);
  return(PHIr)
}
# Function to compute standard error of the PHIr statistic - PHIrSE - #

getPHIrSE <- function (X,pi)
{
  dyadc = X + t(X);
  z = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i < j)
      { 
        if (dyadc[i,j] != 0)
          z = z + 1;
      }
  oddterm = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i < j)
      { 
        if ((dyadc[i,j]%%2) != 0)
        {oddterm = oddterm + (1/(dyadc[i,j]**2))}
        else {oddterm = oddterm}
      }
  oddterm = oddterm/2;
  min = (z/2) + oddterm;
  qqplus2s = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i < j)
      { 
        qqplus2s = qqplus2s + ((4*dyadc[i,j]*((pi[i,j] - 7*(pi[i,j]**2) + 12*(pi[i,j]**3)) - 6*(pi[i,j]**4)) +
                                  + 8*(dyadc[i,j]**2)*(-pi[i,j] + 6*(pi[i,j]**2) - 10*(pi[i,j]**3) + 5*(pi[i,j]**4)) +
                                  + 4*(dyadc[i,j]**3)*(pi[i,j] - 5*(pi[i,j]**2) + 8*(pi[i,j]**3) - 4*(pi[i,j]**4)))/(dyadc[i,j]**4));
      }
  
  PHIrSE = 2*sqrt(qqplus2s)/(2*z-2*min);
  return(PHIrSE)
}
# Function to compute mathematical expectancy of the PHIr statistic - PHIrexpec - #

getPHIrexpec <- function (X,pi)
{
  dyadc = X + t(X);
  z = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i < j)
      {
        if (dyadc[i,j] != 0)
          z = z + 1;
      }
  oddterm = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i < j)
      { 
        if ((dyadc[i,j]%%2) != 0)
        {oddterm = oddterm + (1/(dyadc[i,j]**2))}
        else {oddterm = oddterm}
      }
  oddterm = oddterm/2;
  min = (z/2) + oddterm;
  firstterm = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i < j)
      { 
        firstterm = firstterm + (1/dyadc[i,j])
      }
  secondterm = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      if (dyadc[i,j] == 0.)
      {secondterm = secondterm + 0.}
      if (dyadc[i,j] != 0.)
      {secondterm = secondterm + (((pi[i,j]**2)*(dyadc[i,j]-1))/dyadc[i,j])}
    }
  PHIrexpec = (2*(firstterm + secondterm - min))/(2*z-2*min);
  return(PHIrexpec)
}
# Function to compute the individual contribution to asymmetry as actor - phii - #

getphii <- function (X)
{
  phirmat <- getphirmat(X);
  phii = rowSums(phirmat);
  return(phii)
}
# Function to compute the individual contribution to asymmetry as partner - phij - #

getphij <- function (X)
{
  phirmat <- getphirmat(X);
  phij = colSums(phirmat);
  return(phij)
}
# Function to compute the matrix of directional dyadic asymmetry measures - phirmat - #

getphirmat <- function (X)
{
  dyadc = X + t(X);
  z = 0.;
  for (i in 1:nrow(X))
    for (j in i:ncol(X))
      if (i < j)
      {
        if (dyadc[i,j] != 0)
          z = z + 1;
      }
  max = z;
  oddterm = 0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i < j)
      { 
        if ((dyadc[i,j]%%2) != 0)
        {oddterm = oddterm + (1/(dyadc[i,j]**2))}
        else {oddterm = oddterm}
      }
  oddterm = oddterm/2;
  min = (z/2) + oddterm;
  phirmat = array(dim=c(nrow(X),ncol(X)),0.);
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
      if (i != j)
      {
        if ( ((dyadc[i,j]%%2) > 0) & (dyadc[i,j] != 0) )
        {phirmat[i,j] =(((4*(X[i,j]**2)/(dyadc[i,j]**2)-1))-(1/dyadc[i,j]**2))/(4*max-4*min)}
        if ( ((dyadc[i,j]%%2) == 0) & (dyadc[i,j] != 0) )
        {phirmat[i,j] =((4*(X[i,j]**2)/(dyadc[i,j]**2)-1))/(4*max-4*min)}
      }
  return(phirmat)
}

# Function to estimate p values for asymmetry statistics at different levels of analysis under the null hypothesis #

reciptest2<-function(X,pi,rep=9999,names=NULL,label=FALSE){
  
  if (!squaremat(X))
    return("Error: Matrix X is not square and cannot be analyzed");
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric");
  
  if (!squaremat(pi))
    return("Error: Matrix Pi is not square and cannot be analyzed");
  if ( is.na(pi) || !is.numeric(pi))
    return("Error: Matrix of probabilities must be numeric");
  
  # Limit the number of replications #
  
  if ((rep < 1) | (rep > 1000000))
    return("Error: Number of replications must be between 1 and 1000000");
  
  # Compute the asymmetry statistics at different levels for the original matrix #
  PHIr <- getPHIr(X);
  expec <- getPHIrexpec(X,pi);
  SE <- getPHIrSE(X,pi);
  phirmat <- getphirmat(X);
  PHIij <- getPHIij(X);
  phii <- getphii(X);
  phij <- getphij(X);
  
  # Matrices have to be transformed into vectors in order to pass to a .C routine #
  
  vecX=array(dim=c(nrow(X)*(ncol(X)-1)/2,1)); 
  m=0.;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      m =m+1;
      vecX[m]=X[i,j];}
  vecpi=array(dim=c(nrow(X)*(ncol(X)-1)/2,1));
  m=0.;
  for (i in 1:nrow(pi))
    for (j in 1:ncol(pi))
    {
      m =m+1;
      vecpi[m]=pi[i,j];}
  
  # R wrapper of the .C routine #
  
  out <- .C("recip",
            as.double(vecX),
            as.double(vecpi),
            as.integer(nrow(X)),
            as.integer(rep),
            pvPHIrright=double(1),
            pvPHIrleft=double(1),
            pvphirmatright=double(length(vecX)),
            pvphirmatleft=double(length(vecX)),
            pvphiiright=double(nrow(X)),
            pvphiileft=double(nrow(X)),
            pvphijright=double(nrow(X)),
            pvphijleft=double(nrow(X)),
            pvPHIijright=double(length(vecX)),
            pvPHIijleft=double(length(vecX)),
            PACKAGE="DyaDA")
  
  pvPHIrright = out$pvPHIrright;
  pvPHIrleft = out$pvPHIrleft;
  pvphirmatright = array(dim=c(nrow(X),ncol(X)),0.);
  m = 1;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      pvphirmatright[i,j] = out$pvphirmatright[m];
      m = m+1;
    }
  pvphirmatleft = array(dim=c(nrow(X),ncol(X)),0.);
  m = 1;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      pvphirmatleft[i,j] = out$pvphirmatleft[m];
      m = m+1;
    }
  pvPHIijright = array(dim=c(nrow(X),ncol(X)),0.);
  m = 1;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      pvPHIijright[i,j] = out$pvPHIijright[m];
      m = m+1;
    }
  pvPHIijleft = array(dim=c(nrow(X),ncol(X)),0.);
  m = 1;
  for (i in 1:nrow(X))
    for (j in 1:ncol(X))
    {
      pvPHIijleft[i,j] = out$pvPHIijleft[m];
      m = m+1;
    }
  pvphiiright = out$pvphiiright;
  pvphiileft = out$pvphiileft;
  pvphijright = out$pvphijright;
  pvphijleft = out$pvphijleft;
  
  res <- list(label=label,names=names,PHIr=PHIr,pvPHIrright=pvPHIrright,
              pvPHIrleft=pvPHIrleft,expec=expec,SE=SE,
              phii=phii,pvphiiright=pvphiiright,pvphiileft=pvphiileft,
              phij=phij,pvphijright=pvphijright,pvphijleft=pvphijleft,
              phirmat=phirmat,pvphirmatright=pvphirmatright,pvphirmatleft=pvphirmatleft,
              PHIij=PHIij,pvPHIijright=pvPHIijright,pvPHIijleft=pvPHIijleft)  
  class(res) <- "reciptest2"
  res
}

print.reciptest2 <- function(x,digits=max(4,getOption("digits")-4),...)
{
  # Computation of some summary statistics and print results #
  
  if (x$label == TRUE){
    dimnames(x$phirmat)<- list(c(x$names),c(x$names))
    dimnames(x$pvphirmatright)<- list(c(x$names),c(x$names))
    dimnames(x$pvphirmatleft)<- list(c(x$names),c(x$names))
    dimnames(x$PHIij)<- list(c(x$names),c(x$names))
    dimnames(x$pvPHIijright)<- list(c(x$names),c(x$names))
    dimnames(x$pvPHIijleft)<- list(c(x$names),c(x$names))}  
  
    cat("SOCIAL RECIPROCITY ANALYSIS AT DIFFERENT LEVELS OF ANALYSIS","\n")
    cat("===========================================================","\n")
    cat("    ","\n")
    cat("OVERALL MEASURES OF SOCIAL RECIPROCITY","\n")
    cat("======================================","\n")
    cat("PHIr","\n")
    cat("====","\n")
    cat(round(x$PHIr,digits),"\n");
    cat("Right p value","\n")
    cat("=============","\n")
    cat(round(x$pvPHIrright,digits),"\n")
    cat("Left p value","\n")
    cat("============","\n")
    cat(round(x$pvPHIrleft,digits),"\n")
    cat("    ","\n")
    cat("Mathematical expectancy of PHIr","\n")
    cat("===============================","\n")
    cat(round(x$expec,digits),"\n");
    cat("Standard error of PHIr","\n")
    cat("======================","\n")
    cat(round(x$SE,digits),"\n")
    cat("    ","\n")
    cat("    ","\n")
    cat("INDIVIDUAL MEASURES OF SOCIAL RECIPROCITY","\n")
    cat("=========================================","\n")
    phiimat = array(dim=c(nrow(X),3));
    if (x$label == TRUE){dimnames(phiimat)<- list(c(x$names),c("phii","Right p value","Left p value"))}
    if (x$label == FALSE){colnames(phiimat)=c("phii","Right p value","Left p value")}
    phiimat[,1] = round(x$phii,digits);
    phiimat[,2] = round(x$pvphiiright,digits);
    phiimat[,3] = round(x$pvphiileft,digits);
    cat("Individual contribution to asymmetry as actor","\n")
    cat("=============================================","\n")
    print.table(phiimat,digits);
    phijmat = array(dim=c(nrow(X),3));
    if (x$label == TRUE){dimnames(phijmat)<- list(c(x$names),c("phij","Right p value","Left p value"))};
    if (x$label == FALSE){colnames(phijmat)=c("phij","Right p value","Left p value")};
    phijmat[,1] = round(x$phij,digits);
    phijmat[,2] = round(x$pvphijright,digits);
    phijmat[,3] = round(x$pvphijleft,digits);
    cat("    ","\n")
    cat("Individual contribution to asymmetry as partner","\n")
    cat("===============================================","\n")
    print.table(phijmat);
    cat("    ","\n")
    cat("    ","\n")
    cat("DYADIC MEASURES OF SOCIAL RECIPROCITY","\n")
    cat("=====================================","\n")
    cat("Dyadic directional asymmetry","\n")
    cat("============================","\n")
    print.table(x$phirmat,digits);
    cat("    ","\n")
    cat("Right p value","\n")
    cat("=============","\n")    
    print.table(x$pvphirmatright,digits);
    cat("    ","\n")
    cat("Left p value","\n")
    cat("============","\n")
    print.table(x$pvphirmatleft,digits);
    cat("    ","\n")
    cat("Dyadic contribution to asymmetry","\n")
    cat("================================","\n")
    print.table(x$PHIij,digits)
    cat("    ","\n");
    cat("Right p value","\n")
    cat("=============","\n")
    print.table(x$pvPHIijright,digits);
    cat("    ","\n")
    cat("Left p value","\n")
    cat("============","\n")
    print.table(x$pvPHIijleft,digits);
    invisible(x)
}

