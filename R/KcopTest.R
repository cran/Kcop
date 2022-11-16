#' Nonparametric smooth test for equality of copulas
#'
#' This functions performs the nonparametric smooth test to compare simultaneously
#' K(K>1) copulas. See 'Details' below for further information.
#'
#' @param Kdata  A list of the K dataframe or matrix
#' @param dn Number of copulas coefficients considered
#' @param paired A logical indicating whether to consider the datas as paired
#'
#' @details
#' Recall that we have K multivariate populations of arbitrary sizes, possibly paired
#' with unknow associated copulas \code{C_1,...,C_K} respectively. \code{KcopTest} performs the
#' following hypothesis H0: C_1=C_2=...=C_K against H1: C_l differs from C_m
#' (l different from m and l,m in 1:K). The test is based on copulas
#' cross-moments founded on Legendre polynomials that he called copulas coefficients.
#' See the paper at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03475324v2
#'
#' @return A list with three elements: the p-value of the test, the value of the test statistic and
#' the selected rank of copulas  coefficients (number of terms involved in the test statistic)
#'
#' @author Yves Ismael Ngounou Bakam
#'
#' @export
#'
#' @examples
#' ## simulation of 5 three-dimensional populations of different sizes
#' Packages <- c("copula","gtools","dplyr", "orthopolynom", "stats")
#' lapply(Packages, library, character.only = TRUE) # if necessary
#' set.seed(2022)
#' dat1<-rCopula(50, copula = gumbelCopula(param=6,dim = 2))
#' dat2<-rCopula(60, copula = claytonCopula(param=0.4,dim = 2))
#' dat3<-rCopula(55, copula = claytonCopula(param=0.4,dim = 2))
#' ## Form a list of data
#' Kdata<-list(data1=dat1,data2=dat2,data3=dat3)
#' ## Applying the test
#' KcopTest(Kdata = Kdata)

KcopTest<-function(Kdata,dn=3, paired=FALSE){
  # shifted Legendre polynomials function
  pol<- function(m,x){
    legend=orthopolynom::slegendre.polynomials(m,normalized = TRUE)[m+1]
    unlist(orthopolynom::polynomial.values(legend,x))
  }
  #
  card_S<-function(p,d){
    choose(d+p-1,d)-p
  }
  card_S<-Vectorize(card_S, vectorize.args = "d")
  # function to get the set S(d)
  Spd<-function(p,d){
    # p dimension
    # d total order
    if (p<=1 ||d<=1 || p%%1!=0||d%%1!=0)
      stop(" p and d must be positive integers greater than or equal to 2")
    gril<-as.data.frame(unlist(gtools::permutations(d+1, p, v =0:d,
                                                    repeats.allowed = TRUE)))
    maxi <- som <- NULL
    gril$maxi<-apply(gril,1,max)
    gril$som<-apply(gril[-length(gril)],1,sum)
    gril<-gril[gril$maxi!=gril$som,]
    gril<-subset(gril, select=-c(maxi,som))
    # We restrict now ourselves to the lists whose total order is d
    A<-gril
    A$som<-apply(A,1,sum)
    A<-subset(A, som==d, select=-c(som))
    A<-A[nrow(A):1,]
    rownames(A) <- 1:nrow(A)
    return(A)
  }
  #
  Sdk<-function(p,d,k){
    if (k%%1!=0||k>card_S(p,d))
      stop( "an integer k must be greater
          than 0 and least than or equal to ", card_S(p,d))
    A<-Spd(p,d)
    A<-A[1:k,]
    return(A)
  }
  # discrepency function at j vector.
  rjcarre<-function(j,datX,datY){
    datX <- copula::pobs(datX)
    datY <- copula::pobs(datY)
    n1=nrow(datX)
    n2=nrow(datY)
    p1=ncol(datX) # dimension
    p2=ncol(datY)
    if(p1!=p2) stop("the samples or our data must have
                  the same number of dimension")
    A=matrix(rep(0,n1*p1),ncol = p1)
    B=matrix(rep(0,n2*p1),ncol = p1)
    rm(p2) # release of memories
    for (i in 1:p1) {
      A[,i]<-pol(j[i],datX[,i])
      B[,i]<-pol(j[i],datY[,i])
    }

    A<-sum(apply(A, 1, prod))/n1
    B<-sum(apply(B, 1, prod))/n2
    rm(p1);rm(n1);rm(n2) # release of memories
    return((A-B)^2)
  }
  #
  L1 <-function(x){
    sqrt(3)*(2*x-1)
  }
  L1 <-Vectorize(L1)
  #
  M_is<-function(i,dat){
    dat <- copula::pobs(dat)
    L1(dat[i,1])*L1(dat[i,2])+2*sqrt(3)*
      mean((as.numeric(dat[i,1]<=dat[,1])-dat[,1])*L1(dat[,2]))+2*
      sqrt(3)*mean((as.numeric(dat[i,2]<=dat[,2])-dat[,2])*L1(dat[,1]))
  }
  M_is <- Vectorize(M_is,vectorize.args = "i")
  # variance function of the test
  VarianceTest<-function(datX,datY,paired=paired){
    n1=nrow(datX)
    n2=nrow(datY)
    if (n1!=n2 & paired==TRUE) stop("paired samples must have the same simple size")
    a<-n1/(n1+n2)
    if(paired==0){
      (1-a)*stats::var(M_is(1:n1,datX))*(n1-1)/n1+
        a*stats::var(M_is(1:n2,datY))*(n2-1)/n2
    }
    else{
      stats::var(M_is(1:n1,datX)-M_is(1:n2,datY))*(n1-1)/n1
    }
  }
  #
  T2k<-function(k,datX,datY, paired=paired){
    p=ncol(datX)
    n1=nrow(datX)
    n2=nrow(datY)
    if (paired==0){
      n1*n2*sum(apply(Sdk(p,2,k),1,rjcarre, datX=datX, datY=datY))/(n1+n2)
    }
    else {
      n1*sum(apply(Sdk(p,2,k),1,rjcarre, datX=datX, datY=datY))
    }
  }
  #
  Tdk <- function(d,k,datX,datY,paired=paired){
    p=ncol(datX)
    n1=nrow(datX)
    n2=nrow(datY)
    cd=card_S(p,d-1)
    if (d == 2){
      T2k(k,datX,datY, paired = paired)
    }else if (paired==0){
      Tdk(d-1,cd,datX,datY,paired = paired) + n1*n2*sum(apply(Sdk(p,d,k),1,rjcarre,
                                                              datX=datX, datY=datY))/(n1+n2)
    }else {
      Tdk(d-1,cd,datX,datY,paired = paired) + n1*sum(apply(Sdk(p,d,k),1,
                                                           rjcarre, datX=datX,
                                                           datY=datY))
    }
  }
  # select penalized function
  V_Dn<-function(datX,datY,dn, paired=paired){
    p=ncol(datX)
    n1=nrow(datX)
    n2=nrow(datY)
    K=1:card_S(p,dn)
    Tdk<-Vectorize(Tdk,vectorize.args = "k")
    Tdk(dn,K,datX,datY,paired=paired)
  }

  V_Dn<-Vectorize(V_Dn,vectorize.args = "dn")

  # penalized Test function

  T_penal_Dn<-function(datX,datY, dn,paired=paired){
    px<-ncol(datX)
    py<-ncol(datY)
    nx=nrow(datX)
    ny=nrow(datY)
    if (px!=py) stop("data X and data Y must have the same dimensions")
    if (px<=1 ||dn<=0 || px%%1!=0||dn%%1!=0)
      stop(" p and d must be the positives integer greater than
         or equal to 2")
    if (nx!=ny & paired==TRUE) stop("paired samples must have the same size")
    if(dn==0) stop("dn must be postive integer greater than 0")
    a=nx/(nx+ny)
    coef_penal=log(2*ny*a)  #in paired case, we get n
    long_v=sum(card_S(px,2:(dn+1))) # lenght V
    V<-V_Dn(datX,datY,2:(dn+1),paired = paired)
    V<-unlist(V)
    T_penal<-V-(1:long_v)*coef_penal
    ordre_opt<-which.max(T_penal)
    T_final<-V[ordre_opt]
    return(c(ordre_opt,T_final))
  }
  # Test function in two dimension
  Test2D_Final<-function(datX,datY,dn,paired=paired){
    T_final <- T_penal_Dn(datX,datY,dn,paired)[2]
    ordre_opt <-T_penal_Dn(datX,datY,dn,paired)[1]
    # normalization
    T_final<-T_final/VarianceTest(datX,datY,paired = paired)
    pvalue<-1-stats::pchisq(T_final,1)
    return(list(p_value=pvalue, Stat.Test=T_final,rank_opt=ordre_opt))
  }
  # K sample
  V_k <-function(Kdat,dn,paired=paired){
    # Kdat must be a list of all the K data
    K=length(Kdat)
    A=matrix(rep(NA,K^2),ncol = K)
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        A[i,j]<-T_penal_Dn(Kdat[[i]],Kdat[[j]],dn,paired = paired)[2]
        # [2]  we extract the test statistic
        # [1] we extract the optimal oder
      }
    }
    # Extraction of the index of the smallest value of a matrix
    min_row<-function(x){as.vector(arrayInd(which.min(x), dim(x)))[1]}
    min_col<-function(x){as.vector(arrayInd(which.min(x), dim(x)))[2]}
    xmin<-min_row(A)
    ymin<-min_col(A)
    A<-as.vector(t(A))
    A<-A[which(A!="NA")]
    A<-cumsum(sort(A))
    return(list(MatA=A,xmin=xmin,ymin=ymin))

  }

  K=length(Kdata)
  pn<-as.numeric(unlist(sapply(Kdata, nrow)))
  Kcoef_penal= log(K^(K-1)*prod(pn)/(sum(pn))^(K-1))
  if(K==2){return(Test2D_Final(Kdata[[1]],Kdata[[2]],dn,paired=paired))}
  else{
    vk<-K*(K-1)/2
    Va<-V_k(Kdata,dn,paired = paired)
    V<-Va$MatA
    xind<-Va$xmin
    yind<-Va$ymin
    T_Kpenal <- V-(1:vk)*Kcoef_penal
    Kchoix_opt<-which.max(T_Kpenal)
    T_Kfinal<-V[Kchoix_opt]/VarianceTest(Kdata[[xind]],Kdata[[yind]],paired = paired)
    # normalization of the test by the variance
    pvalue<-1-stats::pchisq(T_Kfinal,1)
    return(list(p_value=pvalue, Stat.Test=T_Kfinal, rank_opt=Kchoix_opt))
  }
}
