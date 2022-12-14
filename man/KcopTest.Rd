% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KcopTest.R
\name{KcopTest}
\alias{KcopTest}
\title{Nonparametric smooth test for equality of copulas}
\usage{
KcopTest(Kdata, dn = 3, paired = FALSE)
}
\arguments{
\item{Kdata}{A list of the K dataframe or matrix}

\item{dn}{Number of copulas coefficients considered}

\item{paired}{A logical indicating whether to consider the datas as paired}
}
\value{
A list with three elements: the p-value of the test, the value of the test statistic and
the selected rank of copulas  coefficients (number of terms involved in the test statistic)
}
\description{
This functions performs the nonparametric smooth test to compare simultaneously
K(K>1) copulas. See 'Details' below for further information.
}
\details{
Recall that we have K multivariate populations of arbitrary sizes, possibly paired
with unknow associated copulas \code{C_1,...,C_K} respectively. \code{KcopTest} performs the
following hypothesis H0: C_1=C_2=...=C_K against H1: C_l differs from C_m
(l different from m and l,m in 1:K). The test is based on copulas
cross-moments founded on Legendre polynomials that he called copulas coefficients.
See the paper at the following HAL weblink: https://hal.archives-ouvertes.fr/hal-03475324v2
}
\examples{
## simulation of 5 three-dimensional populations of different sizes
Packages <- c("copula","gtools","dplyr", "orthopolynom", "stats")
lapply(Packages, library, character.only = TRUE) # if necessary
set.seed(2022)
dat1<-rCopula(50, copula = gumbelCopula(param=6,dim = 2))
dat2<-rCopula(60, copula = claytonCopula(param=0.4,dim = 2))
dat3<-rCopula(55, copula = claytonCopula(param=0.4,dim = 2))
## Form a list of data
Kdata<-list(data1=dat1,data2=dat2,data3=dat3)
## Applying the test
KcopTest(Kdata = Kdata)
}
\author{
Yves Ismael Ngounou Bakam
}
