library(MFPCA);library(tidyverse);library(MASS);require(zipfR);require(fOptions)
require(Matrix);library(parallel);library(foreach);library(doMC);library(rms);library(pec)
library(fOptions);library(fda);library(survAUC);library(abind);library(survival);library(flipTime)
library(Triangulation);library(BPST);library(rARPACK);library(data.table);library(Rfast);library(readr)
library(bigmemory);library(bigalgebra);library(refund)
library(mvtnorm);library(simsurv);library(gtools)


#detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(10)
source("ancillary.R")

source("source_functions.R")

nt = 49 # number of triangles
n = N = 300 # number of women
nboot = 1000 # number of bootstraps
n1 = n2 = ncr = 40 # number of pixels in each direction
npix=n1*n2
u1=seq(0,1,length.out=n1)
v1=seq(0,1,length.out=n2)
uu=rep(u1,each=n2)
vv=rep(v1,times=n1)
Z=as.matrix(cbind(uu,vv))


if(nt == 49){
  load("V1.rda");load("Tr1.rda")
  ind.inside=BPST::inVT(Brain.V1, Brain.Tr1, Z[,1], Z[,2])$ind.inside
  V.est<-Brain.V1; Tr.est<-Brain.Tr1;
}
if(nt == 80){
  load("V2.rda");load("Tr2.rda")
  ind.inside=BPST::inVT(Brain.V2, Brain.Tr2, Z[,1], Z[,2])$ind.inside
  V.est<-Brain.V2; Tr.est<-Brain.Tr2;
}
if(nt == 144){
  load("V3.rda");load("Tr3.rda")
  ind.inside=BPST::inVT(Brain.V3, Brain.Tr3, Z[,1], Z[,2])$ind.inside
  V.est<-Brain.V3; Tr.est<-Brain.Tr3;
}

d.est = 2; r = 1
Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
B <- as.matrix(Bfull.est$B) # DIM: npts * ntriangle*((d+1)(d+2)/2)
ind.inside <- Bfull.est$Ind.inside
Q2 <- Bfull.est$Q2
K <- Bfull.est$K

surv = read_rds('sim_data')

A = as.matrix(cbind(surv$A))
X = as.matrix(cbind(surv$X))
M = surv$M[,ind.inside]

fit.M.res = FDAimage.est.ho(B = B, Q2 = Q2, K = K, X = cbind(1, X), Y = M, lambda = c(0, 1, 10, 100, 1000, 1e4))

beta.res = fit.M.res$beta[[1]]


Yhat.res = cbind(1, X) %*% t(beta.res)
R.res = M - Yhat.res
eta.res = t(tcrossprod(P.band,R.res))
Geta.res = t(eta.res) %*% eta.res/n
e.res <- R.res-eta.res

sig2.res = mean(apply(e.res^2, 2, sum)/n)

I.mat = diag(1, nrow = 921, ncol = 921)

var.res = (Geta.res + I.mat*sig2.res)
Linv.res = solve(t(chol(var.res)))

e.res <- e.res %*% t(Linv.res)
sig2.res = mean(apply(e.res^2, 2, sum)/n)

# insert a personalized test statistic after this step 
M.PW.res  = M %*% t(Linv.res)
B.PW.res =  (Linv.res) %*% B

