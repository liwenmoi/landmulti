rm(list=ls())

library(emdbook)
library(NMOF)
library(survival)
library(prodlim)
library(Matrix)
library(expm)
library(mvtnorm)
library(pec)
library(pracma)
library(landpred)
# for parallel computing
library(doParallel)
library(snow)
# please use the directory of HelperFunctions.R
source("HelperFunctions.R")

# read in the example data set
mydata <- read.csv("simulated.data1.csv",  header = TRUE)[,-1]


grid1 = seq(0.01, 5, length.out=5)
grid2 = seq(0.01, 5, length.out=5)
grid3 = list(seq(0.01, 5, length.out=5),
             seq(0.01, 5, length.out=5))


# 1. the following computes the coefficient estimates for group 1,
# for group 2 at s1 = 3, for group 3 at s2 = 3, 
# and for group 4 at (s1, s2) = (3, 3).
# 2. use SE = T to compute standard error estimates. 
#    Caution: takes longer than 30min.

set.seed(1234)
res <- multipred(mydata, formula = time ~ age + s1(st1) + s2(st2) + delta(outcome),
                  t0 = 5, L = 20, SE = F, gs.method = "loop", gs.cl = 4, SE.gs = T, B = 200,
                  s1_beta1 = c(3), grid1 = seq(0.01, 5, length.out=20),
                  s2_beta2=c(3), grid2 = seq(0.01, 5, length.out=20),
                  s1s2_beta3 = c(3,3), grid3=list(seq(0.01, 5, length.out=20),
                                                  seq(0.01, 5, length.out=20)))
res
# $group1
# $group1$group1_coef
# (Intercept)         age 
# -0.4951346   0.7225254 
# 
# $group1$group1_convergence
# [1] TRUE
# 
# 
# $group2_coef
# Intercept    age
# [1,] 0.4258612 0.7561
# 
# $group2_convergence
# [,1]
# [1,] TRUE
# 
# $group3_coef
# Intercept       age
# [1,] -0.08316788 0.6131926
# 
# $group3_convergence
# [,1]
# [1,] TRUE
# 
# $group4_coef
# Intercept       age
# [1,] 0.8113314 0.5836653
# 
# $group4_convergence
# [,1]
# res4ic TRUE


