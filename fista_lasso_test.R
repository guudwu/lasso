# Testing program for "fista_lasso.R"

DIMENSION <- 10

bb <- seq ( 1 , DIMENSION )

DD <- rep(100,DIMENSION)

AA <- diag(bb)

mu <- svd(AA)$d[1]**2

xx <- numeric(DIMENSION)

tol <- 1e-3

source('fista_lasso.R')
res <- fista_lasso(bb,AA[,-c(1,10)],w=DD[-c(1,10)])
