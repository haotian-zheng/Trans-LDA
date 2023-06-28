load("TransLDA_example.RData")

Niter=1
l1 = T
mse.mat0 <-
  matrix(NA, nrow = Niter, ncol = 8)
pred.err.mat0 <-
  matrix(NA, nrow = Niter, ncol = 8)
pred.err.out.mat0 <-
  matrix(NA, nrow = Niter, ncol = 8)
pred.err.out.YM <-
  matrix(NA, nrow = Niter, ncol = 8)
pred.err.out.mat0.nointcpt <-
  matrix(NA, nrow = Niter, ncol = 8)

iter=1

# One sample LDA
beta.init <- 
  las.kA.tune(X, y, A0 = NULL, n.vec = n.vec, l1=l1, lambda = NA)
mu.0.h = mu.est(X[1:n.vec[1], ], y[1:n.vec[1]], n.vec)
pred.err.out.mat0[iter, 1] = sum((t(t(xtest)-(mu.0.h$mu1.0.h + mu.0.h$mu2.0.h)/2) %*% beta.init + log(mu.0.h$pi2.0.h/mu.0.h$pi1.0.h)) * ytest < 0) / ntest
roc.init = roc(response = as.numeric(ytest<0), predictor = as.numeric(t(t(xtest)-(mu.0.h$mu1.0.h + mu.0.h$mu2.0.h)/2) %*% beta.init + log(mu.0.h$pi2.0.h/mu.0.h$pi1.0.h)) )

# Trans-LDA with adaptive steps for unknown informative studies
set.seed(1)
Itil1 = sample(1:n.vec[1],n.vec[1]/3)
Itil2 = sample(setdiff(1:n.vec[1], Itil1),n.vec[1]/3)
prop.re1.lda <- Trans.lasso(X, y, n.vec, I.til = Itil1, l1 = l1, lasso=F, lda=T, lda.exp = T)
prop.re2.lda <-
  Trans.lasso(X, y, n.vec, I.til = Itil2, l1=l1, lasso=F, lda=T, lda.exp = T)
beta.prop.lda<- (prop.re1.lda$beta.hat + prop.re2.lda$beta.hat) / 2
pred.err.out.mat0[iter, 3] = sum((t(t(xtest)-(mu.0.h$mu1.0.h + mu.0.h$mu2.0.h)/2) %*% beta.prop.lda + log(mu.0.h$pi2.0.h/mu.0.h$pi1.0.h)) * ytest < 0) / ntest
roc.prop = roc(response = as.numeric(ytest<0), predictor = as.numeric(t(t(xtest)-(mu.0.h$mu1.0.h + mu.0.h$mu2.0.h)/2) %*% beta.prop.lda + log(mu.0.h$pi2.0.h/mu.0.h$pi1.0.h)) )

# Naive Trans-LDA
M=3
beta.all <- las.kA.tune(X, y, A0 = 1:M, n.vec = n.vec, l1=l1, lambda = 1.82)
pred.err.out.mat0[iter, 4] = sum((t(t(xtest)-(mu.0.h$mu1.0.h + mu.0.h$mu2.0.h)/2) %*% beta.all + log(mu.0.h$pi2.0.h/mu.0.h$pi1.0.h)) * ytest < 0) / ntest
roc.all = roc(response = as.numeric(ytest<0), predictor = as.numeric(t(t(xtest)-(mu.0.h$mu1.0.h + mu.0.h$mu2.0.h)/2) %*% beta.all + log(mu.0.h$pi2.0.h/mu.0.h$pi1.0.h)) )

mse.mat0[iter, ]
pred.err.out.mat0[iter, ] 
