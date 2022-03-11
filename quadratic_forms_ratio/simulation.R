n = 3000
m = 4000
Z = matrix(nrow=n, ncol=m)
het = rep(0,m)
for(i in 1:m) {
	f = runif(1, min = 0.01, max = 0.99)
	het[i] = 2*f*(1-f)
	Z[,i] = rbinom(n, size=2, prob=f)
}
Z = scale(Z,scale=FALSE)
w = 1/het
G = Z %*% diag(w) %*% t(Z)/m

m = 30
Z = matrix(nrow=n, ncol=m)
for(i in 1:m) {
	f = runif(1, min = 0.01, max = 0.99)
	Z[,i] = rbinom(n, size=2, prob=f)
}
Z = scale(Z)

h2 = 0.3
lambda = h2/(1-h2)
V = G * lambda + diag(n)
r = colSums(Z * solve(V, Z))/n/(1-h2)
summary(r)
