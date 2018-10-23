# Spherical symmetry
n <- 10000
theta <- runif(n, 0, 1)
r <- runif(n, 0 ,1)
x <- sqrt(-2*log(r)) * cos(2*pi*theta)
y <- sqrt(-2*log(r)) * sin(2*pi*theta)

u <- runif(n,0,1)
sqsum <- sqrt(x^2 + y^2)
x1 <- sqrt(u)*x/sqsum
x2 <- sqrt(u)*y/sqsum

scatter.smooth(x1,x2)

# Accept reject
x <- runif(n, -1, 1)
y <- runif(n, -1, 1)
xs <- ifelse(x^2 + y^2 < 1, x, 0)
ys <- ifelse(x^2 + y^2 < 1, y, 0)
scatter.smooth(xs,ys)


# MCMC methods




