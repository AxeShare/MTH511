# Estimation of pi using square
n <- 1e5
xs <- runif(n, min = -1, max = 1)
ys <- runif(n, min = -1, max = 1)
l <- ifelse(xs^2 + ys^2 < 1, 1, 0)

pi1 <- 4*mean(l)
er1 <- 100*abs(pi- pi1)/pi
r1 = 1-mean(l)

# Estimation of pi using equilateral triangle
xt <- vector()
yt <- vector()
u1 <- runif(n, min = 0, max = 1)
u2 <- runif(n, min = 0, max = 1)

xt <- sqrt(3)*u1 + 2*sqrt(3)*u2 
yt <- 3*u1
k <- ifelse((xt - sqrt(3))^2 + (yt - 1)^2 < 1 | (xt - 2*sqrt(3))^2 + (yt - 2)^2 < 1,  1, 0) 

pi2 <- 3*sqrt(3)*mean(k)
er2 <- 100*abs(pi - pi2)/pi
r2 = 1-mean(k)

###########################################################################################################################################################################################################################################

# Estimation of e using monte carlo methods

# 1. Inverse transform method
n <- 1e6
x <- runif(n, min = 0, max = 1)
y <- -log(x)

e1 <- ifelse(y>1,1,0)
exp1 <- 1/mean(e1)

# 2a. Extreme value thoery
n <- 1e5
e2 <- 0

for(j in 1:n){
  u <- runif(10, min = 0, max = 1)
  p = min(which(cumsum(u)>1))
  e2 <- e2 + p
}

exp2 <- e2/n

# 2b. Below given method is an improvement in the above method leading to variance reduction using antithetic variables
n <- 1e5
e3 <- 0
for(j in 1:n){
  u1 <- runif(10, min = 0, max = 1)
  u2 <- 1 - u1
  p <- min(which(cumsum(u1)>1))
  q<-  min(which(cumsum(u2)>1))
  e3 <- e3 + (p+q)/2
}

exp3 <- e3/n

# 3. Bootstrap method
library(boot)
n <- 1e3

dat <- c("A", rep("B", n-1))
indicator <- function(x, ndx)xor("A"%in%x[ndx], TRUE) 

p_hat <- function(dat, m=1000){
  foo <- boot(data=dat, statistic=indicator, R=m) 
  1/mean(foo$t)
} 

e4 <- replicate(100, p_hat(dat))
exp4 <- mean(e4)

# 4. Alternative method
n <- 1e5
U <- sort(runif(n), decreasing = TRUE)
S <- cumsum(1/U)/n
m <- min(which(S >= 1))
e5 = 2/(U[m-1]+U[m])
exp5 <- e5

# 5. estimating from normal distribution using CLT
e6 <- vector()
for(i in 1:20){
  n = 1e6
  h = 0.01
  uc <- runif(n)
  us <- runif(n)
  nc <- sqrt(-2*log(uc))*cos(2*pi*us) 
  ns <- sqrt(-2*log(uc))*sin(2*pi*us)
  nor <- append(nc,ns)

  cdf <- 1-length(nor[nor>sqrt(2)])/length(nor)
  cdfh <- 1-length(nor[nor>sqrt(2)+h])/length(nor)

  pdf <- (cdfh - cdf)/h
  e6[i] <- 1/(pdf*sqrt(2*pi))
}
exp6 <- mean(e6)



# Numerical approximations
n <- 1e5
exp7 <- 1/((1-1/n)^n)

t=0
k=0
while(TRUE){
  m = t
  t = t + (k+1)/factorial(k)
  k = k+1
  if(abs(t-m)<1e-3){
    exp8 = t/2
    break
  }
}

t=0
k=0
while(TRUE){
  m = t
  t = t + (3 - 4*k*k)/factorial(2*k+1)
  k = k+1
  if(abs(t-m)<1e-3){
    exp9 = t
    break
  }
}

cat(sprintf("Estimated value of pi from \n1. Square = %f\n2. Triangle = %f\nPercentage Error in value of pi from\n1. Square = %f\n2. Triangle = %f\nFraction of points rejected\n1. Square = %f\n2. Triangle = %f\n\n\nEstimated value of e from Monte Carlo Techniques\n1. Method 1: %f\n2. Method 2a: %f\n3. Method 2b: %f\n4. Method 3: %f\n5. Method 4: %f\n6. Method 5: %f\nEstimated value of e from Numerical Approximations\n1. Method 1: %f\n2. Method 2: %f\n3. Method 3: %f\n",pi1,pi2,er1,er2,r1,r2,exp1,exp2,exp3,exp4,exp5,exp6,exp7,exp8,exp9))
