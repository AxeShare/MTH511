# Question 1(a)
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
cat(sprintf("Estimated value of pi from \n1. Square = %f\n2. Triangle = %f\nPercentage Error in value of pi from\n1. Square = %f\n2. Triangle = %f\nFraction of points rejected\n1. Square = %f\n2. Triangle = %f\n",pi1,pi2,er1,er2,r1,r2))

# Question 1(b)
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

# 5. estimating from normal distribution using pdf  
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
cat(sprintf("Estimated value of e from Monte Carlo Techniques\n1. Method 1: %f\n2. Method 2a: %f\n3. Method 2b: %f\n4. Method 3: %f\n5. Method 4: %f\n6. Method 5: %f\nEstimated value of e from Numerical Approximations\n7. Method 1: %f\n8. Method 2: %f\n9. Method 3: %f\n",exp1,exp2,exp3,exp4,exp5,exp6,exp7,exp8,exp9))
###################################################################################################################################################################################################################################################################################################################################

#  Question 2(a)
# WLOG,assume that one vertex is at the origin
# Let A = (x1, y1), B = (x2,y2) and C = (x3,y3)
n <- 5*1e4

# Triangle - Direct geometrical method
r1 <- runif(n, min = 0, max = 1)
r2 <- runif(n, min = 0, max = 1)

Px = (1 - sqrt(r1)) * 2 + (sqrt(r1) * (1 - r2)) * 9 
Py = (1 - sqrt(r1)) * 3 + (sqrt(r1) * (1 - r2)) * 9

scatter.smooth(Px, Py)

# Triangle - Accept reject
x1 <- 2
y1 <- 3
x2 <- 9
y2 <- 9
Tx <- vector()
Ty <- vector()

X_2 <- max(x1, x2)
Y_2 <- max(y1, y2)
m1 <- y1/x1
m2 <- y2/x2
m3 <- (y2-y1)/(x2-x1)
u <- runif(n, min = 0, max = X_2)
v <- runif(n, min = 0, max = Y_2)
for(j in 1:n){
  if(v[j] - m1*u[j] < 0 && v[j] - m3*u[j] + m3*x1 - y1  < 0 && v[j] - m2*u[j] > 0){
    Tx[j] = u[j]
    Ty[j] = v[j]
  }
}
scatter.smooth(Tx,Ty)



# Question 2(b)
X <- vector()  
Y <- vector()
# Input X coordinates of the six sided convex polygon
X[1] <- -1
X[2] <- 4
X[3] <- 4
X[4] <- -4
X[5] <- -6
X[6] <- -5
X[7] <- X[1] 

# Input Y coordinates of the six sided convex polygon
Y[6] <- 4
Y[5] <- -1
Y[4] <- -5
Y[3] <- -3
Y[2] <- 2
Y[1] <- 6
Y[7] <- Y[1]

# Polygon - Accept Reject

Jx <- vector()
Jy <- vector()

Xmax <- max(X)
Ymax <- max(Y)
Xmin <- min(X)
Ymin <- min(Y)

u <- runif(n, min = Xmin, max = Xmax)
v <- runif(n, min = Ymin, max = Ymax)

l <- acos(((X[6]-u)*(X[1]-u) + (Y[6]-v)*(Y[1]-v))/(sqrt((X[6]-u)^2 + (Y[6]-v)^2)*sqrt((X[1]-u)^2 + (Y[1]-v)^2)))

for(i in 1:5){
  s <- acos(((X[i]-u)*(X[i+1]-u) + (Y[i]-v)*(Y[i+1]-v))/(sqrt((X[i]-u)^2 + (Y[i]-v)^2)*sqrt((X[i+1]-u)^2 + (Y[i+1]-v)^2)))
  l <- l+s
}
i <- 1
for(j in 1:n){
  if(2*pi - l[j] < 1e-4){
    Jx[i] <- u[j]
    Jy[i] <- v[j]
    i = i+1
  }
}
scatter.smooth(Jx, Jy)

# Polygon - Alias Method
Area <- vector()
Ar = 0
Area[6] <- abs(X[6]*Y[1] - X[1]*Y[6])/2

Ar = Ar + Area[6]
for(i in 1:5){
  Area[i] <- abs(X[i]*Y[i+1] - X[i+1]*Y[i])/2
  Ar <- Ar + Area[i]
}

# Probability vector
p <- Area/Ar

Alias <- array( , dim = c(0,6))
Prob <- array( , dim = c(0,6))
Small <- list()
Large <- list()

k=1
j=1

pmf <- as.vector(6*p)
for(i in 1:6){
  if(pmf[i]<1){
    Small[j] <- i
    j <- j+1
  }
  else{
    Large[k] <- i
    k <- k+1
  }
}

while(length(Small) != 0 & length(Large) != 0){
  l <- as.numeric(Small[[length(Small)]])
  Small[[length(Small)]] <- NULL
  
  g <- as.numeric(Large[[length(Large)]])
  Large[[length(Large)]] <- NULL
  
  Prob[l] <- pmf[l]
  Alias[l] <- g
  
  pmf[g] = (pmf[g]+pmf[l]) - 1
  if((pmf[g]) < 1){
    Small <- append(Small, g, after = 0)
  }
  if((pmf[g]) >= 1){
    Large <- append(Large, g, after = 0)
  }
}

while(length(Large) != 0){
  g <- as.numeric(Large[[length(Large)]])
  Large[[length(Large)]] <- NULL
  Prob[g] <- 1
}

while(length(Small) != 0){
  l <- as.numeric(Small[[length(Small)]])
  Small[[length(Small)]] <- NULL
  Prob[l] <- 1
}

Sx <- vector()
Sy <- vector()
r1 <- runif(n)
r2 <- runif(n)
y  <- runif(n)
x <-  sample(1:6, size = n , replace = TRUE)
t <-  ifelse(y < Prob[x], x, Alias[x])

Sx = (1 - sqrt(r1)) * X[t] + (sqrt(r1) * (1 - r2)) * X[t+1] 
Sy = (1 - sqrt(r1)) * Y[t] + (sqrt(r1) * (1 - r2)) * Y[t+1]

scatter.smooth(Sx,Sy)
########################################################################################################################################################################################

# Question 3 - Please run all the functions given below before executing (a),(b),(c)
Vsp <- function(n){
  return((pi)^(n/2)/gamma(n/2 + 1))
}

norm_v <- function(x){
  return(sqrt(sum(x^2)))
}

pdfn <- function(x){
  if(norm_v(x) < 1){
    return(1/Vsp(length(x)))
  }
  else{return(0)}
}

proposal_vector <- function(x){
  return(runif(1, min = x - 0.2 , max = x + 0.2))
}


# Question 3(a)
# Spherical Symmetry
symm_sph <- data.frame()
for(i in 1:1e4){
  Y <- rnorm(3)               # Alter the dimension here (in place of 3)
  u <- runif(1)
  Y <- Y/norm_v(Y)
  X_ <- Y*(u^(1/3))           # Alter the dimension here (in place of 3)
  symm_sph <- rbind(symm_sph, X_)
}


# Question 3(b)
# Accept Reject
ar_sphere <- data.frame()
effa <- vector()
dim <- vector()
l=1
for(d in c(2,5,10,25,50)){
  ind = 1
  for(i in 1:1e4){
    V <- runif(d, -1, 1)   # Alter the dimension here (in place of 3)
    if(norm_v(V) < 1){
      ar_sphere <- rbind(ar_sphere, V)
      ind = ind+1
    }
  }
  effa[l] = ind/1e4
  dim[l] = d
  l = l+1
  # Acceptance rate = 0.7849 0.1661 0.0021 0.0001 0.0001
}

# Question 3(c)
# MCMC method
# Takes a bit of time to execute
eff <- vector()
dim <- vector()
l = 1

for(d in c(2,5,10,25,50)){
  samp_sph <- data.frame()
  x <- rep(0.1, d)              # Alter the dimension here
  samp_sph = rbind(samp_sph, x)
  ind = 1
  
  for(i in 1:1e4){
    Z = x
    
    prop <- sapply(Z, proposal_vector)
    alpha = min(1, pdfn(prop)/pdfn(Z))
    
    u = runif(1)
    if(u < alpha){
      x = prop
      samp_sph = rbind(samp_sph, x)
      ind = ind+1
    }
    else{
      samp_sph = rbind(samp_sph, x)
    }
  }
  
  eff[l] = ind/1e4
  dim[l] = d
  l = l+1
  # Acceptance rate = 0.9076 0.7726 0.5548 0.1407 0.0041
}

plot(dim,effa,ylab = "Accept reject acceptance")
plot(dim,eff, ylab = "MCMC Acceptance")
#################################################################################################################################################################################################

# Question 4 - Please run all the functions given below before executing individual algorithms.
# Given function definitions
fun <- function(s)
{
  x <- s[1]
  y <- s[2]
  cosh(sin(10*x)*x)*(x*sin(20*y) + y*sin(20*x))^2 + cosh(cos(20*y)*y)*(x*cos(10*y) - y*sin(10*x))^2
}

foo <- function(x,y)
{
  cosh(sin(10*x)*x)*(x*sin(20*y) + y*sin(20*x))^2 + cosh(cos(20*y)*y)*(x*cos(10*y) - y*sin(10*x))^2
}

pdx <- function(x,y){
  2*(-10*y*cos(10*x) + cos(10*y))*cosh(y*cos(20*y))*(x*cos(10*y) - y*sin(10*x)) + 2*cosh(x*sin(10*x))*(20*y*cos(20*x) + sin(20*y))*(y*sin(20*x) + x*sin(20*y)) + (10*x*cos(10*x) + sin(10*x))*((y*sin(20*x) + x*sin(20*y))^2)*(sinh(x*sin(10*x)))
}

pdy <- function(x,y){
  2*cosh(y*cos(20*y))*(x*cos(10*y) - y*sin(10*x))*(-sin(10*x) - 10*x*sin(10*y)) + 2*cosh(x*sin(10*x))*(20*x*cos(20*y) + sin(20*x))*(y*sin(20*x) + x*sin(20*y)) + ((x*cos(10*y) - y*sin(10*x))^2)*(cos(20*y) - 20*y*sin(20*y))*sinh(y*cos(20*y))
}

x <- seq(-3,3,length=100)
y <- x
z <- outer(x,y,foo)
persp(x,y,z,theta=60,phi=30, expand=0.5,col="red",ltheta=100,xlab="x",ticktype="detailed",ylab="y",zlab="z")



#1. Simple grid search
x <- seq(-0.1,0.1,length=100)
y <- x

z <- outer(x,y,foo)
cat(sprintf("Simple grid search gives a minimum of %f with the given sampling of suspected region", min(z)))



#2. Gradient descent algorithm

x = .1
y = .2
alpha = 0.01
s <- vector()
s[1] <- x
s[2] <- y

i = 0
while(TRUE){
  q = x
  r = y
  
  x <- x - alpha*pdx(q,r)
  y <- y - alpha*pdy(q,r)
  i = i+1
  
  if(abs(x) > 10 | (abs(y)) > 10 ){
    cat(sprintf("Steepest descent did not converge"))
    break
  }
  
  if(abs(q-x)<1e-6 & abs(r-y)<1e-6){
    gsolx <- x
    gsoly <- y
    cat(sprintf("Steepest descent converged in %d iterations.\nBest set of parameters - %f and %f.\nFunction Value = %f\n",i,x,y,fun(s)))
    break
  }
}



# 3. Nelder Mead
nelder <- optim(c(-1,1), fun, method = "Nelder-Mead")
q <- nelder$par
v <- nelder$value
it <- unname(nelder$counts)
cat(sprintf("Nelder Mead converged in %d iterations.\nBest set of parameters - %f and %f.\nFunction Value = %f\n",it[1],q[1],q[2],v))


# 4. BFGS Algorithm
bfgs <- optim(c(-1,1), fun, gr = NULL, method = "BFGS")
p <- bfgs$par
v <- bfgs$value
it <- unname(bfgs$counts)
cat(sprintf("BFGS converged in %d iterations.\nBest set of parameters - %f and %f.\nFunction Value = %f\n",it[1],p[1],p[2],v))


# 5. Pure random optimization
X <- vector()
X[1] <- 50
X[2] <- 200

for(i in 1:1e5){
  Xp <- X
  mn <- rnorm(2, mean = 0, sd = 100)
  smn <- sqrt(sum(mn^2))
  mvn <- mn/smn
  if(abs(fun(X)) > abs(fun(mvn))){
    X <- mvn    
  }
}

cat(sprintf("Pure random optimization converged and found\nBest set of parameters = %f and %f.\nFunction Value = %f\n",X[1],X[2],fun(X)))


# 6. Genetic Algorithm
pop <- matrix(runif(16, min = -100, max = 100), nrow = 8, ncol = 2, byrow = TRUE)
cost <- matrix(c(0,0,0,0,0,0,0,0), nrow = 8, ncol = 1)
pair <- vector()
off <- matrix(c(0,0,0,0,0,0,0,0), nrow = 4, ncol = 2)
rank <- c(0.4,0.3,0.2,0.1)
iter = 0

while(TRUE){
  
  for(i in 1:nrow(pop)){
    cost[i] <- fun(pop[i,]) 
  }
  
  pop <- cbind(pop, cost)
  
  pop <- pop[sort.list(pop[,3]), ]
  if(mean(pop[1:4,3])<1e-8){break}
  
  pop <- pop[1:4, ]
  
  u <- runif(4)
  for(i in 1:4){
    pair[i] = min(which(cumsum(rank)>u[i]))
  }
  pop <- pop[,-3]
  s <- runif(1)
  cross <- ifelse(s<0.5, 1, 2)
  alt <- ifelse(s<0.5, 2 , 1)
  beta <- runif(2)
  
  
  off[1,1] <- pop[pair[1],cross] - beta[1]*(pop[pair[1], cross] - pop[pair[2], cross])
  off[1,2] <- pop[pair[2], alt]
  
  off[2,1] <- pop[pair[2],cross] + beta[1]*(pop[pair[1], cross] - pop[pair[2], cross])
  off[2,2] <- pop[pair[1], alt]
  
  off[3,1] <- pop[pair[3],cross] - beta[2]*(pop[pair[3], cross] - pop[pair[4], cross])
  off[3,2] <- pop[pair[4], alt]
  
  off[4,1] <- pop[pair[4],cross] + beta[2]*(pop[pair[3], cross] - pop[pair[4], cross])
  off[4,2] <- pop[pair[3], alt]
  
  pop <- rbind(pop,off)
  
  mrow <- ceiling(10*runif(3, min = 0.1, max = 0.7))
  mcol <- ceiling(10*runif(3, min = 0, max = 0.2))
  mut <- runif(3, min = -100, max = 100)
  
  for(i in 1:3){
    pop[mrow[i], mcol[i]] <- mut[i]
    iter = iter + 1
  }
}

cat(sprintf("Genetic Algorithm converged in %d iterations.\nBest set of parameters - %f and %f.\nFunction Value = %f\n",iter,pop[1,1],pop[1,2],pop[1,3]))


#7. Simulated Annealing
simulated_annealing <- function(f, s0, iter = 1e5, step = 0.001){
  sb <- sc <- sn <- s0
  fb <- fc <- fn <- f(sn)
  
  for(k in 1:iter){
    Temp <- (1-step)^k
    sn <- rnorm(2, sc, 1)
    fn <- f(sn)
    
    if(fn < fc | runif(1) < exp(-(fn - fc)/Temp)){
      sc <- sn
      fc <- fn
    }
    
    if(fn < fb){
      sb <- sn
      fb <- fn
    }
    
  }
  return(list(iterations = iter, best_value = fb, best_state = sb))
}

sol <- simulated_annealing(fun, s0 = c(2,2))
cat(sprintf("Simulated Annealing converged\nBest set of parameters - %f and %f.\nFunction Value = %f\n",sol$best_state[1],sol$best_state[2],sol$best_value))

x <- seq(-1.5,1.5,length=200)
y <- x
z <- outer(x,y,foo)

# Run all the algorithms for getting complete contour plot with estimated global minima
filled.contour(x, y, z, color = terrain.colors, plot.axes = {axis(1); axis(2); 
  points(sol$best_state[1],sol$best_state[2],pch = 3, cex = 2, col = "gray", lwd = 3); 
  points(pop[1,2],pop[2,2], pch = 3, cex = 2, col = "black", lwd = 3); 
  points(X[1],X[2],pch = 3, cex = 2, col = "red", lwd = 3); 
  points(p[1],p[2],pch = 3, cex = 2, col = "blue", lwd = 3);
  points(q[1],q[2],pch = 3, cex = 2, col = "orange", lwd = 3);
  points(gsolx,gsoly,pch = 3, cex = 2, col = "chocolate4", lwd = 3);
  legend(1,1, legend=c("SA", "GA", "Random", "BFGS", "Nelder", "Gradient"), col=c("gray", "black", "red" , "blue", "orange", "chocolate4"), lty = 1, cex=0.8, text.col = "black")
})
#####################################################################################################################################################################################################################

# Question 5
# Given data
na  <- 200
nb  <- 40
no  <- 125
nab <- 10
n   <- na + nb + no + nab

# Intial guess
pin <- 0.7
qin <- 0.1
rin <- 0.2


# EM Algorithm
i = 0
p <- pin
q <- qin
r <- rin
while(TRUE){
  pat <- p
  pbt <- q
  pot <- r
  
  Eaa <- (na*(p)^2)/((p)^2 + 2*p*r)
  Eao <- (na - Eaa)
  Ebb <- (nb*(q)^2)/((q)^2 + 2*q*r)
  Ebo <- (nb - Ebb)
  
  p <- (2*Eaa + Eao + nab)/(2*n)
  q <- (2*Ebb + Ebo + nab)/(2*n)
  r <- (Eao + Ebo + 2*no)/(2*n)
  
  print(c(pat, pbt, pot))
  
  i <- i+1
  
  if(abs(p - pat)< 1e-4 & abs(q - pbt)< 1e-4 & abs(r - pot)< 1e-3 ){
    print(c(p, q, r))
    p <- pin
    q <- qin
    r <- rin
    break
  }
}
sprintf("EM Converged in %d iterations", i)

# Fisher Scoring

# Check if Information matrix becomes singular
f <- function(m) class(try(solve(m),silent=T))=="matrix"

j = 0
while(TRUE){
  
  Qr = c(p,q, 1-p-q)
  print(Qr)
  
  V = matrix(
    c(-2*no/(1-p-q) + na*(1/p - 1/(2-p-2*q)) - 2*nb/(2-q-2*p) + nab/p, -2*no/(1-p-q) + nb*(1/q - 1/(2-q-2*p)) - 2*na/(2-p-2*q) + nab/q),
    nrow = 2,
    ncol = 1)
  
  I = matrix(
    c(-2*no/(1-p-q)^2 - na*(1/p^2 + 1/(2-p-2*q)^2) - 4*nb/(2-q-2*p)^2 - nab/p^2, -2*no/(1-p-q)^2 - 2*na/(2-p-2*q)^2 - 2*nb/(2-q-2*p)^2, -2*no/(1-p-q)^2 - 2*na/(2-p-2*q)^2 - 2*nb/(2-q-2*p)^2, -2*no/(1-p-q)^2 - nb*(1/q^2 + 1/(2-q-2*p)^2) - 4*na/(2-p-2*q)^2 - nab/q^2), 
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
  
  up = solve(I)%*%V
  p <- p - up[1]
  q <- q - up[2]
  
  j <- j+1
  Qt <- c(p, q, 1-p-q)
  
  if(sum(abs(Qt - Qr))  < 1e-3){
    print(Qt)
    break
  }
  
}
sprintf("Fisher Converged in %d iterations, EM converged in %d iterations", j,i)
#############################################################################################################################################################################################

# Question 6
t <- 0
avg <- 0
avf <- 0
i <- 0
deg <- vector()
deg[1] = 5
deg[2] = 10
deg[3] = 25
deg[4] = 50
deg[5] = 100
rootr <- vector()
d = 1

for(x in c("discrete" , "normal" , "cauchy" , "exponential")){
  for(g in c(6, 11, 26, 51, 101)){
    while(TRUE){
      if(x == "discrete"){
        u <- runif(g)
        pol <- ifelse(u<0.5, 1, -1)
      }
      if(x == "normal"){
        u <- runif(g)
        v <- runif(g)
        pol <- sqrt(-2*log(u))*cos(2*pi*v)
      }
      if(x == "cauchy"){
        u <- runif(g)
        v <- runif(g)
        pol <- u/v
      }
      if(x == "exponential"){
        u <- runif(g)
        pol <- -log(u)
      }
      
      
      k <- sum(abs(Im(polyroot(pol))) < 1e-6, na.rm=TRUE)
      
      if(i>2){
        avf <- (t-k)/(i-1)
      }
      i <- i+1
      t <- t + k
      avg <- t/i
      
      if(abs(avg - avf) < 1e-4 & avg != avf){
        rootr[d] <- avg
        d = d+1
        cat(sprintf("Degree of polynomial: %d,  Distribution is %s,   Average number of roots: %f,     Niterations: %d\n", g-1, x, avg, i))
        t <- 0
        avg <- 0
        avf <- 0
        i <- 0
        break
      }
    }
  }
}

plot(deg, rootr[1:5],"b", ylim = c(min(rootr),max(rootr)), col = "black", ylab = "Avg number of real roots", xlab = "degree of polynomial")
lines(deg, rootr[6:10],"b", col = "red")
lines(deg, rootr[11:15],"b", col = "green")
lines(deg, rootr[16:20],"b", col = "blue")
legend(5, 3.1, legend=c("Normal", "Discrete", "Cauchy", "Exponential"), col=c("red", "black", "green" , "blue"), lty = 1, cex=0.8, text.col = "black")
########################################################################################################################################################################################################
