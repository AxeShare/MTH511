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

####################################################################################################################################

#1. Simple grid search
x <- seq(-0.1,0.1,length=100)
y <- x

z <- outer(x,y,foo)
cat(sprintf("Simple grid search gives a minimum of %f with the given sampling of suspected region", min(z)))

####################################################################################################################################

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
    cat(sprintf("Steepest descent converged in %d iterations.\nBest set of parameters - %f and %f.\nFunction Value = %f\n",i,x,y,fun(s)))
    break
  }
}


###############################################################################################################################################################

# 3. Nelder Mead
nelder <- optim(c(-1,1), fun, method = "Nelder-Mead")
p <- nelder$par
v <- nelder$value
it <- unname(nelder$counts)
cat(sprintf("Nelder Mead converged in %d iterations.\nBest set of parameters - %f and %f.\nFunction Value = %f\n",it[1],p[1],p[2],v))

###############################################################################################################################################################

# 4. BFGS Algorithm
bfgs <- optim(c(-1,1), fun, gr = NULL, method = "BFGS")
p <- bfgs$par
v <- bfgs$value
it <- unname(bfgs$counts)
cat(sprintf("BFGS converged in %d iterations.\nBest set of parameters - %f and %f.\nFunction Value = %f\n",it[1],p[1],p[2],v))

################################################################################################################################
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
################################################################################################################################################################

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
###############################################################################################################################################################

#7. Simulated Annealing
simulated_annealing <- function(fn, s0, niter = 1e5, step = 0.01){
  s_b <- s_c <- s_n <- s0
  f_b <- f_c <- f_n <- fn(s_n)
  
  for(k in 1:niter){
    Temp <- (1-step)^k
    s_n <- rnorm(2, s_c, 1)
    f_n <- fn(s_n)
    
    if(f_n < f_c | runif(1) < exp(-(f_n - f_c)/Temp)){
      s_c <- s_n
      f_c <- f_n
    }
    
    if(f_n < f_b){
      s_b <- s_n
      f_b <- f_n
    }
    
  }
  return(list(iterations = niter, best_value = f_b, best_state = s_b))
}

sol <- simulated_annealing(fun, s0 = c(2,2))
cat(sprintf("Simulated Annealing converged\nBest set of parameters - %f and %f.\nFunction Value = %f\n",sol$best_state[1],sol$best_state[2],sol$best_value))
###############################################################################################################################################################