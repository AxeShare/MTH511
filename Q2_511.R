# WLOG,assume that one vertex is at the origin
# Let A = (x1, y1), B = (x2,y2) and C = (x3,y3)


n <- 5*1e4

#################################################################################################################
# Direct geometrical method
r1 <- runif(n, min = 0, max = 1)
r2 <- runif(n, min = 0, max = 1)

Px = (1 - sqrt(r1)) * 2 + (sqrt(r1) * (1 - r2)) * 9 
Py = (1 - sqrt(r1)) * 3 + (sqrt(r1) * (1 - r2)) * 9

scatter.smooth(Px, Py)

###################################################################################################################
# Accept reject for triangle
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

##########################################################################################################################
# Given n points (x1,y1),(x2, y2),.....,(xn,yn), we generate uniformly from the convex polygon formed by the given points

# Accept Reject 
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

######################################################################################################################
# Alias Method
# We assume given order of points form a convex polygon

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

####################################################################################################################



