t <- 0
avg <- 0
avf <- 0
i <- 0

for(g in c(6, 11, 26, 51)){
  for(x in c("discrete" , "normal" , "cauchy" , "exponential")){
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