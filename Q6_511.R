t <- 0
avg <- 0
avf <- 0
i <- 0
deg <- vector()
deg[1] = 5
deg[2] = 10
deg[3] = 25
deg[4] = 50
rootr <- vector()
d = 1

for(x in c("discrete" , "normal" , "cauchy" , "exponential")){
  for(g in c(6, 11, 26, 51)){
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

plot(deg, rootr[1:4],"b", ylim = c(min(rootr),max(rootr)), col = "black", ylab = "Avg number of real roots", xlab = "degree of polynomial")
lines(deg, rootr[5:8],"b", col = "red")
lines(deg, rootr[9:12],"b", col = "green")
lines(deg, rootr[13:16],"b", col = "blue")
legend(5, 3.1, legend=c("Normal", "Discrete", "Cauchy", "Exponential"), col=c("red", "black", "green" , "blue"), lty = 1, cex=0.8, text.col = "black")


