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

################################################################################################################
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

################################################################################################################
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
sprintf("Fisher Converged in %d iterations", j)

