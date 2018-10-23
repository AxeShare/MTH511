# Function plot
x <- seq(-2,2,length=90)
y <- x
foo <- function(x,y)
{
  cosh(sin(10*x)*x)*(x*sin(20*y) + y*sin(20*x))^2 + cosh(cos(20*y)*y)*(x*cos(10*y) - y*sin(10*x))^2
}

z <- outer(x,y,foo)
persp(x,y,z,theta=90,phi=30, expand=0.5,col="red",ltheta=100,xlab="x",ticktype="detailed",ylab="y",zlab="z")