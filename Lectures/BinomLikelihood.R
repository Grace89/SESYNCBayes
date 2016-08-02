N = 3 #number of draws
k = 2 #number of white

C <- choose(N,k)
f <- function(t) C*(t**k)*(1-t)**(N-k)
F <- expression(C*(t**k)*(1-t)**(N-k))
plot(f,0,1, col="red")

#d <- D(f,"t")
g <- function(t) eval(D(F,"t"))

plot(g, 0,1, col="orange", 
     main = "Finding MLE", xlab = "theta",  ylab="Likelihood or Likelihood Derivative")
plot(f,0,1, col="red", add=TRUE) # plot deritive
abline(h=0, col=4, lty=2)
abline(v=0.6666883,col=4,lty=2) 
text(0.87,0.13,"Likelihood ", cex=.7)
text(0.15,0.6,"Derivitive of Likelihood", cex=.7)
text(0.6,1,"Max Likelihood Estimator=0.6666883", cex=.7)

uniroot(g, c(0.0001,.9999)) # finds the a zero in the given interval
