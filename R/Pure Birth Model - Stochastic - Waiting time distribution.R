timeend <- 200		#number of events
p0 <- 10		#initial population size
b <- 1/10		#birth rate

t <- c(0:timeend)		#time vector
p <- c(0:timeend)		#population size at each time step during the simulation

p[1] <- p0


for (i in 1:timeend) {
	R1 = rexp(1,rate=b*p0)		#random waiting time drawn from an exponential distribution with rate b*n(t)
	t[i+1] = t[i]+R1		#updates simulation time to record time of next event
	p1 = p0+1		#population size grows by one after birth event
	p0 = p1
	p[i+1] = p0
}	

mod2 <- data.frame (t,p)		#dataframe containing simulation time and the corresponding population size

plot(t,p,col="Orange",ylim=c(0,100),xlim=c(0,25),main="Stochastic Birth Model (b = 1/10)",ylabel="Population Size",xlabel="Time")

par(new=TRUE)