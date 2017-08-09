timeend <- 2000		#number of events
p0 <- 10		#initial population size
b <- 2/100		#birth rate per person per unit of time
d <- 1/100		#death rate per person per unit of time
E <- b+d		#total rate of events

t <- c(0:timeend)		#time vector
p <- c(0:timeend)		#population size at each time step during the simulation

p[1] <- p0


for (i in 1:timeend) {
	R1 = rexp(1,rate=E*p0)		#random waiting time drawn from an exponential distribution with rate b*n(t)
	t[i+1] = t[i]+R1		#updates simulation time to record time of next event
	R2 = runif(1,0,1)
	if (R2 <= b/E)
	p1 = p0+1		#population size grows by one after birth event
	else
	p1 = p0-1		#population size decreases by one after death event
	p0 = p1
	p[i+1] = p0
}	

mod3 <- data.frame (t,p)		#dataframe containing simulation time and the corresponding population size

plot(t,p,col="Green",ylim=c(0,40),xlim=c(0,100),main="Stochastic Birth-Death Model (net growth rate = 1/10)",ylab="Population Size",xlab="Time")

par(new=TRUE)