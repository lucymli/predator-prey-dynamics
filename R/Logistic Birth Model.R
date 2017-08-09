timeend <- 8000000		#number of events
p0 <- 10		#initial population size
b <- 2/10		#birth rate per person per unit of time (independent of population density)
d <- 1/10		#death rate per person per unit of time
K <- 250			#carrying capacity

t <- c(0:timeend)		#time vector
p <- c(0:timeend)		#population size at each time step during the simulation

p[1] <- p0

for (i in 1:timeend) {
	E <- b*(1-p0/K)+d		#total rate of events
	R1 = rexp(1,rate=E*p0)		#random waiting time drawn from an exponential distribution with rate b*n(t)
	t[i+1] = t[i]+R1		#updates simulation time to record time of next event
	R2 = runif(1,0,1)
	if (R2 <= b*(1-p0/K)/E)
	p1 = p0+1		#population size grows by one after birth event
	else
	p1 = p0-1		#population size decreases by one after death event
	p0 = p1
	p[i+1] = p0
}	

mod4 <- data.frame (t,p)		#dataframe containing simulation time and the corresponding population size

plot(t,p,col="green",ylim=c(0,300),xlim=c(0,1500),main="Logistic birth model (K=250)",ylab="Population Size",xlab="Time",pch=16,cex=.01)
lines(t,p,col="green")

par(new=TRUE)

#50 blue, 100 pink, 150 red, 200 black, 250 green

legend(10,300,c("b=0.2","b=0.4","b=0.6","b=0.8","b=1.0","b=10"),lty=c(1,1),lwd=c(1,1),col=c("hotpink4","ivory4","khaki4","lightseagreen","maroon4","green"))