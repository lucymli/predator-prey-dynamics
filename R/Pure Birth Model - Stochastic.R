timeend <- 100		#length of time that passes before simulation ends
p0 <- 10		#initial population size
b <- 1.05		#number of births expected per unit time per individual

t <- 0:timeend		#discrete time vector
p <- array(0,dim=c(timeend+1,1,1))		#population size at each time step during the simulation

p[1] <- p0


for (i in 1:timeend) {
	p1 = rpois(1,lambda = b*p0)		#general population growth model
	p0 = p1
	p[i+1,1,1] = p0
}	

mod1 <- data.frame (t,p)		#dataframe containing simulation time and the corresponding population size

plot(t,p,ylim=c(0,1000),xlim=c(0,timeend),col="Darkorchid4",main="Stochastic Pure Birth Model (b = 1.05)",ylabel="Population Size",xlabel="Time")

par(new=TRUE)



