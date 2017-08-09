timeend <- 50		#length of time that passes before simulation ends
p0 <- 10		#initial population size
b <- 0.9		#number of births expected per unit time per individual

t <- 0:timeend		#discrete time vector
p <- array(0,dim=c(timeend+1,1,1))		#population size at each time point during the simulation

p[1] <- p0


for (i in 1:timeend) {
	p1 = b*p0		#general population growth model
	p0 = p1
	p[i+1] = p0
}	

mod1 <- data.frame (t,p)		#dataframe containing simulation time and the corresponding population size