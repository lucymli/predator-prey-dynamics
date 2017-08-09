###Deterministic Lotka-Volterra Model with Type 2 Functional Response####
####Ensure deSolve package is installed####

#Specify values below
t <- seq(0,#initial time
		200, #end time; 200 for population sizes~10, and 2000 for population sizes~1000
		by = 0.01 #time interval for solving the ODE; 0.01 for end time~200 and 0.1 for end time~2000
		)
pop <- 10 #assuming equal population sizes of predator and prey, this script substitutes pop for population sizes

#Initial sizes of prey and predator populations
state <- c(Prey=pop, Pred=pop)

#Parameter values
par <- c(preyGrow = 0.5, 		# intrinsic growth rate of prey
		ratePred = 0.2, 		# rate of predation
		assimEff = 0.01,		# assimilation efficiency
		predDeath = 0.5, 		# intrinsic death rate of predator
		handleTime = 0.0005			# handling time
		)


#This function describes the ODEs of the system		
LVDetDD <- function(Time, State, Param){
	with(as.list(c(State, Param)), {
		dXdt <- preyGrow*Prey - ratePred*Prey*Pred/(1+ratePred*handleTime*Prey)			#Change in prey population is dependent on prey birth and prey death due to predation; prey birth rate is density dependent
		dYdt <- assimEff*ratePred*Prey*Pred/(1+ratePred*handleTime*Prey) - predDeath*Pred					#Change in predator population is dependent on predator birth (from assimilation of prey) and predator death
		
	return(list(c(dXdt,dYdt)))
	})
}


#Calls the function that stores the time, initial conditions and parameter values
LVDetDD(t,state,par)

#Solves the ODEs using deSolve Package
out <- ode(state, t, LVDetDD, par)



#Plots phase-portrait
matplot(out[ , 2],out[ , 3],type="l",xlab = "Prey", ylab="Predator",main="c = 1",lty=1,lwd=2,col="Blue")

#Plots time evolution of prey and predator
matplot(out[ , 1], out[ , 2:3], type = "l", xlab = "Time", ylab = "Population size",
main = "c = 1", lty=1,lwd = 2,col = c("Blue","Red"))
legend("topright", c("Prey", "Predator"), col = c("Blue","Red"), lty = 1)



out[20001,]
N <- out[,2]
P <- out[,3]
max(N)
min(N)
max(P)
min(P)
