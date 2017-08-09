###Deterministic Lotka-Volterra Model####
####Ensure deSolve package is installed####

#Specify values below
t <- seq(0,#initial time
		100, #end time
		by = 0.01 #time interval for solving the ODE
		)

#Initial sizes of prey and predator populations
state <- c(Prey=10, Pred=10)

#Parameter values
par <- c(preyGrow = 5, 	# intrinsic growth rate of prey
		ratePred = 0.2, 		# rate of predation
		assimEff = 0.01,		# assimilation efficiency
		predDeath = 0.5, 		# intrinsic death rate of predator
		rateImmPrey = 10,		# rate of immigration of prey
		rateImmPred = 10		# rate of immigration of predator
		)

#This function describes the ODEs of the system		
LVDet <- function(Time, State, Param){
	with(as.list(c(State, Param)), {
		dXdt <- preyGrow*Prey - ratePred*Prey*Pred + rateImmPrey			#Change in prey population is dependent on prey birth and prey death due to predation
		dYdt <- assimEff*ratePred*Prey*Pred - predDeath*Pred + rateImmPred		#Change in predator population is dependent on predator birth (from assimilation of prey) and predator death
		
	return(list(c(dXdt,dYdt)))
	})
}

#Calls the function that stores the time, initial conditions and parametervalues
LVDet(t,state,par)

#Solves the ODEs using deSolve Package
out <- ode(state, t, LVDet, par)

#par(mfrow=c(1,2))

#Plots phase-portrait
matplot(out[ , 2],out[ , 3],type="l",xlab = "Prey", ylab="Predator",main="r = 5",lty=1,lwd=2,col="Blue")

#Plots time evolution of prey and predator
matplot(out[ , 1], out[ , 2:3], type = "l", xlab = "Time", ylab = "Population size",
main = "r = 5", lty=1,lwd = 2,col = c("Blue","Red"))
legend("topright", c("Prey", "Predator"), col = c("Blue","Red"), lty = 1)

max(out[,2])
min(out[,2])
max(out[,3])
min(out[,3])
