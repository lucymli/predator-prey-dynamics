LVPPS <- function(r,a,d,iniPrey,iniPred,numevents){

t <- c(0:(numevents))		#discrete time vector
X <- c(0:(numevents))		#population sizes of prey over time
Y <- c(0:(numevents))		#population sizes of predator over time

X[1] <- iniPrey		#initial prey population size
Y[1] <- iniPred		#initial predator population size

J <- r*X[1]			#rate of prey birth
K <- a*X[1]*Y[1]		#rate of prey death
L <- d*Y[i]			#rate of predator death

for (i in 1:(numevents)) {
	TotRate = J + K + L		#total rate of events
	
	R1 <- rexp(1,rate=TotRate)	#random waiting time drawn from exponential distribution
	t[i+1] = t[i] + *R1			#update waiting time

	R2 <- runif(1,0,1)			#random number to determine which event occurs
	
		if(R2 <= J/TotRate) {
			X[i+1] <- X[i] + 1		#prey birth event
			Y[i+1] <- Y[i]
		}
		if(R2 > J/TotRate & R2 <= (J+K)/TotRate) {
			X[i+1] <- X[i] - 1 		#prey death event
			Y[i+1] <- Y[i] + 1		#predator birth event
		}
		if(R2 > (J+K)/TotRate & R2 <= (J+K+L)/TotRate) {
			X[i+1] <- X[i]
			Y[i+1] <- Y[i] - 1		#predator death event
		}
		
	J <- r*X[i+1]
	K <- a*X[i+1]*Y[i+1]
	L <- d*Y[i+1]
	
}

#===============================================================================
# Lotka predator-prey model (Gillespie, 1977; Kot, 2001)
#===============================================================================

# This version of the Lotka predator-prey model is given by
# dN/dt = r*N - a*N*P
# dP/dt = c*a*N*P - d*P
# consisting of the three reaction channels,,
#      N --r--> N + N 
#	   N + P --a--> P 
#      N + P --c*a--> P + P
#	   P --d--> 0


# Define parameter values
a <- seq(0.0,1.0,0.2)		#assessing the effect of r values between 0 and 2 on extinction risk and mean extinction time


	par <- c(preyGrow = 0.5,		#intrinsic growth rate of prey
			ratePred = a[j],			#rate of predation
			assimEff = 0.01,		#assimilation efficiency
			predDeath = 0.5			#intrinsic death rate of predator
	)

	# Define system
	state <- c(N=10, P=10)                           # Initial state vector
	stateChange <- matrix(c(+1, -1, -1, 0, 0, 0, +1, -1),nrow=2,byrow=T) # State-change matrix
	a  <- c("preyGrow*N", "ratePred*N*P","assimEff*ratePred*N*P","predDeath*P") # Propensity vector  
	tf <- 10                                             # Final time
	simName <- "Lotka Predator-Prey Model"


	# Run the simulations 
	#nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

	# Direct method 
	numRuns <- 1000		#number of runs of the stochastic model
	
for(j in 1:length(a)){
	PreyExtTime <- 0 	#Total time to extinction of prey
	PredExtTime <- 0	#Total time to extinction of predator
	PreyExtinct <- 0	#Number of runs in which the prey population experienced extinction
	PredExtinct <- 0	#Number of runs in which the prey population experienced extinction
	
	#alpha <- rep(0,numRuns)
	
	#faster <- function(o){
		#out <- ssa(state,a,stateChange,par,tf,method="D",simName="",verbose=FALSE)
		#o <- out
		#return("hello")}
			
		
	#mclapply(alpha,faster,mc.preschedule=FALSE)
	
	for (i in 1:numRuns){
		for (i in 1:(tf)) {
	TotRate = J + K + L		#total rate of events
	
	R1 <- rexp(1,rate=TotRate)	#random waiting time drawn from exponential distribution
	t[i+1] = t[i] + *R1			#update waiting time

	R2 <- runif(1,0,1)			#random number to determine which event occurs
	
		if(R2 <= J/TotRate) {
			X[i+1] <- X[i] + 1		#prey birth event
			Y[i+1] <- Y[i]
		}
		if(R2 > J/TotRate & R2 <= (J+K)/TotRate) {
			X[i+1] <- X[i] - 1 		#prey death event
			Y[i+1] <- Y[i] + 1		#predator birth event
		}
		if(R2 > (J+K)/TotRate & R2 <= (J+K+L)/TotRate) {
			X[i+1] <- X[i]
			Y[i+1] <- Y[i] - 1		#predator death event
		}
		
	J <- r*X[i+1]
	K <- a*X[i+1]*Y[i+1]
	L <- d*Y[i+1]
	
	}
		out <- ssa(state,a,stateChange,par,tf,method="D",simName="",verbose=FALSE)
		TimeSeries <- out$data[,1]
		PreyData <- out$data[,2]
		PredData <- out$data[,3]
		
		if(isTRUE(any(PreyData==0))){	#If the prey population reaches 0 at any point during the simulation
			x <- which.min(PreyData)	#what is the indices of the earliest appearance of 0 pop size
			y <- TimeSeries[x]				#at what time did the above population size appear
			PreyExtTime <- PreyExtTime + y  #what is the total extinction time accumulated in the simulations so far
			PreyExtinct <- PreyExtinct + 1  #adds 1 to the number of times the prey has gone extinct
		}
		
		if(isTRUE(any(PredData==0))){	#If the predator population reaches 0 at any point during the simulation
			p <- which.min(PredData)	#what is the indices of the earliest appearance of 0 pop size
			q <- TimeSeries[p]				#at what time did the above population size appear
			PredExtTime <- PredExtTime + q  #what is the total extinction time accumulated in the simulations so far
			PredExtinct <- PredExtinct + 1  #adds 1 to the number of times the predator has gone extinct

		}

	}
	
	
	PreyExtProb <- PreyExtinct/numRuns				#probability that the prey population goes extinct in tf time
	alpha <- PreyExtProb						#records the prey extinction probability
	meanPreyExtTime <- PreyExtTime/PreyExtinct		#of the extinct runs, what is the mean time to prey extinction
	beta <- meanPreyExtTime					#records the mean prey extinction time
	PredExtProb <- PredExtinct/numRuns				#probability that the predator population goes extinct in tf time
	gamma <- PredExtProb						#records the predator extinction probability
	meanPredExtTime <- PredExtTime/PredExtinct		#of the extinct runs, what is the mean time to predator extinction
	delta <- meanPredExtTime					#records the mean predator extinction time
	print(PreyExtProb)
	print(meanPreyExtTime)
	print(PredExtProb)
	print(meanPredExtTime)
}
ionext <- as.data.frame(c(a,alpha,beta,gamma,delta),colnames("a","Prey Extinction Prob.","Prey Mean Extinction Time","Predator Extinction Prob.","Predator Mean Extinction Time"))