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
#r <- seq(0.0,2.0,0.2)		#assessing the effect of r values between 0 and 2 on extinction risk and mean extinction time
#ronext <- data.frame(matrix(ncol=5,nrow=length(r)))
#names <- c("Run No.","Prey Extinction Prob.","Prey Mean Extinction Time","Predator Extinction Prob.","Predator Mean Extinction Time")
#colnames(ronext) <- names
#ronext[,1] <- r
#for(j in 1:length(r)){
	par <- c(preyGrow = 0.5,		#intrinsic growth rate of prey
			ratePred = 0.2,			#rate of predation
			assimEff = 0.01,		#assimilation efficiency
			predDeath = 0.5,			#intrinsic death rate of predator
			immPrey = 0.2,			#immigration constant for prey
			immPred = 0.2			#immigration constant for pred
	)

	# Define system
	state <- c(N=10, P=10)                           # Initial state vector
	stateChange <- matrix(c(+1, -1, +1, -1, 0, 0, 0, 0, 0, +1, -1, +1),nrow=2,byrow=T) # State-change matrix
	a  <- c("preyGrow*N", "ratePred*N*P","immPrey","assimEff*ratePred*N*P","predDeath*P","immPred") # Propensity vector  
	tf <- 10                                             # Final time
	simName <- "Lotka Predator-Prey Model"


	# Run the simulations 
	#nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

	# Direct method 
	numRuns <- 1000		#number of runs of the stochastic model
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
	#ronext[j,2] <- PreyExtProb						#records the prey extinction probability
	meanPreyExtTime <- PreyExtTime/PreyExtinct		#of the extinct runs, what is the mean time to prey extinction
	#ronext[j,3] <- meanPreyExtTime					#records the mean prey extinction time
	PredExtProb <- PredExtinct/numRuns				#probability that the predator population goes extinct in tf time
	#ronext[j,4] <- PredExtProb						#records the predator extinction probability
	meanPredExtTime <- PredExtTime/PredExtinct		#of the extinct runs, what is the mean time to predator extinction
	#ronext[j,5] <- meanPredExtTime					#records the mean predator extinction time
	print(PreyExtProb)
	print(meanPreyExtTime)
	print(PredExtProb)
	print(meanPredExtTime)
#}