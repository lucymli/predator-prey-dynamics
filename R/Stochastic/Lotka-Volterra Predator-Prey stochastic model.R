#===============================================================================
# Lotka predator-prey model (Gillespie, 1977; Kot, 2001)
#===============================================================================

# This version of the Lotka predator-prey model is given by
# dN/dt = r*N - a*N*P
# dP/dt = c*a*N*P - d*P
# consisting of the three reaction channels,
#      N --r--> N + N 
#	   N + P --a--> P 
#      N + P --c*a--> P + P
#	   P --d--> 0


# 1a. Define parameter values
par <- c(preyGrow = 0.5,		#intrinsic growth rate of prey
		ratePred = 0.2,			#rate of predation
		assimEff = 0.01,		#assimilation efficiency
		predDeath = 0.5			#intrinsic death rate of predator
)

# 1b. Parameter to be changed
parChange <- seq(1.0,2,0.2)		#vector containing the parameter values to be tested

# 2. Define system
state <- c(N=10, P=10)        # Initial state vector
stateChange <- matrix(c(+1, -1, -1, 0, 0, 0, +1, -1),nrow=2,byrow=T) # State-change matrix
a  <- c("preyGrow*N", "ratePred*N*P","assimEff*ratePred*N*P","predDeath*P") # Propensity vector  
tf <- 10         # Final time
numRuns <- 1000		#number of runs of the stochastic model for each parameter value

# 3. Result vectors - empty matrices to be filled by results from simulations - each row represents the results for one parameter value
alpha <- matrix(data = 0, nrow=length(parChange),ncol=1)		#Prey extinction probability matrix
beta <- matrix(data = 0, nrow=length(parChange),ncol=1)		#Mean prey extinction time
gamma <- matrix(data = 0, nrow=length(parChange),ncol=1)		#Predator extinction probability
delta <- matrix(data = 0, nrow=length(parChange),ncol=1)		#Mean predator extinction time

# 4. 'for' loop for calculating prey/predator extinction probabilities and mean prey/predator extinction times for each parameter value
for(j in 1:length(parChange)){
	PreyExtTime <- 0 	#Total time to extinction of prey
	PredExtTime <- 0	#Total time to extinction of predator
	PreyExtinct <- 0	#Number of runs in which the prey population experienced extinction
	PredExtinct <- 0	#Number of runs in which the prey population experienced extinction
	
	#parameter to be varied
	par[1] <- parChange[j]	#specifies which parameter is altered
	
	#'for' loop for running numRuns number of simulations for each parameter value
	for (i in 1:numRuns){
		
		#GillespieSSA step
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
	alpha[j] <- PreyExtProb							#records the prey extinction probability
	
	meanPreyExtTime <- PreyExtTime/PreyExtinct		#of the extinct runs, what is the mean time to prey extinction
	beta[j] <- meanPreyExtTime							#records the mean prey extinction time
	
	PredExtProb <- PredExtinct/numRuns				#probability that the predator population goes extinct in tf time
	gamma[j] <- PredExtProb							#records the predator extinction probability
	
	meanPredExtTime <- PredExtTime/PredExtinct		#of the extinct runs, what is the mean time to predator extinction
	delta[j] <- meanPredExtTime						#records the mean predator extinction time
	print(par)
	print(PreyExtProb)
	print(meanPreyExtTime)
	print(PredExtProb)
	print(meanPredExtTime)
}

# 5. Handling results
result <- matrix(data=0, nrow=length(parChange),ncol=5)		#empty matrix to store results
result[,1] <- parChange		#the first column contains the parameter values
result[,2] <- alpha[1:length(parChange)]		#the second column contains prey extinction probabilities for each parameter value
result[,3] <- beta[1:length(parChange)] 		#the third column contains mean prey extinction time for each parameter value
result[,4] <- gamma[1:length(parChange)]		#the fourth column contains predator extinction probability for each parameter value
result[,5] <- delta[1:length(parChange)]		#the fifth column contains mean predator extinction time for each parameter value
result <- as.data.frame(result)		#changes result matrix to a data frame
colnames(result) <- c("c","Prey Extinction Probability","Mean Prey Extinction Time","Predator Extinction Probability","Mean Predator Extinction Time")		#set dataframe column headings


# 6. Plots results
par(mfrow=c(1,2))
plot(result[,1],result[,2],type="n",ylim=c(0,1),xlab="Prey Growth Rate, r",ylab="Probability of Extinction")
lines(result[,1],result[,2],col="blue")
lines(result[,1],result[,4],col="red")
legend("topright",c("Prey","Predator"),col=c("blue","red"),lty=1)
plot(result[,1],result[,5],type="n",xlab="Prey Growth Rate, r",ylab="Mean Time to Extinction",ylim=c(0,6.5))
lines(result[,1],result[,3],col="blue")
lines(result[,1],result[,5],col="red")
legend("topright",c("Prey","Predator"),col=c("blue","red"),lty=1)

# 7. Exports results as csv file
write.csv(result,"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/2. Stochastic Predator Prey Models/Data2a.csv")