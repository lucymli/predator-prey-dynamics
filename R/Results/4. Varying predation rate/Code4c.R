#===============================================================================
# 2 consumer - 1 resource stochastic Lotka Volterra with Density Dependence
#===============================================================================

# This version of the Lotka predator-prey model is given by
# dP1/dt = P1 * (c1 * a1 * N - d1)			First predator species
# dP2/dt = P2 * (c2 * a2 * N - d2)			Second predator species
# dN/dt = N * (r * (1 - N/K) - a1 * P1 - a2 * P2)			Prey species
# 
# consisting of the seven events,
#      N --r(1-N/K)--> N + N 			Prey birth
#	   N + P1 --a1--> P1				Prey death due to predation by Predator 1
#      N + P2 --a2--> P2 				Prey death due to predation by Predator 2
#      N + P1 --c1*a1--> P1 + P1		Prey death due to predation by Predator 1 leading to birth of a new Predator 1
#      N + P2 --c2*a2--> P2 + P2		Prey death due to predation by Predator 2 leading to birth of a new Predator 2
#	   P1 --d1--> 0						Predator 1 death
#      P2 --d2--> 0						Predator 2 death


# 1a. Define parameter values
par <- c(r = 1, 	# intrinsic growth rate of prey
			
		a1 = 0.2, 		# rate of predation of Predator 1
		c1 = 0.1,		# assimilation efficiency of Predator 1
		d1 = 0.1,		# death rate of Predator 1
		
		a2 = 0.2,		# rate of predation of Predator 2
		c2 = 0.1,		# assimilation efficiency of Predator 2
		d2 = 0.1, 		# death rate of Predator 2
			
		K = 1000		# carrying capacity of prey
		)

RstarRatio <- seq(0.1,4,0.1)				#R*1 to R*2 ratios probed
Rstar2 <- par[7]/par[5]/par[6]			#R* of predator 2 (d2/a2/c2)
Rstar1 <- RstarRatio * Rstar2			#R* of predator 1
parChange <- par[4]/Rstar1/par[3]		#Corresponding parameter values for given R*1:R*2 ratios


# 2. Define system
state <- c(N = 1000, P1 = 500, P2 = 500)        # Initial state vector
stateChange <- matrix(c(+1, -1, -1, -1, -1,  0,  0,
						0,   0,  0, +1,  0, -1,  0,
						0,   0,  0,  0, +1,  0, -1),nrow=3,byrow=T) # State-change matrix
propVec  <- c("r*N*(1-N/K)", "a1*N*P1","a2*N*P2","c1*a1*N*P1","c2*a2*N*P2","d1*P1","d2*P2") # Propensity vector  
tf <- 100         # Final time
numRuns <- 1000		#number of runs of the stochastic model for each parameter value



#Empty vectors collecting overall results for all parameter values
Pred1ExtResults <- rep(0,length(parChange))		#Results vector recording the probability that predator 1 goes extinct first
Pred2ExtResults <- rep(0,length(parChange))		#Results vector recording the probability that predator 2 goes extinct first
Pred0ExtResults <- rep(0,length(parChange))		#Results vector recording the probability that neither predators goes extinct

PreyExtProbResults <- rep(0,length(parChange))		#Results vector recording the probability that the prey goes extinct before tf
Pred1ExtProbResults <- rep(0,length(parChange))		#Results vector recording the probability that predator 1 goes extinct before tf
Pred2ExtProbResults <- rep(0,length(parChange))		#Results vector recording the probability that predator 2 goes extinct before tf

meanPreyExtTimeResults <- rep(0,length(parChange))	#Results vector recording the mean time to Prey Extinction
meanPred1ExtTimeResults <- rep(0,length(parChange)) #Results vector recording the mean time to Predator 1 Extinction
meanPred2ExtTimeResults <- rep(0,length(parChange)) #Results vector recording the mean time to Predator 2 Extinction


#'for' loop which repeats numRuns number of simulations for each parameter value
for(j in 1:length(parChange)){
	
	par[2] <- parChange[j]					#The parameter to be varied
	#Empty vectors collecting results for each parameter value
	PreyExtTime <- 0 	#Total time to extinction (sum of all extinction times in #numRuns) of Prey
	Pred1ExtTime <- 0	#Total time to extinction (sum of all extinction times in #numRuns) of Predator 1
	Pred2ExtTime <- 0	#Total time to extinction (sum of all extinction times in #numRuns) of Predator 2
	numExtinctPred1 <- 0  #Number of times the predator 1 population has gone extinct
	numExtinctPred2 <- 0  #Number of times the predator 2 population has gone extinct

	PreyExtinct <- 0	#Number of runs in which the Prey population experienced extinction
	Pred1Extinct <- 0	#Number of runs in which Predator 1 goes extinct before Predator 2
	Pred2Extinct <- 0	#Number of runs in which Predator 2 goes extinct before Predator 1
	Pred0Extinct <- 0   #Number of runs in which neither Predator 1 nor Predator 2 populations experienced extinction
	
	#'for' loop for running numRuns number of simulations
	for (i in 1:numRuns){
		
		#GillespieSSA step
		out <- ssa(state,propVec,stateChange,par,tf,method="D",simName="",verbose=FALSE)
		
		TimeSeries <- out$data[,1]
		PreyData <- out$data[,2]
		Pred1Data <- out$data[,3]
		Pred2Data <- out$data[,4]
		
		if(isTRUE(any(PreyData==0))){		#If the prey population reaches 0 at any point during the simulation
			x <- which.min(PreyData)		#what is the indices of the earliest appearance of 0 pop size
			y <- TimeSeries[x]				#at what time did the above population size appear
			PreyExtTime <- PreyExtTime + y  #what is the total extinction time accumulated in the simulations so far
			PreyExtinct <- PreyExtinct + 1  #adds 1 to the number of times the prey has gone extinct
		}
		
		if(isTRUE(any(Pred1Data==0)) && isTRUE(any(Pred2Data==0))){	#If the predator 1 population reaches 0 at any point during the simulation
			if(which.min(Pred1Data) < which.min(Pred2Data)){		#If the minimum Predator 1 population size occured earlier than minimum Predator 2 population size,
				Pred1Extinct <- Pred1Extinct + 1					#add 1 to the number of times Predator 1 has gone extinct
			}
			if(which.min(Pred2Data) < which.min(Pred1Data)){		#If the minimum Predator 2 population size occured earlier than minimum Predator 1 population size,
				Pred2Extinct <- Pred2Extinct + 1					#add 1 to the number of times Predator 2 has gone extinct
			}
		}
		
		if(isTRUE(any(Pred1Data==0)) && isTRUE(any(Pred2Data==0))==FALSE){		#If Predator 2 population still persists at tf but Predator 1 has gone extinct, 
			Pred1Extinct <- Pred1Extinct + 1									#add 1 to the number of times Predator 1 has gone extinct
		}
		
		if(isTRUE(any(Pred2Data==0)) && isTRUE(any(Pred1Data==0))==FALSE){		#If Predator 1 population still persists at tf but Predator 2 has gone extinct,
			Pred2Extinct <- Pred2Extinct + 1									#add 1 to the number of times Predator 2 has gone extinct
		}
				
		if(isTRUE(any(Pred1Data==0))==FALSE && isTRUE(any(Pred2Data==0))==FALSE){	#If neither predator populations has gone extinct,
			Pred0Extinct <- Pred0Extinct + 1										#add 1 to the number of times 'no extinctions' has occured
		}

		if(isTRUE(any(Pred1Data==0))){												#If Predator 1 population has gone extinct
			Pred1ExtTime <- Pred1ExtTime + TimeSeries[which.min(Pred1Data)]			#Add the extinction time to the total extinction time for Predator 1
			numExtinctPred1 <- numExtinctPred1 + 1									#Add 1 to the total number of times Predator 1 has gone extinct
		}
	
		if(isTRUE(any(Pred2Data==0))){												#If Predator 2 population has gone extinct
			Pred2ExtTime <- Pred2ExtTime + TimeSeries[which.min(Pred2Data)]			#Add the extinction time to the total extinction time for Predator 2
			numExtinctPred2 <- numExtinctPred2 + 1									#Add 1 to the total number of times Predator 2 has gone extinct
		}

		#print(i) #real-time output to show which run the simulation is at
	}
	currentTime <- Sys.time()
	print(currentTime)
	print(j)

	#Result analysis
	
	Pred1ExtResults[j] <- Pred1Extinct/numRuns				#probability that Predator 1 population goes extinct in tf time before Predator 2
	Pred2ExtResults[j] <- Pred2Extinct/numRuns				#probability that Predator 2 population goes extinct in tf time before Predator 1
	Pred0ExtResults[j] <- Pred0Extinct/numRuns

	PreyExtProbResults[j] <- PreyExtinct/numRuns					#probability that the prey population goes extinct in tf time
	Pred1ExtProbResults[j] <- numExtinctPred1						#probability that the predator 1 population goes extinct in tf time
	Pred2ExtProbResults[j] <- numExtinctPred2						#probability that the predator 2 population goes extinct in tf time

	meanPreyExtTimeResults[j] <- PreyExtTime/PreyExtinct			#of the extinct runs, what is the mean time to prey extinction
	meanPred1ExtTimeResults[j] <- Pred1ExtTime/numExtinctPred1		#of the extinct runs, what is the mean time to Predator 1 extinction
	meanPred2ExtTimeResults[j] <- Pred2ExtTime/numExtinctPred2		#of the extinct runs, what is the mean time to Predator 2 extinction
	
}

#Combine results vectors into results dataframe
results <- data.frame(RstarRatio,Pred1ExtResults,Pred2ExtResults,Pred0ExtResults,PreyExtProbResults,Pred1ExtProbResults,Pred2ExtProbResults,meanPreyExtTimeResults,meanPred1ExtTimeResults,meanPred2ExtTimeResults)
colnames(results) <- c('R*1:R*2 Ratio','Prob Pred1 extinct first','Prob Pred2 extinct first','Prob neither goes extinct','Prey Extinction Prob','Pred 1 Extinction Prob','Pred 2 Extinction Prob','Mean Time to Prey Extinction','Mean Time to Pred 1 Extinction','Mean Time to Pred2 Extionction')

#Exports results as csv file
write.csv(results,"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/4. Varying predation rate/Data4c.csv")

#Plots results
plot(RstarRatio,Pred1ExtResults,type="n",xlab="R*1:R*2 Ratio",ylab="Probability of Earlier Extinction",ylim=c(0,1))
lines(RstarRatio,Pred1ExtResults,col="Red")
lines(RstarRatio,Pred2ExtResults,col="Purple")
lines(RstarRatio,Pred0ExtResults,col="Green")
legend("topright",c("Pred 1","Pred 2","Neither extinct"),col=c("Red","Purple","green"),lty=1)