###Deterministic Lotka-Volterra Predator-Prey model with Density Dependence####
####Ensure deSolve package is installed####
	
#Specify values below
t <- seq(0,#initial time
		1000, #end time
		by = 0.01 #time interval for solving the ODE
		)
	
#Initial sizes of prey and predator populations
state <- c(N=10,	#initial prey population size
		   P1=5)	#initial Predator 1 population size

#Parameter values
par <- c(r = 1, 		# intrinsic growth rate of prey
			
		a1 = 0.2, 		# rate of predation of Predator 1
		c1 = 0.1,		# assimilation efficiency of Predator 1
		d1 = 0.1,		# death rate of Predator 1
			
		K = 10,		# carrying capacity of prey
		
		iN = 0.1,		#immigration rate of prey
		iP = 0.2)		#immigration rate of predator
		
parChange <- seq(0,10,0.1)

#empty results vectors
Rstarlist <- rep(0,length(parChange))
eqmP <- rep(0,length(parChange))
minN <- rep(0,length(parChange))
maxN <- rep(0,length(parChange))
minP <- rep(0,length(parChange))
maxP <- rep(0,length(parChange))
Rstarlist_i <- rep(0,length(parChange))
eqmP_i <- rep(0,length(parChange))
minN_i <- rep(0,length(parChange))
maxN_i <- rep(0,length(parChange))
minP_i <- rep(0,length(parChange))
maxP_i <- rep(0,length(parChange))
for(i in 1:length(parChange)){
	
	par[4] <- parChange[i]
	
	#When just prey and predator 1 are present		
	LVDet <- function(Time, State, Param){
			with(as.list(c(State, Param)), {
				dNdt <- N*(r*(1-N/K) - a1*P1)			#Change in prey population is dependent on prey birth and prey death due to predation
				dP1dt <- P1*(c1*a1*N - d1)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			
				return(list(c(dNdt,dP1dt)))
				})
			}
	LVDet2 <- function(Time, State, Param){
			with(as.list(c(State, Param)), {
				dNdt <- N*(r*(1-N/K) - a1*P1) + iN			#Change in prey population is dependent on prey birth and prey death due to predation
				dP1dt <- P1*(c1*a1*N - d1) + iP		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			
				return(list(c(dNdt,dP1dt)))
				})
			}

	#Calls the LVDet function using the end time, initial conditions and parameter values
	LVDet(t,state,par[c(1:5)])
	LVDet2(t,state,par)
	
	#Solves the ODEs using deSolve Package
	out <- ode(state, t, LVDet, par[c(1:5)])
	out2 <- ode(state, t, LVDet, par)
	
	#scaling of text
	scaleText <- 1.5
	
	#maximum and minimum values of prey, predator 1 and predator 2	
	maxN[i] <- max(out[,2])
	minN[i] <- min(out[,2])
	maxP[i] <- max(out[,3])
	minP[i] <- min(out[,3])
	
	maxN_i[i] <- max(out2[,2])
	minN_i[i] <- min(out2[,2])
	maxP_i[i] <- max(out2[,3])
	minP_i[i] <- min(out2[,3])
	
	#final population sizes at the end of the run
	Rstarlist[i] <- out[length(out[1,]),2]
	eqmP[i] <- out[length(out[1,]),3]
	
	Rstarlist_i[i] <- out2[length(out2[1,]),2]
	eqmP_i[i] <- out2[length(out2[1,]),3]
}	

results <- data.frame(parChange,Rstarlist,Rstarlist_i,eqmP,eqmP_i,minN,minN_i,maxN,maxN_i,minP,minP_i,maxP,maxP_i)
colnames(results) <- c('d','R*','R* with imm','eqm P','eqm P with imm','min prey','min prey with imm','max prey','max prey with imm','min pred','min pred with imm','max pred','max pred with imm')
write.csv(results,"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/8. Deterministic immigration/Data8a.csv")

pdf("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/8. Deterministic immigration/Fig8a.pdf")
#Plots phase-portrait of prey against predator 1 and against predator 2
matplot(results[,1],results[ , 2:3],type="l",xlab = "Predator Death Rate", ylab="R*",lty=1,lwd=2,col=c("Red","Purple"),cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText)
legend("topright",c("Without Immigration","With Immigration"),col=c("Red","Purple"),lty = 1,cex=scaleText)
dev.off()
#exports results into a csv file
#write.csv(out,"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/2. Underlying parameters of R star in 2-consumer deterministic model/Data2b.csv")
