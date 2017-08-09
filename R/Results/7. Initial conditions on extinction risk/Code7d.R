###Deterministic 2-consumer Model####
####Ensure deSolve package is installed####
	
#Specify values below
t <- seq(0,#initial time
		100, #end time
		by = 0.01 #time interval for solving the ODE
		)
	
#Initial sizes of prey and predator populations
state <- c(N=10,	#initial prey population size
		   P1=5,	#initial Predator 1 population size
		   P2=5)	#initial Predator 2 population size
	
filenames <- c("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7d1.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7d2.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7d3.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7d4.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7d5.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7d6.pdf")

csvnames <- c("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7d1.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7d2.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7d3.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7d4.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7d5.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7d6.csv")

for(i in 1:6){
#Parameter values
par <- c(r = 1, 		# intrinsic growth rate of prey
			
		a1 = 0.2, 		# rate of predation of Predator 1
		c1 = 0.1,		# assimilation efficiency of Predator 1
		d1 = 0.1,		# death rate of Predator 1
			
		a2 = 0.2,		# rate of predation of Predator 2
		c2 = 0.1,		# assimilation efficiency of Predator 2
		d2 = 0.1, 		# death rate of Predator 2
			
		K = 10)		# carrying capacity of prey


	if(i==1){
		par[4] <- 0.09
	}
	if(i==4){
		par[4] <- 0.11
	}
	if(i==2){
		par[2] <- 2/9
	}
	if(i==5){
		par[2] <- 18/99
	}
	if(i==3){
		par[3] <- 1/9
	}
	if(i==6){
		par[3] <- 9/99
	}
	
#This function describes the ODEs of the system	with 1 prey and 2 predators	
LVDet <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1 - a2*P2)			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			dP2dt <- P2*(c2*a2*N - d2)
			
			return(list(c(dNdt,dP1dt,dP2dt)))
			})
		}

#When just prey and predator 1 are present		
LVDet2 <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1)			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			
			return(list(c(dNdt,dP1dt)))
			})
		}

#When just prey and predator 2 are present		
LVDet3 <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a2*P2)			#Change in prey population is dependent on prey birth and prey death due to predation
			dP2dt <- P2*(c2*a2*N - d2)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
		
			return(list(c(dNdt,dP2dt)))
			})
		}
	
#Calls the LVDet function using the end time, initial conditions and parameter values
LVDet(t,state,par)
LVDet2(t,state[1:2],par[c(1:4,8)])
LVDet3(t,state[c(1,3)],par[c(1,5:8)])
	
#Solves the ODEs using deSolve Package
out <- ode(state, t, LVDet, par)
out2 <- ode(state[1:2], t, LVDet2, par[c(1:4,8)])
out3 <- ode(state[c(1,3)],t,LVDet3, par[c(1,5:8)])

pdf(filenames[i])	

#scaling of text
scaleText <- 1.5
	
#Plots phase-portrait of prey against predator 1 and against predator 2

#Plots time evolution of prey, predator 1 and predator 2
matplot(out[ , 1], out[ , 2:4],type = "l", xlab = "Time", ylab = "Population size", main = "b", lty=1,lwd = 2,col = c("Blue","Red","Purple"),cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText)
legend("topright", c("Prey", "Predator 1","Predator 2"), col = c("Blue","Red","Purple"), lty = 1,cex=scaleText)
dev.off()

#maximum and minimum values of prey, predator 1 and predator 2	
print("Max and min values of prey")
max(out[,2])
min(out[,2])
print("Max and min values of Predator 1")
max(out[,3])
min(out[,3])
print("Max and min values of Predator 2")
max(out[,4])
min(out[,4])

#final population sizes at the end of the run
print("final values of simluation with 2 predators and 1 prey")
out[length(t),]
print("final values of simulation with Predator 1 and prey")
out2[length(t),]
print("final values of simulation with Predator 2 and prey")
out3[length(t),]

write.csv(out,csvnames[i])
}	


###Deterministic 2-consumer Model####
####Ensure deSolve package is installed####
	
#Specify values below
t <- seq(0,#initial time
		100, #end time
		by = 0.01 #time interval for solving the ODE
		)
	
#Initial sizes of prey and predator populations
state <- c(N=10,	#initial prey population size
		   P1=5,	#initial Predator 1 population size
		   P2=5)	#initial Predator 2 population size
	
filenames <- c("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7e1.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7e2.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7e3.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7e4.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7e5.pdf","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Only top right graph/Fig7e6.pdf")

csvnames <- c("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7e1.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7e2.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7e3.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7e4.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7e5.csv","/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/Data7e6.csv")

for(i in 1:6){
#Parameter values
par <- c(r = 1, 		# intrinsic growth rate of prey
			
		a1 = 0.2, 		# rate of predation of Predator 1
		c1 = 0.1,		# assimilation efficiency of Predator 1
		d1 = 0.1,		# death rate of Predator 1
			
		a2 = 0.2,		# rate of predation of Predator 2
		c2 = 0.1,		# assimilation efficiency of Predator 2
		d2 = 0.1, 		# death rate of Predator 2
			
		K = 10)		# carrying capacity of prey


	if(i==1){
		par[4] <- 0.05
	}
	if(i==2){
		par[2] <- 0.4
	}
	if(i==3){
		par[3] <- 0.2
	}
	if(i==4){
		par[4] <- 0.15
	}
	if(i==5){
		par[2] <- 0.1+1/30
	}
	if(i==6){
		par[3] <- 2/30
	}
	
#This function describes the ODEs of the system	with 1 prey and 2 predators	
LVDet <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1 - a2*P2)			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			dP2dt <- P2*(c2*a2*N - d2)
			
			return(list(c(dNdt,dP1dt,dP2dt)))
			})
		}

#When just prey and predator 1 are present		
LVDet2 <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1)			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			
			return(list(c(dNdt,dP1dt)))
			})
		}

#When just prey and predator 2 are present		
LVDet3 <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a2*P2)			#Change in prey population is dependent on prey birth and prey death due to predation
			dP2dt <- P2*(c2*a2*N - d2)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
		
			return(list(c(dNdt,dP2dt)))
			})
		}
	
#Calls the LVDet function using the end time, initial conditions and parameter values
LVDet(t,state,par)
LVDet2(t,state[1:2],par[c(1:4,8)])
LVDet3(t,state[c(1,3)],par[c(1,5:8)])
	
#Solves the ODEs using deSolve Package
out <- ode(state, t, LVDet, par)
out2 <- ode(state[1:2], t, LVDet2, par[c(1:4,8)])
out3 <- ode(state[c(1,3)],t,LVDet3, par[c(1,5:8)])

pdf(filenames[i])	

#scaling of text
scaleText <- 1.5
	
	
#Plots time evolution of prey, predator 1 and predator 2
matplot(out[ , 1], out[ , 2:4],type = "l", xlab = "Time", ylab = "Population size", main = "b", lty=1,lwd = 2,col = c("Blue","Red","Purple"),cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText)
legend("topright", c("Prey", "Predator 1","Predator 2"), col = c("Blue","Red","Purple"), lty = 1,cex=scaleText)
dev.off()

#maximum and minimum values of prey, predator 1 and predator 2	
print("Max and min values of prey")
max(out[,2])
min(out[,2])
print("Max and min values of Predator 1")
max(out[,3])
min(out[,3])
print("Max and min values of Predator 2")
max(out[,4])
min(out[,4])

#final population sizes at the end of the run
print("final values of simluation with 2 predators and 1 prey")
out[length(t),]
print("final values of simulation with Predator 1 and prey")
out2[length(t),]
print("final values of simulation with Predator 2 and prey")
out3[length(t),]

write.csv(out,csvnames[i])
}	
#exports results into a csv file
#write.csv(out,"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/2. Underlying parameters of R star in 2-consumer deterministic model/Data2b.csv")

