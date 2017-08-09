###Deterministic 2-consumer Model####
####Ensure deSolve package is installed####
	
#Specify values below
t <- seq(0,#initial time
		100, #end time
		by = 0.001 #time interval for solving the ODE
		)
	
#Initial sizes of prey and predator populations
state <- c(N=100,	#initial prey population size
		   P1=5,	#initial Predator 1 population size
		   P2=5)	#initial Predator 2 population size
	
#Parameter values
par <- c(r = 1, 		# intrinsic growth rate of prey
			
		a1 = 0.05, 		# rate of predation of Predator 1
		c1 = 0.1,		# assimilation efficiency of Predator 1
		d1 = 0.1,		# death rate of Predator 1
			
		a2 = 0.2,		# rate of predation of Predator 2
		c2 = 0.1,		# assimilation efficiency of Predator 2
		d2 = 0.1, 		# death rate of Predator 2
			
		K = 100)		# carrying capacity of prey
	
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
#out2 <- ode(state[1:2], t, LVDet2, par[c(1:4,8)])
#out3 <- ode(state[c(1,3)],t,LVDet3, par[c(1,5:8)])
	
#par(mfrow=c(2,2),mar=c(4,5,4,4))

#scaling of text
scaleText <- 1.5
	
#Plots phase-portrait of prey against predator 1 and against predator 2
#matplot(out[ , 2],out[ , 3:4],type="l",xlab = "Prey", ylab="Predator",main="a",lty=1,lwd=2,col=c("Red","Purple"),cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText)
#legend("topright",c("Predator 1","Predator 2"),col=c("Red","Purple"),lty = 1,cex=scaleText)
pdf("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/9. Extras/Fig9f.pdf")	
#Plots time evolution of prey, predator 1 and predator 2
matplot(out[ , 1], out[ , 2:4],type = "l", xlab = "Time", ylab = "Population size", lty=1,lwd = 2,col = c("Blue","Red","Purple"),cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText)
legend("topright", c("Prey", "Predator 1","Predator 2"), col = c("Blue","Red","Purple"), lty = 1,cex=scaleText)
dev.off()
#Plots time evolution of prey and predator 1 in the absence of predator 2
#matplot(out2[ , 1], out2[ , 2:3],type = "l", xlab = "Time", ylab = "Population size", main = "c", lty=1,lwd = 2,col = c("Blue","Red"),cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText)
#legend("topright", c("Prey", "Predator 1"), col = c("Blue","Red"), lty = 1,cex=scaleText)

#Plots time evolution of prey and predator 2 in the absence of predator 1
#matplot(out3[ , 1], out3[ , 2:3],type = "l", xlab = "Time", ylab = "Population size", main = "d", lty=1,lwd = 2,col = c("Blue","Purple"),cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText)
#legend("topright", c("Prey","Predator 2"), col = c("Blue","Purple"), lty = 1,cex=scaleText)


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
	
#exports results into a csv file
#write.csv(out,"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/2. Underlying parameters of R star in 2-consumer deterministic model/Data2b.csv")
