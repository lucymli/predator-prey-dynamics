###Deterministic 2-consumer Model####
####Ensure deSolve package is installed####
	
#Specify values below
t <- seq(0,#initial time
		1000, #end time
		by = 0.01 #time interval for solving the ODE
		)
	
#Initial sizes of prey and predator populations
state <- c(N=10,	#initial prey population size
		   P1=5,	#initial Predator 1 population size
		   P2=5)	#initial Predator 2 population size
	
#Parameter values
par <- c(r = 1, 		# intrinsic growth rate of prey
			
		a1 = 0.2, 		# rate of predation of Predator 1
		c1 = 0.1,		# assimilation efficiency of Predator 1
		d1 = 0.1,		# death rate of Predator 1
			
		a2 = 0.2,		# rate of predation of Predator 2
		c2 = 0.1,		# assimilation efficiency of Predator 2
		d2 = 0.1, 		# death rate of Predator 2
			
		K = 10,		# carrying capacity of prey
		
		i_N = 0.1,
		i_1 = 0.1,
		i_2 = 0.1
		)

RstarRatio <- seq(0,4,0.1)
Rstar2 <- par[7]/par[5]/par[6]
Rstar1 <- RstarRatio*Rstar2
parChange <- Rstar1*par[2]*par[3]
	
#This function describes the ODEs of the system	with 1 prey and 2 predators	
LVDet <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1 - a2*P2)			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1)		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			dP2dt <- P2*(c2*a2*N - d2)
			
			return(list(c(dNdt,dP1dt,dP2dt)))
			})
	}
	
LVDet2 <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1 - a2*P2) + i_N			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1) + i_1		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			dP2dt <- P2*(c2*a2*N - d2) + i_2
			
			return(list(c(dNdt,dP1dt,dP2dt)))
			})
	}
	
LVDet3 <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1 - a2*P2) + i_N			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1) + i_1		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			dP2dt <- P2*(c2*a2*N - d2) + i_2
			
			return(list(c(dNdt,dP1dt,dP2dt)))
			})
	}

LVDet4 <- function(Time, State, Param){
		with(as.list(c(State, Param)), {
			dNdt <- N*(r*(1-N/K) - a1*P1 - a2*P2) + i_N			#Change in prey population is dependent on prey birth and prey death due to predation
			dP1dt <- P1*(c1*a1*N - d1) + i_1		#Change in Predator 1 population is dependent on predator birth (from assimilation of prey) and predator death
			dP2dt <- P2*(c2*a2*N - d2) + i_2
			
			return(list(c(dNdt,dP1dt,dP2dt)))
			})
	}
	
	
#Calls the LVDet function using the end time, initial conditions and parameter values
LVDet(t,state,par[1:8])
LVDet2(t,state,par)
LVDet3(t,state,par)
LVDet4(t,state,par)


LVDetPred1Ext <- rep(0,length(parChange))
LVDetPred2Ext <- rep(0,length(parChange))
LVDetPred0Ext <- rep(0,length(parChange))

LVDet2Pred1Ext <- rep(0,length(parChange))
LVDet2Pred2Ext <- rep(0,length(parChange))
LVDet2Pred0Ext <- rep(0,length(parChange))

LVDet3Pred1Ext <- rep(0,length(parChange))
LVDet3Pred2Ext <- rep(0,length(parChange))
LVDet3Pred0Ext <- rep(0,length(parChange))

LVDet4Pred1Ext <- rep(0,length(parChange))
LVDet4Pred2Ext <- rep(0,length(parChange))
LVDet4Pred0Ext <- rep(0,length(parChange))


for(i in 1:length(parChange)){	
par[4] <- parChange[i]

#Solves the ODEs using deSolve Package
out <- ode(state, t, LVDet, par)
if(out[length(out[,1]),3]<1 && out[length(out[,1]),4]>1){
	LVDetPred1Ext[i] <- 1
}
if(out[length(out[,1]),4]<1 && out[length(out[,1]),3]>1){
	LVDetPred2Ext[i] <- 1
}
if(out[length(out[,1]),3]<1 && out[length(out[,1]),4]<1){
	if(which.min(out[,3])<which.min(out[,4])){
		LVDetPred1Ext[i] <- 1
	}
	if(which.min(out[,4])<which.min(out[,3])){
		LVDetPred2Ext[i] <- 1
	}
}
if(out[length(out[,1]),3]>1 && out[length(out[,1]),4]>1){
	Pred0Ext[i] <- 1
}

par[10] <- 0.001
out2 <- ode(state, t, LVDet2, par)

if(out2[length(out2[,1]),3]<1 && out2[length(out2[,1]),4]>1){
	LVDet2Pred1Ext[i] <- 1
}
if(out2[length(out2[,1]),4]<1 && out2[length(out2[,1]),3]>1){
	LVDet2Pred2Ext[i] <- 1
}
if(out2[length(out2[,1]),3]<1 && out2[length(out2[,1]),4]<1){
	if(which.min(out2[,3])<which.min(out2[,4])){
		LVDet2Pred1Ext[i] <- 1
	}
	if(which.min(out2[,4])<which.min(out2[,3])){
		LVDet2Pred2Ext[i] <- 1
	}
}
if(out2[length(out2[,1]),3]>1 && out2[length(out2[,1]),4]>1){
	LVDet2Pred0Ext[i] <- 1
}


par[10] <- 0.1
out3 <- ode(state, t, LVDet3, par)

if(out3[length(out3[,1]),3]<1 && out3[length(out3[,1]),4]>1){
	LVDet3Pred1Ext[i] <- 1
}
if(out3[length(out3[,1]),4]<1 && out3[length(out3[,1]),3]>1){
	LVDet3Pred2Ext[i] <- 1
}
if(out3[length(out3[,1]),3]<1 && out3[length(out3[,1]),4]<1){
	if(which.min(out3[,3])<which.min(out3[,4])){
		LVDet3Pred1Ext[i] <- 1
	}
	if(which.min(out3[,4])<which.min(out3[,3])){
		LVDet3Pred2Ext[i] <- 1
	}
}
if(out3[length(out3[,1]),3]>1 && out3[length(out3[,1]),4]>1){
	LVDet3Pred0Ext[i] <- 1
}

par[10] <- 10
out4 <- ode(state, t, LVDet4, par)

if(out4[length(out4[,1]),3]<1 && out4[length(out4[,1]),4]>1){
	LVDet4Pred1Ext[i] <- 1
}
if(out4[length(out4[,1]),4]<1 && out4[length(out4[,1]),3]>1){
	LVDet4Pred2Ext[i] <- 1
}
if(out4[length(out4[,1]),3]<1 && out4[length(out4[,1]),4]<1){
	if(which.min(out4[,3])<which.min(out4[,4])){
		LVDet4Pred1Ext[i] <- 1
	}
	if(which.min(out4[,4])<which.min(out4[,3])){
		LVDet4Pred2Ext[i] <- 1
	}
}
if(out4[length(out4[,1]),3]>1 && out4[length(out4[,1]),4]>1){
	LVDet4Pred0Ext[i] <- 1
}
}	


scaleText <- 1.5

pdf("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/9. Extras/Fig9h.pdf")
plot(RstarRatio,LVDetPred1Ext,type="l",xlab="R*1:R*2 Ratio",ylab="Probability of Earlier Extinction",lty=1,lwd=2,col="Red",cex.main=scaleText,cex.lab=scaleText,cex.axis=scaleText,xlim=c(0,2))
lines(RstarRatio,LVDet2Pred1Ext,col="Blue",lty=1,lwd=2)
lines(RstarRatio,LVDet3Pred1Ext,col="Green",lty=1,lwd=2)
lines(RstarRatio,LVDet4Pred1Ext,col="Orange",lty=1,lwd=2)
#points(RstarRatio,Pred2Ext,type="l",lty=1,lwd=2,col="purple")
#points(RstarRatio,Pred0Ext,type="l",lty=1,lwd=2,col="green")
legend("topright",c("Without imm","Imm Ratio=0.01","Imm Ratio=1","Imm Ratio=100"),col=c("Red","Blue","Green","Orange"),lty=1,cex=scaleText)
dev.off()
#exports results into a csv file
#write.csv(out,"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/2. Underlying parameters of R star in 2-consumer deterministic model/Data2b.csv")
