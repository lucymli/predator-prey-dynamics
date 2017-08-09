#===============================================================================
# Lotka predator-prey model (Gillespie, 1977; Kot, 2001)
#===============================================================================

# This version of the Lotka predator-prey model is given by
# dY1/dt = c1*Y1 - c2*Y1*Y2
# dY2/dt = c2*Y1*Y2 - c3*Y2
# consisting of the three reaction channels,,
#      Y1 --c1--> Y1 + Y1 
# Y1 + Y2 --c2--> Y2 + Y2 
#      Y1 --c3--> 0

# Define parameters
parms <- c(c1=0.5, c2=0.2, c3=0.01, c4=0.5, c5=0.0005)

# Define system
x0 <- c(Y1=10, Y2=10)                           # Initial state vector
nu <- matrix(c(+1, -1, -1, 0, 0, 0, 1, -1),nrow=2,byrow=T) # State-change matrix
a  <- c("c1*Y1", "c2*Y1*Y2/(1+c2*c5*Y1)","c3*c2*Y1*Y2/(1+c2*c5*Y1)","c4*Y2")                # Propensity vector  
tf <- 2                                             # Final time
simName <- "Lotka predator-prey model"

# Run the simulations 
nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

# Direct method 
#set.seed(1)
out <- ssa(x0,a,nu,parms,tf,method="D",simName="K=0",verbose=FALSE,consoleInterval=1)
ssa.plot(out,show.title=TRUE,show.legend=TRUE)
