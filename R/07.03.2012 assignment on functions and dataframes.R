#Function that takes two values and add them together
Add <- function(x,y) {
	return(x + y)
}


#Function that takes two vectors of the same length and produce a dataframe with the two vectors as columns
DF <- function (x,y) {
	n<-data.frame(x,y)
	return(n)
}


#Function that takes a dataframe and returns every other row of it
DF2 <- function(datfra) {
	x <- datfra
	n <- nrow(datfra)
	y <- seq(1,n,by=2)
	print(y)
	return(x[y,])
}