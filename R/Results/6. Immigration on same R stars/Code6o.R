scaleText <- 1.5
results <- read.csv("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/6. Immigration on same R stars/Data6o.csv")
myline.fit <- lm(results[,3] ~ results[,2])
summary(myline.fit)
pdf("/Users/lmq/Documents/IC/Final Year Project/Reports/Results/6. Immigration on same R stars/Fig6o.pdf")
plot(results[,2],results[,3],xlab="R*1:R*2 Ratio",ylab="Critical Immigration Ratio",ylim=c(0,3.5),xlim=c(0,2),cex.lab=scaleText,cex.axis=scaleText)
abline(myline.fit)
dev.off()