par <- c(r=1,d=0.1,a=0.2,c=0.1,K=10)
Rstar <- seq(0,10,0.1)
dlist <- Rstar*par[3]*par[4]
alist <- par[2]/Rstar/par[4]
clist <- par[2]/Rstar/par[3]
eqmP1d <- par[1]/par[3]*(1-dlist/par[4]/par[3]/par[5])
eqmP1a <- par[1]/alist*(1-par[2]/par[4]/alist/par[5])
eqmP1c <- par[1]/par[3]*(1-par[2]/clist/par[3]/par[5])
write.csv(data.frame(Rstar,eqmP1d,eqmP1a,eqmP1c),"/Users/lmq/Documents/IC/Final Year Project/Reports/Results/7. Initial conditions on extinction risk/eqmRstar.csv")

plot(Rstar,eqmP1d,col="red")
points(Rstar,eqmP1a,col="green")
points(Rstar,eqmP1c,col="blue")