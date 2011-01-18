######################################
######################################
#######BM PROCESS#####################
######################################
######################################
#postscript(file="/home/djhwueng/Desktop/bm1bm5.ps",
postscript(file="C:\\Users\\Tony\\Desktop\\bm1.ps",
paper="special",
width=10,
height=8,
horizontal=FALSE)

#sim.ou <- function(alpha,sigma){
#  x <- rep(0,path)
#  for (i in 2:path)
#   x[i] <- x[i-1]-alpha * (x[i - 1]-mu)+ sigma*rnorm(1)
#    x # returns the value of x
#  }

sigma<-0.5
#alpha<-0.2
mu<-20
path<-15
repath<-3

sigma<-1
#BM
X1 <- replicate(repath, cumsum(c(0, sigma*rnorm(path-1))))
#OU

sigma<-5
X2 <- replicate(repath, cumsum(c(0, sigma*rnorm(path-1))))
#X2 <- replicate(repath, sim.ou(alpha,sigma))

#layout(matrix(1:2, c(1,2)))
#yl <- range(X2)
yl <- range(X1)
matplot(X1, ylim = yl, type = "l", col = 1:repath, cex.lab=1.5 ,xlab = "Time t", ylab = "Trait Value y(t)" ,cex.main=2,main = expression(paste(sigma,"=1"))) 

#matplot(X2, ylim = yl, type = "l", col = 1:repath, cex.lab=1.5 ,xlab = "Time t", ylab = "Trait Value y(t)" ,cex.main=2,main = expression(paste(sigma,"=5"))) 

dev.off()


#plot(density(X1[path ,]))
#plot(density(X2[path ,]))

var(X1[path ,])
var(X2[path ,])

######################################
######################################
#######OU PROCESS
######################################
######################################
postscript(file="/home/djhwueng/Desktop/ou1ou5.ps",
paper="special",
width=10,
height=8,
horizontal=FALSE)

sim.ou <- function(alpha,sigma){
  x <- rep(0,path)
  for (i in 2:path)
    x[i] <- x[i-1]-alpha * (x[i - 1]-mu)+ sigma*rnorm(1)
    x # returns the value of x
  }


path<-100
repath<-5

sigma<-1
mu<-20

alpha<-0.05
X1 <- replicate(repath, sim.ou(alpha,sigma))
alpha<-0.15
X2 <- replicate(repath, sim.ou(alpha,sigma))

layout(matrix(1:2, c(1,2)))
yl <- range(X1)
matplot(X1, ylim = yl, type = "l", col = 1:repath, cex.lab=1.5 ,xlab = "Time t", ylab = "Trait Value y(t)" ,cex.main=2,main = expression(paste(alpha,"=0.05",", " ,mu,"=20"))) 

matplot(X2, ylim = yl, type = "l", col = 1:repath, cex.lab=1.5 ,xlab = "Time t", ylab = "Trait Value y(t)" ,cex.main=2,main = expression(paste(alpha,"=0.15",", " ,mu,"=20"))) 

dev.off()
