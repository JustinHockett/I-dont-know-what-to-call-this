# analytical workflows : coding demo
# Oct 12 2021


# Clear workspace and close all graphics ----------------------------------

rm(list=ls())
graphics.off()


# generate some data ------------------------------------------------------

n<- 30
m<- 2.2
b<- 1.1

x<- rnorm(n=n,mean= 0,sd=1)   #indep. var
z<- rnorm(n=n,mean= 0,sd=1)   #noise
y<- m * x + b + z             #dep. var

plot(x,y) #debug plot



# Do linear regression ----------------------------------------------------

fit<- glm(y ~ x)



# Plot the results --------------------------------------------------------

par(mar=c(6,7,4,1))
plot(x,y,
     xlab="Normalized velocity",
     ylab="Normalized\nenergy consumption"
     )
abline(fit)


# save data ---------------------------------------------------------------

#save.image("./data/my_expensive_calculations.Rdata")



# load data ---------------------------------------------------------------

#load(file="./data/my_expensive_calculations.Rdata")


#getwd()
#setwd("C:/GIT/I-dont-know-what-to-call-this/code/")
