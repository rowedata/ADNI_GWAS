#one million iterations analysis

#install coda if not installed
library(coda)

#set your working directory for R
setwd("C:/Users/Benazir/Documents")

# path and prefix to mcmc output files for 100million iterations run

folder ="dissertation/convergence/adni_binary_10M_run1"

#read the  output files

path1 <- read.table(paste0(folder,".path.txt"), header= T)
snp1 <- read.table(paste0(folder,".snp.txt"), header= T)
mcmc1 <- read.table(paste0(folder,".mcmc.txt"), header= T)
gamma1 <- read.table(paste0(folder,".gamma.txt"), header= T)

# exploratory visuals
head(path1)
hist1info <- hist(path1$hh,freq =FALSE, labels=TRUE,main="Density plot for one million steps")
hist1info$density


########################################################################

#gelman rubin (MCMC diagnostics based on multiple runs)

# convert data to mcmc object
chain1 = as.mcmc(path1$hh)
summary(chain1)
plot(chain1)

# put path to the output of the second chain including the prefix

folder_run2 ="dissertation/convergence/adni_binary_10M_run1"
run2<- read.table(paste0(folder_run2,".path.txt"), header= T)

# convert data to mcmc object
chain2 = as.mcmc(run2$hh[1:100000])
summary(chain2)
plot(chain2)

# combine two (or more) chains 
combinedchains = mcmc.list(chain1, chain2)

#use the output of the next three lines as a main info to base convergence decision on
plot(combinedchains, main="Combined chains")
gelman.diag(combinedchains,confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)
gelman.plot(combinedchains, main= "Scale reduction factor over time")

#additional geweke  diagnostics
geweke.diag(chain1, frac1=0.1, frac2=0.5)
geweke.plot(chain1, frac1 = 0.1, frac2 = 0.5)


####################################################################################
#analysis of subchains (useful but not necessary step) add zeros for larger chains
hhone1 =path1$hh[1:10000]
hhone2 =path1$hh[10001:20000]
hhone3 =path1$hh[20001:30000]
hhone4 =path1$hh[30001:40000]
hhone5 =path1$hh[40001:50000]
hhone6 =path1$hh[50001:60000]
hhone7 =path1$hh[60001:70000]
hhone8 =path1$hh[70001:80000]
hhone9 =path1$hh[80001:90000]
hhone10 =path1$hh[90001:100000]

means <- matrix(,10,4)
for (i in 1:10)
{
  means[i,1] = mean(get(paste("hhone",i,sep="")))
  means[i,2] = median(get(paste("hhone",i,sep="")))
  means[i,3] = sd(get(paste("hhone",i,sep="")))
  
}
means
colnames(means)= c("mean","median", "sd", "gr")
write.table(means, file=paste0(folder,"onemillion_binary.txt"), quote = F, row.names = F)
mean(path1$hh)
median(path1$hh)
sd(path1$hh)

#plots
par(mfrow=c(1,1))
plot(path1$hh, type="l")

par(mfrow=c(2,2))
plot(hhone1, type="l", main = "1:10000")
plot(hhone2, type="l", main = "10001:20000")
plot(hhone3, type="l", main = "20001:30000")
plot(hhone4, type="l", main = "30001:40000")

plot(hhone5, type="l", main = "40001:50000")
plot(hhone6, type="l", main = "50001:60000")
plot(hhone7, type="l", main = "60001:70000")
plot(hhone8, type="l", main = "70001:80000")
par(mfrow=c(1,2))
plot(hhone9, type="l", main = "80001:90000")
plot(hhone10, type="l", main = "90001:100000")


#zoom in 
hist(path1$hh, labels=TRUE,main="Density plot for one million steps", breaks=c(0,0.0001, 0.0011, 0.0021, 0.0031, 0.0041, 0.0051, 0.0061, 0.0071, 0.0081, 0.0091, 0.021))
