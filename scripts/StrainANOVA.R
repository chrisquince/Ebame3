#!/usr/bin/Rscript

#load libraries
library(getopt)

spec = matrix(c('gfile','c',1,"character",'samplefile','s',1,"character",'stub','t',1,"character"),byrow=TRUE,ncol=4)

opt=getopt(spec)



gammaFile <- opt$gfile
sampleFile <- opt$samplefile
stub <- opt$stub

Gamma <- read.csv(gammaFile,header=TRUE,row.names=1)

rownames(Gamma) <- gsub("m","M",rownames(Gamma))

rownames(Gamma) <- gsub("_5M","_05M",rownames(Gamma))

Samples <- read.csv(sampleFile,header=TRUE,row.names=1)

SamplesR <- Samples[rownames(Gamma),]

NC <- ncol(Gamma) 
pvalues <- rep(-1,NC)
tvalues <- rep(-1,NC)  

gammaMean <- colMeans(Gamma)

for (i in 1:NC){
  Test.df <- data.frame(Region=SamplesR$region,Gamma=Gamma[,i])
  
  test <- kruskal.test(Gamma ~ Region, data = Test.df) 
  
  cat(sprintf("%s,%d,%f,%f,%f\n",stub,i,gammaMean[i],test$'p.value',test$statistic))
}
