#### 3. Fit the models       -----------------------------
library(ggplot2)
library(tidyverse)
library(Hmsc)
library(corrplot)
library(colorspace)
library(vioplot)
library(writexl)




#A stand for automated, can be ran a an R job

load(file = "models/unfitted_BR4")

# Set up the amount of chains, iterations and samples
#samples_list = c(5,250,250, 250, 500, 5000) #500
#thin_list    = c(1,1,  10, 100, 100, 100)  #100

samples_list = c(5000) 
thin_list    = c(100)  

nst = length(thin_list)
nChains = 4


# Run the models
# This is a loop which will run the different models with the different amount of iterations and samples

for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  nm = length(models)
  for (model in 1:nm) {
    print(paste0("model = ",modelnames[model]))
    m = models[[model]]
    
    m = sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nParallel = nChains,
                   nChains = nChains)
    
    models[[model]] = m
  }
  filename = paste("models/models_BR4_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  
  save(models,modelnames,file=filename)
}