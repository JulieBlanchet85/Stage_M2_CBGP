#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)


#You may wish to loop over samples_list and thinning_list as done in Script S3 
nChains = 4
thin = 100
samples = 5000
transient = round(0.5*samples*thin)
nCores=12 #often set to nChains, but with the new pcomputepredict function can be set higher for more efficiency depending on your amount of cpu's

print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))

filename_in = paste("models/models_BR4_", as.character(thin),
                    "_samples_", as.character(samples),
                    "_chains_",as.character(nChains),
                    ".Rdata",sep = "")

load(file = filename_in) #models, modelnames

nm = length(models)

MF = list()
MFCV = list()
WAIC = list()

#zie code overfitting.R en boek p163-164
#partition.spat = createPartition(m, nfolds = 2, column = 3)

for(model in 1:nm){
  print(paste0("model = ", as.character(model)))
  m = models[[model]]
  preds = computePredictedValues(m)
  MF[[model]] = evaluateModelFit(hM=m, predY=preds)
  partition = createPartition(m, nfolds = 3) 
  preds = pcomputePredictedValues(m,partition=partition, nParallel = nChains) #this is the new function for better cpu disribution
  MFCV[[model]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[model]] = computeWAIC(m)
}

filename_out = paste("models/BR4_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")

save(MF,MFCV,WAIC,modelnames,file = filename_out)
