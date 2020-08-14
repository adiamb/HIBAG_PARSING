## author : ambati@stanford.edu
arg=commandArgs(trailingOnly = T)
#arg1 = model
#arg2 = plink bed files in hg19 coordinates
#arg3 = output file
#'''load the libraries and setup a cluster'''
sink(paste0(arg[3], '.log'))
PackList = rownames(installed.packages())
if('data.table' %in%  PackList & 'parallel' %in% PackList){
  require(data.table)
  require(ggplot2)
  require(stringr)
} else {
  install.packages('data.table')
}
library(parallel)
require(data.table)
if('HIBAG' %in% PackList){
  require(HIBAG)
} else {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("HIBAG")
  require(HIBAG)
}
### make n-1 clusters
n.cores = detectCores() -1
cl <- makeCluster(n.cores)
model.list <- get(load(arg[1])) 
yourgeno = hlaBED2Geno(bed.fn = paste0(arg[2], ".bed"),
                       bim.fn = paste0(arg[2], ".bim"),
                       fam.fn = paste0(arg[2], ".fam"))
summary(yourgeno)
output_pred_QC=list()
hlas_id = names(model.list)
#'''predict function to wrap in a lapply call'''
pred_hibag=function(hla, yourgeno, model.list){
  model = hlaModelFromObj(model.list[[hla]])
  print(model)
  pred.guess = predict(model, yourgeno, type = "response+prob", cl=cl, verbose = T, match.type="Position")
  cat('PREDICTED ', hla, timestamp())
  return(pred.guess$value)
}

output_pred_QC=lapply(hlas_id, pred_hibag,  yourgeno = yourgeno, model.list = model.list)
names(output_pred_QC) = hlas_id
## stop the parallel cluster
stopCluster(cl = cl)
## combined the outputs and write out
combinedPreds =data.frame(sample.id=output_pred_QC[[1]]$sample.id)
for(i in seq(1:length(model.list))){
  hla.names = c(paste0(hlas_id[i], '.1'), paste0(hlas_id[i], '.2'), paste0(hlas_id[i], '.prob'))
  temp = output_pred_QC[[i]]
  names(temp)[2:4] = hla.names
  if (identical(combinedPreds$sample.id, temp$sample.id)){
    combinedPreds=cbind.data.frame(combinedPreds, temp[, 2:4])
  } else {
    stop('SAMPLE IDS are not in same order check !!')
  }
}
print(dim(combinedPreds))
## write out the file
fwrite(combinedPreds, file=arg[3])
sink()
