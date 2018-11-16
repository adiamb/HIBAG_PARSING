#setwd('/media/labcomp/HDD2/GWAS_EARLY_RUNS/')
arg=commandArgs(trailingOnly = T)
require(HIBAG)
library(parallel)
require(data.table)
cl <- makeCluster(8)
## mark dups and reconvert to plink in shell
#tped = fread('Plate41_take2_nodup_CHR6.tped')
#tped=tped[!duplicated(V2)]
#write.table(tped, sep = " ", file = 'Plate41_take2_nodup2_CHR6.tped', row.names = F, quote = F, col.names = F)
model.list <- get(load(arg[1]))
#model.list <- get(load('~/Downloads/InfiniumMEGA-European-HLA4-hg19.RData'))
# plink='/media/labcomp/HDD2/NMDA_GWAS_2018/NMDAR_and_controls'
# yourgeno = hlaBED2Geno(bed.fn = paste0(plink, ".bed"),
#                        bim.fn = paste0(plink, ".bim"),
#                        fam.fn = paste0(plink, ".fam"), verbose = T)

yourgeno = hlaBED2Geno(bed.fn = paste0(arg[2], ".bed"),
                       bim.fn = paste0(arg[2], ".bim"),
                       fam.fn = paste0(arg[2], ".fam"))
summary(yourgeno)

output_pred_QC=list()
hlas_id = names(model.list)

pred_hibag=function(hla, yourgeno, model.list){
  model = hlaModelFromObj(model.list[[hla]])
  print(model)
  pred.guess = predict(model, yourgeno, type = "response+prob", cl=cl, verbose = T, match.type="Position")
  return(pred.guess)
}

output_pred_QC=lapply(hlas_id, pred_hibag,  yourgeno = yourgeno, model.list = model.list)
names(output_pred_QC) = hlas_id[4]


hibag_extract=function(hla.num, output_pred_QC, filename){
  fwrite(output_pred_QC[[hla.num]]$value, file = paste0('IMPUTED_', output_pred_QC[[hla.num]]$locus, filename, ".txt"), row.names = T)
  fwrite(as.data.frame(output_pred_QC[[hla.num]]$postprob), file= paste0('POSTERIORS_', output_pred_QC[[hla.num]]$locus, filename, ".txt"), row.names = T)
  
}

lapply(1:length(output_pred_QC),  hibag_extract, output_pred_QC = output_pred_QC, filename=arg[3])

# for (i in 1:length(output_pred_QC)) {
#   print(output_pred_QC[[i]]$value)
#   fwrite(output_pred_QC[[i]]$value, file = paste0('IMPUTED_', output_pred_QC[[i]]$locus, arg[3], ".txt"), row.names = T)
#   fwrite(as.data.frame(output_pred_QC[[i]]$postprob), file= paste0('POSTERIORS_', output_pred_QC[[i]]$locus, arg[3], ".txt"), row.names = T)
# }

