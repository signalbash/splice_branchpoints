args <- commandArgs(trailingOnly = TRUE)

k1=as.character(args[1])

setwd(k1)
options(stringsAsFactors = F)

chrom=c("chr1","chr2","chr3", "chr4","chr5",
        "chr6", "chr7", "chr8", "chr9","chr10",
        "chr11", "chr12", "chr13", "chr14","chr15",
        "chr16", "chr17", "chr18", "chr19","chr20",
        "chr21","chr22", "chrX","chrY")

library(caret)
library(kernlab)

#load in info files for each chromosome (!21) 
for(c in 1:length(chrom)){
branchpoint_df_s=read.csv(file=paste0("data/outputs/branchpoint_df_with_seq_r2_", chrom[c],".csv"))
branchpoint_df=branchpoint_df_s[,c(31,2,19,30,32,33,35:50)]

if(exists("mega_branchpoint_df")){
  mega_branchpoint_df=rbind(mega_branchpoint_df, branchpoint_df)
}else{
  mega_branchpoint_df= branchpoint_df
  
}
message(paste(chrom[c], "loaded"))
}
branchpoint_df=mega_branchpoint_df
rm(mega_branchpoint_df)

#only HC or NEG
branchpoint_df_HCN=branchpoint_df[branchpoint_df$set=="HC" | branchpoint_df$set=="NEG",]

#only variables
branchpoint_df_HCN_vars=branchpoint_df_HCN[,-c(1,2)]

#remove N values if lacking nt seqs
TF=apply(branchpoint_df_HCN_vars, 2, function(x) x=="N")
TF_v=apply(TF,1, any)
rm=which(TF_v==T)
if(length(rm) >0){
  branchpoint_df_HCN_vars=branchpoint_df_HCN_vars[-rm,]
  branchpoint_df_HCN=branchpoint_df_HCN[-rm,]
}

message("NA removed")

#make dummy vars for nucleotides A/T/C/G –> 1/0
dummies <- dummyVars(set ~ ., data=branchpoint_df_HCN[,-2])
branchpoint_df_HCN_vars_withDum <- predict(dummies, newdata=branchpoint_df_HCN[,-2])
branchpoint_df_HCN_vars_withDum <- cbind(set=branchpoint_df_HCN$set, branchpoint_df_HCN_vars_withDum)

message("Dummies made")

#make sure values are numeric...
for(i in 2:length(colnames(branchpoint_df_HCN_vars_withDum))){
	branchpoint_df_HCN_vars_withDum[,i] <- as.numeric(branchpoint_df_HCN_vars_withDum[,i])
}
message("as numeric")

#rename
filteredDescr <- branchpoint_df_HCN_vars_withDum

message("preprocessed")

save.image("data/Preprocessed_data.RData")
save(filteredDescr, branchpoint_df_HCN,dummies,TF_v, rm, file="data/Preprocessed_data_small.RData")
