################################################################################
#                      Process Genome-wide predictions                         #
################################################################################

options(scipen = 999)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(plyr)

#load information from training/testing sets to match + excude from predictions
load("data/Preprocessed_data.RData")
branchpoint_df <- branchpoint_df[branchpoint_df$set != "LC",]

#process files in parts
split_files = TRUE

if(split_files){
  files <- list.files("data/genome_predictions/", pattern = "gencode.v12_split")
  
  files_pred <- files[grep("pred", files)]
  files_attr <- files[grep("attr", files)]
  
  for(f in seq_along(files_attr)){
    gencode_v12_part <- as.data.frame(fread(paste0("data/genome_predictions/", files_pred[f])))
    gencode_v12_attributtes_part <- as.data.frame(fread(paste0("data/genome_predictions/", files_attr[f])))
    
    if(exists("gencode_v12")){
      gencode_v12 <- rbind(gencode_v12, gencode_v12_part)
      gencode_v12_attributtes <- rbind(gencode_v12_attributtes, gencode_v12_attributtes_part)
    }else{
      gencode_v12 <- gencode_v12_part
      gencode_v12_attributtes <- gencode_v12_attributtes_part
    }
  }
  #gencodev24
  files <- list.files("data/genome_predictions/", pattern = "gencode.v24_split")
  
  files_pred <- files[grep("pred", files)]
  files_attr <- files[grep("attr", files)]
  
  for(f in seq_along(files_attr)){
    gencode_v24_part <- as.data.frame(fread(paste0("data/genome_predictions/", files_pred[f])))
    gencode_v24_attributtes_part <- as.data.frame(fread(paste0("data/genome_predictions/", files_attr[f])))
    
    if(exists("gencode_v24")){
      gencode_v24 <- rbind(gencode_v24, gencode_v24_part)
      gencode_v24_attributtes <- rbind(gencode_v24_attributtes, gencode_v24_attributtes_part)
    }else{
      gencode_v24 <- gencode_v24_part
      gencode_v24_attributtes <- gencode_v24_attributtes_part
    }
  }
  
}else{
  
  gencode_v12 <- as.data.frame(fread("data/genome_predictions/gencode.v12_predictions.txt"))
  gencode_v12_attributtes <- as.data.frame(fread("data/genome_predictions/gencode.v12_attributes.txt"))
  
  gencode_v24 <- as.data.frame(fread("data/genome_predictions/gencode.v24_predictions.txt"))
  gencode_v24_attributtes <- as.data.frame(fread("data/genome_predictions/gencode.v24_attributes.txt"))
  
  gencode_v19 <- as.data.frame(fread("data/genome_predictions/gencode.v19_predictions.txt"))
  gencode_v19_attributtes <- as.data.frame(fread("data/genome_predictions/gencode.v19_attributes.txt"))
  
}

model_ids <- with(branchpoint_df, paste(seq_neg5,seq_neg4,seq_neg3,seq_neg2,seq_neg1, seq_pos0,
                                     seq_pos1, seq_pos2, seq_pos3, seq_pos4, seq_pos5,dist.2,
                                     ppt_start,ppt_run_length,canon_hit1,canon_hit2,canon_hit3,canon_hit4,canon_hit5, sep="_"))
model_chrom <- grep("chr",unlist(str_split(branchpoint_df$new_ID,"_")), value=TRUE)
model_ids <- paste0(model_chrom, model_ids)
sitep <- sapply(str_split(branchpoint_df$new_ID,c("[+]")),"[[",1)
siten <- sapply(str_split(branchpoint_df$new_ID,c("-")),"[[",1)

site <- sitep
site[grep("-", branchpoint_df$new_ID)] <- siten[grep("-", branchpoint_df$new_ID)]
site <- unlist(str_split(site,"_"))
site <- as.numeric(site[-grep("chr",site)]) +1

model_ids <- paste0(site,model_ids)

gencode_ids <- with(gencode_v12_attributtes, paste(seq_neg5,seq_neg4,seq_neg3,seq_neg2,seq_neg1, seq_pos0,
                                          seq_pos1, seq_pos2, seq_pos3, seq_pos4, seq_pos5,to_3prime,
                                          ppt_start,ppt_run_length,canon_hit1,canon_hit2,canon_hit3,canon_hit4,canon_hit5, sep="_"))
gencode_ids <- paste0(gencode_v12_attributtes$chromosome, gencode_ids)
gencode_ids <- paste0(gencode_v12_attributtes$end, gencode_ids)
m <- match(gencode_ids, model_ids)
gencode_v12_attributtes$in_testtrain <- 0
gencode_v12_attributtes$in_testtrain[which(!is.na(m))] <-1
gencode_v12_attributtes$in_testtrain[which(!is.na(m))] <- branchpoint_df$set[m[which(!is.na(m))]]

#combine attributes and predictions
att_name <- with(gencode_v12, paste(id,distance, allele_status, sep="_"))
m <- match(att_name, gencode_v12_attributtes$id)
gencode_v12 <- cbind(gencode_v12, gencode_v12_attributtes[m,])

gencode_v12 <- gencode_v12[,c(1,20,19,6,8,7,9,10,2,11,39,20:38)]
gencode_v12_seqs <- gencode_v12_attributtes$seq[m]

save(gencode_v12, gencode_v12_seqs, file="data/genome_predictions/gencode.v12.RData")

gencode_v24_seqs <- gencode_v24_attributtes$seq

m <- match(gencode_v24_seqs, gencode_v12_seqs)
gencode_v24_attributtes$in_testtrain <- 0
gencode_v24_attributtes$in_testtrain[which(!is.na(m))] <- gencode_v12$in_testtrain[m[which(!is.na(m))]]

#combine attributes and predictions
att_name <- with(gencode_v24, paste(id,distance, allele_status, sep="_"))
m <- match(att_name, gencode_v24_attributtes$id)
gencode_v24 <- cbind(gencode_v24, gencode_v24_attributtes[m,])
gencode_v24 <- gencode_v24[,c(1,20,19,6,8,7,9,10,2,11,39,20:38)]

save(gencode_v24, gencode_v24_seqs, file="data/genome_predictions/gencode.v24.RData")

gencode_v19_seqs <- gencode_v19_attributtes$seq

m <- match(gencode_v19_seqs, gencode_v12_seqs)
gencode_v19_attributtes$in_testtrain <- 0
gencode_v19_attributtes$in_testtrain[which(!is.na(m))] <- gencode_v19$in_testtrain[m[which(!is.na(m))]]

#combine attributes and predictions
att_name <- with(gencode_v19, paste(id,distance, allele_status, sep="_"))
m <- match(att_name, gencode_v19_attributtes$id)
gencode_v19 <- cbind(gencode_v19, gencode_v19_attributtes[m,])
gencode_v19 <- gencode_v24[,c(1,20,19,6,8,7,9,10,2,11,39,20:38)]

save(gencode_v19, gencode_v19_seqs, file="data/genome_predictions/gencode.v24.RData")

