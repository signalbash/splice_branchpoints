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
  files <- list.files("data/genome_predictions/", pattern = "gencode_v19_split")
  
  files_pred <- files[grep("pred", files)]

  for(f in seq_along(files_pred)){
    gencode_v19_part <- fread(paste0("data/genome_predictions/", files_pred[f]), data.table = FALSE)
    #gencode_v19_attributtes_part <- as.data.frame(fread(paste0("data/genome_predictions/", files_attr[f])))
    message(f)
    if(exists("gencode_v19")){
      gencode_v19 <- rbind(gencode_v19, gencode_v19_part)
      #gencode_v19_attributtes <- rbind(gencode_v19_attributtes, gencode_v19_attributtes_part)
    }else{
      gencode_v19 <- gencode_v19_part
      #gencode_v19_attributtes <- gencode_v19_attributtes_part
    }
  }
  #gencodev26
  files <- list.files("data/genome_predictions/", pattern = "gencode_v26_r2_split")
  
  files_pred <- files[grep("pred", files)]
  files_attr <- files[grep("attr", files)]
  
  for(f in seq_along(files_pred)){
    gencode_v26_part <- fread(paste0("data/genome_predictions/", files_pred[f]), data.table = FALSE)
    #gencode_v26_attributtes_part <- as.data.frame(fread(paste0("data/genome_predictions/", files_attr[f])))
    message(f)
    if(exists("gencode_v26")){
      gencode_v26 <- rbind(gencode_v26, gencode_v26_part)
      #gencode_v26_attributtes <- rbind(gencode_v26_attributtes, gencode_v26_attributtes_part)
    }else{
      gencode_v26 <- gencode_v26_part
      #gencode_v26_attributtes <- gencode_v26_attributtes_part
    }
  }
  
}else{
  
  #gencode_v19 <- as.data.frame(fread("data/genome_predictions/gencode_v19_predictions.txt"))
  #gencode_v19_attributtes <- as.data.frame(fread("data/genome_predictions/gencode.v19_attributes.txt"))
  
  gencode_v26 <- fread("data/genome_predictions/gencode_v26_predictions.txt", data.table = FALSE)
  gencode_v26_attributtes <- fread("data/genome_predictions/gencode.v26_attributes.txt", data.table = FALSE)
  
  gencode_v19 <- fread("data/genome_predictions/gencode_v19_predictions.txt", data.table = FALSE)
  #gencode_v19_attributtes <- fread("data/genome_predictions/gencode.v19_attributes.txt", data.table = FALSE)
  
}
branchpoint_df$seq_motif <- with(branchpoint_df, paste0(seq_neg5,seq_neg4,seq_neg3,seq_neg2,seq_neg1, seq_pos0,
                                        seq_pos1, seq_pos2, seq_pos3, seq_pos4, seq_pos5))
branchpoint_df$chrom <- unlist(lapply(str_split(branchpoint_df$new_ID, "_"), "[[", 1))
branchpoint_df$strand <- "-"
branchpoint_df$strand[grepl("[+]", branchpoint_df$new_ID)] <- "+"

model_ids <- with(branchpoint_df, paste(chrom, strand, seq_motif,dist.2, sep="_"))
gencode_ids <- with(gencode_v19, paste(chromosome, strand, seq_motif, to_3prime,sep="_"))

m <- match(gencode_ids, model_ids)
gencode_v19$in_testtrain <- 0
gencode_v19$in_testtrain[which(!is.na(m))] <-1
gencode_v19$in_testtrain[which(!is.na(m))] <- branchpoint_df$set[m[which(!is.na(m))]]

save(gencode_v19, file="data/genome_predictions/gencode_v19.RData")

gencode_ids_v26 <- with(gencode_v26, paste(chromosome, strand, seq_motif, to_3prime,sep="_"))

m <- match(gencode_ids_v26, model_ids)
gencode_v26$in_testtrain <- 0
gencode_v26$in_testtrain[which(!is.na(m))] <-1
gencode_v26$in_testtrain[which(!is.na(m))] <- branchpoint_df$set[m[which(!is.na(m))]]

save(gencode_v26, file="data/genome_predictions/gencode_v26.RData")
