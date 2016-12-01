################################################################################
#         Human genetic variation at branchpoints - common variants            #
################################################################################

options(stringsAsFactors = F)

library(stringr)
library(branchpointer)
library(data.table)
library(biomaRt)
options(scipen=999)

###### Load in GTEx SNPs ######
gtex_vars <- fread("data/GTEX/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt")
gtex_vars <- as.data.frame(gtex_vars)

gtex_vars$Chr <- paste0("chr", gtex_vars$Chr)

norm <- nchar(gtex_vars$Ref_b37)
alt <- nchar(gtex_vars$Alt)

keep <- which(norm == 1 & alt == 1)
gtex_vars <- gtex_vars[keep,]

snp_info_query <- data.frame(id=gtex_vars[,3], chromosome=gtex_vars[,1], 
                          chrom_start=gtex_vars[,2], strand=2,
                          allele_ref=gtex_vars[,4],
                          allele_alt=gtex_vars[,5])

reps <- ceiling(length(gtex_vars$VariantID) / 500000)

for(r in 1:reps){
  to_pos=min(r*500000)
  from_pos=to_pos-(500000 +1)
  to_pos=min(length(gtex_vars$VariantID), to_pos)

  snp_info_query_small=snp_info_query[from_pos:to_pos,]
  write.csv(snp_info_query_small,file=paste0("data/GTEX/gtex_snps_part_",r,".csv"), row.names=FALSE, quote=FALSE)
}

###### branchpointer ######

source("scripts/analysis/variants/predictions_to_stats.R")

exons <- readExonAnnotation("data/genome_predictions/gencode.v19.exons.txt")
for(r in 1:reps){
   
  query <- readQueryFile(paste0("data/GTeX/gtex_snps_part_",r,".csv"),query_type = "SNP")
   
  queryLoc <- getQueryLoc(query, query_type = "SNP", exons=exons, filter=TRUE, max_dist=50)

  query_attributes <- getBranchpointSequence(queryLoc,
                                          query_type = "SNP",
                                          genome = "data/genome_predictions/GRCh37.p13.genome.fa",
                                          bedtools_location="/Applications/apps/bedtools2/bin/bedtools")

  write.table(query_attributes, file=paste0("data/GTEX/GTEX_part_", r,"_attributes.txt"), sep="\t",
              row.names=FALSE, quote=FALSE)

  branchpoint_predictions <- predictBranchpoints(query_attributes) 

  write.table(branchpoint_predictions, file=paste0("data/GTEX/GTEX_part_", r,"_predictions.txt"), sep="\t",
              row.names=FALSE, quote=FALSE)

  pred_id <- with(branchpoint_predictions, paste0(id,"_",distance,"_",allele_status))
  m <- match(pred_id, query_attributes$id)

  GTEX_predictions <- cbind(branchpoint_predictions, query_attributes[m,-c(1:7)])

  snp_ids <- unique(GTEX_predictions$id)

  source("scripts/analysis/variation/predictions_to_stats.R")

  GTEX_processed <- predictions_to_stats(GTEX_predictions, snp_ids, query)

  write.table(GTEX_processed , file=paste0("data/GTEX/GTEX_part_", r,"_processed.txt"), sep="\t",
              row.names=FALSE, quote=FALSE)
}

nt_cols=c("#359646","#4D7ABE","#FAA859","#CB3634")

files <- list.files("data/GTeX/")
processed_files <- files[grep("processed", files)]

rm(processed)

for(f in processed_files){
  processed_part <- as.data.frame(fread(paste0("data/GTeX/", f)))
  
  if(exists("processed")){
    processed <- rbind(processed, processed_part)
  }else{
    processed <- processed_part
  }
}

processed$multi_BP[processed$BP_REF_num ==0] <-0
processed$multi_BP[processed$BP_REF_num >1 ] <- "2+"

processed_filtered <- processed[processed$created_n >0 | processed$deleted_n > 0,]
processed_GTEX <- processed
processed_GTEX_filtered <- processed_GTEX[processed_GTEX$created_n > 0 | processed_GTEX$deleted_n > 0,]

#load in sQTL calls from GTEX
sQTL_ALL <- read.csv("data/GTeX/sQTL_ALL.csv", row.names=1)
eQTL_ALL <- read.csv("data/GTeX/eQTL_ALL.csv", row.names=1)
GTEX_variants <- as.data.frame(fread("data/GTeX/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt"))

GTEX_ID <- gsub("_neg","",gsub("_pos", "",processed_GTEX_filtered$id))
GTEX_ID[which(str_sub(GTEX_ID,1,1) == "X")] <- gsub("X","", GTEX_ID[which(str_sub(GTEX_ID,1,1) == "X")])

GTEX_rsID <- GTEX_variants$RS_ID_dbSNP142_CHG37p13[match(GTEX_ID,GTEX_variants$VariantID)]
m <- match(sQTL_ALL$snpId,GTEX_rsID)
sQTL_ALL_BP <- cbind(sQTL_ALL[which(!is.na(m)),], processed_GTEX_filtered[m[which(!is.na(m))],])

m <- match(eQTL_ALL$snp, GTEX_ID)
eQTL_ALL_BP <- cbind(eQTL_ALL[which(!is.na(m)),], processed_GTEX_filtered[m[which(!is.na(m))],])

m <- match(GTEX_rsID,sQTL_ALL$snpId)
processed_GTEX_filtered$is_sQTL <- "no"
processed_GTEX_filtered$is_sQTL[which(!is.na(m))] <- "yes"

m <- match(GTEX_ID,eQTL_ALL$snp)
processed_GTEX_filtered$is_eQTL <- "no"
processed_GTEX_filtered$is_eQTL[which(!is.na(m))] <- "yes"

processed_GTEX$dist_to_BP <- processed_GTEX$dist_to_BP_REF
processed_GTEX$dist_to_BP[processed_GTEX$BP_REF_num == 0] <- processed_GTEX$dist_to_BP_ALT[processed_GTEX$BP_REF_num == 0]
processed_GTEX$multi_BP[processed_GTEX$BP_REF_num == 0] <- "0"

processed_GTEX$snp_in <- NA
processed_GTEX$snp_in[processed_GTEX$dist_to_BP <= -1] <- "3'SS to BP"
processed_GTEX$snp_in[processed_GTEX$dist_to_BP >= 3] <- "BP to 5'SS"
processed_GTEX$snp_in[processed_GTEX$dist_to_BP == 0 | processed_GTEX$dist_to_BP == 0 | processed_GTEX$dist_to_BP == 0] <- "BP"
processed_GTEX$snp_in[processed_GTEX$dist_to_exon < 3] <- "3'SS"

snp_locs <- as.data.frame(table(processed_GTEX$snp_in[processed_GTEX$multi_BP != 0], processed_GTEX$multi_BP[processed_GTEX$multi_BP != 0]))

snp_locs$percent_Freq <- (snp_locs$Freq / rep(as.data.frame(table(processed_GTEX$snp_in[processed_GTEX$multi_BP != 0]))$Freq,2)) 

snp_locs$percent_Freq <- snp_locs$percent_Freq / c(rep(length(which(processed_GTEX$multi_BP == 1))/length(which(processed_GTEX$multi_BP != 0)),4),
rep(length(which(processed_GTEX$multi_BP == "2+"))/length(which(processed_GTEX$multi_BP != 0)),4)) - 1

snp_locs$Var2 <- as.character(snp_locs$Var2)
snp_locs$Var2[snp_locs$Var2 == "1"] <- "single BP"
snp_locs$Var2[snp_locs$Var2 == "2+"] <- "multiple BPs"

snp_locs$Variant_location <- paste(snp_locs$Var1,snp_locs$Var2, sep=" | ")

save(snp_locs, processed_GTEX, file="data/commonVariants.Rdata")

write.csv(processed_GTEX_filtered, file="Tables/TableS3.csv", row.names=FALSE)
