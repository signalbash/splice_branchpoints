################################################################################
#         Human genetic variation at branchpoints - common variants            #
################################################################################

options(stringsAsFactors = F)

library(stringr)
library(branchpointer)
library(data.table)
library(biomaRt)
options(scipen=999)
library(GenomicRanges)

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



exons <- gtfToExons("data/gencode.v19.annotation.gtf")

q_ranges <- GRanges(seqnames=Rle(snp_info_query$chromosome),
                    ranges=IRanges(start=snp_info_query$chrom_start, end=snp_info_query$chrom_start),
                    strand=Rle("+"),
                    id=snp_info_query$id)

introns <- exons
posInd <- which(strand(exons) == "+")
start(introns[posInd]) <- start(introns[posInd]) - 51
end(introns[posInd]) <- start(introns[posInd]) + 50
negInd <- which(strand(exons) == "-")
end(introns[negInd]) <- end(introns[negInd]) + 51
start(introns[negInd]) <- end(introns[negInd]) - 50

o <- GenomicRanges::findOverlaps(q_ranges, introns)
keep <- unique(o@from)

snp_info_query <- snp_info_query[keep,]


reps <- ceiling(dim(snp_info_query)[1] / 1000)

for(r in 1:reps){
  to_pos <- r*1000
  from_pos <- to_pos-1000 +1
  to_pos <- min(nrow(snp_info_query), to_pos)

  snp_info_query_small <- snp_info_query[from_pos:to_pos,]
  write.csv(snp_info_query_small,file=paste0("data/GTEX/gtex_snps_part_",r,".csv"), row.names=FALSE, quote=FALSE)
}

###### branchpointer ######

for(r in 1:reps){
   
  query <- readQueryFile(paste0("data/GTEX/gtex_snps_part_",r,".csv"),queryType = "SNP", exons=exons, filter=FALSE)
   
  query_pred <- predictBranchpoints(query,queryType = "SNP",
                                    BSgenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  
  
  write.table(mcols(query_pred), file=paste0("data/GTEX/gtex_snps_part_", r,"_predictions.txt"), sep="\t",
              row.names=FALSE, quote=FALSE)
  message(r)

}

for(r in 1:reps){
  
  query <- readQueryFile(paste0("data/GTEX/gtex_snps_part_",r,".csv"),queryType = "SNP", exons=exons, filter=TRUE)
  
  query_pred <- data.table::fread(paste0("data/GTEX/gtex_snps_part_", r,"_predictions.txt"), data.table=FALSE)
  m <- match(query_pred$id, query$id)
  
  predictions <- query[m]
  mcols(predictions) <- query_pred
  
  summary <- branchpointer::predictionsToSummary(query, predictions, probabilityCutoff = 0.48, probabilityChange = 0.15)
  
  if(exists("summary_all")){
    summary_all <- c(summary_all, summary)
  }else{
    summary_all <- summary
  }
  
  message(r)
  
}

gtex_summary <- summary_all

gtex_summary$multi_BP[gtex_summary$BP_num_REF == 0] <-0
gtex_summary$multi_BP[gtex_summary$BP_num_REF == 1] <-1
gtex_summary$multi_BP[gtex_summary$BP_num_REF >1 ] <- "2+"

gtex_summary$dist_to_BP <- gtex_summary$dist_to_BP_REF
gtex_summary$dist_to_BP[gtex_summary$BP_num_REF == 0] <- gtex_summary$dist_to_BP_ALT[gtex_summary$BP_num_REF == 0]

#load in sQTL calls from GTEX
sQTL_ALL <- read.csv("data/GTEX/sQTL_ALL.csv", row.names=1)
eQTL_ALL <- read.csv("data/GTEX/eQTL_ALL.csv", row.names=1)
GTEX_variants <- as.data.frame(fread("data/GTEX/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt"))

GTEX_ID <- gsub("_neg","",gsub("_pos", "",gtex_summary$id))
GTEX_ID[which(str_sub(GTEX_ID,1,1) == "X")] <- gsub("X","", GTEX_ID[which(str_sub(GTEX_ID,1,1) == "X")])

GTEX_rsID <- GTEX_variants$RS_ID_dbSNP142_CHG37p13[match(GTEX_ID,GTEX_variants$VariantID)]
m <- match(sQTL_ALL$snpId,GTEX_rsID)
sQTL_ALL_BP <- cbind(sQTL_ALL[which(!is.na(m)),], gtex_summary[m[which(!is.na(m))]])

m <- match(eQTL_ALL$snp, GTEX_ID)
eQTL_ALL_BP <- cbind(eQTL_ALL[which(!is.na(m)),], gtex_summary[m[which(!is.na(m))]])

m <- match(GTEX_rsID,sQTL_ALL$snpId)
gtex_summary$is_sQTL <- "no"
gtex_summary$is_sQTL[which(!is.na(m))] <- "yes"

m <- match(GTEX_ID,eQTL_ALL$snp)
gtex_summary$is_eQTL <- "no"
gtex_summary$is_eQTL[which(!is.na(m))] <- "yes"

gtex_summary$snp_in <- NA
gtex_summary$snp_in[gtex_summary$dist_to_BP <= -1] <- "3'SS to BP"
gtex_summary$snp_in[gtex_summary$dist_to_BP >= 3] <- "BP to 5'SS"
gtex_summary$snp_in[gtex_summary$dist_to_BP == 0 | gtex_summary$dist_to_BP == 2] <- "BP"
gtex_summary$snp_in[gtex_summary$to_3prime < 3] <- "3'SS"

gtf <- rtracklayer::import("data/gencode.v19.annotation.gtf")
gtf <- gtf[gtf$type=="exon"]
o <- GenomicRanges::findOverlaps(gtex_summary, gtf)
rm <- unique(o@from)
gtex_summary_filtered <- gtex_summary[-rm]

gtex_summary_filtered <- as.data.frame(gtex_summary_filtered)[(gtex_summary_filtered$created_n > 0 | gtex_summary_filtered$deleted_n > 0) &
                                                                      gtex_summary_filtered$to_3prime > 2,]
keep <- which((gtex_summary_filtered$max_prob_REF != gtex_summary_filtered$max_prob_ALT) | 
                (gtex_summary_filtered$max_U2_REF != gtex_summary_filtered$max_U2_ALT) | 
              is.na(gtex_summary_filtered$max_prob_ALT) | is.na(gtex_summary_filtered$max_prob_REF))
gtex_summary_filtered <- gtex_summary_filtered[keep,]

snp_locs <- as.data.frame(table(gtex_summary$snp_in[gtex_summary$multi_BP != 0], 
                                gtex_summary$multi_BP[gtex_summary$multi_BP != 0]))

snp_locs$percent_Freq <- (snp_locs$Freq / rep(as.data.frame(table(gtex_summary$snp_in[gtex_summary$multi_BP != 0]))$Freq,2)) 

snp_locs$percent_Freq <- snp_locs$percent_Freq / c(rep(length(which(gtex_summary$multi_BP == 1))/length(which(gtex_summary$multi_BP != 0)),4),
rep(length(which(gtex_summary$multi_BP == "2+"))/length(which(gtex_summary$multi_BP != 0)),4)) - 1

snp_locs$Var2 <- as.character(snp_locs$Var2)
snp_locs$Var2[snp_locs$Var2 == "1"] <- "single BP"
snp_locs$Var2[snp_locs$Var2 == "2+"] <- "multiple BPs"
snp_locs$Variant_location <- paste(snp_locs$Var1,snp_locs$Var2, sep=" | ")

save(snp_locs, gtex_summary, gtex_summary_filtered, file="data/commonVariants.Rdata")

write.csv(gtex_summary_filtered, file="Tables/TableS4.csv", row.names=FALSE)
