################################################################################
#                   Human genetic variation at branchpoints                    #
################################################################################

options(stringsAsFactors = F)

library(biomaRt)
library(stringr)
library(branchpointer)
library(data.table)
nt_cols=c("#359646","#4D7ABE","#FAA859","#CB3634")

#load clinVar annotation
#downloaded 2016/07/26
ClinVar_variant_summary <- read.delim("data/ClinVar/variant_summary_clinVar_downloaded20161106.txt")
ClinVar_variant_summary_hg38 <- ClinVar_variant_summary[ClinVar_variant_summary$Assembly=="GRCh38",]
ClinVar_variant_summary_hg38 <- ClinVar_variant_summary_hg38[ClinVar_variant_summary_hg38$Type=="single nucleotide variant",]

#convert ambigous codes
altallele <- ClinVar_variant_summary_hg38$AlternateAllele
ClinVar_variant_summary_hg38$AlternateAllele[altallele=="Y"] <- "C/T"
ClinVar_variant_summary_hg38$AlternateAllele[altallele=="H"] <- "A/C/T"
ClinVar_variant_summary_hg38$AlternateAllele[altallele=="R"] <- "A/G"
ClinVar_variant_summary_hg38$AlternateAllele[altallele=="D"] <- "A/G/T"

ClinVar_variant_summary_hg38 <- ClinVar_variant_summary_hg38[-which(altallele=="na"),]

ClinVar_variant_summary_hg38$new_id <- with(ClinVar_variant_summary_hg38, paste0("clinVar:",X.AlleleID,"+rs",RS...dbSNP.))

ClinVar_variant_summary_hg38_small <- ClinVar_variant_summary_hg38[,c("new_id","Chromosome",
                                                                      "Start","X.AlleleID",
                                                                      "ReferenceAllele","AlternateAllele")]

colnames(ClinVar_variant_summary_hg38_small)<- c("id","chromosome","chrom_start","strand", "ref_allele","alt_allele")
ClinVar_variant_summary_hg38_small$strand <- 2

snp_info_clinVar_single <- ClinVar_variant_summary_hg38_small[which(ClinVar_variant_summary_hg38_small$alt_allele %in% c("A","T","C",'G')),]
snp_info_clinVar_multi <- ClinVar_variant_summary_hg38_small[-which(ClinVar_variant_summary_hg38_small$alt_allele %in% c("A","T","C",'G')),]

#convert multiple alternatives to single
alts <- str_split(snp_info_clinVar_multi$alt_allele,"/")
for(i in seq_along(snp_info_clinVar_multi$alt_allele)){
  snp_info_clinVar_multi$alt_allele[i] <- alts[[i]][1] 
  
  new_line <- snp_info_clinVar_multi[rep(i, length(alts[[i]])-1),]
  new_line$alt_allele <- alts[[i]][-1]
  
  snp_info_clinVar_multi <- rbind(snp_info_clinVar_multi, 
                                new_line)
}
snp_info_clinVar <- rbind(snp_info_clinVar_single,
                        snp_info_clinVar_multi)
snp_info_clinVar$chromosome <- paste0("chr",snp_info_clinVar$chromosome)
snp_info_clinVar$id <- with(snp_info_clinVar, paste0(id, "_",ref_allele,"/",alt_allele))
write.csv(snp_info_clinVar, file="data/ClinVar/clinvar_info_24.csv", row.names = FALSE)

###### submit to branchpointer ######
exons <- readExonAnnotation("data/genome_annotations/gencode.v24.annotation.exons.txt")
 
query <- readQueryFile("data/ClinVar/clinvar_info_24.csv",query_type = "SNP")
 
queryLoc <- getQueryLoc(query, query_type = "SNP", exons=exons, filter=TRUE, max_dist=50)

query_attributes <- getBranchpointSequence(queryLoc,
                                        query_type = "SNP",
                                        genome = "data/genome_annotations/GRCh38.p5.genome.fa",
                                        bedtools_location="/Applications/apps/bedtools2/bin/bedtools")

write.table(query_attributes, file=paste0("data/ClinVar/clinvar_attributes.txt"), sep="\t",
            row.names=FALSE, quote=FALSE)

branchpoint_predictions <- predictBranchpoints(query_attributes) 

write.table(branchpoint_predictions, file=paste0("data/ClinVar/clinvar_predictions.txt"), sep="\t",
            row.names=FALSE, quote=FALSE)

##### process SNP effects #######

clinvar_attributes <- as.data.frame(fread("data/ClinVar/clinvar_attributes.txt"))
clinvar_predictions <- as.data.frame(fread("data/ClinVar/clinvar_predictions.txt"))

pred_id <- with(clinvar_predictions, paste0(id,"_",distance,"_",allele_status))
m <- match(pred_id, clinvar_attributes$id)

clinvar_predictions <- cbind(clinvar_predictions, clinvar_attributes[m,-c(1:7)])

snp_ids <- unique(clinvar_predictions$id)

source("scripts/analysis/variation/predictions_to_stats.R")

clinVar_processed <- predictions_to_stats(clinvar_predictions, snp_ids, query)

clinVar_processed$allele_ref_feat_strand <- clinVar_processed$ref_allele
clinVar_processed$allele_ref_feat_strand[clinVar_processed$strand=="-" & clinVar_processed$ref_allele=="A"] <- "T"
clinVar_processed$allele_ref_feat_strand[clinVar_processed$strand=="-" & clinVar_processed$ref_allele=="T"] <- "A"
clinVar_processed$allele_ref_feat_strand[clinVar_processed$strand=="-" & clinVar_processed$ref_allele=="G"] <- "C"
clinVar_processed$allele_ref_feat_strand[clinVar_processed$strand=="-" & clinVar_processed$ref_allele=="C"] <- "G"

clinVar_processed$multi_BP <- "1"
clinVar_processed$multi_BP[clinVar_processed$BP_REF_num >1] <- "2+"

###### ClinVar Attributes ######

clinVar_disease_names <- read.delim("data/ClinVar/clinVar_disease_names.txt")

#Manually retrieved OMIM "branchpoint" SNP entries
OMIM_branchpoint_SNPs_clinVar <- read.csv("data/ClinVar/OMIM_branchpoint_SNPs_clinVar.csv")
OMIM_branchpoint_SNPs_clinVar <- OMIM_branchpoint_SNPs_clinVar[!is.na(OMIM_branchpoint_SNPs_clinVar$ClinVar.Allele.ID),]

m <- match(OMIM_branchpoint_SNPs_clinVar$ClinVar.Allele.ID, ClinVar_variant_summary_hg38$X.AlleleID)
OMIM_branchpoint_SNPs_clinVar <- cbind(OMIM_branchpoint_SNPs_clinVar[which(!is.na(m)),], ClinVar_variant_summary_hg38[m[!is.na(m)],])

#some loss due to incomplete annotation in the clinVar variant summary file
OMIM_branchpoint_SNPs_clinVar <- OMIM_branchpoint_SNPs_clinVar[which(!duplicated(OMIM_branchpoint_SNPs_clinVar$dbSNP.ID)),]
m <- match(OMIM_branchpoint_SNPs_clinVar$new_id, clinVar_processed$id)

index <- vector()
for(i in seq_along(OMIM_branchpoint_SNPs_clinVar$new_id)){
  index[i] <- grep(gsub("-",".",gsub("[+]",".",gsub(":",".",OMIM_branchpoint_SNPs_clinVar$new_id[i]))), clinVar_processed$id)
}

OMIM_branchpoint_SNPs_clinVar <- cbind(OMIM_branchpoint_SNPs_clinVar, clinVar_processed[index,])
#we are not using 15562 and 15299 as their descriptions confirm that the affected branchpoint is outside the branchpoint window.

min_diff <- min((OMIM_branchpoint_SNPs_clinVar$max_prob_REF - OMIM_branchpoint_SNPs_clinVar$max_prob_ALT)[OMIM_branchpoint_SNPs_clinVar$dist_to_exon > 17 & OMIM_branchpoint_SNPs_clinVar$deleted_n>0])

clinVar_processed$clinVar_AlleleID <- matrix(unlist(str_split(str_sub(clinVar_processed$id,9,16), "[.]")), ncol=2, byrow = T)[,1]

phenotypes <-  ClinVar_variant_summary_hg38$PhenotypeIDS[match(clinVar_processed$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]

#keep those with a phenotype
keep <- grepl("MedGen:", phenotypes)
clinVar_processed <- clinVar_processed[keep,]
phenotypes <- phenotypes[keep]

MedGen <- list()
for(p in seq_along(phenotypes)){
  split_p <- unlist(str_split(phenotypes[p],";"))
  split_p <- unlist(str_split(split_p,","))
  
  split_p <- split_p[grep("MedGen", split_p)]
  split_p <- unlist(str_split(split_p,":"))
  split_p <- split_p[-grep("MedGen", split_p)]
  
  MedGen[[p]] <- split_p
}


clinVar_processed = data.frame(cbind(clinVar_processed, clinVar_disease_names[match(unlist(lapply(MedGen, `[[`, 1))
, clinVar_disease_names$ConceptID),]))

clinVar_processed$ClinicalSignificance = ClinVar_variant_summary_hg38$ClinicalSignificance[match(clinVar_processed$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]
clinVar_processed$Origin = ClinVar_variant_summary_hg38$Origin[match(clinVar_processed$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]
clinVar_processed$GeneSymbol = ClinVar_variant_summary_hg38$GeneSymbol[match(clinVar_processed$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]

clinVar_processed$dist_to_BP <- clinVar_processed$dist_to_BP_REF
clinVar_processed$dist_to_BP[clinVar_processed$BP_REF_num == 0] <- clinVar_processed$dist_to_BP_ALT[clinVar_processed$BP_REF_num == 0]
clinVar_processed$multi_BP[clinVar_processed$BP_REF_num == 0] <- "0"

clinVars_filtered <- clinVar_processed[grep("athog", clinVar_processed$ClinicalSignificance),]
clinVars_filtered <- clinVars_filtered[clinVars_filtered$deleted_n > 0 | clinVars_filtered$created_n >0,]

write.csv(clinVars_filtered, file="Tables/TableS2.csv", row.names=FALSE)
write.csv(OMIM_branchpoint_SNPs_clinVar, file="Tables/Table2.csv", row.names=FALSE)

save.image(file="data/diseaseVariants.Rdata")