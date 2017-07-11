################################################################################
#                   Human genetic variation at branchpoints                    #
################################################################################

options(stringsAsFactors = F)

library(biomaRt)
library(stringr)
library(branchpointer)
library(data.table)
library(GenomicRanges)

#load clinVar annotation
#downloaded 2017/06/02
ClinVar_variant_summary <- read.delim("data/ClinVar/variant_summary_clinVar_downloaded20170602.txt")
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
write.csv(snp_info_clinVar, file="data/ClinVar/clinvar_info.csv", row.names = FALSE)

###### run branchpointer ######

exons <- gtfToExons("data/gencode.v26.annotation.gtf")
 
query <- readQueryFile("data/ClinVar/clinvar_info.csv",queryType = "SNP",
                       exons, filter=TRUE)
 
reps <- ceiling(length(query) / 1000)

preds <- list()
for(i in 1:reps){
index_end <- i*1000
index_start <- index_end - 999  
index_end <- min(index_end, length(query))

query_pred <- predictBranchpoints(query[index_start:index_end],queryType = "SNP",
                                        BSgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

preds[[i]] <- query_pred

write.table(mcols(query_pred), file=paste0("data/ClinVar/clinvar_attributes_",i,".txt"), sep="\t",
            row.names=FALSE, quote=FALSE)
message(i)
}

clinvar_predictions <- do.call("c", preds)

clinvar_summary <- branchpointer::predictionsToSummary(query, clinvar_predictions, 
                                                     probabilityCutoff = 0.52, 
                                                     probabilityChange = 0.15)

clinvar_summary$allele_ref_feat_strand <- clinvar_summary$ref_allele
clinvar_summary$allele_ref_feat_strand[clinvar_summary$strand=="-" & clinvar_summary$ref_allele=="A"] <- "T"
clinvar_summary$allele_ref_feat_strand[clinvar_summary$strand=="-" & clinvar_summary$ref_allele=="T"] <- "A"
clinvar_summary$allele_ref_feat_strand[clinvar_summary$strand=="-" & clinvar_summary$ref_allele=="G"] <- "C"
clinvar_summary$allele_ref_feat_strand[clinvar_summary$strand=="-" & clinvar_summary$ref_allele=="C"] <- "G"

clinvar_summary$multi_BP <- 0
clinvar_summary$multi_BP[clinvar_summary$BP_num_REF == 1] <- 1
clinvar_summary$multi_BP[clinvar_summary$BP_num_REF > 1] <- "2+"

###### ClinVar Branchpoints Summary ######

clinVar_disease_names <- read.delim("data/ClinVar/clinVar_disease_names.txt")

#Manually retrieved OMIM "branchpoint" SNP entries
OMIM_branchpoint_SNPs_clinVar <- read.csv("data/ClinVar/OMIM_branchpoint_SNPs_clinVar.csv")
OMIM_branchpoint_SNPs_clinVar <- OMIM_branchpoint_SNPs_clinVar[!is.na(OMIM_branchpoint_SNPs_clinVar$ClinVar.Allele.ID),]

m <- match(OMIM_branchpoint_SNPs_clinVar$ClinVar.Allele.ID, ClinVar_variant_summary_hg38$X.AlleleID)
OMIM_branchpoint_SNPs_clinVar <- cbind(OMIM_branchpoint_SNPs_clinVar[which(!is.na(m)),], ClinVar_variant_summary_hg38[m[!is.na(m)],])

#some loss due to incomplete annotation in the clinVar variant summary file
OMIM_branchpoint_SNPs_clinVar <- OMIM_branchpoint_SNPs_clinVar[which(!duplicated(OMIM_branchpoint_SNPs_clinVar$ClinVar.Allele.ID)),]
m <- match(OMIM_branchpoint_SNPs_clinVar$new_id, clinvar_summary$id)

index <- vector()
for(i in seq_along(OMIM_branchpoint_SNPs_clinVar$new_id)){
  index[i] <- grep(gsub("-",".",gsub("[+]",".",gsub(":",".",OMIM_branchpoint_SNPs_clinVar$new_id[i]))), clinvar_summary$id)
}

OMIM_branchpoint_SNPs_clinVar <- cbind(OMIM_branchpoint_SNPs_clinVar, clinvar_summary[index])
#we are not using 15562 and 15299 as their descriptions confirm that the affected branchpoint is outside the branchpoint window.
OMIM_branchpoint_SNPs_clinVar <- OMIM_branchpoint_SNPs_clinVar[OMIM_branchpoint_SNPs_clinVar$to_3prime %in% 18:44,]

clinvar_summary$clinVar_AlleleID <- matrix(unlist(str_split(str_sub(clinvar_summary$id,9,16), "[.]")), ncol=2, byrow = T)[,1]
phenotypes <-  ClinVar_variant_summary_hg38$PhenotypeIDS[match(clinvar_summary$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]

#keep those with a phenotype
keep <- which(grepl("MedGen:", phenotypes))
clinvar_summary <- clinvar_summary[keep]
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


clinvar_summary$Disease_names <- clinVar_disease_names$X.DiseaseName[match(unlist(lapply(MedGen, `[[`, 1)), clinVar_disease_names$ConceptID)]
clinvar_summary$ClinicalSignificance = ClinVar_variant_summary_hg38$ClinicalSignificance[match(clinvar_summary$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]
clinvar_summary$Origin = ClinVar_variant_summary_hg38$Origin[match(clinvar_summary$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]
clinvar_summary$GeneSymbol = ClinVar_variant_summary_hg38$GeneSymbol[match(clinvar_summary$clinVar_AlleleID, ClinVar_variant_summary_hg38$X.AlleleID)]

clinvar_summary$dist_to_BP <- clinvar_summary$dist_to_BP_REF
clinvar_summary$dist_to_BP[clinvar_summary$BP_REF_num == 0] <- clinvar_summary$dist_to_BP_ALT[clinvar_summary$BP_REF_num == 0]

clinvar_summary$snp_in <- NA
clinvar_summary$snp_in[clinvar_summary$dist_to_BP <= -1] <- "3'SS to BP"
clinvar_summary$snp_in[clinvar_summary$dist_to_BP >= 3] <- "BP to 5'SS"
clinvar_summary$snp_in[clinvar_summary$dist_to_BP == 0 | clinvar_summary$dist_to_BP == 2] <- "BP"
clinvar_summary$snp_in[clinvar_summary$to_3prime < 3] <- "3'SS"

gtf <- rtracklayer::import("data/gencode.v26.annotation.gtf")
gtf <- gtf[gtf$type=="exon"]
o <- GenomicRanges::findOverlaps(clinvar_summary, gtf)
rm <- unique(o@from)
clinvar_summary_filtered <- clinvar_summary[-rm]

clinvar_summary_filtered <- as.data.frame(clinvar_summary_filtered)[(clinvar_summary_filtered$created_n > 0 | clinvar_summary_filtered$deleted_n > 0) &
                                                       clinvar_summary_filtered$to_3prime > 2,]
keep <- which((clinvar_summary_filtered$max_prob_REF != clinvar_summary_filtered$max_prob_ALT) | 
                (clinvar_summary_filtered$max_U2_REF != clinvar_summary_filtered$max_U2_ALT) | 
                is.na(clinvar_summary_filtered$max_prob_ALT) | is.na(clinvar_summary_filtered$max_prob_REF))
clinvar_summary_filtered <- clinvar_summary_filtered[keep,]


TableS3 <- as.data.frame(clinvar_summary_filtered)

write.csv(TableS3, file="Tables/TableS3.csv", row.names=FALSE)

m <- match(clinvar_predictions$id,OMIM_branchpoint_SNPs_clinVar$id)
n <- which(!is.na(m))
OMIM_preds <- clinvar_predictions[n]

ref_BP_pos <- aggregate(to_3prime_point ~ id,
          data=as.data.frame(OMIM_preds[OMIM_preds$branchpoint_prob >= 0.52 & OMIM_preds$status=="REF"]),
          function(x) paste(as.character(x),collapse=","))

alt_BP_pos <- aggregate(to_3prime_point ~ id,
          data=as.data.frame(OMIM_preds[OMIM_preds$branchpoint_prob >= 0.52 & OMIM_preds$status=="ALT"]),
          function(x) paste(as.character(x),collapse=","))

OMIM_branchpoint_SNPs_clinVar$BPs_ref <- NA
OMIM_branchpoint_SNPs_clinVar$BPs_ref[match(ref_BP_pos$id, OMIM_branchpoint_SNPs_clinVar$id)] <- ref_BP_pos$to_3prime_point

OMIM_branchpoint_SNPs_clinVar$BPs_alt <- NA
OMIM_branchpoint_SNPs_clinVar$BPs_alt[match(alt_BP_pos$id, OMIM_branchpoint_SNPs_clinVar$id)] <- alt_BP_pos$to_3prime_point

OMIM_branchpoint_SNPs_clinVar$Allele <- paste0(OMIM_branchpoint_SNPs_clinVar$ReferenceAllele, "/", OMIM_branchpoint_SNPs_clinVar$AlternateAllele)

Table2 <- OMIM_branchpoint_SNPs_clinVar[,c('GeneSymbol','ClinVar.Allele.ID',"Allele",'BPs_ref',"BPs_alt")]
colnames(Table2) <- c("Gene","ClinVar ID","Allele","BPs (ref.)","BPs (alt.)")
write.csv(Table2, file="Tables/Table2.csv", row.names=FALSE)

query_clinvar <- query

save(clinVars_filtered, OMIM_branchpoint_SNPs_clinVar, clinvar_summary, 
     query_clinvar,clinvar_predictions, ClinVar_variant_summary_hg38,
     snp_info_clinVar, file="data/diseaseVariants.Rdata")
