################################################################################
#            Genome wide prediction of splicing branchpoints                   #
################################################################################

options(stringsAsFactors = F)

library(data.table)
library(ggplot2)
library(stringr)
library(plyr)
library(parallel)
library(entropy)
library(DEXSeq)
library(branchpointer)

#Variables
cutoff = 0.5

#from process_predictions.R
load("data/genome_predictions/gencode.v12.RData")

exons <- readExonAnnotation("data/genome_predictions/gencode.v12.annotation.gtf")
exons$exon_id <- with(exons, paste0(exon_id,"_",exon_number))
G12_all <- gencode_v12

rm(gencode_v12_seqs, gencode_v12)
###### Make intron centric annotation and annotate structural features ######

#make intron-centric annotation
introns <- unique(G12_all$id)
branchpoints_introns <- data.frame(introns)
branchpoints_introns$status <- "predicted"

#testing/training introns all contained at least one HC
BPs <- which(G12_all$in_testtrain == "HC")
branchpoints_introns$status[which(!is.na(match(branchpoints_introns$introns,
                                               G12_all$id[BPs])))] <- "test/train"

#number of known (HC) BPs
branchpoints_introns$known_BPs <- 0
known_BPs <- as.data.frame(table(G12_all$id[BPs]))
m <- match(branchpoints_introns$introns, known_BPs$Var1)
branchpoints_introns$known_BPs[which(!is.na(m))] <- known_BPs$Freq[m[!is.na(m)]]

#number of predicted BPs
branchpoints_introns$predicted_BPs <- 0
predicted_BPs <- as.data.frame(table(G12_all$id[G12_all$branchpoint_prob >= cutoff]))
m <- match(branchpoints_introns$introns, predicted_BPs$Var1)
branchpoints_introns$predicted_BPs[which(!is.na(m))] <-
  predicted_BPs$Freq[m[which(!is.na(m))]]

#number of annotated BPs
#If a BP is in HC dataset, use that, else use predicted scores
branchpoints_introns$annotated_BPs <- 0
annotated_BPs <- as.data.frame(table(G12_all$id[G12_all$branchpoint_prob >= cutoff |
                                          G12_all$in_testtrain == "HC"]))
m <- match(branchpoints_introns$introns, annotated_BPs$Var1)
branchpoints_introns$annotated_BPs[which(!is.na(m))] <-
  annotated_BPs$Freq[m[which(!is.na(m))]]

branchpoints_introns$status[branchpoints_introns$annotated_BPs == 0] <- "unknown"

#get intron size and biotype
m <- match(branchpoints_introns$introns, G12_all$id)
branchpoints_introns$intron_size <- G12_all$to_3prime[m] + G12_all$to_5prime[m]
n <- match(G12_all$id, exons$exon_id)
G12_all$gene_biotype <- exons$gene_type[n]
branchpoints_introns$gene_biotype <- G12_all$gene_biotype[m]

#fix biotype to be broader
biotype <- branchpoints_introns$gene_biotype
lncRNAs <- which(
  biotype == "3prime_overlapping_ncrna" |
    biotype == "antisense" |
    biotype == "bidirectional_promoter_lncrna" |
    biotype == "lincRNA" |
    biotype == "processed_transcript" |
    biotype == "sense_intronic" |
    biotype == "sense_overlapping"
)
branchpoints_introns$gene_biotype_broad <- biotype
branchpoints_introns$gene_biotype_broad[lncRNAs] <- "lncRNA"
pseudo <- grep("pseudo", biotype)
branchpoints_introns$gene_biotype_broad[pseudo] <- "pseudogene"
other <- which(!(branchpoints_introns$gene_biotype_broad %in%
    c("lncRNA", "protein_coding","pseudogene")))
branchpoints_introns$gene_biotype_broad[other] <- "other"

#convert number of branchpoints to a factor variable
branchpoints_introns$annotated_BPs_factor <- branchpoints_introns$annotated_BPs
branchpoints_introns$annotated_BPs_factor[
  which(branchpoints_introns$annotated_BPs > 3)] <- "4+"
branchpoints_introns$predicted_BPs_factor <- branchpoints_introns$predicted_BPs
branchpoints_introns$predicted_BPs_factor[
  which(branchpoints_introns$predicted_BPs > 3)] <- "4+"

#number of alternative 5'exons by gencode annotation

exons_pos <-
  arrange(exons, transcript_id, chromosome, start)
exons_pos <-
  exons_pos[exons_pos$strand == "+",]

exons_neg <-
  arrange(exons, transcript_id, chromosome, plyr::desc(end))
exons_neg <-
  exons_neg[exons_neg$strand == "-",]

exons_pos_introns <-
  data.frame(start = exons_pos$end[-length(exons_pos$end)],
    end = exons_pos$start[-1],
    start_transcript_id = exons_pos$transcript_id[-length(exons_pos$end)],
    end_transcript_id = exons_pos$transcript_id[-1])
exons_pos_introns <-
  exons_pos_introns[(exons_pos_introns$start_transcript_id ==
      exons_pos_introns$end_transcript_id),]

exons_neg_introns <-
  data.frame(
    start = exons_neg$end[-1],
    end = exons_neg$start[-length(exons_neg$end)],
    start_transcript_id = exons_neg$transcript_id[-1],
    end_transcript_id = exons_neg$transcript_id[-length(exons_neg$end)]
  )
exons_neg_introns <-
  exons_neg_introns[(exons_neg_introns$start_transcript_id ==
      exons_neg_introns$end_transcript_id),]

m <- match(exons_pos_introns$start_transcript_id, 
    exons_pos$transcript_id)
exons_pos_introns <- cbind(exons_pos_introns, 
                                        exons_pos[m,c(1,5,6,7,9)])
m <- match(exons_neg_introns$start_transcript_id,
           exons_neg$transcript_id)
exons_neg_introns <- cbind(exons_neg_introns,
                                        exons_neg[m,c(1,5,6,7,9)])

rm(exons_neg,exons_pos)

exons_pos_introns$chr_end_gencode <-
  with(exons_pos_introns, 
       paste(chromosome, end,sep = "_"))
exons_pos_introns$chr_startend_gencode <-
  with(exons_pos_introns, 
       paste(chromosome,start, end,sep = "_"))
exons_pos_introns <- exons_pos_introns[
  !duplicated(exons_pos_introns$chr_startend_gencode),]
exon5count_pos <-
  as.data.frame(table(exons_pos_introns$chr_end_gencode))

exons_neg_introns$chr_start_gencode <-
  with(exons_neg_introns, 
       paste(chromosome, start,sep = "_"))
exons_neg_introns$chr_startend_gencode <-
  with(exons_neg_introns, 
       paste(chromosome,start, end,sep = "_"))
exons_neg_introns <-
  exons_neg_introns[
    !duplicated(exons_neg_introns$chr_startend_gencode),]
exon5count_neg <-
  as.data.frame(table(exons_neg_introns$chr_start_gencode))

G12_all$exon2_start <- exons$start[match(G12_all$exon_3prime, exons$exon_id)]
G12_all$exon2_end <- exons$end[match(G12_all$exon_3prime, exons$exon_id)]

chr_start_BP <- with(G12_all, paste(chromosome,exon2_start,sep = "_"))
chr_start_BP <- gsub(" ","", chr_start_BP)
chr_end_BP <- with(G12_all, paste(chromosome,exon2_end,sep = "_"))
chr_end_BP <- gsub(" ","", chr_end_BP)

m_p <- match(chr_start_BP, exon5count_pos$Var1)
m_n <- match(chr_end_BP, exon5count_neg$Var1)

gencode_alternative5 <- exon5count_pos$Freq[m_p]
gencode_alternative5[G12_all$strand == "-"] <-
  (exon5count_neg$Freq[m_n])[G12_all$strand == "-"]
G12_all$gencode_alternative5 <- gencode_alternative5

#match up site-wise and intron-wise attributes
m <- match(branchpoints_introns$introns, G12_all$id)
branchpoints_introns$exon1_contains <- G12_all$exon1_contains[m]
branchpoints_introns$exon2_contains <- G12_all$exon2_contains[m]
branchpoints_introns$gencode_alternative5 <-
  G12_all$gencode_alternative5[m]

m <- match(G12_all$id,branchpoints_introns$introns)
G12_all$exon1_contains <- branchpoints_introns$exon1_contains[m]
G12_all$exon2_contains <- branchpoints_introns$exon2_contains[m]
G12_all$known_BPs <- branchpoints_introns$known_BPs[m]
G12_all$predicted_BPs <- branchpoints_introns$predicted_BPs[m]
G12_all$annotated_BPs <- branchpoints_introns$annotated_BPs[m]
G12_all$predicted_BPs_factor <-
  branchpoints_introns$predicted_BPs_factor[m]
G12_all$annotated_BPs_factor <-
  branchpoints_introns$annotated_BPs_factor[m]
G12_all$intron_size <- branchpoints_introns$intron_size[m]

rm(annotated_BPs, known_BPs, predicted_BPs, G12_pc, gencode_with_cdsutr,
  gencode_with_cdsutr_chr, gencode_with_cds,gencode_with_utr)

rm(annotation,biotype,BPs,c,cds_ind,cds_locations_n,cds_locations_p,chroms,
  introns,lncRNAs,m,next_entry_same,other,prev_entry_same,
  pseudo,u3_ind,u5_ind,utr3_locations_n,utr3_locations_p,
  utr5_locations_n,utr5_locations_p)
rm(exon5count_pos,exon5count_neg,exons_pos_introns,
  exons_neg_introns,chr_end_BP,chr_start_BP,
  m_n,m_p,gencode_alternative5)

###### Information on "probability scores" and U2 binding energy ######

#get fivemer frequencies & mean probability scores
fivemers <- with(G12_all, paste(seq_neg3,seq_neg2,seq_neg1,
                                seq_pos0,seq_pos1, sep = ""))
G12_all$fivemers <- fivemers
fivemer_summary <- as.data.frame(table(G12_all$fivemers))

fivemer_summary <- cbind(fivemer_summary,
                         aggregate(branchpoint_prob ~ fivemers, 
                                   data = G12_all, min)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers,data = G12_all, quantile,0.25)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G12_all, median)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G12_all, mean)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G12_all, quantile,0.75)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G12_all, max)$branchpoint_prob)
colnames(fivemer_summary) <- c("fivemer", "count","min_branchpoint_prob","Q1_branchpoint_prob",
                               "median_branchpoint_prob","mean_branchpoint_prob","Q3_branchpoint_prob","max_branchpoint_prob")

#plot 5mer freq by BP/Neg ratios
fivemer_freq_BP <-
  as.data.frame(table(G12_all$fivemers[G12_all$in_testtrain == "HC"]))
m <- match(fivemer_summary$fivemer,fivemer_freq_BP$Var1)
fivemer_summary$num_BP <- fivemer_freq_BP$Freq[m]
fivemer_summary$num_BP[is.na(fivemer_summary$num_BP)] <- 0
fivemer_freq_N <-
  as.data.frame(table(G12_all$fivemers[G12_all$in_testtrain == "NEG"]))
m <- match(fivemer_summary$fivemer,fivemer_freq_N$Var1)
fivemer_summary$num_NEG <- fivemer_freq_N$Freq[m]
fivemer_summary$num_NEG[is.na(fivemer_summary$num_NEG)] <- 0

fivemer_summary$percent_BP <- fivemer_summary$num_BP / 
  (fivemer_summary$num_NEG + fivemer_summary$num_BP)

rm(fivemers, m, fivemer_freq_N,fivemer_freq_BP)

#### U2 binding energy###
seqs <- c("A","T","C","G")
eightmers <- apply(expand.grid(seqs,seqs,seqs,seqs,seqs,seqs,seqs,seqs),
                   1, paste, collapse = "")
U2_binding_df <- data.frame(U2 = "GTGTAGTA", eightmers)

U2_to_eightmers_output <-
  read.csv("data/U2_to_eightmers_output.txt", header = FALSE, sep = ":")

m1 <- matrix(unlist(str_split((U2_to_eightmers_output$V2), "[(]")),
             ncol = 2, byrow = T)[,2]
m2 <- as.numeric(str_sub(m1, 2,5))

U2_binding_df <- cbind(U2_binding_df, energy = m2)
eightmer_df <- expand.grid(seqs,seqs,seqs,seqs,seqs,seqs,seqs,seqs)
colnames(eightmer_df) <- c("seq_neg5", "seq_neg4","seq_neg3","seq_neg2",
                           "seq_neg1","seq_pos1","seq_pos2","seq_pos3")
U2_binding_df <- cbind(U2_binding_df, eightmer_df)

u2_eightmers <- with(G12_all, paste(seq_neg5,seq_neg4,seq_neg3,
                                    seq_neg2,seq_neg1,seq_pos1,
                                    seq_pos2,seq_pos3, sep =""))
m <- match(u2_eightmers, U2_binding_df$eightmers)
G12_all$U2_binding_energy <- U2_binding_df$energy[m]
G12_all$U2_eightmers <- u2_eightmers

U2_summary <- aggregate(U2_binding_energy ~ fivemers, 
                        data = G12_all, min)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G12_all, quantile,0.25)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G12_all, median)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G12_all, mean)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G12_all, quantile,0.75)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G12_all, max)$U2_binding_energy)
colnames(U2_summary) <- c("fivemer", "min_U2","Q1_U2","median_U2",
                          "mean_U2","Q3_U2","max_U2")
fivemer_summary <- cbind(fivemer_summary, 
                         U2_summary[match(fivemer_summary$fivemer, 
                                          U2_summary$fivemer),-1])

fivemer_summary$BP_nt <- str_sub(fivemer_summary$fivemer, 4,4)
fivemer_summary <- fivemer_summary[fivemer_summary$BP_nt!="N",]


rm(eightmer_df,eightmers,m,m1,m2,seqs,U2_binding_df,
   u2_eightmers,U2_summary,U2_to_eightmers_output)

save(G12_all, branchpoints_introns, fivemer_summary,
     file = "data/Figure_files_G12.Rdata")

