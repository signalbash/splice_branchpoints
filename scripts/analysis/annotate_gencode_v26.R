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
cutoff = 0.48

#from process_predictions.R
load("data/gencode_v26.RData")
gtf <- gtfToExons("data/gencode.v26.annotation.gtf")

G26_all <- gencode_v26
rm(gencode_v26)

G26_all$exon3prime_start <- start(ranges(gtf[match(G26_all$exon_id, gtf$exon_id)]))
G26_all$exon3prime_end <- end(ranges(gtf[match(G26_all$exon_id, gtf$exon_id)]))
G26_all$exon5prime_start <- start(ranges(gtf[match(G26_all$exon_5prime, gtf$exon_id)]))
G26_all$exon5prime_end <- end(ranges(gtf[match(G26_all$exon_5prime, gtf$exon_id)]))

###### Make intron centric annotation and annotate structural features ######
#make intron-centric annotation
introns <- unique(G26_all$exon_id)
branchpoints_introns <- data.frame(introns)
branchpoints_introns$status <- "predicted"

#testing/training introns all contained at least one HC
BPs <- which(G26_all$in_testtrain == "HC")
branchpoints_introns$status[which(!is.na(match(branchpoints_introns$introns,
                                               G26_all$exon_id[BPs])))] <- "test/train"

#number of known (HC) BPs
branchpoints_introns$known_BPs <- 0
known_BPs <- as.data.frame(table(G26_all$exon_id[BPs]))
m <- match(branchpoints_introns$introns, known_BPs$Var1)
branchpoints_introns$known_BPs[which(!is.na(m))] <- known_BPs$Freq[m[!is.na(m)]]

#number of predicted BPs
branchpoints_introns$predicted_BPs <- 0
predicted_BPs <- as.data.frame(table(G26_all$exon_id[G26_all$branchpoint_prob >= cutoff]))
m <- match(branchpoints_introns$introns, predicted_BPs$Var1)
branchpoints_introns$predicted_BPs[which(!is.na(m))] <-
  predicted_BPs$Freq[m[which(!is.na(m))]]

#number of annotated BPs
#If a BP is in HC dataset, use that, else use predicted scores
branchpoints_introns$annotated_BPs <- 0
annotated_BPs <- as.data.frame(table(G26_all$exon_id[G26_all$branchpoint_prob >= cutoff |
                                          G26_all$in_testtrain == "HC"]))
m <- match(branchpoints_introns$introns, annotated_BPs$Var1)
branchpoints_introns$annotated_BPs[which(!is.na(m))] <-
  annotated_BPs$Freq[m[which(!is.na(m))]]

branchpoints_introns$status[branchpoints_introns$annotated_BPs == 0] <- "unknown"

#get intron size and biotype
G26_all$intron_size <- G26_all$to_3prime + G26_all$to_5prime
m <- match(branchpoints_introns$introns, G26_all$exon_id)
branchpoints_introns$intron_size <- G26_all$intron_size[m]
n <- match(G26_all$exon_id, gtf$exon_id)
G26_all$gene_biotype <- gtf$gene_type[n]
branchpoints_introns$gene_biotype <- G26_all$gene_biotype[m]

#fix biotype to be broader
biotype <- branchpoints_introns$gene_biotype
lncRNAs <- which(
  biotype %in% c("3prime_overlapping_ncrna","antisense",
                 "bidirectional_promoter_lncrna","lincRNA",
                 "processed_transcript","sense_intronic","sense_overlapping"))

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
  which(branchpoints_introns$annotated_BPs > 1)] <- "2+"
branchpoints_introns$predicted_BPs_factor <- branchpoints_introns$predicted_BPs
branchpoints_introns$predicted_BPs_factor[
  which(branchpoints_introns$predicted_BPs > 1)] <- "2+"

#match up site-wise and intron-wise attributes

m <- match(G26_all$exon_id,branchpoints_introns$introns)
G26_all$known_BPs <- branchpoints_introns$known_BPs[m]
G26_all$predicted_BPs <- branchpoints_introns$predicted_BPs[m]
G26_all$annotated_BPs <- branchpoints_introns$annotated_BPs[m]
G26_all$predicted_BPs_factor <-
  branchpoints_introns$predicted_BPs_factor[m]
G26_all$annotated_BPs_factor <-
  branchpoints_introns$annotated_BPs_factor[m]

###### Information on "probability scores" and U2 binding energy ######

#get fivemer frequencies & mean probability scores
G26_all$fivemers <- str_sub(G26_all$seq_motif, start=3, end=7)

fivemer_summary <- as.data.frame(table(G26_all$fivemers))

fivemer_summary <- cbind(fivemer_summary,
                         aggregate(branchpoint_prob ~ fivemers, 
                                   data = G26_all, min)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers,data = G26_all, quantile,0.25)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G26_all, median)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G26_all, mean)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G26_all, quantile,0.75)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G26_all, max)$branchpoint_prob)
colnames(fivemer_summary) <- c("fivemer", "count","min_branchpoint_prob","Q1_branchpoint_prob",
                               "median_branchpoint_prob","mean_branchpoint_prob","Q3_branchpoint_prob","max_branchpoint_prob")

#plot 5mer freq by BP/Neg ratios
fivemer_freq_BP <-
  as.data.frame(table(G26_all$fivemers[G26_all$in_testtrain == "HC"]))
m <- match(fivemer_summary$fivemer,fivemer_freq_BP$Var1)
fivemer_summary$num_BP <- fivemer_freq_BP$Freq[m]
fivemer_summary$num_BP[is.na(fivemer_summary$num_BP)] <- 0
fivemer_freq_N <-
  as.data.frame(table(G26_all$fivemers[G26_all$in_testtrain == "NEG"]))
m <- match(fivemer_summary$fivemer,fivemer_freq_N$Var1)
fivemer_summary$num_NEG <- fivemer_freq_N$Freq[m]
fivemer_summary$num_NEG[is.na(fivemer_summary$num_NEG)] <- 0

fivemer_summary$percent_BP <- fivemer_summary$num_BP / 
  (fivemer_summary$num_NEG + fivemer_summary$num_BP)

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

u2_eightmers <- str_sub(G26_all$seq_motif, start=1, end=8)

m <- match(u2_eightmers, U2_binding_df$eightmers)
#G26_all$U2_binding_energy <- U2_binding_df$energy[m]
G26_all$U2_eightmers <- u2_eightmers

U2_summary <- aggregate(U2_binding_energy ~ fivemers, 
                        data = G26_all, min)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G26_all, quantile,0.25)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G26_all, median)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G26_all, mean)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G26_all, quantile,0.75)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G26_all, max)$U2_binding_energy)
colnames(U2_summary) <- c("fivemer", "min_U2","Q1_U2","median_U2",
                          "mean_U2","Q3_U2","max_U2")
fivemer_summary <- cbind(fivemer_summary, 
                         U2_summary[match(fivemer_summary$fivemer, 
                                          U2_summary$fivemer),-1])

fivemer_summary$BP_nt <- str_sub(fivemer_summary$fivemer, 4,4)
fivemer_summary <- fivemer_summary[fivemer_summary$BP_nt!="N",]

###### DEXSeq exon expression ######

#DEXseq specific gtf
gencode.dexseq <- rtracklayer::import("data/genome_annotations/gencode.v26.dexseq.gtf")
gencode.dexseq <- gencode.dexseq[gencode.dexseq$type=="exonic_part"]
gencode.dexseq$exon_names <- with(mcols(gencode.dexseq), paste(gene_id, exonic_part_number, sep=":"))

#combine exonic parts into whole exons
gencode.dexseq$exon_group <- gencode.dexseq$exon_names
same_exon <- which((end(gencode.dexseq)+1)[-length(gencode.dexseq)] == start(gencode.dexseq[-1]) &
                     (gencode.dexseq$gene_id[-length(gencode.dexseq)] == gencode.dexseq$gene_id[-1]) &
                     (gencode.dexseq$exon_group[-length(gencode.dexseq)] != gencode.dexseq$exon_group[-1]))

while (length(same_exon) > 0) {
  gencode.dexseq$exon_group[same_exon + 1] <-
    gencode.dexseq$exon_group[same_exon]
  same_exon <- which((end(gencode.dexseq)+1)[-length(gencode.dexseq)] == start(gencode.dexseq[-1]) &
                       (gencode.dexseq$gene_id[-length(gencode.dexseq)] == gencode.dexseq$gene_id[-1]) &
                       (gencode.dexseq$exon_group[-length(gencode.dexseq)] != gencode.dexseq$exon_group[-1]))
  message(length(same_exon))
}

#DEXSeq expression
rep1.dexseq <-
  read.delim(
    "data/ENCODE_K562/dexseq/K562_rep1.dexseq.txt", header =
      FALSE
  )
rep2.dexseq <-
  read.delim(
    "data/ENCODE_K562/dexseq/K562_rep2.dexseq.txt", header =
      FALSE
  )

m <- match(rep1.dexseq$V1, rep2.dexseq$V1)
dexseq <- cbind(rep1.dexseq, rep2.dexseq[m,2])
colnames(dexseq) <- c("exon_name", "rep1","rep2")
rownames(dexseq) <- dexseq$exon_name
dexseq$exon_name <- NULL
size_factors <- estimateSizeFactorsForMatrix(dexseq)
dex_seq_norm <- dexseq / size_factors

exon_groups <- unique(gencode.dexseq$exon_group)

gencode.dexseq$mean_DEX_count <- NA
gencode.dexseq$dex_rep1 <- dex_seq_norm$rep1[match(gencode.dexseq$exon_names,rownames(dex_seq_norm))]
gencode.dexseq$dex_rep2 <- dex_seq_norm$rep2[match(gencode.dexseq$exon_names,rownames(dex_seq_norm))]

#get counts per exons groups, not exon parts
rep1_total_counts <- aggregate(dex_rep1 ~ exon_group, mcols(gencode.dexseq), sum)
rep2_total_counts <- aggregate(dex_rep2 ~ exon_group, mcols(gencode.dexseq), sum)
total_counts <- cbind(rep1_total_counts,rep2_total_counts[,-1])
total_counts$mean <- rowMeans(total_counts[,-1])
m <- match(gencode.dexseq$exon_group, total_counts$exon_group)
gencode.dexseq$group_mean_count <- total_counts$mean[m]
gencode.dexseq$mean_DEX_count <- rowMeans(as.data.frame(mcols(gencode.dexseq)[,c("dex_rep1", "dex_rep2")]))

# copy over to bp summary
# match exons
G26_all.names <- with(G26_all, paste0(chromosome, "_",exon3prime_start,"_",strand))
G26_all.names[G26_all$strand == "-"] <- with(G26_all[G26_all$strand == "-",], paste0(chromosome, "_",exon3prime_end,"_",strand))
dex.names <- paste0(seqnames(gencode.dexseq), "_", start(gencode.dexseq), "_", strand(gencode.dexseq))
negInd <- which(strand(gencode.dexseq) == "-")
dex.names[negInd] <- paste0(seqnames(gencode.dexseq)[negInd], "_", end(gencode.dexseq)[negInd], "_", strand(gencode.dexseq)[negInd])

m <- match(G26_all.names, dex.names)

G26_all$exon_count <- gencode.dexseq$mean_DEX_count[m]
G26_all$exon_group_count <- gencode.dexseq$group_mean_count[m]
G26_all$dex_exon_name <- gencode.dexseq$exon_names[m]
G26_all$dex_exon_group <- gencode.dexseq$exon_group[m]
G26_all$exon_group_log10count <- log10(G26_all$exon_group_count + 0.1)

#add to intron annotation
m <- match(branchpoints_introns$introns, G26_all$exon_id)
branchpoints_introns$dex_exon_group <- G26_all$dex_exon_group[m]
branchpoints_introns$exon_group_count <- G26_all$exon_group_count[m]
branchpoints_introns$exon_group_log10count <- G26_all$exon_group_log10count[m]

##### Conservation ######
files <- list.files("data/conservation", pattern="conservation_df_", full.names = TRUE)

for(f in seq_along(files)){
  
  conservation_df_orginal <- fread(files[f], data.table = FALSE)
  
  if(exists("conservation_df_orginal_all")){
    conservation_df_orginal_all <- rbind(conservation_df_orginal_all,conservation_df_orginal)
  }else{
    conservation_df_orginal_all <- conservation_df_orginal
  }
}

conservation_df <- conservation_df_orginal_all
conservation_df$exon_id <- unlist(lapply(str_split(conservation_df$id, "_"), "[[", 1))
m <- match(conservation_df$exon_id, branchpoints_introns$introns)
conservation_df$status <- branchpoints_introns$status[m]
conservation_df$annotated_BPs <- branchpoints_introns$annotated_BPs_factor[m]
conservation_df$gene_biotype_broad <- branchpoints_introns$gene_biotype_broad[m]

melted <- melt(conservation_df[,c(1:12,21:23)], id.vars=c("id","status","annotated_BPs","gene_biotype_broad"))
melted$variable <- gsub("cons_neg", "-", melted$variable)
melted$variable <- gsub("cons_pos", "", melted$variable)
melted$variable <- as.numeric(melted$variable)

conservation_summary <- rbind(
  cbind(set="known",multiplicity="single_BP",
        as.data.frame(aggregate(value ~ variable, 
                                data=melted[melted$annotated_BPs == 1 & 
                                              melted$gene_biotype_broad=="protein_coding" &
                                              melted$status == "test/train",], median))),
  cbind(set="predicted",multiplicity="single_BP",
        as.data.frame(aggregate(value ~ variable, 
                                data=melted[melted$annotated_BPs == 1 & 
                                              melted$gene_biotype_broad=="protein_coding" & 
                                              melted$status == "predicted",], median))),
  cbind(set="known",multiplicity="multi_BP",
        as.data.frame(aggregate(value ~ variable, 
                                data=melted[melted$annotated_BPs == "2+" & 
                                              melted$gene_biotype_broad=="protein_coding" & 
                                              melted$status == "test/train",], median))),
  cbind(set="predicted",multiplicity="multi_BP",
        as.data.frame(aggregate(value ~ variable, 
                                data=melted[melted$annotated_BPs == "2+" & 
                                              melted$gene_biotype_broad=="protein_coding" & 
                                              melted$status == "predicted",], median))))

conservation_summary$se_lower <- rbind(as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == 1 & 
                                                                             melted$gene_biotype_broad=="protein_coding" &
                                                                             melted$status == "test/train",], function(x) (median(x) - sd(x)/sqrt(length(x))))),
                                       as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == 1 & 
                                                                             melted$gene_biotype_broad=="protein_coding" & 
                                                                             melted$status == "predicted",], function(x) (median(x) - sd(x)/sqrt(length(x))))),
                                       as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == "2+" & 
                                                                             melted$gene_biotype_broad=="protein_coding" & 
                                                                             melted$status == "test/train",], function(x) (median(x) - sd(x)/sqrt(length(x))))),
                                       as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == "2+" & 
                                                                             melted$gene_biotype_broad=="protein_coding" & 
                                                                             melted$status == "predicted",], function(x) (median(x) - sd(x)/sqrt(length(x))))))[,2]
conservation_summary$se_upper <- rbind(as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == 1 & 
                                                                             melted$gene_biotype_broad=="protein_coding" &
                                                                             melted$status == "test/train",], function(x) (median(x) + sd(x)/sqrt(length(x))))),
                                       as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == 1 & 
                                                                             melted$gene_biotype_broad=="protein_coding" & 
                                                                             melted$status == "predicted",], function(x) (median(x) + sd(x)/sqrt(length(x))))),
                                       as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == "2+" & 
                                                                             melted$gene_biotype_broad=="protein_coding" & 
                                                                             melted$status == "test/train",], function(x) (median(x) + sd(x)/sqrt(length(x))))),
                                       as.data.frame(aggregate(value ~ variable, 
                                                               data=melted[melted$annotated_BPs == "2+" & 
                                                                             melted$gene_biotype_broad=="protein_coding" & 
                                                                             melted$status == "predicted",], function(x) (median(x) + sd(x)/sqrt(length(x))))))[,2]

conservation_summary$multiplicity <- factor(conservation_summary$multiplicity, levels=c("single_BP","multi_BP"))

G26_all$id <- with(G26_all, paste0(exon_id, "_", to_3prime, "_",branchpoint_nt))

m <- match(G26_all$id, conservation_df$id)
G26_all <- cbind(G26_all, conservation_df[m,c(1:18)])

m <- match(branchpoints_introns$introns, G26_all$exon_id[which(G26_all$branchpoint_prob > cutoff)])
branchpoints_introns <-
  cbind(branchpoints_introns, G26_all[which(G26_all$branchpoint_prob > cutoff)[m],c("cons_SS_p0","cons_SS_n1",
                                                                                    "cons_exon_mean","cons_exon_med",
                                                                                    "cons_intron_mean","cons_intron_med")])

###### multi-branchpoint distances ######

w <- which(branchpoints_introns$predicted_BPs > 1)

m <- match(G26_all$exon_id[G26_all$branchpoint_prob > cutoff], branchpoints_introns$introns[w])
mt <- as.data.frame(table(m))
mt$name <- branchpoints_introns$introns[w[mt$m]]

small_G26 <- G26_all[G26_all$branchpoint_prob > cutoff ,c("exon_id","to_3prime")]

for(i in 1:max(mt$Freq)){
  n <- match(mt$name, small_G26$exon_id)
  mt <- cbind(mt, small_G26$to_3prime[n])
  small_G26 <- small_G26[-n[which(!is.na(n))],]
}

todim <- dim(mt)[2]

mt$min_dist_diff <- apply(mt[,c(4:todim)],1,function(x) min(diff(sort(as.numeric(x)))) )
mt$max_dist_diff <- apply(mt[,c(4:todim)],1,function(x) max(diff(sort(as.numeric(x)))) )

mt$max_dist <- apply(mt[,c(4:todim)],1,max, na.rm=T)
mt$min_dist <- apply(mt[,c(4:todim)],1,min, na.rm=T)
# mt$max_dist_diff <- mt$max_dist - mt$min_dist

m <- match(branchpoints_introns$introns, mt$name) 
branchpoints_introns <- cbind(branchpoints_introns, mt[m, c("max_dist","min_dist","max_dist_diff","min_dist_diff")])    

save(G26_all, branchpoints_introns, fivemer_summary,conservation_summary,
     file = "data/Figure_files.Rdata")