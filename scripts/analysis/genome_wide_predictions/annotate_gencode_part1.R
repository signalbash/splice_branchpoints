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
load("data/genome_predictions/gencode.v24.RData")

exons <- readExonAnnotation("data/genome_annotations/gencode.v24.annotation.exons.txt")

G24_all <- gencode_v24

rm(gencode_v24)
###### Make intron centric annotation and annotate structural features ######
colnames(G24_all)[1] <- "id"
#make intron-centric annotation
introns <- unique(G24_all$id)
branchpoints_introns <- data.frame(introns)
branchpoints_introns$status <- "predicted"

#testing/training introns all contained at least one HC
BPs <- which(G24_all$in_testtrain == "HC")
branchpoints_introns$status[which(!is.na(match(branchpoints_introns$introns,
                                               G24_all$id[BPs])))] <- "test/train"

#number of known (HC) BPs
branchpoints_introns$known_BPs <- 0
known_BPs <- as.data.frame(table(G24_all$id[BPs]))
m <- match(branchpoints_introns$introns, known_BPs$Var1)
branchpoints_introns$known_BPs[which(!is.na(m))] <- known_BPs$Freq[m[!is.na(m)]]

#number of predicted BPs
branchpoints_introns$predicted_BPs <- 0
predicted_BPs <- as.data.frame(table(G24_all$id[G24_all$branchpoint_prob >= cutoff]))
m <- match(branchpoints_introns$introns, predicted_BPs$Var1)
branchpoints_introns$predicted_BPs[which(!is.na(m))] <-
  predicted_BPs$Freq[m[which(!is.na(m))]]

#number of annotated BPs
#If a BP is in HC dataset, use that, else use predicted scores
branchpoints_introns$annotated_BPs <- 0
annotated_BPs <- as.data.frame(table(G24_all$id[G24_all$branchpoint_prob >= cutoff |
                                          G24_all$in_testtrain == "HC"]))
m <- match(branchpoints_introns$introns, annotated_BPs$Var1)
branchpoints_introns$annotated_BPs[which(!is.na(m))] <-
  annotated_BPs$Freq[m[which(!is.na(m))]]

branchpoints_introns$status[branchpoints_introns$annotated_BPs == 0] <- "unknown"

#get intron size and biotype
m <- match(branchpoints_introns$introns, G24_all$id)
branchpoints_introns$intron_size <- G24_all$to_3prime[m] + G24_all$to_5prime[m]
n <- match(G24_all$id, exons$exon_id)
G24_all$gene_biotype <- exons$gene_type[n]
branchpoints_introns$gene_biotype <- G24_all$gene_biotype[m]

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
  which(branchpoints_introns$annotated_BPs > 1)] <- "2+"
branchpoints_introns$predicted_BPs_factor <- branchpoints_introns$predicted_BPs
branchpoints_introns$predicted_BPs_factor[
  which(branchpoints_introns$predicted_BPs > 1)] <- "2+"

#Get exon types for protein coding genes CDS/UTR
#gtf annotation with CDS and UTR types only
#from format_cdsutr_annotations.sh on gencode 24
gencode_with_cdsutr <- as.data.frame(fread("data/genome_annotations/gencode.v24.condensed.cdsutr.txt"))
#arrange by transcript (not gene), then location - to order in +ve strand fashion
gencode_with_cdsutr <- arrange(gencode_with_cdsutr, V8,V1,V3,V4)

#find transcript boundries to annotate 1st and last UTRs for each transcript
next_entry_same <-
  gencode_with_cdsutr$V8[1:(length(gencode_with_cdsutr$V1) - 1)] ==
  gencode_with_cdsutr$V8[2:(length(gencode_with_cdsutr$V1))]
gencode_with_cdsutr$next_entry_same <- c(next_entry_same, F)
prev_entry_same <-
  gencode_with_cdsutr$V8[2:(length(gencode_with_cdsutr$V1))] ==
  gencode_with_cdsutr$V8[1:(length(gencode_with_cdsutr$V1) - 1)]
gencode_with_cdsutr$prev_entry_same <- c(F,prev_entry_same)
gencode_with_cdsutr$V2[which(gencode_with_cdsutr$prev_entry_same == F &
                               gencode_with_cdsutr$V2 == "UTR")] <- "UTR1"
gencode_with_cdsutr$V2[which(gencode_with_cdsutr$next_entry_same == F &
                               gencode_with_cdsutr$V2 == "UTR")] <- "UTR2"

#extend UTR annotations if over multiple exons
gencode_with_cdsutr$next_entry <-
  c(gencode_with_cdsutr$V2[2:(length(gencode_with_cdsutr$V1))],"CDS")
gencode_with_cdsutr$prev_entry <-
  c("CDS",gencode_with_cdsutr$V2[1:(length(gencode_with_cdsutr$V1) - 1)])

while (any(gencode_with_cdsutr$V2 == "UTR")) {
  gencode_with_cdsutr$V2[which(gencode_with_cdsutr$V2 == "UTR" &
                                 gencode_with_cdsutr$prev_entry == "UTR1")] <- "UTR1"
  gencode_with_cdsutr$V2[which(gencode_with_cdsutr$V2 == "UTR" &
                                 gencode_with_cdsutr$next_entry == "UTR2")] <- "UTR2"
  #reannotate
  gencode_with_cdsutr$next_entry <-
    c(gencode_with_cdsutr$V2[2:(length(gencode_with_cdsutr$V1))],"CDS")
  gencode_with_cdsutr$prev_entry <-
    c("CDS",gencode_with_cdsutr$V2[1:(length(gencode_with_cdsutr$V1) - 1)])
}

#convert to 3/5' UTRs based on strand
gencode_with_cdsutr$V2[which(gencode_with_cdsutr$V5 == "+" &
                               gencode_with_cdsutr$V2 == "UTR1")] <- "UTR5"
gencode_with_cdsutr$V2[which(gencode_with_cdsutr$V5 == "-" &
                               gencode_with_cdsutr$V2 == "UTR2")] <- "UTR5"
gencode_with_cdsutr$V2[which(gencode_with_cdsutr$V5 == "+" &
                               gencode_with_cdsutr$V2 == "UTR2")] <- "UTR3"
gencode_with_cdsutr$V2[which(gencode_with_cdsutr$V5 == "-" &
                               gencode_with_cdsutr$V2 == "UTR1")] <- "UTR3"

#rearrange by position
gencode_with_cdsutr <- arrange(gencode_with_cdsutr, V1,V3,V4,V8)

#annotate main dataset
G24_all$exon1_contains <- NA
G24_all$exon2_contains <- NA

chroms <- unique(gencode_with_cdsutr$V1)
for (c in seq_along(chroms)) {
  gencode_with_cdsutr_chr <-
    gencode_with_cdsutr[gencode_with_cdsutr$V1 == chroms[c],]
  
  #get CDS locations as vectors for each strand
  gencode_with_cds <-
    as.data.frame(gencode_with_cdsutr_chr[gencode_with_cdsutr_chr$V2 == "CDS",])
  cds_locations_n <-
    apply(gencode_with_cds[which(gencode_with_cds$V5 == "-"),c(3,4)], 1, 
          function(x) {seq(from = x[1], to = x[2])})
  cds_locations_p <-
    apply(gencode_with_cds[which(gencode_with_cds$V5 == "+"),c(3,4)], 1, 
          function(x) {seq(from = x[1], to = x[2])})
  cds_locations_p <- unlist(cds_locations_p)
  cds_locations_n <- unlist(cds_locations_n)
  
  #get 5'UTR locations as vectors for each strand
  gencode_with_utr <-
    as.data.frame(gencode_with_cdsutr_chr[gencode_with_cdsutr_chr$V2 == "UTR5",])
  utr5_locations_n <-
    apply(gencode_with_utr[which(gencode_with_utr$V5 == "-"),c(3,4)], 1, 
          function(x) {seq(from = x[1], to = x[2])})
  utr5_locations_p <-
    apply(gencode_with_utr[which(gencode_with_utr$V5 == "+"),c(3,4)], 1, 
          function(x) {seq(from = x[1], to = x[2])})
  utr5_locations_p <- unlist(utr5_locations_p)
  utr5_locations_n <- unlist(utr5_locations_n)
  
  #get 3'UTR locations as vectors for each strand
  gencode_with_utr <-
    as.data.frame(gencode_with_cdsutr_chr[gencode_with_cdsutr_chr$V2 == "UTR3",])
  utr3_locations_n <-
    apply(gencode_with_utr[which(gencode_with_utr$V5 == "-"),c(3,4)], 1, 
          function(x) {seq(from = x[1], to = x[2])})
  utr3_locations_p <-
    apply(gencode_with_utr[which(gencode_with_utr$V5 == "+"),c(3,4)], 1, 
          function(x) {seq(from = x[1], to = x[2])})
  utr3_locations_p <- unlist(utr3_locations_p)
  utr3_locations_n <- unlist(utr3_locations_n)
  
  G24_pc <- G24_all[which(G24_all$chromosome == chroms[c] &
                    G24_all$gene_biotype == "protein_coding"),]
  
  m3 <- match(G24_pc$exon_3prime, exons$exon_id)
  m5 <- match(G24_pc$exon_5prime, exons$exon_id)
  
  #find the end of the 5'exon in each of the location vectors
  u5_ind <-
    match(exons$end[m5[G24_pc$strand == "+"]], utr5_locations_p)
  u3_ind <-
    match(exons$end[m5[G24_pc$strand == "+"]], utr3_locations_p)
  cds_ind <-
    match(exons$end[m5[G24_pc$strand == "+"]], cds_locations_p)
  
  
  annotation <- rep(NA, length(G24_pc$id))
  #some exons annotated as both 5' and 3' utrs
  annotation[which(!is.na(u5_ind) & !is.na(u3_ind))] <- "UTR3+5"
  annotation[which(!is.na(u5_ind) & is.na(u3_ind))] <- "UTR5"
  annotation[which(is.na(u5_ind) & !is.na(u3_ind))] <- "UTR3"
  #CDS annotation trumps all (d/t isoforms)
  annotation[which(!is.na(cds_ind))] <- "CDS"
  G24_pc$exon1_contains <- annotation
  
  #find the start of the 3'exon in each of the location vectors
  u5_ind <-
    match(exons$start[m3][G24_pc$strand == "-"], utr5_locations_n)
  u3_ind <-
    match(exons$start[m3][G24_pc$strand == "-"], utr3_locations_n)
  cds_ind <-
    match(exons$start[m3][G24_pc$strand == "-"], cds_locations_n)
  
  annotation <- rep(NA, length(G24_pc$id))
  #some exons annotated as both 5' and 3' utrs
  annotation[which(!is.na(u5_ind) & !is.na(u3_ind))] <- "UTR3+5"
  annotation[which(!is.na(u5_ind) & is.na(u3_ind))] <- "UTR5"
  annotation[which(is.na(u5_ind) & !is.na(u3_ind))] <- "UTR3"
  #CDS annotation trumps all (d/t isoforms)
  annotation[which(!is.na(cds_ind))] <- "CDS"
  G24_pc$exon2_contains <- annotation
  
  m <- match(G24_pc$id, G24_all$id)
  
  G24_all$exon1_contains[m] <- G24_pc$exon1_contains
  G24_all$exon2_contains[m] <- G24_pc$exon2_contains
  
  message(c)
}

#number of alternative 5'exons by gencode annotation (not SJ counts)
#and maxentscan score

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

G24_all$exon2_start <- exons$start[match(G24_all$exon_3prime, exons$exon_id)]
G24_all$exon2_end <- exons$end[match(G24_all$exon_3prime, exons$exon_id)]

chr_start_BP <- with(G24_all, paste(chromosome,exon2_start,sep = "_"))
chr_start_BP <- gsub(" ","", chr_start_BP)
chr_end_BP <- with(G24_all, paste(chromosome,exon2_end,sep = "_"))
chr_end_BP <- gsub(" ","", chr_end_BP)

m_p <- match(chr_start_BP, exon5count_pos$Var1)
m_n <- match(chr_end_BP, exon5count_neg$Var1)

gencode_alternative5 <- exon5count_pos$Freq[m_p]
gencode_alternative5[G24_all$strand == "-"] <-
  (exon5count_neg$Freq[m_n])[G24_all$strand == "-"]
G24_all$gencode_alternative5 <- gencode_alternative5
branchpoints_introns$gencode_alternative5 <- 
  G24_all$gencode_alternative5[match(branchpoints_introns$introns, G24_all$id)]

#make bed file for getting .fa sequences
bed <- exons_pos_introns[,c(5,1,1,11,5,6)]
bed$chromosome.1 <- 0
head(bed)
bed$start <- bed$start-3
bed$start.1 <- bed$start.1 + 6

write.table(bed, sep = "\t", file = "data/5prime_seqs_pos.bed",
  row.names = F,col.names = F,quote = F)

bed <- exons_neg_introns[,c(5,2,2,11,5,6)]
bed$chromosome.1 <- 0
head(bed)
bed$end <- bed$end-7
bed$end.1 <- bed$end.1+2

write.table(bed, sep = "\t", file = "data/5prime_seqs_neg.bed",
            row.names = F,col.names = F,quote = F)

cmd <- paste0("/Applications/apps/bedtools2/bin/bedtools"," getfasta -fi ", 
              "data/genome_annotations/GRCh38.p5.genome.fa",
              " -bed data/5prime_seqs_pos.bed -fo ",
              "data/5prime_seqs_pos.fa -name -s")
#system(cmd)
cmd <- paste0("/Applications/apps/bedtools2/bin/bedtools"," getfasta -fi ", 
              "data/genome_annotations/GRCh38.p5.genome.fa",
              " -bed data/5prime_seqs_neg.bed -fo ",
              "data/5prime_seqs_neg.fa -name -s")
#system(cmd)
#Run MaxEntScan 

fiveprime_seqs_neg_MaxEntScan <- read.csv("data/5prime_seqs_neg_MaxEntScan_5scoresplice_Output.txt", header=FALSE, sep=";")
s <- fiveprime_seqs_neg_MaxEntScan[seq(2,dim(fiveprime_seqs_neg_MaxEntScan)[1],by = 2),1]
df <- as.data.frame(str_split(s, "\t"))
df <- as.data.frame(t(df))
rownames(df) <- NULL
df[,2] <- as.numeric(gsub("MAXENT: ","", df[,2]))
df[,3] <- as.numeric(gsub("MDD: ","", df[,3]))
df[,4] <- as.numeric(gsub("MM: ","", df[,4]))
df[,5] <- as.numeric(gsub("WMM: ","", df[,5]))
df[,6] <- NULL
colnames(df) <- c("sequence", "MAXENT","MDD","MM","WMM")
df$id <- gsub(">","",fiveprime_seqs_neg_MaxEntScan[seq(1,dim(fiveprime_seqs_neg_MaxEntScan)[1],by = 2),1])
splicesite5_neg <- df

exons_neg_introns$MAXENTSCAN <- splicesite5_neg$MAXENT[match(exons_neg_introns$chr_startend_gencode, splicesite5_neg$id)]
fiveprime_seqs_pos_MaxEntScan <- read.csv("data/5prime_seqs_pos_MaxEntScan_5scoresplice_Output.txt", header=FALSE, sep=";")
s <- fiveprime_seqs_pos_MaxEntScan[seq(2,dim(fiveprime_seqs_pos_MaxEntScan)[1],by = 2),1]
df <- as.data.frame(str_split(s, "\t"))
df <- as.data.frame(t(df))
rownames(df) <- NULL
df[,2] <- as.numeric(gsub("MAXENT: ","", df[,2]))
df[,3] <- as.numeric(gsub("MDD: ","", df[,3]))
df[,4] <- as.numeric(gsub("MM: ","", df[,4]))
df[,5] <- as.numeric(gsub("WMM: ","", df[,5]))
df[,6] <- NULL
colnames(df) <- c("sequence", "MAXENT","MDD","MM","WMM")
df$id <- gsub(">","",fiveprime_seqs_pos_MaxEntScan[seq(1,dim(fiveprime_seqs_pos_MaxEntScan)[1],by = 2),1])

splicesite5_pos <- df
exons_pos_introns$MAXENTSCAN <- splicesite5_pos$MAXENT[match(exons_pos_introns$chr_startend_gencode, splicesite5_pos$id)]

#summary for multiple 5' exons
max_maxedntscan <- aggregate(MAXENTSCAN ~ chr_end_gencode, exons_pos_introns, max)
var_maxedntscan <- aggregate(MAXENTSCAN ~ chr_end_gencode, exons_pos_introns, var)
min_maxedntscan <- aggregate(MAXENTSCAN ~ chr_end_gencode, exons_pos_introns, min)
med_maxedntscan <- aggregate(MAXENTSCAN ~ chr_end_gencode, exons_pos_introns, median)
m_p <- match(chr_start_BP,max_maxedntscan$chr_end_gencode)
G24_all$max_maxentscan <- max_maxedntscan$MAXENTSCAN[m_p]
G24_all$var_maxentscan <- var_maxedntscan$MAXENTSCAN[m_p]
G24_all$min_maxentscan <- min_maxedntscan$MAXENTSCAN[m_p]
G24_all$med_maxentscan <- med_maxedntscan$MAXENTSCAN[m_p]
max_maxedntscan <- aggregate(MAXENTSCAN ~ chr_start_gencode, exons_neg_introns, max)
var_maxedntscan <- aggregate(MAXENTSCAN ~ chr_start_gencode, exons_neg_introns, var)
min_maxedntscan <- aggregate(MAXENTSCAN ~ chr_start_gencode, exons_neg_introns, min)
med_maxedntscan <- aggregate(MAXENTSCAN ~ chr_start_gencode, exons_neg_introns, median)
m_p <- match(chr_end_BP,max_maxedntscan$chr_start_gencode)
G24_all$max_maxentscan[G24_all$strand=="-"] <- max_maxedntscan$MAXENTSCAN[m_p][G24_all$strand=="-"] 
G24_all$var_maxentscan[G24_all$strand=="-"] <- var_maxedntscan$MAXENTSCAN[m_p][G24_all$strand=="-"] 
G24_all$min_maxentscan[G24_all$strand=="-"] <- min_maxedntscan$MAXENTSCAN[m_p][G24_all$strand=="-"] 
G24_all$med_maxentscan[G24_all$strand=="-"] <- med_maxedntscan$MAXENTSCAN[m_p][G24_all$strand=="-"] 

branchpoints_introns$max_maxentscan <- G24_all$max_maxentscan[match(branchpoints_introns$introns, G24_all$id)]
branchpoints_introns$var_maxentscan <- G24_all$var_maxentscan[match(branchpoints_introns$introns, G24_all$id)]
branchpoints_introns$min_maxentscan <- G24_all$min_maxentscan[match(branchpoints_introns$introns, G24_all$id)]
branchpoints_introns$med_maxentscan <- G24_all$med_maxentscan[match(branchpoints_introns$introns, G24_all$id)]

#match up site-wise and intron-wise attributes
m <- match(branchpoints_introns$introns, G24_all$id)
branchpoints_introns$exon1_contains <- G24_all$exon1_contains[m]
branchpoints_introns$exon2_contains <- G24_all$exon2_contains[m]
branchpoints_introns$gencode_alternative5 <-
  G24_all$gencode_alternative5[m]

m <- match(G24_all$id,branchpoints_introns$introns)
G24_all$exon1_contains <- branchpoints_introns$exon1_contains[m]
G24_all$exon2_contains <- branchpoints_introns$exon2_contains[m]
G24_all$known_BPs <- branchpoints_introns$known_BPs[m]
G24_all$predicted_BPs <- branchpoints_introns$predicted_BPs[m]
G24_all$annotated_BPs <- branchpoints_introns$annotated_BPs[m]
G24_all$predicted_BPs_factor <-
  branchpoints_introns$predicted_BPs_factor[m]
G24_all$annotated_BPs_factor <-
  branchpoints_introns$annotated_BPs_factor[m]
G24_all$intron_size <- branchpoints_introns$intron_size[m]

rm(annotated_BPs, known_BPs, predicted_BPs, G24_pc, gencode_with_cdsutr,
  gencode_with_cdsutr_chr, gencode_with_cds,gencode_with_utr)

rm(annotation,biotype,BPs,c,cds_ind,cds_locations_n,cds_locations_p,chroms,
  introns,lncRNAs,m,next_entry_same,other,prev_entry_same,
  pseudo,u3_ind,u5_ind,utr3_locations_n,utr3_locations_p,
  utr5_locations_n,utr5_locations_p)
rm(exon5count_pos,exon5count_neg,exons_pos_introns,
  exons_neg_introns,chr_end_BP,chr_start_BP,
  m_n,m_p,gencode_alternative5)
rm(bed, df,fiveprime_seqs_pos_MaxEntScan,
   fiveprime_seqs_neg_MaxEntScan,
   max_maxedntscan,med_maxedntscan, min_maxedntscan, 
   var_maxedntscan, splicesite5_pos,splicesite5_neg)
###### Information on "probability scores" and U2 binding energy ######

#get fivemer frequencies & mean probability scores
fivemers <- with(G24_all, paste(seq_neg3,seq_neg2,seq_neg1,
                                seq_pos0,seq_pos1, sep = ""))
G24_all$fivemers <- fivemers
fivemer_summary <- as.data.frame(table(G24_all$fivemers))

fivemer_summary <- cbind(fivemer_summary,
                         aggregate(branchpoint_prob ~ fivemers, 
                                   data = G24_all, min)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers,data = G24_all, quantile,0.25)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G24_all, median)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G24_all, mean)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G24_all, quantile,0.75)$branchpoint_prob)
fivemer_summary <- cbind(fivemer_summary,
        aggregate(branchpoint_prob ~ fivemers, data = G24_all, max)$branchpoint_prob)
colnames(fivemer_summary) <- c("fivemer", "count","min_branchpoint_prob","Q1_branchpoint_prob",
                               "median_branchpoint_prob","mean_branchpoint_prob","Q3_branchpoint_prob","max_branchpoint_prob")

#plot 5mer freq by BP/Neg ratios
fivemer_freq_BP <-
  as.data.frame(table(G24_all$fivemers[G24_all$in_testtrain == "HC"]))
m <- match(fivemer_summary$fivemer,fivemer_freq_BP$Var1)
fivemer_summary$num_BP <- fivemer_freq_BP$Freq[m]
fivemer_summary$num_BP[is.na(fivemer_summary$num_BP)] <- 0
fivemer_freq_N <-
  as.data.frame(table(G24_all$fivemers[G24_all$in_testtrain == "NEG"]))
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

u2_eightmers <- with(G24_all, paste(seq_neg5,seq_neg4,seq_neg3,
                                    seq_neg2,seq_neg1,seq_pos1,
                                    seq_pos2,seq_pos3, sep =""))
m <- match(u2_eightmers, U2_binding_df$eightmers)
G24_all$U2_binding_energy <- U2_binding_df$energy[m]
G24_all$U2_eightmers <- u2_eightmers

U2_summary <- aggregate(U2_binding_energy ~ fivemers, 
                        data = G24_all, min)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G24_all, quantile,0.25)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G24_all, median)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G24_all, mean)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G24_all, quantile,0.75)$U2_binding_energy)
U2_summary <- cbind(U2_summary,aggregate(U2_binding_energy ~ fivemers, 
                                         data = G24_all, max)$U2_binding_energy)
colnames(U2_summary) <- c("fivemer", "min_U2","Q1_U2","median_U2",
                          "mean_U2","Q3_U2","max_U2")
fivemer_summary <- cbind(fivemer_summary, 
                         U2_summary[match(fivemer_summary$fivemer, 
                                          U2_summary$fivemer),-1])

fivemer_summary$BP_nt <- str_sub(fivemer_summary$fivemer, 4,4)
fivemer_summary <- fivemer_summary[fivemer_summary$BP_nt!="N",]


rm(eightmer_df,eightmers,m,m1,m2,seqs,U2_binding_df,
   u2_eightmers,U2_summary,U2_to_eightmers_output)

###### alternative splicing ######

###### K562 STAR ######
rep1_SJ <- read.delim("data/ENCODE_K562/K562_rep1_other_files/SJ.out.tab", header = FALSE)
rep2_SJ <- read.delim("data/ENCODE_K562/K562_rep2_other_files/SJ.out.tab", header = FALSE)

colnames(rep1_SJ) <- c("chromosome", "intron_start","intron_end","strand",
                       "intron_motif","annotated","num_unique_reads",
                       "num_multi_reads", "max_overhang")
colnames(rep2_SJ) <- c("chromosome", "intron_start","intron_end","strand",
                       "intron_motif","annotated","num_unique_reads",
                       "num_multi_reads", "max_overhang")

rep1_SJ$reads <- rep1_SJ$num_unique_reads + rep1_SJ$num_multi_reads
rep2_SJ$reads <- rep2_SJ$num_unique_reads + rep2_SJ$num_multi_reads

get_SJ_exons <-
  function(SJ_data_frame, read_cutoff = 10, value = 0,message = FALSE) {
    SJ_pos <-
      as.data.frame(table(SJ_data_frame$intron_end[SJ_data_frame$strand == 1],SJ_data_frame$chromosome[SJ_data_frame$strand ==
                                                                                                         1]))
    SJ_neg <-
      as.data.frame(table(SJ_data_frame$intron_start[SJ_data_frame$strand == 2],SJ_data_frame$chromosome[SJ_data_frame$strand ==
                                                                                                           2]))
    SJ_pos <- SJ_pos[SJ_pos$Freq > 0,]
    SJ_neg <- SJ_neg[SJ_neg$Freq > 0,]
    SJ_pos$strand <- "+"
    SJ_neg$strand <- "-"
    colnames(SJ_pos) <-
      c("exon_start","chromosome","alternative5'annotated","strand")
    colnames(SJ_neg) <-
      c("exon_start","chromosome","alternative5'annotated","strand")
    SJ <- rbind(SJ_pos,SJ_neg)
    
    SJ$number_alt_exons <- NA
    
    exon_names <- with(SJ, paste(chromosome, exon_start,sep = "_"))
    SJ_names_pos <-
      with(SJ_data_frame, paste(chromosome, intron_end,sep = "_"))
    SJ_names_neg <-
      with(SJ_data_frame, paste(chromosome, intron_start,sep = "_"))
    
    w <-  which(SJ$`alternative5'annotated` > 1)
    
    x_pos <- match(SJ_names_pos,exon_names[SJ$strand == "+"])
    x_neg <- match(SJ_names_neg,exon_names[SJ$strand == "-"])
    
    SJ$exon_names <- exon_names
    
    SJ_pos <- SJ[which(SJ$strand == "+")[x_pos[which(!is.na(x_pos))]],]
    SJ_pos$match <- which(!is.na(x_pos))
    SJ_pos$reads <- SJ_data_frame$reads[SJ_pos$match]
    
    number_alt_exons <-
      aggregate(
        reads ~ exon_names, data = SJ_pos, FUN = function(reads) {
          length(reads >= read_cutoff)
        }
      )
    m <- match(SJ$exon_names,number_alt_exons$exon_names)
    SJ$number_alt_exons[which(SJ$strand == "+")] <-
      number_alt_exons$reads[m[which(SJ$strand == "+")]]
   
    
    SJ_neg <- SJ[which(SJ$strand == "-")[x_neg[which(!is.na(x_neg))]],]
    SJ_neg$match <- which(!is.na(x_neg))
    SJ_neg$reads <- SJ_data_frame$reads[SJ_neg$match]
    
    number_alt_exons <-
      aggregate(
        reads ~ exon_names, data = SJ_neg, FUN = function(reads) {
          length(reads >= read_cutoff)
        }
      )
    m <- match(SJ$exon_names,number_alt_exons$exon_names)
    SJ$number_alt_exons[which(SJ$strand == "-")] <-
      number_alt_exons$reads[m[which(SJ$strand == "-")]]
   
    SJ$exon_start <- as.numeric(levels(SJ$exon_start))[SJ$exon_start]
    SJ$chromosome <- as.character(SJ$chromosome)
    
    return(SJ)
  }

rep1 <- get_SJ_exons(rep1_SJ, value=NA)
rep2 <- get_SJ_exons(rep2_SJ, value=NA)

rep1_names <- with(rep1, paste(exon_names, strand,sep = "_"))
rep2_names <- with(rep2, paste(exon_names, strand,sep = "_"))

SJ_d <- data.frame(names = unique(c(rep1_names,rep2_names)))
m1 <- match(SJ_d$names, rep1_names)
m2 <- match(SJ_d$names, rep2_names)
SJ_d$`alternative5'annotated` <-
  apply(
    cbind(
      rep1$`alternative5'annotated`[m1], rep2$`alternative5'annotated`[m2]
    ), 1, max, na.rm = T
  )
SJ_d$number_alt_exons <-
  apply(cbind(rep1$number_alt_exons[m1], rep2$number_alt_exons[m2]), 1, mean, na.rm =
          T)

m <- match(SJ_d$names, rep1_names)
SJ_d <-
  cbind(SJ_d, rep1[m,c("chromosome", "exon_start","strand")])
m <- match(SJ_d$names, rep2_names)
SJ_d[which(is.na(SJ_d$chromosome)),c("chromosome", "exon_start","strand")] <-
  rep1[m[which(is.na(SJ_d$chromosome))],c("chromosome", "exon_start","strand")]

chroms = unique(SJ_d$chromosome)
index_m = rep(NA, length(G24_all$id))
index_n = rep(NA, length(G24_all$id))

for (c in chroms) {
  m <- match(G24_all$exon2_start,SJ_d$exon_start + 1)
  n <- match(G24_all$exon2_end,SJ_d$exon_start - 1)
  
  index_m[which(G24_all$chromosome == c)] <-
    m[which(G24_all$chromosome == c)]
  index_n[which(G24_all$chromosome == c)] <-
    n[which(G24_all$chromosome == c)]
  
}

index <- index_m
index[G24_all$strand == "-"] <- index_n[G24_all$strand == "-"]
G24_all <- cbind(G24_all, SJ_d[index,c("alternative5'annotated", "number_alt_exons")])

m <- match(branchpoints_introns$introns, G24_all$id)
branchpoints_introns <-
  cbind(branchpoints_introns, G24_all[m,c(
    "alternative5'annotated",
    "number_alt_exons"
  )])

rm(
  c,chroms,get_SJ_exons,index,index_m,index_n,m,
  m1,m2,n,rep1,rep1_names,
  rep1_SJ,rep2,rep2_names,rep2_SJ,SJ_d
)

###### DEXSeq exon expression ######

#DEXseq specific gtf
gencode.v24.annotation_DEXseq <-
  read.delim("data/genome_annotations/gencode.v24.annotation.dexseq2.gtf", header = FALSE)

dex_exons <-
  gencode.v24.annotation_DEXseq[gencode.v24.annotation_DEXseq$V3 == "exonic_part",]
dex_exons$gene_id <-
  gsub(" ","",gsub("gene_id","",grep("gene_id",unlist(
    str_split(dex_exons$V9,";")
  ),value = T)))
dex_exons$transcripts <-
  gsub(" ","",gsub("transcripts","",grep(
    "transcripts",unlist(str_split(dex_exons$V9,";")),value = T
  )))
dex_exons$exonic_part_number <-
  gsub(" ","",gsub(
    "exonic_part_number","",grep("exonic_part_number",unlist(str_split(dex_exons$V9,";")),value = T)
  ))

exon_names <-
  with(dex_exons, paste(gene_id,exonic_part_number,sep = ":"))
dex_exons$exon_name <- exon_names

#combine exonic parts into whole exons
dex_exons$exon_group <- dex_exons$exon_name
same_exon <-
  which((dex_exons$V5[-(length(dex_exons$V5))] == (dex_exons$V4[-1] - 1)) &
          (dex_exons$exon_group[-(length(dex_exons$V5))] != dex_exons$exon_group[-1])
        &
          (dex_exons$gene_id[-(length(dex_exons$V5))] == dex_exons$gene_id[-1])
  )

while (length(same_exon) > 0) {
  dex_exons$exon_group[same_exon + 1] <-
    dex_exons$exon_group[same_exon]
  same_exon <-
    which((dex_exons$V5[-(length(dex_exons$V5))] == (dex_exons$V4[-1] - 1)) &
            (dex_exons$exon_group[-(length(dex_exons$V5))] != dex_exons$exon_group[-1])
          &
            (dex_exons$gene_id[-(length(dex_exons$V5))] == dex_exons$gene_id[-1])
    )
}

rm(gencode.v24.annotation_DEXseq,same_exon,exon_names)

#DEXSeq expression
rep1.dexseq <-
  read.delim(
    "data/ENCODE_K562/DEXseq/K562_rep1.dexseq.txt", header =
      FALSE
  )
rep2.dexseq <-
  read.delim(
    "data/ENCODE_K562/DEXseq/K562_rep2.dexseq.txt", header =
      FALSE
  )

m <- match(rep1.dexseq$V1, rep2.dexseq$V1)
dexseq <- cbind(rep1.dexseq, rep2.dexseq[m,2])
colnames(dexseq) <- c("exon_name", "rep1","rep2")
rownames(dexseq) <- dexseq$exon_name
dexseq$exon_name <- NULL
size_factors <- estimateSizeFactorsForMatrix(dexseq)
dex_seq_norm <- dexseq / size_factors

exon_groups <- unique(dex_exons$exon_group)

dex_exons$mean_DEX_count <- NA
dex_exons <-
  cbind(dex_exons,dex_seq_norm[match(dex_exons$exon_name,
                                     rownames(dex_seq_norm)),])
#get counts per exons groups, not exon parts
rep1_total_counts <- aggregate(rep1 ~ exon_group, dex_exons, sum)
rep2_total_counts <- aggregate(rep2 ~ exon_group, dex_exons, sum)

total_counts <- cbind(rep1_total_counts,rep2_total_counts[,-1])
total_counts$mean <- rowMeans(total_counts[,-1])
m <- match(dex_exons$exon_group, total_counts$exon_group)
dex_exons$group_mean_count <- total_counts$mean[m]
dex_exons$mean_DEX_count <- rowMeans(dex_exons[,c("rep1", "rep2")])

#add to site wide annotation
m1 <- merge(
  G24_all[G24_all$strand == "+",],
  dex_exons[which(dex_exons$V7 == "+"),c(1,4,13,14)],
  by.x = c("chromosome","exon2_start"), by.y = c("V1","V4")
)
m2 <- merge(
  G24_all[G24_all$strand == "-",],
  dex_exons[which(dex_exons$V7 == "-"),c(1,5,13,14)],
  by.x = c("chromosome","exon2_end"), by.y = c("V1","V5")
)
####### check ######
m1_m <- match(G24_all$id, m1$id)
m2_m <- match(G24_all$id, m2$id)

G24_all$dex_exon_name <- m1$exon_name[m1_m]
G24_all$dex_exon_group <- m1$exon_group[m1_m]
G24_all$dex_exon_name[which(is.na(m1_m))] <-
  m2$exon_name[m2_m][which(is.na(m1_m))]
G24_all$dex_exon_group[which(is.na(m1_m))] <-
  m2$exon_group[m2_m][which(is.na(m1_m))]
m <- match(G24_all$dex_exon_name, dex_exons$exon_name)
G24_all$dex_mean_count_group <- dex_exons$group_mean_count[m]
G24_all$dex_mean_count <- dex_exons$mean_DEX_count[m]
G24_all$dex_mean_log <- log10(G24_all$dex_mean_count_group + 0.1)

#add to intron annotation
m <- match(branchpoints_introns$introns, G24_all$id)
branchpoints_introns <-
  cbind(branchpoints_introns, G24_all[m,c("dex_exon_group",
                                          "dex_mean_count_group",
                                          "dex_mean_log")])

rm(
  m1,m2, m1_m,m2_m, dex_exons,dex_seq_norm,dexseq,
  rep1_total_counts,rep1.dexseq,rep2_total_counts,rep2.dexseq,
  size_factors,total_counts
)

save(G24_all, branchpoints_introns, fivemer_summary,
     file = "data/Figure_files.Rdata")

