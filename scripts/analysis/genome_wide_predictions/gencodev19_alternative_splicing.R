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
nt_cols=c("#359646","#4D7ABE","#FAA859","#CB3634")

#from process_predictions.R
load("data/genome_predictions/gencode.v19.RData")

exons <- readExonAnnotation("data/genome_annotations/gencode.v19.annotation.exons.txt")
G19_all <- gencode_v19

rm(gencode_v19)
###### Make intron centric annotation and annotate structural features ######

#make intron-centric annotation
introns <- unique(G19_all$id)
branchpoints_introns <- data.frame(introns)
branchpoints_introns$status <- "predicted"

#testing/training introns all contained at least one HC
BPs <- which(G19_all$in_testtrain == "HC")
branchpoints_introns$status[which(!is.na(match(branchpoints_introns$introns,
                                               G19_all$id[BPs])))] <- "test/train"

#number of known (HC) BPs
branchpoints_introns$known_BPs <- 0
known_BPs <- as.data.frame(table(G19_all$id[BPs]))
m <- match(branchpoints_introns$introns, known_BPs$Var1)
branchpoints_introns$known_BPs[which(!is.na(m))] <- known_BPs$Freq[m[!is.na(m)]]

#number of predicted BPs
branchpoints_introns$predicted_BPs <- 0
predicted_BPs <- as.data.frame(table(G19_all$id[G19_all$branchpoint_prob >= cutoff]))
m <- match(branchpoints_introns$introns, predicted_BPs$Var1)
branchpoints_introns$predicted_BPs[which(!is.na(m))] <-
  predicted_BPs$Freq[m[which(!is.na(m))]]

#number of annotated BPs
#If a BP is in HC dataset, use that, else use predicted scores
branchpoints_introns$annotated_BPs <- 0
annotated_BPs <- as.data.frame(table(G19_all$id[G19_all$branchpoint_prob >= cutoff |
                                          G19_all$in_testtrain == "HC"]))
m <- match(branchpoints_introns$introns, annotated_BPs$Var1)
branchpoints_introns$annotated_BPs[which(!is.na(m))] <-
  annotated_BPs$Freq[m[which(!is.na(m))]]

branchpoints_introns$status[branchpoints_introns$annotated_BPs == 0] <- "unknown"

#get intron size and biotype
m <- match(branchpoints_introns$introns, G19_all$id)
branchpoints_introns$intron_size <- G19_all$to_3prime[m] + G19_all$to_5prime[m]
n <- match(G19_all$id, exons$exon_id)
G19_all$gene_biotype <- exons$gene_type[n]
branchpoints_introns$gene_biotype <- G19_all$gene_biotype[m]

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

G19_all$exon2_start <- exons$start[match(G19_all$exon_3prime, exons$exon_id)]
G19_all$exon2_end <- exons$end[match(G19_all$exon_3prime, exons$exon_id)]

chr_start_BP <- with(G19_all, paste(chromosome,exon2_start,sep = "_"))
chr_start_BP <- gsub(" ","", chr_start_BP)
chr_end_BP <- with(G19_all, paste(chromosome,exon2_end,sep = "_"))
chr_end_BP <- gsub(" ","", chr_end_BP)

m_p <- match(chr_start_BP, exon5count_pos$Var1)
m_n <- match(chr_end_BP, exon5count_neg$Var1)

gencode_alternative5 <- exon5count_pos$Freq[m_p]
gencode_alternative5[G19_all$strand == "-"] <-
  (exon5count_neg$Freq[m_n])[G19_all$strand == "-"]
G19_all$gencode_alternative5 <- gencode_alternative5

#match up site-wise and intron-wise attributes
m <- match(branchpoints_introns$introns, G19_all$id)
branchpoints_introns$exon1_contains <- G19_all$exon1_contains[m]
branchpoints_introns$exon2_contains <- G19_all$exon2_contains[m]
branchpoints_introns$gencode_alternative5 <-
  G19_all$gencode_alternative5[m]

m <- match(G19_all$id,branchpoints_introns$introns)
G19_all$exon1_contains <- branchpoints_introns$exon1_contains[m]
G19_all$exon2_contains <- branchpoints_introns$exon2_contains[m]
G19_all$known_BPs <- branchpoints_introns$known_BPs[m]
G19_all$predicted_BPs <- branchpoints_introns$predicted_BPs[m]
G19_all$annotated_BPs <- branchpoints_introns$annotated_BPs[m]
G19_all$predicted_BPs_factor <-
  branchpoints_introns$predicted_BPs_factor[m]
G19_all$annotated_BPs_factor <-
  branchpoints_introns$annotated_BPs_factor[m]
G19_all$intron_size <- branchpoints_introns$intron_size[m]

rm(annotated_BPs, known_BPs, predicted_BPs)

rm(biotype,BPs,introns,lncRNAs,m,pseudo)
rm(exon5count_pos,exon5count_neg,exons_pos_introns,
  exons_neg_introns,chr_end_BP,chr_start_BP,
  m_n,m_p,gencode_alternative5)

#### U2 binding energy ###
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

u2_eightmers <- with(G19_all, paste(seq_neg5,seq_neg4,seq_neg3,
                                    seq_neg2,seq_neg1,seq_pos1,
                                    seq_pos2,seq_pos3, sep =""))
m <- match(u2_eightmers, U2_binding_df$eightmers)
G19_all$U2_binding_energy <- U2_binding_df$energy[m]
G19_all$U2_eightmers <- u2_eightmers

top_U2 <- aggregate(U2_binding_energy ~ id, G19_all[G19_all$branchpoint_prob >= cutoff,], max)
branchpoints_introns$top_U2 <- top_U2$U2_binding_energy[match(branchpoints_introns$introns, top_U2$id)]

###### exon skipping ######

skipped_annotations <- read.delim("data/skipped_exons_and_control_annotations_gencode19.txt")

m <- match(branchpoints_introns$introns, G19_all$id)
branchpoints_introns$start <- G19_all$exon2_start[m]
branchpoints_introns$end <- G19_all$exon2_end[m]
branchpoints_introns$chromosome <- G19_all$chromosome[m]
branchpoints_introns$strand <- G19_all$strand[m]
branchpoints_introns$AS_label <- NA

bi_newName <- with(branchpoints_introns, paste(chromosome, start,strand,sep="_"))
sa_newName <- with(skipped_annotations, paste(chromosome, start,strand,sep="_"))

m <- match(sa_newName, bi_newName)

branchpoints_introns_sa <- branchpoints_introns[m[which(!is.na(m))],]
branchpoints_introns_sa$sa1 <- skipped_annotations$annotation[which(!is.na(m))]
branchpoints_introns_sa$sa2 <- skipped_annotations$annotation_2[which(!is.na(m))]

table(branchpoints_introns_sa$sa1, branchpoints_introns_sa$sa2)

rm <- which((branchpoints_introns_sa$sa2 == "E2" & branchpoints_introns_sa$sa1 == "skipped_exon" & branchpoints_introns_sa$gencode_alternative5 < 2) |
                (branchpoints_introns_sa$sa2 == "E2" & branchpoints_introns_sa$sa1 == "control" & branchpoints_introns_sa$gencode_alternative5 > 1))
branchpoints_introns_sa <- branchpoints_introns_sa[-rm,]

table(branchpoints_introns_sa$gencode_alternative5[which(branchpoints_introns_sa$sa2 == "E2" & branchpoints_introns_sa$sa1 == "skipped_exon")])
table(branchpoints_introns_sa$gencode_alternative5[which(branchpoints_introns_sa$sa2 == "E2" & branchpoints_introns_sa$sa1 == "control")])

#no difference between number of branchpoints in skipped/constitutive exons
tbl <- (table(branchpoints_introns_sa$annotated_BPs_factor[branchpoints_introns_sa$sa2=="SE" & branchpoints_introns_sa$annotated_BPs != 0], 
              branchpoints_introns_sa$sa1[branchpoints_introns_sa$sa2=="SE"& branchpoints_introns_sa$annotated_BPs != 0]))
chisq.test(tbl)
#no significant difference in branchpoint strength
wilcox.test(branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="E2" & branchpoints_introns_sa$annotated_BPs != 0],
            branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="E2"& branchpoints_introns_sa$annotated_BPs != 0] )

wilcox.test(branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="SE" & branchpoints_introns_sa$annotated_BPs != 0],
            branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="SE"& branchpoints_introns_sa$annotated_BPs != 0] )

wilcox.test(branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="E1" & branchpoints_introns_sa$annotated_BPs != 0],
            branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="E1"& branchpoints_introns_sa$annotated_BPs != 0] )


#get 5'SS strength
bedp <- branchpoints_introns_sa[branchpoints_introns_sa$strand=="+",c("chromosome","end","end","introns","status","strand")]
bedp[,5] <- 0
head(bedp)
bedp[,2] <- bedp[,2]-3
bedp[,3] <- bedp[,3] + 6
bedp[,1] <- gsub("chr", "", bedp[,1])

bedn <- branchpoints_introns_sa[branchpoints_introns_sa$strand=="-",c("chromosome","start","start","introns","status","strand")]
bedn[,5] <- 0
head(bedn)
bedn[,2] <- bedn[,2]-7
bedn[,3] <- bedn[,3] + 2
bedn[,1] <- gsub("chr", "", bedn[,1])

colnames(bedp) <- c("chromosome","start","end","name","score","strand")
colnames(bedn) <- c("chromosome","start","end","name","score","strand")

bed <- rbind(bedp, bedn)

write.table(bed, sep = "\t", file = "data/5prime_seqs_se.bed",
            row.names = F,col.names = F,quote = F)
cmd <- paste0("/Applications/apps/bedtools2/bin/bedtools"," getfasta -fi ", 
              "data/genome_annotations/GRCh37.p13.genome.fa",
              " -bed data/5prime_seqs_se.bed -fo ",
              "data/5prime_seqs_se.fa -name -s")
system(cmd)

fiveprime_seqs_se_MaxEntScan <- read.csv("data/5prime_seqs_se_MaxEntScan_5scoresplice_Output.txt", header=FALSE, sep=";")
s <- fiveprime_seqs_se_MaxEntScan[seq(2,dim(fiveprime_seqs_se_MaxEntScan)[1],by = 2),1]
df <- as.data.frame(str_split(s, "\t"))
df <- as.data.frame(t(df))
rownames(df) <- NULL
df[,2] <- as.numeric(gsub("MAXENT: ","", df[,2]))
df[,3] <- as.numeric(gsub("MDD: ","", df[,3]))
df[,4] <- as.numeric(gsub("MM: ","", df[,4]))
df[,5] <- as.numeric(gsub("WMM: ","", df[,5]))
df[,6] <- NULL
colnames(df) <- c("sequence", "MAXENT","MDD","MM","WMM")
df$id <- gsub(">","",fiveprime_seqs_se_MaxEntScan[seq(1,dim(fiveprime_seqs_se_MaxEntScan)[1],by = 2),1])

splicesite5_se <- df

#get 3'SS strength
bedp <- branchpoints_introns_sa[branchpoints_introns_sa$strand=="+",c("chromosome","start","start","introns","status","strand")]
bedp[,5] <- 0
head(bedp)
bedp[,2] <- bedp[,2] - 21
bedp[,3] <- bedp[,3] + 2
bedp[,1] <- gsub("chr", "", bedp[,1])

bedn <- branchpoints_introns_sa[branchpoints_introns_sa$strand=="-",c("chromosome","end","end","introns","status","strand")]
bedn[,5] <- 0
head(bedn)
bedn[,2] <- bedn[,2] - 3
bedn[,3] <- bedn[,3] + 20
bedn[,1] <- gsub("chr", "", bedn[,1])

colnames(bedp) <- c("chromosome","start","end","name","score","strand")
colnames(bedn) <- c("chromosome","start","end","name","score","strand")

bed <- rbind(bedp, bedn)

write.table(bed, sep = "\t", file = "data/3prime_seqs_se.bed",
            row.names = F,col.names = F,quote = F)
cmd <- paste0("/Applications/apps/bedtools2/bin/bedtools"," getfasta -fi ", 
              "data/genome_annotations/GRCh37.p13.genome.fa",
              " -bed data/3prime_seqs_se.bed -fo ",
              "data/3prime_seqs_se.fa -name -s")
system(cmd)

threeprime_seqs_se_MaxEntScan <- read.csv("data/3prime_seqs_se_MaxEntScan_3scoresplice_Output.txt", header=FALSE, sep=";")
s <- threeprime_seqs_se_MaxEntScan[seq(2,dim(threeprime_seqs_se_MaxEntScan)[1],by = 2),1]
df <- as.data.frame(str_split(s, "\t"))
df <- as.data.frame(t(df))
rownames(df) <- NULL
df[,c(3,5,7,8)] <- NULL
df[,2] <- as.numeric(gsub("MAXENT: ","", df[,2]))
df[,3] <- as.numeric(gsub("MM: ","", df[,3]))
df[,4] <- as.numeric(gsub("WMM: ","", df[,4]))
colnames(df) <- c("sequence", "MAXENT","MM","WMM")
df$id <- gsub(">","",threeprime_seqs_se_MaxEntScan[seq(1,dim(threeprime_seqs_se_MaxEntScan)[1],by = 2),1])

splicesite3_se <- df

m <- match(branchpoints_introns_sa$introns, splicesite5_se$id)
branchpoints_introns_sa$splicesite_5strength <- splicesite5_se$MAXENT[m]
m <- match(branchpoints_introns_sa$introns, splicesite3_se$id)
branchpoints_introns_sa$splicesite_3strength <- splicesite3_se$MAXENT[m]

make_boxPlotDf <- function(variable, condition1, condition2, branchpoints_introns_sa){
    stats <- boxplot.stats(branchpoints_introns_sa[branchpoints_introns_sa$sa1==condition1 & branchpoints_introns_sa$sa2==condition2, match(variable, colnames(branchpoints_introns_sa))])$stats
    df <- data.frame(variable=variable,condition1=condition1,condition2=condition2, 
                     ymin=stats[1], lower=stats[2], middle=stats[3], upper=stats[4], ymax=stats[5])
    return(df)
}
splicing_strength <- make_boxPlotDf("splicesite_5strength", "skipped_exon","E1", branchpoints_introns_sa)
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("splicesite_5strength", "control","E1", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("top_U2", "skipped_exon","SE", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("top_U2", "control","SE", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("splicesite_3strength", "skipped_exon","SE", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("splicesite_3strength", "control","SE", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("splicesite_5strength", "skipped_exon","SE", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("splicesite_5strength", "control","SE", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("top_U2", "skipped_exon","E2", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("top_U2", "control","E2", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("splicesite_3strength", "skipped_exon","E2", branchpoints_introns_sa))
splicing_strength <- rbind(splicing_strength, make_boxPlotDf("splicesite_3strength", "control","E2", branchpoints_introns_sa))


splicing_strength$variable <- gsub("top_U2", "branchpoint strength",splicing_strength$variable)
splicing_strength$variable <- gsub("splicesite_3strength", "3'SS strength",splicing_strength$variable)
splicing_strength$variable <- gsub("splicesite_5strength", "5'SS strength",splicing_strength$variable)

splicing_strength$condition2 <- gsub("E1", "exon 1",splicing_strength$condition2)
splicing_strength$condition2 <- gsub("SE", "skipped exon",splicing_strength$condition2)
splicing_strength$condition2 <- gsub("E2", "exon 2",splicing_strength$condition2)

splicing_strength$variable <- factor(splicing_strength$variable, levels=c("branchpoint strength", "3'SS strength","5'SS strength"))
splicing_strength$condition2 <- factor(splicing_strength$condition2, levels=c("exon 1", "skipped exon","exon 2"))

save(branchpoints_introns_sa, splicing_strength, file="data/exon_skipping.Rdata")
