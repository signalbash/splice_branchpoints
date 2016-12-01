################################################################################
#                  Branchpoints in alternative splicing                        #
################################################################################

options(stringsAsFactors = F)
library(data.table)
library(ggplot2)
library(stringr)
library(plyr)
library(parallel)
library(entropy)
library(equivalence)
library(cowplot)

load("data/Figure_files_G12.Rdata")
cutoff = 0.5

#Miso annotations from: http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_hg19_v2.zip

###### read and reconfigure MISO Annotations ######
#skipped exons
SE.hg19 <- read.delim("data/MISO/hg19/SE.hg19.gff3", header=FALSE)
#Mutaually exclusive
MXE.hg19 <- read.delim("data/MISO/hg19/MXE.hg19.gff3", header=FALSE)
#retained intron
RI.hg19 <- read.delim("data/MISO/hg19/RI.hg19.gff3", header=FALSE)
#alt ss
A5SS.hg19 <- read.delim("data/MISO/hg19/A5SS.hg19.gff3", header=FALSE)
A3SS.hg19 <- read.delim("data/MISO/hg19/A3SS.hg19.gff3", header=FALSE)

#get "exons" only - not genes/transcripts
SE.hg19 <- SE.hg19[SE.hg19$V3=="exon",]
MXE.hg19 <- MXE.hg19[MXE.hg19$V3=="exon",]
RI.hg19 <- RI.hg19[RI.hg19$V3=="exon",]
A3SS.hg19 <- A3SS.hg19[A3SS.hg19$V3=="exon",]
A5SS.hg19 <- A5SS.hg19[A5SS.hg19$V3=="exon",]

#fix ids
id_fix <- function(id_list){
  ids <- unlist(str_split(id_list, ";"))
  ids <- grep("ID", ids, value=T)                   
  ids2 <- unlist(str_split(ids, "[.]"))
  ids2 <- (ids2)[seq(3,length(ids2), 3)]
}
SE.hg19$type <- id_fix(SE.hg19$V9)
SE.hg19$event <- "skipped"
MXE.hg19$type <- id_fix(MXE.hg19$V9)
MXE.hg19$event <- "mut_exl"
RI.hg19$type <- id_fix(RI.hg19$V9)
RI.hg19$event <- "retained_intron"
A3SS.hg19$type <- id_fix(A3SS.hg19$V9)
A3SS.hg19$event <- "alt3"
A5SS.hg19$type <- id_fix(A5SS.hg19$V9)
A5SS.hg19$event <- "alt5"

#combine and name events
miso_alts <- rbind(RI.hg19[(RI.hg19$type=="dn" | RI.hg19$type=="withRI"),],
                MXE.hg19[MXE.hg19$type=="mxe1"| MXE.hg19$type=="mxe2",],
                SE.hg19[(SE.hg19$type=="se" | SE.hg19$type=="dn"),],
                A3SS.hg19[A3SS.hg19$type!="up",],
                A5SS.hg19[A5SS.hg19$type!="dn",])

miso_alts$event[miso_alts$event=="alt3"] <-  paste0("alt3_", miso_alts$type[miso_alts$event=="alt3"])
miso_alts$event[miso_alts$event=="alt5"] <-  paste0("alt5_", miso_alts$type[miso_alts$event=="alt5"])
miso_alts$event[miso_alts$event=="retained_intron"] <-  paste0("retained_intron_", miso_alts$type[miso_alts$event=="retained_intron"])
miso_alts$event[miso_alts$event=="skipped"] <-  paste0("skipped_", miso_alts$type[miso_alts$event=="skipped"])

miso_alts$type <- miso_alts$event

#combine annotations where an exon is involved in multiple events
chroms <- unique(miso_alts$V1)

for(c in 1:length(chroms)){
  chrom <- chroms[c]
  miso_ind <- which(miso_alts$V1==chrom & miso_alts$V7=="+")
  starts <- unique(miso_alts$V4[miso_ind])
  m <- match(miso_alts$V4[miso_ind], starts)
  
  for(n in seq_along(starts)){
    line <- miso_alts[miso_ind[which(m==n)[1]],-c(2)]
    line$type <- paste(miso_alts$type[miso_ind][which(m==n)], collapse = ";")
    if(exists("miso_annotation")){
      miso_annotation <- rbind(miso_annotation,line)
    }else{
      miso_annotation <- line
    }
  }
  
  miso_ind <- which(miso_alts$V1==chrom & miso_alts$V7=="-")
  starts <- unique(miso_alts$V5[miso_ind])
  m <- match(miso_alts$V5[miso_ind], starts)
  
  for(n in seq_along(starts)){
    line <- miso_alts[miso_ind[which(m==n)[1]],-c(2)]
    line$type <- paste(miso_alts$type[miso_ind][which(m==n)], collapse = ";")
    if(exists("miso_annotation")){
      miso_annotation <- rbind(miso_annotation,line)
    }else{
      miso_annotation <- line
    }
  }
  message(chrom)
}

write.table(miso_annotation, "data/miso_annotation12.txt", sep="\t", row.names = F)
miso_annotation <- read.delim("miso_annotation12.txt", header=T)

Homo_sapiens.GRCh37 <- read.delim("data/MISO/hg19/Homo_sapiens.GRCh37.65.gff", header=FALSE)
Homo_sapiens.GRCh37 <- Homo_sapiens.GRCh37[Homo_sapiens.GRCh37$V3=="exon",]
Homo_sapiens.GRCh37$V1 <- paste0("chr", Homo_sapiens.GRCh37$V1)

m <- match(branchpoints_introns$introns, G12_all$id)
branchpoints_introns$start <- G12_all$exon2_start[m]
branchpoints_introns$end <- G12_all$exon2_end[m]
branchpoints_introns$chromosome <- G12_all$chromosome[m]
branchpoints_introns$strand <- G12_all$strand[m]
branchpoints_introns$AS_label <- NA

chroms <- unique(branchpoints_introns$chromosome)

for(c in seq_along(chroms)){
  chrom <- chroms[c]

  m <- match(branchpoints_introns$start[branchpoints_introns$chromosome==chrom & branchpoints_introns$strand=="+"], 
             miso_annotation$V4[miso_annotation$V1==chrom & miso_annotation$V7=="+"])
  n <- match(branchpoints_introns$start[branchpoints_introns$chromosome==chrom & branchpoints_introns$strand=="+"], 
          Homo_sapiens.GRCh37$V4[Homo_sapiens.GRCh37$V1==chrom & Homo_sapiens.GRCh37$V7=="+"])
  
  branchpoints_introns$AS_label[branchpoints_introns$chromosome == chrom & branchpoints_introns$strand == "+"] <- 
    miso_annotation$type[miso_annotation$V1 == chrom & miso_annotation$V7 == "+"][m]
  
  branchpoints_introns$AS_label[branchpoints_introns$chromosome==chrom & branchpoints_introns$strand=="+"][which(is.na(m) & !is.na(n))] <- 
    "normal"
  
  m <- match(branchpoints_introns$end[branchpoints_introns$chromosome==chrom & branchpoints_introns$strand=="-"], 
             miso_annotation$V5[miso_annotation$V1==chrom & miso_annotation$V7=="-"])
  n <- match(branchpoints_introns$end[branchpoints_introns$chromosome==chrom & branchpoints_introns$strand=="-"], 
          Homo_sapiens.GRCh37$V5[Homo_sapiens.GRCh37$V1==chrom & Homo_sapiens.GRCh37$V7=="-"])
  branchpoints_introns$AS_label[branchpoints_introns$chromosome==chrom & branchpoints_introns$strand=="-"] <-  
    miso_annotation$type[miso_annotation$V1==chrom & miso_annotation$V7=="-"][m]
  branchpoints_introns$AS_label[branchpoints_introns$chromosome==chrom & branchpoints_introns$strand=="-"][which(is.na(m) & !is.na(n))] <-  
    "normal"
  
  message(chrom)
}

#get top probability score and U2 binding energy
G12_all <- arrange(G12_all, plyr::desc(branchpoint_prob))
m <- match(branchpoints_introns$introns, G12_all$id)
branchpoints_introns$top_branchpoint_prob <- G12_all$branchpoint_prob[m]
G12_all <- arrange(G12_all, plyr::desc(U2_binding_energy))
m <- match(branchpoints_introns$introns, G12_all$id)
branchpoints_introns$top_U2 <- G12_all$U2_binding_energy[m]

#set single comparison group to NA
branchpoints_introns$AS_group <- NA

branchpoints_introns$multiplicity <- 0
branchpoints_introns$multiplicity[branchpoints_introns$predicted_BPs ==1 ] <- 1
branchpoints_introns$multiplicity[branchpoints_introns$predicted_BPs > 1 ] <- "2+"

#comparisions
alt_splicing_comparison <- function(condition1="normal", condition2, condition3, branchpoint_introns,G12_all, gene_filter=FALSE, filter_multiplicity=0,ep=0.31){

  branchpoints_introns$AS_group <- NA

  branchpoints_introns$AS_group[unlist(lapply(str_split(branchpoints_introns$AS_label, ';'), 
                                              function(x){any(x==condition3 | x==condition2)}))] <- 
    paste(condition2,condition3,sep="+")
  branchpoints_introns$AS_group[unlist(lapply(str_split(branchpoints_introns$AS_label, ';'), 
                                              function(x){any(x==condition2) & all(x!=condition3)}))] <- 
    condition2
  branchpoints_introns$AS_group[unlist(lapply(str_split(branchpoints_introns$AS_label, ';'), 
                                              function(x){any(x==condition3) & all(x!=condition2)}))] <- 
    condition3

  branchpoints_introns$AS_group[which(branchpoints_introns$AS_label == 
                                          condition1)] <- condition1

  branchpoints_introns$AS_group <- factor(branchpoints_introns$AS_group, levels=c(condition3,paste0(condition2,"+", condition3),condition2,condition1))

  dodge <- position_dodge(width = 0.85)
  splicingcols <- c('#2b83ba','#abdda4','#ffffbf','#d7191c')

  Figure7_part <<- ggplot(branchpoints_introns[branchpoints_introns$predicted_BPs >0 & !is.na(branchpoints_introns$AS_group),], 
                          aes(fill=AS_group, y=top_U2, x=AS_group)) + 
    geom_violin(scale="width",position=dodge, width=0.75,lwd=0.25) + 
    geom_boxplot(alpha=0.1, width=0.5,position=dodge, outlier.size = 0.25, lwd=0.25)+
    scale_y_continuous(name="Top U2 binding energy", limits=c(0,12)) + 
    scale_fill_manual(values = splicingcols, name="MISO annotation") + 
    theme_bw() + 
    theme(legend.position = "none",text=element_text(size=8))+
     scale_x_discrete(name="Annotated branchpoints") 

  cond1_u2 <- branchpoints_introns$top_U2[(branchpoints_introns$AS_group==condition1 & 
    !is.na(branchpoints_introns$AS_group) & branchpoints_introns$annotated_BPs!=0 &
                                               branchpoints_introns$predicted_BPs >=1)]
  cond2_u2 <- branchpoints_introns$top_U2[(branchpoints_introns$AS_group==condition2 & 
    !is.na(branchpoints_introns$AS_group)& branchpoints_introns$annotated_BPs!=0 & 
                                               branchpoints_introns$predicted_BPs>=1)]
  cond3_u2 <- branchpoints_introns$top_U2[(branchpoints_introns$AS_group==condition3 & 
    !is.na(branchpoints_introns$AS_group)& branchpoints_introns$annotated_BPs!=0 & 
                                               branchpoints_introns$predicted_BPs >=1)]
  cond23_u2 <- branchpoints_introns$top_U2[(branchpoints_introns$AS_group==paste(condition2,condition3,sep="+") & 
    !is.na(branchpoints_introns$AS_group)& branchpoints_introns$annotated_BPs!=0 & 
                                                branchpoints_introns$predicted_BPs >=1)]

  w1 <- wilcox.test(cond1_u2,cond2_u2)
  w2 <- wilcox.test(cond1_u2,c(cond2_u2, cond23_u2))
  w3 <- wilcox.test(cond1_u2,cond3_u2)
  w4 <- wilcox.test(cond1_u2,c(cond3_u2, cond23_u2))
  w5 <- wilcox.test(cond2_u2,cond3_u2)
  w6 <- wilcox.test(cond1_u2,c(cond2_u2,cond3_u2, cond23_u2))

  wilcox_df <- data.frame(condition1=c(condition1,condition1,condition1,condition1,condition2,condition1),
                          condition2=c(condition2,paste0(condition2,"/",condition2,"+",condition3), condition3,
                                       paste0(condition3,"/",condition2,"+",condition3), condition3,
                                       paste0(condition2,"/",condition3,"/",condition2,"+",condition3)),
                          cond1_mean=c(mean(cond1_u2),mean(cond1_u2),mean(cond1_u2),
                                       mean(cond1_u2),mean(cond2_u2),mean(cond1_u2)),
                          cond2_mean=c(mean(cond2_u2),mean(c(cond2_u2,cond23_u2)),mean(cond3_u2),
                                       mean(c(cond3_u2,cond23_u2)),mean(cond3_u2),
                                       mean(c(cond3_u2,cond23_u2, cond2_u2))),
                          cond1_median=c(median(cond1_u2),median(cond1_u2),median(cond1_u2),
                                       median(cond1_u2),median(cond2_u2),median(cond1_u2)),
                          cond2_median=c(median(cond2_u2),median(c(cond2_u2,cond23_u2)),median(cond3_u2),
                                       median(c(cond3_u2,cond23_u2)),median(cond3_u2),
                                       median(c(cond3_u2,cond23_u2, cond2_u2))),
                          wilcox.pvalue = c(w1$p.value, w2$p.value, w3$p.value,
                                     w4$p.value, w5$p.value, w6$p.value),
                          wilcox.W=c(w1$statistic, w2$statistic, w3$statistic,
                              w4$statistic, w5$statistic, w6$statistic))

  return(wilcox_df)
}


stats_frame <- alt_splicing_comparison("normal", "skipped_dn","skipped_se", branchpoint_introns, G12_all,ep = 0.05)
Figure7A <- Figure7_part

stats_frame <- rbind(stats_frame,alt_splicing_comparison("normal", "retained_intron_dn","retained_intron_withRI", branchpoint_introns, G12_all,ep = 0.05))
Figure7B <- Figure7_part

stats_frame <- rbind(stats_frame,alt_splicing_comparison("normal", "alt3_core","alt3_coreAndExt", branchpoint_introns, G12_all,ep = 0.05))
Figure7C <- Figure7_part

Figure7=ggdraw() + draw_plot(Figure7A, 0,0,0.33,1) + 
  draw_plot(Figure7B, 0.33,0,0.33,1) + 
  draw_plot(Figure7C, 0.66,0,0.33,1) + 
  draw_plot_label(c("A","B","C"), c(0,0.33,0.66), c(1,1,1), size=10)

Figure7
pdf("Figures/Figure7_cowplot.pdf", useDingbats = F, height=2.5, width=6.69)
Figure7
dev.off()
