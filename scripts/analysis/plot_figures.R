################################################################################
#           Genome wide prediction of splicing branchpoints - plots            #
################################################################################

# Produces Figures 2 - 4
options(stringsAsFactors = F)

library(branchpointer)
library(ggplot2)
library(cowplot)
library(stringr)
library(pheatmap)
library(data.table)
source("scripts/analysis/genome_wide_predictions/quick_wilcox.R")

theme_figure <- theme_bw()+ theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_blank(),
                                  axis.line.x = element_line(), 
                                  axis.line.y = element_line(), 
                                  panel.background = element_rect(colour = "black", size=1, fill=NA),
                                  plot.title = element_text(hjust = 0.5))

number_introns_pal=rev(c('#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'))
nt_cols=c("#359646","#4D7ABE","#FAA859","#CB3634")
BP_multi_colors=c("#bababa", "#1a9850", "#a6d96a", "#ffffbf","#fdae61")
BP_multi_colors=c("#bababa", "#deebf7", "#9ecae1", "#4292c6","#08519c")

### Figure 1A - locations of Mercer Branchpoints ###

files <- list.files("data/outputs/branchpoints_for_training/by_type" ,
                    pattern = "hc",full.names = TRUE)
for(f in files){
    Mercer_branchpoints_hc_chrom <- read.csv(f, row.names=1)
    if(exists("Mercer_branchpoints_hc")){
        Mercer_branchpoints_hc <- rbind(Mercer_branchpoints_hc, Mercer_branchpoints_hc_chrom)
    }else{
        Mercer_branchpoints_hc <- Mercer_branchpoints_hc_chrom
    }
}

files <- list.files("data/outputs/branchpoints_for_training/by_type" ,
                    pattern = "lc",full.names = TRUE)
for(f in files){
    Mercer_branchpoints_lc_chrom <- read.csv(f, row.names=1)
    if(exists("Mercer_branchpoints_lc")){
        Mercer_branchpoints_lc <- rbind(Mercer_branchpoints_lc, Mercer_branchpoints_lc_chrom)
    }else{
        Mercer_branchpoints_lc <- Mercer_branchpoints_lc_chrom
    }
}
Mercer_branchpoints_hc$set <- "HC"
Mercer_branchpoints_lc$set <- "LC"

Mercer_branchpoints <- rbind(Mercer_branchpoints_hc,Mercer_branchpoints_lc)

load("data/Figure_files.Rdata")

Mercer_branchpoints_pluspred <- Mercer_branchpoints[,c("dist.2","set")]
G24_BPs <- G24_all[G24_all$branchpoint_prob >= 0.5 & G24_all$in_testtrain == 0, c("to_3prime", "in_testtrain")]
G24_BPs$in_testtrain <- "pred"
colnames(G24_BPs) <- colnames(Mercer_branchpoints_pluspred)
Mercer_branchpoints_pluspred <- rbind(Mercer_branchpoints_pluspred, G24_BPs)
Mercer_branchpoints_pluspred$set <- factor(Mercer_branchpoints_pluspred$set, levels=c("pred","LC","HC"))

Figure3A <- ggplot(Mercer_branchpoints_pluspred[Mercer_branchpoints_pluspred$set != "LC",], 
                   aes(x=dist.2*-1, fill=set)) + 
  geom_vline(xintercept = quantile(Mercer_branchpoints_hc$dist.2, c(0.05,.95))*-1, 
             color="grey",linetype = "dashed") +
  geom_bar() + 
  scale_fill_manual(values=nt_cols[c(2,4)] , name=element_blank(), 
                    labels=c("Predicted",
                             "Mercer et al.")) + 
  scale_x_continuous(limits=c(-100,0), 
                     labels=c("100","75","50","25","0"), 
                     name="Distance to 3' exon") +
  theme_figure + 
  theme(legend.position=c(0.25,0.8))

###### Figure 2C&D - Branchpoint motifs ######

Figure2C=ggplot(fivemer_summary, aes(x=percent_BP, y=median_branchpoint_prob, size=num_BP, col=BP_nt)) + 
  geom_point() +scale_color_manual(values = nt_cols, name="BP nt")+
  scale_x_continuous(name="Relative motif frequency\n(BP/Negative)") + 
  scale_y_continuous(name="Median BP Probability Score") + 
  theme_figure +
  scale_size_continuous(name="BP count", range=c(0.5,3))

Figure2D=ggplot(fivemer_summary, aes(y=median_branchpoint_prob, x=mean_U2,size=num_BP, col=BP_nt)) + geom_point() +
  scale_color_manual(values = nt_cols, name="BP nt")+
  scale_y_continuous(name="Median BP Probability Score") + 
  scale_x_continuous(name="Mean U2 binding energy") +
  theme_figure +
  scale_size_continuous(name="BP count", range=c(0.5,3))

Figure2_bot = ggdraw() + 
  draw_plot(Figure2C, 0,0,0.6,1) + 
  draw_plot(Figure2D + theme(legend.position = "none"), 0.6,0,0.4,1) + 
  draw_plot_label(c("C","D"), c(0,0.6), c(1,1), size=18)

pdf("Figures/Figure2_bot.pdf", useDingbats = F, height=2.7, width=6.69)
Figure2_bot
dev.off()


branchpoints_introns$gene_biotype_broad2 <- branchpoints_introns$gene_biotype_broad
branchpoints_introns$gene_biotype_broad2[grep("unprocessed_pseudogene", branchpoints_introns$gene_biotype)] <- "unprocessed_pseudogene"
branchpoints_introns$gene_biotype_broad2[grepl("processed_pseudogene", branchpoints_introns$gene_biotype) & 
                                           branchpoints_introns$gene_biotype_broad2 !="unprocessed_pseudogene"] <- "processed_pseudogene"

number_bp_by_gene_type=as.data.frame(table(branchpoints_introns$gene_biotype_broad2, branchpoints_introns$predicted_BPs_factor)/rowSums(table(branchpoints_introns$gene_biotype_broad2, branchpoints_introns$predicted_BPs_factor)))

Figure3C <- ggplot(number_bp_by_gene_type[number_bp_by_gene_type$Var1 %in% names(which(table(branchpoints_introns$gene_biotype_broad2) > 2000)),], 
                   aes(x=Var1, y=Freq, fill=Var2)) + 
  scale_fill_manual(values = BP_multi_colors, name="Branchpoint\ncount") + 
  geom_bar(stat="identity") + theme_figure +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(name="Gene Biotype") + 
  scale_y_continuous(name="Fraction of Total Introns") + ggtitle(label = "Branchpoint Count")

##### Conservation
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
m <- match(conservation_df$id, branchpoints_introns$introns)
conservation_df$status <- branchpoints_introns$status[m]
conservation_df$annotated_BPs <- branchpoints_introns$annotated_BPs_factor[m]
conservation_df$gene_biotype_broad <- branchpoints_introns$gene_biotype_broad[m]

melted <- melt(conservation_df[,c(1:12,20:22)], id.vars=c("id","status","annotated_BPs","gene_biotype_broad"))
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
Figure3D = ggplot(conservation_summary, 
       aes(x=as.factor(variable), col=set, group=set, y=value, ymin=se_lower, ymax=se_upper)) + 
  geom_point() + geom_line() +
  geom_errorbar(width=0.25) +
  facet_wrap(~multiplicity) +
  scale_color_manual(values=nt_cols[c(4,2)] , name=element_blank(), 
                    labels=c("Mercer et al.","Predicted")) +
  theme_figure + theme(legend.position = c(0.8, 0.7),strip.text.x = element_blank()) + 
  scale_y_continuous(name="phyloP conservation") + 
  scale_x_discrete(name="Distance to Branchpoint") + 
  ggtitle("Single BP Multi BP")

###### Figure 2C&D - Prediction of splicing branchpoints in GENCODE introns ######

#want number of introns annotated by mercer capture for the trest/train data
branchpoints_introns$bps <- branchpoints_introns$annotated_BPs
branchpoints_introns$bps[branchpoints_introns$status == 'test/train'] <- 
  branchpoints_introns$known_BPs[branchpoints_introns$status == 'test/train']
branchpoints_introns$bps[branchpoints_introns$bps >= 4] ="4+"

FigureS4=ggplot(branchpoints_introns, aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black") + 
  scale_fill_manual(values = BP_multi_colors, name="Branchpoint\ncount") + 
  theme_figure +
  theme(legend.position = c(0.8,0.7)) + 
  scale_y_continuous(name="Introns") +
  scale_x_discrete(labels=c("Predicted","Mercer et al.","No Annotation"), name=element_blank())

pdf("Figures/FigureS4.pdf", width=3.31, height=3.3)
FigureS4
dev.off()

Figure3B=ggplot(branchpoints_introns[branchpoints_introns$status!="unknown",], 
                aes(x=dex_mean_log, fill=status)) + 
  geom_density(alpha=0.5) +
  theme_figure +
  scale_fill_manual(values=nt_cols[c(2,4)] , name=element_blank(), 
                    labels=c("Predicted",
                             "Mercer et al.")) + 
  theme(legend.position = c(0.8,0.8)) + 
  scale_x_continuous(name="Exon expression (log10)")

Figure3=
    ggdraw() + 
    draw_plot(Figure3A, 0,0.5,0.6,.5) + 
    draw_plot(Figure3B, 0.6,0.5,0.4,.5) + 
    draw_plot(Figure3C, 0,0,0.33,0.5) + 
    draw_plot(Figure3D, 0.33,0,0.64,0.5) + 
    draw_plot_label(c("A","B","C","D"), c(0,0.6,0,0.33), c(1,1,.5,.5), size=18)

pdf("Figures/Figure3.pdf", useDingbats = F, height=5, width=6.69)
Figure3
dev.off()

FigureS7=ggplot(branchpoints_introns[branchpoints_introns$gene_biotype_broad != "other",], aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black") + 
  facet_wrap(~gene_biotype_broad, scales = "free")+ 
  scale_fill_manual(values = BP_multi_colors, name="Branchpoint\ncount") +
  theme_figure + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(drop=FALSE,labels=c("Predicted","Mercer et al.","No Annotation"), name=element_blank())

pdf("Figures/FigureS7.pdf", 6.69,3,useDingbats = F)
FigureS7
dev.off()


###### Figure S7,  3 - Features of introns with branchpoint multiplicity ######

quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "intron_size")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted" & branchpoints_introns$gene_biotype == "protein_coding",], "dex_mean_log")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "cons_intron_mean")

max_U2 <- aggregate(U2_binding_energy ~ id, G24_all, max)
branchpoints_introns$max_U2 <- max_U2$U2_binding_energy[match(branchpoints_introns$introns, max_U2$id)]

quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "max_U2")

aggregate(max_U2 ~ gene_biotype_broad, branchpoints_introns[branchpoints_introns$predicted_BPs==2,], median)

wilcox.test(branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs!=0 & 
                                          branchpoints_introns$gene_biotype_broad == "lncRNA"],
            branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs!=0 & 
                                          branchpoints_introns$gene_biotype_broad == "protein_coding"])

FigureS7A=ggplot(branchpoints_introns[branchpoints_introns$status=="predicted",], 
                       aes(predicted_BPs_factor,intron_size, fill=factor(predicted_BPs_factor))) + 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) +
  scale_y_log10(name="Intron size (nt)")+ 
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none",text=element_text(size=10))+ scale_x_discrete(name="Annotated branchpoints")

Figure4A=ggplot(branchpoints_introns[branchpoints_introns$status=="predicted",], 
       aes(x=intron_size, col=factor(predicted_BPs_factor))) + 
  stat_ecdf(geom="step") +
  scale_x_log10(name="Intron size (nt)", limits=c(50, 500000))+ 
  scale_color_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure + 
  theme(legend.position = "none") +
  scale_y_continuous(name="Cumulative Fraction") + ggtitle("Shorter Introns")

summary <- Rmisc::summarySE(branchpoints_introns[branchpoints_introns$status=="predicted",], 
                            measurevar="intron_size", 
                            groupvars=c("predicted_BPs_factor"))

summary$median <- aggregate(intron_size ~ predicted_BPs_factor, branchpoints_introns[branchpoints_introns$status=="predicted",], median)[,2]

Figure4A_subset=ggplot(summary, aes(x=predicted_BPs_factor, y=median,fill=predicted_BPs_factor)) +
  geom_bar(stat="identity") + geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2) +
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure +
  scale_y_continuous(name="Intron size") + 
  scale_x_discrete(name="BP per intron") + 
  theme(legend.position = "none")

Figure4A_all <- ggdraw() + draw_plot(Figure4A, 0,0,1,1) + 
  draw_plot(Figure4A_subset, 0.5,0.1,0.45,.7)

FigureS7A=ggplot(branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & 
                                      branchpoints_introns$status=="predicted",], 
               aes(x=dex_mean_log, col=factor(predicted_BPs_factor))) + 
  stat_ecdf(geom="step") +
  scale_x_continuous(name="Exon Expression (log10)")+ 
  scale_color_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure +
  theme(legend.position = "none") +
  scale_y_continuous(name="Cumulative Fraction") + ggtitle("Higher Expression")

summary <- Rmisc::summarySE(branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & 
                                                   branchpoints_introns$status=="predicted",], 
                            measurevar="dex_mean_log", 
                            groupvars=c("predicted_BPs_factor"))

summary$median <- aggregate(dex_mean_log ~ predicted_BPs_factor, branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & 
                                                                                       branchpoints_introns$status=="predicted",], median)[,2]

FigureS7A_subset = ggplot(summary, aes(x=predicted_BPs_factor, y=median,fill=predicted_BPs_factor)) +
  geom_bar(stat="identity", color="black") + geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2) +
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure +
  scale_y_continuous(name="Exon Expression (log10)") + 
  scale_x_discrete(name="BP per intron") + 
  theme(legend.position = "none")

FigureS7B=ggplot(branchpoints_introns[which(branchpoints_introns$status=="predicted" & !
                                            is.na(branchpoints_introns$cons_intron_mean)),], 
               aes(x=cons_intron_mean, col=factor(predicted_BPs_factor))) + 
  stat_ecdf(geom="step") +
  scale_x_continuous(name="Mean Intron Conservation")+ 
  scale_color_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure + 
  theme(legend.position = "none") +
  scale_y_continuous(name="Cumulative Fraction") + ggtitle("Higher Conservation")

summary <- Rmisc::summarySE(branchpoints_introns[which(branchpoints_introns$status=="predicted" & 
                                                         !is.na(branchpoints_introns$cons_intron_mean)),], 
                            measurevar="cons_intron_mean", 
                            groupvars=c("predicted_BPs_factor"))

summary$median <- aggregate(cons_intron_mean ~ predicted_BPs_factor, branchpoints_introns[which(branchpoints_introns$status=="predicted" & 
                                                                                              !is.na(branchpoints_introns$cons_intron_mean)),], median)[,2]

FigureS7B_subset =ggplot(summary, aes(x=predicted_BPs_factor, y=median,fill=predicted_BPs_factor)) +
  geom_bar(stat="identity", color="black") + geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2) +
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure + 
  scale_y_continuous(name="Mean Intron Conservation") + 
  scale_x_discrete(name="BP per intron") + 
  theme(legend.position = "none")

FigureS7A_all <- ggdraw() + draw_plot(FigureS7A, 0,0,1,1) + 
  draw_plot(FigureS7A_subset, 0.6,0.15,0.35,.6)

FigureS7B_all <- ggdraw() + draw_plot(FigureS7B, 0,0,1,1) + 
  draw_plot(FigureS7B_subset, 0.6,0.15,0.35,.6)

FigureS7 <- ggdraw() + 
  draw_plot(FigureS7A_all, 0.00,0, 0.5,1) + 
  draw_plot(FigureS7B_all, 0.5,0, 0.5,1) + 
  draw_plot_label(c("A","B"), c(0,0.5), c(1,1), size=18)

pdf("Figures/FigureS7.pdf", height=3,width=6.69, useDingbats = FALSE)
FigureS7
dev.off()

#no significant difference in multiplicity between gene biotype class
tbl <- table(branchpoints_introns$gene_biotype_broad[branchpoints_introns$gene_biotype_broad !="other" & branchpoints_introns$status == "predicted"], 
      branchpoints_introns$bps[branchpoints_introns$gene_biotype_broad !="other" & branchpoints_introns$status == "predicted"])
chisq.test(tbl)

#no significant difference in multiplicity in CDS/UTR flanked introns
#to use chi.sq assumption that expected counts > 5 need to filter out intron types and use "multiplicity", not number of BPs.
branchpoints_introns$intron_type <- paste(branchpoints_introns$exon1_contains,branchpoints_introns$exon2_contains, sep="-")
branchpoints_introns$multiplicity <- branchpoints_introns$bps
branchpoints_introns$multiplicity[branchpoints_introns$bps %in% c(2,3,"4+")] <- "2+"
tbl <- table(branchpoints_introns$intron_type[!is.na(branchpoints_introns$exon1_contains) & !is.na(branchpoints_introns$exon2_contains) & 
                                                branchpoints_introns$exon1_contains != "UTR3+5" & branchpoints_introns$exon2_contains != "UTR3+5" &
                                                ((branchpoints_introns$exon1_contains == "CDS" | branchpoints_introns$exon2_contains == "CDS") |
                                                branchpoints_introns$exon1_contains == branchpoints_introns$exon2_contains) &
                                                branchpoints_introns$status == "predicted"], 
             branchpoints_introns$multiplicity[!is.na(branchpoints_introns$exon1_contains) & !is.na(branchpoints_introns$exon2_contains) & 
                                        branchpoints_introns$exon1_contains != "UTR3+5" & branchpoints_introns$exon2_contains != "UTR3+5" &
                                        ((branchpoints_introns$exon1_contains == "CDS" | branchpoints_introns$exon2_contains == "CDS") |
                                        branchpoints_introns$exon1_contains == branchpoints_introns$exon2_contains) &
                                        branchpoints_introns$status == "predicted"])
chisq.test(tbl)

#difference in multiplicity for extensively alternatively spliced introns
#by gencode annot
branchpoints_introns$gencode_alternative5_factor=branchpoints_introns$gencode_alternative5
branchpoints_introns$gencode_alternative5_factor[branchpoints_introns$gencode_alternative5 > 6] ="7+"
tbl <- table(branchpoints_introns$gencode_alternative5_factor[branchpoints_introns$status == "predicted"], 
             branchpoints_introns$multiplicity[branchpoints_introns$status == "predicted"])
chisq.test(tbl)

norm_tabl <- tbl
norm_tabl[,1] <- tbl[,1]/sum(tbl[,1]) *10000
norm_tabl[,2] <- tbl[,2]/sum(tbl[,2]) *10000

tbl_pheat <- norm_tabl
tbl_pheat[,1] <- norm_tabl[,1]/rowSums(norm_tabl)
tbl_pheat[,2] <- norm_tabl[,2]/rowSums(norm_tabl)

pdf("Figures/Figure4C.pdf", height=2,width=2, useDingbats = FALSE)
pheatmap(tbl_pheat, cluster_cols = F,cluster_rows = F,annotation_names_row=T, display_numbers = tbl)
dev.off()

#and SJ reads
branchpoints_introns$`alternative5'annotated_factor`=branchpoints_introns$`alternative5'annotated`
branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$`alternative5'annotated` > 4] ="5+"
tbl <- table(branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$status == "predicted"], 
             branchpoints_introns$multiplicity[branchpoints_introns$status == "predicted"])
chisq.test(tbl)

wilcox.test((branchpoints_introns$var_maxentscan[branchpoints_introns$predicted_BPs > 1 & branchpoints_introns$gencode_alternative5 > 1]),
            (branchpoints_introns$var_maxentscan[branchpoints_introns$predicted_BPs == 1 & branchpoints_introns$gencode_alternative5 > 1]))

#no association between number of alternative 5' exons and spacing of multiple branchpoints

branchpoints_introns$min_dist_diff_factor <- branchpoints_introns$min_dist_diff
branchpoints_introns$min_dist_diff_factor[which(branchpoints_introns$min_dist_diff < 10 & !is.na(branchpoints_introns$min_dist_diff))] <- "close"
branchpoints_introns$min_dist_diff_factor[which(branchpoints_introns$min_dist_diff >= 10 & !is.na(branchpoints_introns$min_dist_diff))] <- "far"

tbl <- table(as.character(branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$predicted_BPs==2]), 
             branchpoints_introns$min_dist_diff_factor[branchpoints_introns$predicted_BPs==2])
chisq.test(tbl)

#locations of common variants (GTEX)
load(file="data/commonVariants.Rdata")

snp_locs$Variant_location <- gsub("single BP", "1 ", snp_locs$Variant_location)
snp_locs$Variant_location <- gsub("multiple BPs", "2+", snp_locs$Variant_location)
snp_locs$Var2 <- gsub("multiple BPs", "2+", snp_locs$Var2)
snp_locs$Var2 <- gsub("single BP", "1 ", snp_locs$Var2)

tbl <- table(processed_GTEX$snp_in[processed_GTEX$multi_BP!="0" & !is.na(processed_GTEX$snp_in)],
      processed_GTEX$multi_BP[processed_GTEX$multi_BP!="0" & !is.na(processed_GTEX$snp_in)])

c <- chisq.test(tbl)
c$p.value

snp_locs$Variant_location_reorder <- factor(snp_locs$Variant_location,levels=c("BP to 5'SS | 1 ",
                                                                                  "BP to 5'SS | 2+",
                                                                                  "BP | 1 ",
                                                                                  "BP | 2+",
                                                                                  "3'SS to BP | 1 ",
                                                                                  "3'SS to BP | 2+",
                                                                                  "3'SS | 1 ",
                                                                                  "3'SS | 2+"))

Figure4B=ggplot(snp_locs, aes(x=Variant_location_reorder , fill=Var2, y=percent_Freq)) + 
    geom_bar(stat="identity",position="dodge") +  
    geom_text(aes(x=Variant_location_reorder , y=percent_Freq + 0.05 * sign(percent_Freq), label=Freq), hjust=0.5, size=2) +
    theme_figure +
    theme(legend.position="none") + 
    scale_fill_manual(values = BP_multi_colors[c(2,5,4,3)],name="Annotated branchpoints",labels=c("1","2+")) + 
    scale_x_discrete(name="GTEx Variant Location", labels=rep(c("5'SS to BP", "BP","3'SS to BP","3'SS"), each=2)) + 
    scale_y_continuous(name="Relative Common\nSNP Frequency") + ggtitle("Enrichment of Common Variants at BPs")

Figure4_top <- ggdraw() + 
    draw_plot(Figure4A_all, 0,0.2, 0.4,0.8) + 
    draw_plot(Figure4B, 0.4,0, 0.6,1) + 
    draw_plot_label(c("A","B"), c(0,0.4), c(1,1), size=18)

pdf("Figures/Figure4_top.pdf", height=3,width=6.69, useDingbats = FALSE)
Figure4_top
dev.off()

####### Figure 4 Supplement #######

#Intron size trend present in Mercer annotation
#non-annotated introns are longer

Mercer_branchpoints$intron_size <- Mercer_branchpoints$dist.1 + Mercer_branchpoints$dist.2
#all
Mercer_introns <- as.data.frame(table(Mercer_branchpoints$exon_id.2[Mercer_branchpoints$set=="HC"]))
Mercer_introns$size <- Mercer_branchpoints$intron_size[match(Mercer_introns$Var1, Mercer_branchpoints$exon_id.2)]
Mercer_introns$Freq[Mercer_introns$Freq >= 4] <- "4+"

FigureS6=ggplot(Mercer_introns, aes(x=size, col=Freq)) + 
    stat_ecdf(geom="step") +
    scale_x_log10(name="Intron size (nt)") + 
    scale_color_manual(values = BP_multi_colors[-1], name="Annotated branchpoints") + 
    theme_figure +
    theme(legend.position = "none") +
    scale_y_continuous(name="Cumulative Fraction")

summary <- Rmisc::summarySE(Mercer_introns, 
                            measurevar="size", 
                            groupvars=c("Freq"))

summary$median <- aggregate(size ~ Freq, Mercer_introns, median)[,2]

FigureS6_subset=ggplot(summary, aes(x=Freq, y=median,fill=Freq)) +
  geom_bar(stat="identity") + geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2) +
  scale_fill_manual(values = BP_multi_colors[-1], name="Annotated branchpoints") + 
  theme_figure +
  scale_y_continuous(name="Intron size") + 
  scale_x_discrete(name="BP per intron") + 
  theme(legend.position = "none")

FigureS6_all <- ggdraw() + draw_plot(FigureS5, 0,0,1,1) + 
  draw_plot(FigureS6_subset, 0.6,0.15,0.35,.6)
pdf("Figures/FigureS6.pdf", height=3, width=3.31)
FigureS6_all
dev.off()


wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 1], 
            Mercer_introns$size[Mercer_introns$Freq == 2])

wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 2], 
            Mercer_introns$size[Mercer_introns$Freq == 3])

wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 3], 
            Mercer_introns$size[Mercer_introns$Freq == "4+"])

####### Supplementary Figure 8 -- exon skipping #######
load("data/exon_skipping.Rdata")

FigureS8 <- ggplot(splicing_strength, aes(x=variable, lower=lower, upper=upper, 
                                middle=middle, ymin=ymin, ymax=ymax, 
                                fill=condition1)) + 
    geom_boxplot(stat="identity") + 
    facet_grid(.~condition2) + 
    scale_fill_manual(values = nt_cols[c(2,4)], name="exon set") + 
    scale_y_continuous(name="element strength") + 
    #scale_fill_manual(values = splicingcols, name="MISO annotation") + 
    theme_bw() + 
    theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches"),axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_x_discrete(name="splicing element")

pdf("Figures/FigureS8.pdf",useDingbats = F, height=3.5, width=6.69)
FigureS8
dev.off()

wilcox.test(branchpoints_introns_sa$splicesite_5strength[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="E1"],
            branchpoints_introns_sa$splicesite_5strength[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="E1"] )
wilcox.test(branchpoints_introns_sa$splicesite_5strength[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="SE"],
            branchpoints_introns_sa$splicesite_5strength[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="SE"] )

wilcox.test(branchpoints_introns_sa$splicesite_3strength[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="E2"],
            branchpoints_introns_sa$splicesite_3strength[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="E2"] )
wilcox.test(branchpoints_introns_sa$splicesite_3strength[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="SE"],
            branchpoints_introns_sa$splicesite_3strength[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="SE"] )

wilcox.test(branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="E2"],
            branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="E2"] )
wilcox.test(branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="skipped_exon" & branchpoints_introns_sa$sa2=="SE"],
            branchpoints_introns_sa$top_U2[branchpoints_introns_sa$sa1=="control" & branchpoints_introns_sa$sa2=="SE"] )


#no significant differnece in number of branchpoints
tbl <- table(branchpoints_introns_sa$sa1[branchpoints_introns_sa$status == "predicted" & branchpoints_introns_sa$sa2=="SE"], 
             branchpoints_introns_sa$annotated_BPs_factor[branchpoints_introns_sa$status == "predicted"& branchpoints_introns_sa$sa2=="SE"])
chisq.test(tbl)


###### Figure S9, 4 - Disease causing variants affecting BPs #######
#locations of clinVar intronic variants and GTEX variants

load("data/diseaseVariants.Rdata")
load("data/commonVariants.Rdata")

clinVar_processed$snp_in <- NA
clinVar_processed$snp_in[clinVar_processed$dist_to_BP <= -1] <- "3'SS to BP"
clinVar_processed$snp_in[clinVar_processed$dist_to_BP >= 3] <- "BP to 5'SS"
clinVar_processed$snp_in[clinVar_processed$dist_to_BP == 0 | clinVar_processed$dist_to_BP == 0 | clinVar_processed$dist_to_BP == 0] <- "BP"
clinVar_processed$snp_in[clinVar_processed$dist_to_exon < 3] <- "3'SS"

clinVar_processed$set <- "ClinVar"
processed_GTEX$set <- "GTEx"
processed_GTEX$exon_3prime=NULL

snp_locs_all <- rbind(processed_GTEX, clinVar_processed[,match(colnames(processed_GTEX), colnames(clinVar_processed))])

Figure5A=ggplot(snp_locs_all[snp_locs_all$multi_BP!="0" & !is.na(snp_locs_all$snp_in),], 
       aes(x=snp_in, fill=multi_BP)) + 
    geom_bar() +  
    theme_figure +
    theme(strip.text.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("ClinVar GTEx") +
    scale_fill_manual(name="Branchpoints\nper Intron",values = BP_multi_colors[c(2,5)]) +
    facet_wrap(~set, scales="free") + 
      scale_x_discrete(name=element_blank(), labels=c("3'SS", "3'SS to BP", "BP","5'SS to BP"))

pdf("Figures/Figure5A.pdf", width=6.69, height=3)
Figure5A
dev.off()

chisq.test(table(clinVar_processed$snp_in[clinVar_processed$multi_BP!=0], clinVar_processed$multi_BP[clinVar_processed$multi_BP!=0]))

chisq.test(table(processed_GTEX$snp_in[processed_GTEX$multi_BP!=0], processed_GTEX$multi_BP[processed_GTEX$multi_BP!=0]))

source("scripts/analysis/variation/plotBranchpointWindow_Figure.R")

#FECH branchpoint mutation
pdf("Figures/Figure5_bot.pdf", width=6.69, height=6)
plotBranchpointWindow_Figure(clinVars_filtered$id[clinVars_filtered$GeneSymbol=="FECH"], clinvar_predictions,clinvar_attributes,exons)
dev.off()
