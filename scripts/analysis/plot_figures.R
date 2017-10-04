################################################################################
#           Genome wide prediction of splicing branchpoints - plots            #
################################################################################

# Produces Figures 1 - 4
options(stringsAsFactors = F)

library(branchpointer)
library(ggplot2)
library(cowplot)
library(stringr)
library(pheatmap)
library(data.table)
library(grid)

source("~/Documents/Projects/splice_branchpoints/scripts/analysis/quick_wilcox.R")
#source("scripts/analysis/quick_wilcox.R")

theme_figure <- theme_bw()+ theme(text=element_text(size=6),legend.key.size=unit(0.1, "inches"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_rect(colour = "black", size=0.5, fill=NA),
                                  axis.line.x = element_line(), 
                                  axis.line.y = element_line(),
                                  axis.ticks = element_line(size=0.25),
                                  #legend.position = "none",
                                  axis.ticks.length = unit(0.05, "cm"),
                                  panel.background = element_blank(),
                                  plot.title = element_text(hjust = 0.5))

number_introns_pal=rev(c('#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'))
nt_cols=c("#359646","#4D7ABE","#FAA859","#CB3634")
BP_multi_colors=c("#bababa", "#1a9850", "#a6d96a", "#ffffbf","#fdae61")
BP_multi_colors=c("#bababa", "#deebf7", "#9ecae1", "#4292c6","#08519c")

ggplot_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

cutoff = 0.48

# Figure 1. Development and performance of the branchpointer model

load("data/performance_objects.Rdata")
load("data/Figure_files.Rdata")

# variable importance in GBM model
Figure1B <- ggplot(gbm_importance[!duplicated(gbm_importance$variableLO),], 
                   aes(x=Overall, y=meanF1Change*-1, col=color, shape=grouped, label=variableName)) + 
  geom_point(size=0.5) + geom_text(hjust = 0, nudge_x = 0.05, size=2) +
  scale_color_manual(values=c(nt_cols[c(4,1,2,3)])) +
  scale_x_log10(name="GBM variable importance") + scale_y_log10(name="Change in F1") +
  theme_figure + theme(legend.position = c(0.8,0.2),legend.title = element_blank()) + 
  ggtitle(label="Variable Importance")

Figure1C <- ggplot(U2_df[abs(U2_df$pos) < 4,], aes(x=pos, y=importanceScaled, fill=nt)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_x_continuous(name="Distance to Branchpoint", breaks = seq(-3,3,1), labels=c("-3","-2","-1","BP","1","2","3")) + 
  scale_y_continuous(name="Feature Importance\nin Model", breaks=0:5, labels=c("0.00", "1.00","2.00","3.00","4.00","5.00")) + 
  scale_fill_manual(values=c(nt_cols)) + 
  theme_figure + theme(legend.position = "none")

Figure1 <- ggdraw() +
  draw_plot(Figure1B, 0,0.4,1,0.35) +
  draw_plot(Figure1C, 0.15,0,0.85,0.25) +
  draw_plot_label(c("A","B","C"), c(0,0,0), c(1,0.75,0.4), size=12)

pdf("Figures/Figure1.pdf", useDingbats = F, height=6, width=3.25)
Figure1
dev.off()

#ROC curves
Figure2A <- ggplot(roc_curves, aes(x = FPR,y = TPR, col = method)) +
  geom_line() +
  scale_color_manual(values=nt_cols) + theme_figure +
  theme(legend.position = c(0.8,0.2))  + ggtitle(label="ROC Curve")
#PR curves
Figure2B <- ggplot(pr_curves, aes(x = X1,y = X2, col = method)) +
  geom_line() +
  scale_color_manual(values=nt_cols) +
  labs(x="Recall",y="Precision") +
  theme_figure + theme(legend.position = c(0.75,0.75)) + ggtitle(label="Precision Recall")

Figure2 <- ggdraw() +
  draw_plot(Figure2A, 0,0.15,0.5,0.85) +
  draw_plot(Figure2B, 0.5,0.15,0.5,0.85) +
  draw_plot(ggplot_legend(Figure2A), 0.5,0,0.5,0.2) +
  draw_plot_label(c("A","B"), c(0,0.5), c(1,1), size=12)

pdf("Figures/Figure2.pdf", useDingbats = F, height=2, width=3.25)
Figure2
dev.off()

Figure3 <- ggplot(conservation_summary, 
                   aes(x=as.factor(variable), col=set, group=set, y=value)) + 
  geom_point(size=0.5) + geom_line(size=0.5) +
  #geom_errorbar(width=0.25) +
  facet_wrap(~multiplicity, ncol=2) +
  scale_color_manual(values=c("grey60","grey30") , name=element_blank(), 
                     labels=c("Mercer et al.","Predicted")) +
  theme_figure + theme(legend.position = c(0.8, 0.7)) + 
  scale_y_continuous(name="phyloP Conservation") + 
  scale_x_discrete(name="Distance to Branchpoint",breaks = seq(-5,5,1), labels=c("-5","-4","-3","-2","-1","BP","1","2","3","4","5"))


pdf("Figures/Figure3.pdf", useDingbats = F, height=1.5, width=3.25)
Figure3
dev.off()

Figure4A <- ggplot(fivemer_summary, aes(x=percent_BP, y=median_branchpoint_prob, size=num_BP, col=BP_nt)) + 
  geom_point(alpha=0.5) +
  scale_color_manual(values = nt_cols, name="BP nt")+
  scale_x_continuous(name="Relative Motif Frequency\n(BP/Negative)") + 
  scale_y_continuous(name="Median BP\nProbability Score") + 
  theme_figure +
  scale_size_continuous(name="BP count", range=c(0.5,3))

Figure4B <- ggplot(fivemer_summary, aes(y=median_branchpoint_prob, x=mean_U2,size=num_BP, col=BP_nt)) + 
  geom_point(alpha=0.5) +
  scale_color_manual(values = nt_cols, name="BP nt")+
  scale_y_continuous(name="Median BP\nProbability Score") + 
  scale_x_continuous(name="Mean U2 binding energy\n") +
  theme_figure +
  scale_size_continuous(name="BP count", range=c(0.5,3))


#F1 cutoffs
Figure4C <- ggplot(cutoff_performance, aes(x=vals, y=F1))  +
  #geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=0.5, shape=3)+ 
  geom_vline(xintercept = cutoff, linetype=2, color="grey") +
  scale_x_continuous(name="Branchpointer Probability Score") +
  scale_y_continuous(name="F1") + 
  theme_figure


G26_all$short_motif <- factor(paste0(stringr::str_sub(G26_all$seq_motif, start=4, end=4),
                                     "N", stringr::str_sub(G26_all$seq_motif, start=6, end=6)), levels=c(
                                       "ANA","CNA","GNA","TNA",
                                       "ANC","CNC","GNC","TNC", 
                                       "ANG","CNG","GNG","TNG", 
                                       "ANT","CNT","GNT","TNT"))

max(G26_all$branchpoint_prob[!(G26_all$short_motif %in% c("ANA","CNA","GNA","TNA"))])

Figure4D <- ggplot(G26_all[G26_all$branchpoint_prob >= 0.505,], aes(x=branchpoint_prob, fill = short_motif)) +
  geom_histogram(binwidth=0.01) +
  guides(fill = guide_legend(ncol = 2)) + 
  scale_fill_manual(values = c("#00441b","#238b45","#74c476","#c7e9c0",
                               "#08306b","#2171b5",'#6baed6',"#c6dbef",
                               "#7f2704","#d94801",'#fd8d3c',"#fdd0a2",
                               "#67000d","#cb181d","#fb6a4a","#fcbba1"), drop=FALSE, name="Motif") +
  scale_x_continuous(name="Branchpointer Probability Score") +
  scale_y_continuous(name="Count") +
  theme_figure

Figure4 <- ggdraw() +
  # draw_plot(ggplot_legend(Figure2B), 0.4,0.9,0.6,0.1) +
  draw_plot(Figure4C, 0,0,0.35,0.5) + 
  draw_plot(Figure4D, 0.35,0,0.65,0.5) +
  draw_plot(Figure4A+ theme(legend.position = "none"),0,0.5,0.45,0.5) +
  draw_plot(ggplot_legend(Figure4A), 0.45,0.5,0.1,0.5) +
  draw_plot(Figure4B+ theme(legend.position = "none"), 0.55,0.5,0.45,0.5) +
  draw_plot_label(c("A","B","C","D"), c(0,0.55,0,0.35), c(1,1,0.5,0.5), size=12)


pdf("Figures/Figure4.pdf", useDingbats = F, height=3.5, width=5)
Figure4
dev.off()

# Figure S2. Branchpointer classification performance metrics for probability cutoffs 0.01-0.99

SuppFigure2A = ggplot(cutoff_performance, aes(x=vals, y=PPV))  +
  #geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=0.5, shape=3)+ 
  scale_x_continuous(name="Branchpointer Probability Score") +
  scale_y_continuous(name="Positive Predictive Value") + 
  theme_figure

SuppFigure2B = ggplot(cutoff_performance, aes(x=vals, y=Sensitivity))  +
  #geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=0.5, shape=3)+ 
  scale_x_continuous(name="Branchpointer Probability Score") +
  scale_y_continuous(name="Sensitivity") + 
  theme_figure

SuppFigure2C=ggplot(cutoff_performance, aes(x=vals, y=Accuracy))  +
  #geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=0.5, shape=3)+ 
  scale_x_continuous(name="Branchpointer Probability Score") +
  scale_y_continuous(name="Accuracy") + 
  theme_figure

SuppFigure2D=ggplot(cutoff_performance, aes(x=vals, y=Balanced_Accuracy))  +
  #geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=0.5, shape=3)+ 
  scale_x_continuous(name="Branchpointer Probability Score") +
  scale_y_continuous(name="Balanced Accuracy") + 
  theme_figure

FigureS2 = ggdraw() + draw_plot(SuppFigure2A, 0,0.5,0.5,0.5) + 
  draw_plot(SuppFigure2B, 0.5,0.5,0.5,0.5) + 
  draw_plot(SuppFigure2C, 0,0,0.5,0.5) + 
  draw_plot(SuppFigure2D, 0.5,0,0.5,0.5) + 
  draw_plot_label(c("A","B","C","D"), c(0,0.5,0,0.5), c(1,1,0.5,0.5), size=12)

pdf("Figures/FigureS2.pdf", useDingbats = F, height=4, width=4.5)
FigureS2
dev.off()

# Figure 2. Nucleotide sequence motif importance in branchpointer model development
### Figure 3. Prediction of splicing branchpoints in GENCODE introns

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

Mercer_branchpoints_pluspred <- Mercer_branchpoints[,c("dist.2","set")]
G26_BPs <- G26_all[G26_all$branchpoint_prob >= cutoff & G26_all$in_testtrain == 0, c("to_3prime", "in_testtrain")]
G26_BPs$in_testtrain <- "pred"
colnames(G26_BPs) <- colnames(Mercer_branchpoints_pluspred)
Mercer_branchpoints_pluspred <- rbind(Mercer_branchpoints_pluspred, G26_BPs)
Mercer_branchpoints_pluspred$set <- factor(Mercer_branchpoints_pluspred$set, levels=c("pred","LC","HC"))

branchpoints_introns$gene_biotype_broad2 <- branchpoints_introns$gene_biotype_broad
branchpoints_introns$gene_biotype_broad2[grep("unprocessed_pseudogene", branchpoints_introns$gene_biotype)] <- "unprocessed_pseudogene"
branchpoints_introns$gene_biotype_broad2[grepl("processed_pseudogene", branchpoints_introns$gene_biotype) & 
                                           branchpoints_introns$gene_biotype_broad2 !="unprocessed_pseudogene"] <- "processed_pseudogene"

number_bp_by_gene_type=as.data.frame(table(branchpoints_introns$gene_biotype_broad2[branchpoints_introns$status %in% c("predicted", "unknown")], 
                                           branchpoints_introns$predicted_BPs_factor[branchpoints_introns$status %in% c("predicted", "unknown")])/
                                       rowSums(table(branchpoints_introns$gene_biotype_broad2[branchpoints_introns$status %in% c("predicted", "unknown")], 
                                                     branchpoints_introns$predicted_BPs_factor[branchpoints_introns$status %in% c("predicted", "unknown")])))

number_bp_by_gene_type$Var1 <- factor(number_bp_by_gene_type$Var1, levels = rev(c("protein_coding","lncRNA", "unprocessed_pseudogene",
                                                                              "processed_pseudogene", "pseudogene","other")))


Figure6A <- ggplot(Mercer_branchpoints_pluspred[Mercer_branchpoints_pluspred$set != "LC",], 
                   aes(x=dist.2*-1, fill=set)) + 
  geom_vline(xintercept = quantile(Mercer_branchpoints_hc$dist.2, c(0.05,.95))*-1, 
             color="grey",linetype = "dashed") +
  geom_bar(col="grey60", size=0.25) + 
  scale_fill_manual(values=c("black", "white") , name=element_blank(), 
                    labels=c("Predicted",
                             "Mercer et al.")) + 
  scale_x_continuous(limits=c(-100,0), 
                     labels=c("100","75","50","25","0"), 
                     name="Distance to 3' Exon") +
  scale_y_continuous(name="Count") +
  theme_figure + 
  theme(legend.position=c(0.25,0.8))

Figure6B <- ggplot(branchpoints_introns[branchpoints_introns$status!="unknown",], 
                aes(x=exon_group_log10count, fill=status, linetype=status)) + 
  geom_density(alpha=0.5, size=0.25) +
  theme_figure +
  scale_fill_manual(values=c("grey60", NA) , name=element_blank(), 
                    labels=c("Predicted",
                             "Mercer et al.")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  #theme(legend.position = c(0.7,0.8)) + 
  scale_x_continuous(name="Exon Expression (log10)") + 
  scale_y_continuous(name="Density")


Figure6C <- ggplot(number_bp_by_gene_type[number_bp_by_gene_type$Var1 %in% names(which(table(branchpoints_introns$gene_biotype_broad2) > 2000)),], 
                   aes(x=Var1, y=Freq, fill=Var2, linetype=Var2)) + 
  scale_fill_manual(values = c("white","grey90", "grey40"), name="Branchpoint\ncount") + 
  scale_linetype_manual(values=c("solid", "dashed","solid")) +
  geom_bar(stat="identity", col="black", size=0.25) + 
  coord_flip() + 
  theme_figure +
  theme(legend.position = "none") +
  scale_x_discrete(name="Gene Biotype") + 
  scale_y_continuous(name="Fraction of Total Introns")

#dashed linetype here goes **super weird**
Figure6D <- ggplot(branchpoints_introns[branchpoints_introns$status=="predicted",], 
                     aes(x=intron_size, linetype=factor(predicted_BPs_factor))) + 
  stat_ecdf(geom="step") +
  scale_x_log10(name="Intron Size (nt)", limits=c(50, 2000000))+ 
  scale_linetype_manual(values = c("solid", "solid"), name="Annotated branchpoints") + 
  theme_figure + 
  theme(legend.position = "none") +
  scale_y_continuous(name="Cumulative Fraction")

summary <- Rmisc::summarySE(branchpoints_introns[branchpoints_introns$status=="predicted",], 
                            measurevar="intron_size", 
                            groupvars=c("predicted_BPs_factor"))

summary$median <- aggregate(intron_size ~ predicted_BPs_factor, branchpoints_introns[branchpoints_introns$status=="predicted",], median)[,2]

Figure6D_subset=ggplot(summary, aes(x=predicted_BPs_factor, y=median,fill=predicted_BPs_factor, linetype=predicted_BPs_factor)) +
  geom_bar(stat="identity", col="black", size=0.25) + geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2, size=0.25) +
  scale_fill_manual(values = c("grey90", "grey40"), name="Annotated branchpoints") + 
  theme_figure +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_y_continuous(name="Intron Size") + 
  scale_x_discrete(name="BP per intron") + 
  theme(legend.position = "none", axis.title.x = element_blank())

Figure6D_all <- ggdraw() + draw_plot(Figure6D, 0,0,1,1) + 
  draw_plot(Figure6D_subset, 0.5,0.2,0.45,.6)

Figure6D <- Figure6D_all

Figure5 <- ggdraw() + 
  draw_plot(Figure6A, 0.0,0.5,0.5,0.5) + 
  draw_plot(Figure6B, 0.5,0.5,0.5,0.5) + 
  draw_plot(Figure6C, 0.0,0.0,0.4,0.5) + 
  draw_plot(Figure6D, 0.4,0.0,0.35,0.5) +
  draw_plot_label(c("A","B","C","D"), c(0,0.5,0,0.4), c(1,1,.5,.5), size=12)

pdf("Figures/Figure6.pdf", useDingbats = F, height=3, width=6.5)
Figure6
dev.off()


# Figure 5

load("data/branchpointer_example_figures.Rdata")

plot.prob.ref <- ggplot(rs587776767_attributes_df[rs587776767_attributes_df$status=="REF",], aes(x = to_3prime_point*-1, y = branchpoint_prob, fill = seq_pos0,
                                                                                                 col=seq_pos0, alpha=U2_binding_energy)) +
    geom_bar(stat = "identity", size=0.25) +
    scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                       breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                          "0.60","0.70","0.80","0.90","1.00")) +
    geom_hline(yintercept = 0.52, col="grey", linetype="dashed") +
    scale_fill_manual(values = nt_cols, drop = TRUE,
                      limits = c("A","C","G","T")) +
    scale_color_manual(values = nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) + theme_figure +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),axis.ticks.y = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

#U2 binding energy by sequence position
plot.U2.ref <- ggplot(rs587776767_attributes_df[rs587776767_attributes_df$status=="REF" & rs587776767_attributes_df$branchpoint_prob > 0.52,],
                      aes(x =  to_3prime_point*-1, y = U2_binding_energy)) +
    geom_bar(stat = "identity",width = 1, size=0.25) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
    geom_hline(yintercept = 5, col="grey", linetype="dashed") +
    theme_figure +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),axis.ticks.y = element_blank()) +
    scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                       breaks = seq(0,9,3), labels = c("0","3","6","9"))

#Sequence identity
plot.seq.ref <- ggplot(rs587776767_attributes_df[rs587776767_attributes_df$status=="REF",], aes(x = to_3prime_point*-1, y = 1, col = seq_pos0,label = seq_pos0)) +
    geom_text(size = 3) +
    scale_color_manual(values = nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) +
    theme_figure +
    theme(legend.position = "none",axis.text.y = element_blank(),
          axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

plot.prob.alt <- ggplot(rs587776767_attributes_df[rs587776767_attributes_df$status=="ALT",], aes(x = to_3prime_point*-1, y = branchpoint_prob, fill = seq_pos0,
                                                                                                 col=seq_pos0, alpha=U2_binding_energy)) +
    geom_bar(stat = "identity", size=0.25) +
    scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                       breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                          "0.60","0.70","0.80","0.90","1.00")) +
    geom_hline(yintercept = 0.52, col="grey", linetype="dashed") +
    scale_fill_manual(values = nt_cols, drop = TRUE,
                      limits = c("A","C","G","T")) +
    scale_color_manual(values = nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) + theme_figure +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),axis.ticks.y = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

#U2 binding energy by sequence position
plot.U2.alt <- ggplot(rs587776767_attributes_df[rs587776767_attributes_df$status=="ALT" & rs587776767_attributes_df$branchpoint_prob > 0.52,],
                      aes(x =  to_3prime_point*-1, y = U2_binding_energy)) +
    geom_bar(stat = "identity",width = 1, size=0.25) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
    geom_hline(yintercept = 5, col="grey", linetype="dashed") +
    theme_figure +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank(),axis.ticks.y = element_blank()) +
    scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                       breaks = seq(0,9,3), labels = c("0","3","6","9"))

#Sequence identity
plot.seq.alt <- ggplot(rs587776767_attributes_df[rs587776767_attributes_df$status=="ALT",], aes(x = to_3prime_point*-1, y = 1, col = seq_pos0,label = seq_pos0)) +
    geom_text(size = 3) +
    scale_color_manual(values = nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) +
    theme_figure +
    theme(legend.position = "none",axis.text.y = element_blank(),
          axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))


plot.prob.brca <- ggplot(brca2_attributes_df, aes(x = to_3prime_point*-1, y = branchpoint_prob, fill = seq_pos0,
                                                  col=seq_pos0, alpha=U2_binding_energy)) +
    geom_bar(stat = "identity", size=0.25) +
    scale_y_continuous(name = "Branchpointer\nProbability Score",limits = c(0,1),
                       breaks = seq(0,1,0.5), labels = c("0.0","0.5","1.0")) +
    geom_hline(yintercept = 0.52, col="grey", linetype="dashed") +
    scale_fill_manual(values = nt_cols, drop = TRUE,
                      limits = c("A","C","G","T")) +
    scale_color_manual(values = nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) + theme_figure +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

#U2 binding energy by sequence position
plot.U2.brca2 <- ggplot(brca2_attributes_df[brca2_attributes_df$branchpoint_prob > 0.52,],
                        aes(x =  to_3prime_point*-1, y = U2_binding_energy)) +
    geom_bar(stat = "identity",width = 1, size=0.25) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
    geom_hline(yintercept = 5, col="grey", linetype="dashed") +
    theme_figure +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
    scale_y_continuous(name = "U2 Binding\nEnergy",limits = c(0,10),
                       breaks = seq(0,9,3), labels = c("0.0","3.0","6.0","9.0"))

#Sequence identity
plot.seq.brca2 <- ggplot(brca2_attributes_df, aes(x = to_3prime_point*-1, y = 1, col = seq_pos0,label = seq_pos0)) +
    geom_text(size = 3) +
    scale_color_manual(values = nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) +
    theme_figure +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_y_continuous(name = "U2 Binding\nEnergy",limits = c(-1,3),
                       breaks = seq(-1,3,1), labels = c("0.0","0.0","3.0","6.0","9.0")) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)")


Figure5 <- ggdraw() + 
    draw_plot(plot.seq.brca2, 0.05,0.05,0.3,0.1) + 
    draw_plot(plot.prob.brca, 0.05,0.1,0.3,0.3) + 
    draw_plot(plot.U2.brca2, 0.05,0.0,0.3,0.1) +
    draw_plot(plot.seq.ref, 0.5,0.05,0.25,0.1) + 
    draw_plot(plot.prob.ref, 0.5,0.1,0.25,0.3) + 
    draw_plot(plot.U2.ref, 0.5,0.0,0.25,0.1) +
    draw_plot(plot.seq.alt, 0.75,0.05,0.25,0.1) + 
    draw_plot(plot.prob.alt, 0.75,0.1,0.25,0.3) + 
    draw_plot(plot.U2.alt, 0.75,0.0,0.25,0.1) 

pdf("Figures/Figure5.pdf", height=4, width=6.5)
Figure5
dev.off()
library(pheatmap)
pheatmap(pheat_attributes_rs587776767[,c(28:55)], cluster_rows = F,cluster_cols = F, border_color = NA, filename="Figures/pheat_top_r_b.pdf",height=1,width=1.8, fontsize=5,legend=F)
pheatmap(pheat_attributes_rs587776767[,c(1:27, 55)], cluster_rows = F,cluster_cols = F, border_color = NA, filename="Figures/pheat_top_r_a.pdf",height=1,width=1.8, fontsize=5,legend=F)
pheatmap(pheat_attributes_brca2, cluster_rows = F,cluster_cols = F, border_color = NA, filename="Figures/pheat_top_l.pdf",height=1,width=1.8, fontsize=5,legend=F)


# Features of introns with branchpoint multiplicity

quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "intron_size")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted" & branchpoints_introns$gene_biotype == "protein_coding",], "exon_group_log10count")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "cons_intron_mean")

max_U2 <- aggregate(U2_binding_energy ~ exon_id, G26_all, max)
branchpoints_introns$max_U2 <- max_U2$U2_binding_energy[match(branchpoints_introns$introns, max_U2$exon_id)]

quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "max_U2")


# lncRNAs have weaker BPs in general (not a huge difference)

aggregate(max_U2 ~ gene_biotype_broad, branchpoints_introns[branchpoints_introns$predicted_BPs==1,], summary)

wilcox.test(branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs==1 & 
                                          branchpoints_introns$gene_biotype_broad == "lncRNA"],
            branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs==1 & 
                                          branchpoints_introns$gene_biotype_broad == "protein_coding"])

aggregate(max_U2 ~ gene_biotype_broad, branchpoints_introns[branchpoints_introns$predicted_BPs>1,], summary)

wilcox.test(branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs>1 & 
                                          branchpoints_introns$gene_biotype_broad == "lncRNA"],
            branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs>1 & 
                                          branchpoints_introns$gene_biotype_broad == "protein_coding"])

# Figure S3. Human introns with annotated branchpoints from the high confidence Mercer annotation, 
# the branchpointer model (Predicted) and those with no branchpoint (No Annotation). 

#want number of introns annotated by mercer capture for the trest/train data
branchpoints_introns$bps <- branchpoints_introns$annotated_BPs
branchpoints_introns$bps[branchpoints_introns$status == 'test/train'] <- 
  branchpoints_introns$known_BPs[branchpoints_introns$status == 'test/train']
branchpoints_introns$bps[branchpoints_introns$bps >= 4] ="4+"

FigureS3=ggplot(branchpoints_introns, aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black", size=0.35) + 
  scale_fill_manual(values = BP_multi_colors, name="BP per\nIntron") + 
  theme_figure +
  theme(legend.position = c(0.8,0.7)) + 
  scale_y_continuous(name="Introns") +
  scale_x_discrete(labels=c("Predicted","Mercer et al.","No Annotation"), name=element_blank())

pdf("Figures/FigureS3.pdf", width=2, height=2)
FigureS3
dev.off()

# Figure S4. Introns with annotated branchpoints for gene biotypes.

branchpoints_introns$gene_biotype_broad2 <- factor(branchpoints_introns$gene_biotype_broad2, 
                                                   levels = c("protein_coding","lncRNA", "unprocessed_pseudogene",
                                                              "processed_pseudogene", "pseudogene","other"))

FigureS4=ggplot(branchpoints_introns[branchpoints_introns$gene_biotype_broad2 %in% 
                                       c("protein_coding","lncRNA", "unprocessed_pseudogene",
                                         "processed_pseudogene"),], aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black", size=0.25) + 
  facet_wrap(~gene_biotype_broad2, scales = "free")+ 
  scale_fill_manual(values = BP_multi_colors, name="BP per\nIntron") +
  theme_figure + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous(name="Count") + 
  scale_x_discrete(drop=FALSE,labels=c("Predicted","Mercer et al.","No Annotation"), name=element_blank())

pdf("Figures/FigureS4.pdf", width=5,height=5,useDingbats = F)
FigureS4
dev.off()

#Figure S5. Size of introns with single or multiple annotated branchpoints from the Mercer et al. branchpoint annotation

Mercer_branchpoints$intron_size <- Mercer_branchpoints$dist.1 + Mercer_branchpoints$dist.2
#all
Mercer_introns <- as.data.frame(table(Mercer_branchpoints$exon_id.2[Mercer_branchpoints$set=="HC"]))
Mercer_introns$size <- Mercer_branchpoints$intron_size[match(Mercer_introns$Var1, Mercer_branchpoints$exon_id.2)]
Mercer_introns$Freq[Mercer_introns$Freq >= 4] <- "4+"

FigureS5=ggplot(Mercer_introns, aes(x=size, col=Freq)) + 
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

FigureS5_subset=ggplot(summary, aes(x=Freq, y=median,fill=Freq)) +
  geom_bar(stat="identity", col="black", size=0.25) + 
  geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2, size=0.25) +
  scale_fill_manual(values = BP_multi_colors[-1], name="Annotated branchpoints") + 
  theme_figure +
  scale_y_continuous(name="Intron size") + 
  scale_x_discrete(name="BP per Intron") + 
  theme(legend.position = "none")

FigureS5_all <- ggdraw() + draw_plot(FigureS5, 0,0,1,1) + 
  draw_plot(FigureS5_subset, 0.6,0.2,0.35,.65)
pdf("Figures/FigureS5.pdf", height=2, width=3)
FigureS5_all
dev.off()

# size in mercer introns
wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 1], 
            Mercer_introns$size[Mercer_introns$Freq == 2])

wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 2], 
            Mercer_introns$size[Mercer_introns$Freq == 3])

wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 3], 
            Mercer_introns$size[Mercer_introns$Freq == "4+"])


# Figure S6. Features of intron usage associated with branchpoint multiplicity. 
FigureS6A=ggplot(branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & 
                                      branchpoints_introns$status=="predicted",], 
               aes(x=exon_group_log10count, col=factor(predicted_BPs_factor))) + 
  stat_ecdf(geom="step") +
  scale_x_continuous(name="Exon Expression (log10)")+ 
  scale_color_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure +
  theme(legend.position = "none") +
  scale_y_continuous(name="Cumulative Fraction") + ggtitle("Higher Expression")

summary <- Rmisc::summarySE(branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & 
                                                   branchpoints_introns$status=="predicted",], 
                            measurevar="exon_group_log10count", 
                            groupvars=c("predicted_BPs_factor"))

summary$median <- aggregate(exon_group_log10count ~ predicted_BPs_factor, branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & 
                                                                                       branchpoints_introns$status=="predicted",], median)[,2]

FigureS6A_subset = ggplot(summary, aes(x=predicted_BPs_factor, y=median,fill=predicted_BPs_factor)) +
  geom_bar(stat="identity", color="black", size=0.25) + 
  geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2, size=0.25) +
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure +
  scale_y_continuous(name="Exon Expression (log10)") + 
  scale_x_discrete(name="BP per Intron") + 
  theme(legend.position = "none")

FigureS6B=ggplot(branchpoints_introns[which(branchpoints_introns$status=="predicted" & !
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

FigureS6B_subset =ggplot(summary, aes(x=predicted_BPs_factor, y=median,fill=predicted_BPs_factor)) +
  geom_bar(stat="identity", color="black", size=0.25) + geom_errorbar(aes(ymin=median-se, ymax=median+se), width=0.2, size=0.25) +
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_figure + 
  scale_y_continuous(name="Mean Intron Conservation") + 
  scale_x_discrete(name="BP per Intron") + 
  theme(legend.position = "none")

FigureS6A_all <- ggdraw() + draw_plot(FigureS6A, 0,0,1,1) + 
  draw_plot(FigureS6A_subset, 0.6,0.15,0.35,.6)

FigureS6B_all <- ggdraw() + draw_plot(FigureS6B, 0,0,1,1) + 
  draw_plot(FigureS6B_subset, 0.6,0.15,0.35,.6)

FigureS6 <- ggdraw() + 
  draw_plot(FigureS6A_all, 0.00,0, 0.5,1) + 
  draw_plot(FigureS6B_all, 0.5,0, 0.5,1) + 
  draw_plot_label(c("A","B"), c(0,0.5), c(1,1), size=12)

pdf("Figures/FigureS6.pdf", height=2.5,width=5, useDingbats = FALSE)
FigureS6
dev.off()

#lncrnas have less bps per intron -- can be explained by expression
branchpoints_introns$multiplicity <- branchpoints_introns$bps
branchpoints_introns$multiplicity[branchpoints_introns$bps %in% c(2,3,"4+")] <- "2+"

tbl <- table(branchpoints_introns$gene_biotype_broad[branchpoints_introns$gene_biotype_broad !="other" & branchpoints_introns$status == "predicted"], 
      branchpoints_introns$multiplicity[branchpoints_introns$gene_biotype_broad !="other" & branchpoints_introns$status == "predicted"])

chisq.test(tbl)
tbl/rowSums(tbl)

tbl <- table(branchpoints_introns$gene_biotype_broad[branchpoints_introns$exon_group_log10count > 2 & branchpoints_introns$gene_biotype_broad !="other" & branchpoints_introns$status == "predicted"], 
             branchpoints_introns$multiplicity[branchpoints_introns$exon_group_log10count > 2 & branchpoints_introns$gene_biotype_broad !="other" & branchpoints_introns$status == "predicted"])
chisq.test(tbl)
tbl/rowSums(tbl)

# Figure S7. Frequencies and locations of intronic GTEx and ClinVar SNPs

# locations of common snps within BPs
load(file="data/commonVariants.Rdata")
snp_locs$Variant_location <- gsub("single BP", "1 ", snp_locs$Variant_location)
snp_locs$Variant_location <- gsub("multiple BPs", "2+", snp_locs$Variant_location)
snp_locs$Var2 <- gsub("multiple BPs", "2+", snp_locs$Var2)
snp_locs$Var2 <- gsub("single BP", "1 ", snp_locs$Var2)

tbl <- table(gtex_summary$snp_in[gtex_summary$multi_BP!="0" & !is.na(gtex_summary$snp_in)],
             gtex_summary$multi_BP[gtex_summary$multi_BP!="0" & !is.na(gtex_summary$snp_in)])

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

FigureS7 <- ggplot(snp_locs, aes(x=Variant_location_reorder , fill=Var2, y=percent_Freq)) + 
  geom_bar(stat="identity",position="dodge", col="black", size=0.25) +  
  geom_text(aes(x=Variant_location_reorder , y=percent_Freq + 0.05 * sign(percent_Freq), label=Freq), hjust=0.5, size=2) +
  theme_figure +
  theme(legend.position=c(0.9,0.9)) + 
  scale_fill_manual(values = BP_multi_colors[c(2,5,4,3)],name="BP per Intron",labels=c("1","2+")) + 
  scale_x_discrete(name="GTEx Variant Location", labels=rep(c("5'SS to BP", "BP","3'SS to BP","3'SS"), each=2)) + 
  scale_y_continuous(name="Relative Common SNP Frequency") + ggtitle("Enrichment of Common Variants at BPs")

pdf("Figures/FigureS7.pdf", height=3, width=5)
FigureS7
dev.off()

load("data/diseaseVariants.Rdata")

clinvar_summary$set <- "ClinVar"
gtex_summary$set <- "GTEX"

#gtex_summary$exon_3prime=NULL
gtex_summary <- as.data.frame(gtex_summary)
clinvar_summary <- as.data.frame(clinvar_summary)
cols <- colnames(clinvar_summary)[which(colnames(clinvar_summary) %in% colnames(gtex_summary))]
snp_locs_all <- rbind(clinvar_summary[,cols],gtex_summary[,cols])

FigureS8=ggplot(snp_locs_all[snp_locs_all$multi_BP!="0" & !is.na(snp_locs_all$snp_in),], 
       aes(x=snp_in, fill=multi_BP)) + 
    geom_bar(col="black", size=0.25) +  
    theme_figure +
    theme(strip.text.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("ClinVar GTEx") +
    scale_fill_manual(name="BP\nper Intron",values = BP_multi_colors[c(2,5)]) +
    facet_wrap(~set, scales="free") + 
      scale_x_discrete(name=element_blank(), labels=c("3'SS", "3'SS to BP", "BP","5'SS to BP")) +
  scale_y_continuous(name="Count")

pdf("Figures/FigureS8.pdf", width=5, height=2.5)
FigureS8
dev.off()

# different proportions are at multi/single BP introns
chisq.test(table(clinvar_summary$snp_in[clinvar_summary$multi_BP!=0], clinvar_summary$multi_BP[clinvar_summary$multi_BP!=0]))
chisq.test(table(gtex_summary$snp_in[gtex_summary$multi_BP!=0], gtex_summary$multi_BP[gtex_summary$multi_BP!=0]))

# different distrubution of SNP locations (relative to BPs) in clinvar vs. Gtex
snp_in <- c(clinvar_summary$snp_in[clinvar_summary$multi_BP=="1"],
            clinvar_summary$snp_in[clinvar_summary$multi_BP=="2+"],
            gtex_summary$snp_in[gtex_summary$multi_BP=="1"],
            gtex_summary$snp_in[gtex_summary$multi_BP=="2+"])
chisq.test(table(snp_in, c(rep("clinvar - single", length(clinvar_summary$snp_in[clinvar_summary$multi_BP=="1"])),
                           rep("clinvar - multi", length(clinvar_summary$snp_in[clinvar_summary$multi_BP=="2+"])),
                rep("gtex - single", length(gtex_summary$snp_in[gtex_summary$multi_BP=="1"])),
                rep("gtex - multi", length(gtex_summary$snp_in[gtex_summary$multi_BP=="2+"])))))

prop.table(table(gtex_summary$snp_in[gtex_summary$multi_BP!=0]))
prop.table(table(clinvar_summary$snp_in[clinvar_summary$multi_BP!=0]))

# Figure S9. effect of rs2269219 of FECH
query_name <- clinvar_summary[grep("rs2269219", clinvar_summary$id),]$id

m <- which(!is.na(match(clinvar_predictions$id, query_name)))

prediction_subset <- as.data.frame(clinvar_predictions[m])

prediction_subset <- prediction_subset[prediction_subset$branchpoint_prob > 0.52,]

dummy_line <- prediction_subset[3,]
dummy_line$U2_binding_energy <- 0
dummy_line$status <- "REF"

prediction_subset <- rbind(prediction_subset, dummy_line)

FigureS9_right <- ggplot(prediction_subset, 
       aes(x=to_3prime_point*-1, y=U2_binding_energy, fill=status, group=status)) + 
  geom_hline(yintercept = max(prediction_subset$U2_binding_energy[prediction_subset$status == "REF"]), 
             col=nt_cols[2],linetype = "dashed", size=0.5) +
  geom_hline(yintercept = max(prediction_subset$U2_binding_energy[prediction_subset$status == "ALT"]), 
             col=nt_cols[4],linetype = "dashed", size=0.5) +
  geom_bar(stat="identity",position="dodge", width=1,col="black", size=0.25) + theme_figure + 
  theme(legend.position = c(0.2,0.8)) +
  scale_y_continuous(name="Branchpoint Strength\n(U2 Binding Energy)", limits=c(0,5)) +
  scale_x_continuous(name="Distance to 3' Exon (nt)", breaks=c(-26,-24,-22,-20), labels=c(26,24,22,20)) +
  scale_fill_manual(values=nt_cols[c(4,2)], name="Allele")
  

FigureS9 <- ggdraw() + 
  draw_plot(FigureS9_right, 0.6,0, 0.4,1)

pdf("Figures/FigureS9.pdf", width=5, height=2)
FigureS9
dev.off()
