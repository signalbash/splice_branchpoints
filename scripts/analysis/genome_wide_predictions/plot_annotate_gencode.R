################################################################################
#           Genome wide prediction of splicing branchpoints - plots            #
################################################################################

# Produces Figures 4 - 6
options(stringsAsFactors = F)

library(ggplot2)
library(cowplot)
library(stringr)

number_introns_pal=rev(c('#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'))
nt_cols=c("#359646","#4D7ABE","#FAA859","#CB3634")

load("data/Figure_files.Rdata")

###### Figure 4 - Branchpoint detection model scores for fivemer motifs ######

Figure4A=ggplot(fivemer_summary, aes(x=percent_BP, y=median_branchpoint_prob, size=num_BP, col=BP_nt)) + 
  geom_point() +scale_color_manual(values = nt_cols, name="BP nucleotide")+theme_bw()  +
  scale_x_continuous(name="Relative motif frequency (BPs/Negatives)") + scale_y_continuous(name="Median branchpointer probability score") + 
  theme(text=element_text(size=8),legend.key.size=unit(0.2, "inches")) + scale_size_continuous(name="Number of BPs\nin training dataset", range=c(0.5,3))

Figure4B=ggplot(fivemer_summary, aes(y=median_branchpoint_prob, x=mean_U2,size=num_BP, col=BP_nt)) + geom_point() +
  scale_color_manual(values = nt_cols, name="BP nucleotide")+theme_bw()+
  scale_y_continuous(name="Median branchpointer probability score") + 
  scale_x_continuous(name="Mean U2 binding energy") +
  theme(text=element_text(size=8),legend.key.size=unit(0.2, "inches"))+ scale_size_continuous(name="Number of BPs\nin training dataset", range=c(0.5,3))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend_Figure4=g_legend(Figure4A)

Figure4=ggdraw() + draw_plot(Figure4A + theme(legend.position="none"), 0,0,0.4,1) + 
  draw_plot(Figure4B+ theme(legend.position="none"), 0.4,0,0.4,1) + 
  draw_grob(legend_Figure4,0.8,0,0.2,1)+
  draw_plot_label(c("A","B"), c(0,0.4), c(1,1), size=10)

Figure4
pdf("Figures/Figure4.pdf", useDingbats = F, height=3, width=6.69)
Figure4
dev.off()

###### Figure 5 - Prediction of splicing branchpoints in GENCODE introns ######

#want number of introns annotated by mercer capture for the trest/train data
branchpoints_introns$bps <- branchpoints_introns$annotated_BPs
branchpoints_introns$bps[branchpoints_introns$status == 'test/train'] <- 
  branchpoints_introns$known_BPs[branchpoints_introns$status == 'test/train']
branchpoints_introns$bps[branchpoints_introns$bps >= 4] ="4+"

Figure5A=ggplot(branchpoints_introns, aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black") + 
  scale_fill_manual(values = number_introns_pal, name="Annotated\nbranchpoints") + 
  theme_bw() + 
  theme(legend.position = c(0.8,0.7),text=element_text(size=8),legend.key.size=unit(0.2, "inches")) + 
  scale_y_continuous(name="introns")

Figure5B=ggplot(branchpoints_introns, aes(x=dex_mean_log, fill=status)) + 
  geom_density(alpha=0.5) +
  scale_fill_brewer(palette = "Paired") + 
  theme_bw() + 
  theme(legend.position = c(0.8,0.8),text=element_text(size=8),legend.key.size=unit(0.2, "inches")) + 
  scale_x_continuous(name="log10(exon normalised counts + 0.1)")

Figure5=ggdraw() + draw_plot(Figure5A, 0,0,0.5,1) + 
  draw_plot(Figure5B, 0.5,0,0.5,1) + 
  draw_plot_label(c("A","B"), c(0,0.5), c(1,1), size=10)

pdf("Figures/Figure5_cowplot.pdf", useDingbats = F, height=3, width=6.69)
Figure5
dev.off()

#Supplementary
#Figure5A by gene biotype
FigureS3=ggplot(branchpoints_introns, aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black") + 
  facet_wrap(~gene_biotype_broad, scales = "free")+ 
  scale_fill_manual(values = number_introns_pal, name="Annotated\nbranchpoints") + 
  theme_bw() +
  theme(text=element_text(size=8),legend.key.size=unit(0.2, "inches")) +
  scale_x_discrete(drop=FALSE)

pdf("Figures/FigureS3.pdf", 6.69,6,useDingbats = F)
FigureS3
dev.off()


###### Figure 6 - Features of introns with branchpoint multiplicity ######

Figure6A_top=ggplot(branchpoints_introns[branchpoints_introns$status=="predicted",], 
                    aes(fill=factor(predicted_BPs_factor), x=predicted_BPs_factor)) + 
  geom_bar(color="black", lwd=0.25) + 
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(),text=element_text(size=8)) + 
  scale_x_discrete(name="Annotated branchpoints")+ scale_y_continuous(name="introns")

Figure6A_bottom=ggplot(branchpoints_introns[branchpoints_introns$status=="predicted",], 
                       aes(predicted_BPs_factor,intron_size, fill=factor(predicted_BPs_factor))) + 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) +
  scale_y_log10(name="Intron size (nt)")+ 
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none",text=element_text(size=8))+ scale_x_discrete(name="Annotated branchpoints")

Figure6B_top=ggplot(branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & branchpoints_introns$status=="predicted",], 
                    aes(fill=factor(predicted_BPs_factor), x=predicted_BPs_factor)) + 
  geom_bar(color="black", lwd=0.25) + 
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(),text=element_text(size=8)) + 
  scale_x_discrete(name="Annotated branchpoints") + 
  scale_y_continuous(name="introns")

Figure6B_bottom=ggplot(branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & branchpoints_introns$status=="predicted",], 
                       aes(y=dex_mean_log, fill=factor(predicted_BPs_factor), x=predicted_BPs_factor)) + 
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw()+ 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) + 
  scale_x_discrete(name="Annotated branchpoints") + 
  scale_y_continuous(name="log10(exon normalised counts + 0.1)") + 
  theme(legend.position = "none",text=element_text(size=8)) 

Figure6C_top=ggplot(branchpoints_introns[which(branchpoints_introns$status=="predicted" & !is.na(branchpoints_introns$cons_intron_mean)),], 
                    aes(fill=factor(predicted_BPs_factor), x=predicted_BPs_factor)) + 
  geom_bar(color="black", lwd=0.25) + 
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(),text=element_text(size=8)) +
  scale_x_discrete(name="Annotated branchpoints") + 
  scale_y_continuous(name="introns")

Figure6C_bottom=ggplot(branchpoints_introns[which(branchpoints_introns$status=="predicted" & !is.na(branchpoints_introns$cons_intron_mean)),], 
  aes(y=cons_intron_mean, fill=factor(predicted_BPs_factor), x=predicted_BPs_factor)) + 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) +
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw() +
  scale_x_discrete(name="Annotated branchpoints") + 
  scale_y_continuous(name="mean intron conservation") +
  theme(legend.position = "none",text=element_text(size=8))

Figure6A=plot_grid(Figure6A_top,Figure6A_bottom, ncol=1, rel_heights = c(1,3), align="v")
Figure6B=plot_grid(Figure6B_top,Figure6B_bottom, ncol=1, rel_heights = c(1,3), align="v")
Figure6C=plot_grid(Figure6C_top,Figure6C_bottom, ncol=1, rel_heights = c(1,3), align="v")

Figure6=ggdraw() + draw_plot(Figure6A, 0,0,0.33,1) + 
  draw_plot(Figure6B, 0.33,0,0.33,1) + 
  draw_plot(Figure6C, 0.66,0,0.33,1) + 
  draw_plot_label(c("A","B","C"), c(0,0.33,0.66), c(1,1,1), size=10)

pdf("Figures/Figure6.pdf", height=3,width=6.69, useDingbats = FALSE)
Figure6
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

tbl_percent <- as.data.frame(cbind((tbl[,1]/rowSums(tbl)),tbl[,2]/rowSums(tbl)))
colnames(tbl_percent) <- c("single_BP", "multiple_BP")
tbl_percent <- cbind(tbl_percent, number_alt_exons=rownames(tbl_percent))

pheatmap(tbl_percent[,-3], cluster_rows = F, cluster_cols = F)

#and SJ reads
branchpoints_introns$`alternative5'annotated_factor`=branchpoints_introns$`alternative5'annotated`
branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$`alternative5'annotated` > 4] ="5+"
tbl <- table(branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$status == "predicted"], 
             branchpoints_introns$multiplicity[branchpoints_introns$status == "predicted"])
chisq.test(tbl)
tbl[,1]/rowSums(tbl)


branchpoints_introns$max_dist_diff_factor <- branchpoints_introns$max_dist_diff
branchpoints_introns$max_dist_diff_factor[which(branchpoints_introns$max_dist_diff < 10 & !is.na(branchpoints_introns$max_dist_diff))] <- "close"
branchpoints_introns$max_dist_diff_factor[which(branchpoints_introns$max_dist_diff >= 10 & !is.na(branchpoints_introns$max_dist_diff))] <- "far"

tbl <- table(as.character(branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$predicted_BPs==3]), 
             branchpoints_introns$max_dist_diff_factor[branchpoints_introns$predicted_BPs==3])
chisq.test(tbl)

tbl <- table(as.character(branchpoints_introns$gencode_alternative5_factor[branchpoints_introns$predicted_BPs==3]), 
             branchpoints_introns$max_dist_diff_factor[branchpoints_introns$predicted_BPs==3])
chisq.test(tbl)
tbl[,1]/rowSums(tbl)


ggplot(branchpoints_introns[which(!is.na(branchpoints_introns$max_dist_diff)),], 
       aes(x=factor(`alternative5'annotated`), y=min_dist_diff)) + 
    geom_violin() +
    facet_wrap(~predicted_BPs)

#Supplementary
#Figure 6A SUPP
gene_types=as.data.frame(table(branchpoints_introns$gene_biotype))
keep_gene_types=gene_types$Var1[gene_types$Freq > 1000]
FigureS4=ggplot(branchpoints_introns[(branchpoints_introns$gene_biotype %in% keep_gene_types & branchpoints_introns$status=="predicted"),], 
       aes(y=intron_size, fill=factor(annotated_BPs_factor), x=annotated_BPs_factor)) + 
  geom_violin(scale="width") + 
  geom_boxplot(alpha=0.1, width=0.5) +
  scale_y_log10(name="Intron size (nt)") + 
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_x_discrete(name="Annotated branchpoints") + 
  facet_wrap(~gene_biotype)

#Intron size trend present in Mercer annotation
#non-annotated introns are longer
FigureS5=ggplot(branchpoints_introns[(branchpoints_introns$status!="predicted" & branchpoints_introns$bps !=0),], 
       aes(y=intron_size, fill=factor(bps), x=bps)) + 
  geom_violin(scale="width") + 
  geom_boxplot(alpha=0.1, width=0.5) +
  scale_y_log10(name="Intron size (nt)") + 
  scale_fill_manual(values = number_introns_pal[-1], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_x_discrete(name="Annotated branchpoints") 

pdf("Figures/FigureS4.pdf", 6.69,6,useDingbats = F)
FigureS4
dev.off()
pdf("Figures/FigureS5.pdf", 3.35,3,useDingbats = F)
FigureS5
dev.off()

#variable_class == bps
quick_wilcox <- function(data_frame, variable_test,variable_class="bps"){
  class1 <- 1
  class2 <- 2
  class3 <- 3
  class4 <- "4+"
  
  #1 v 2
  c1 <- which(data_frame[,variable_class] %in% class1)
  c2 <- which(data_frame[,variable_class] %in% class2)
  
  w1 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  #1 v all
  c1 <- which(data_frame[,variable_class] %in% class1)
  c2 <- which(data_frame[,variable_class] %in% c(class2,class3,class4))
  
  w2 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  
  #2 v 3
  c1 <- which(data_frame[,variable_class] %in% class2)
  c2 <- which(data_frame[,variable_class] %in% class3)
  
  w3 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  #3 v 4
  c1 <- which(data_frame[,variable_class] %in% class3)
  c2 <- which(data_frame[,variable_class] %in% class4)
  
  w4 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  df <- data.frame(variable=variable_test,class1=c(1,1,2,3), class2=c(2,"2,3,4+", 3,"4+"), pval=c(w1$p.value, w2$p.value,w3$p.value,w4$p.value))
  
  return(df)
}

quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "intron_size")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted" & branchpoints_introns$gene_biotype == "protein_coding",], "dex_mean_log")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "cons_intron_mean")


