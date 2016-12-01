################################################################################
#           Genome wide prediction of splicing branchpoints - plots            #
################################################################################

# Produces Figures 2 - 4
options(stringsAsFactors = F)

library(ggplot2)
library(cowplot)
library(stringr)
library(pheatmap)
library(data.table)
source("scripts/analysis/genome_wide_predictions/quick_wilcox.R")

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

Figure1A=ggplot(Mercer_branchpoints, aes(x=dist.2*-1, fill=set)) + 
    geom_bar() + 
    theme_bw() + 
    scale_fill_manual(values=nt_cols[c(4,3)] , name=element_blank(), 
                      labels=c("High confidence\nMercer branchpoints",
                               "Low confidence\nMercer branchpoints")) + 
    scale_x_continuous(limits=c(-100,0), 
                       labels=c("100","75","50","25","0"), 
                       name="Distance to 3' exon") +
    geom_vline(xintercept = quantile(Mercer_branchpoints_hc$dist.2, c(0.05,.95))*-1, 
               color="grey",linetype = "dashed") +
    theme(text=element_text(size=6),legend.key.size=unit(0.2, "inches"),
          legend.position=c(0.25,0.8), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.border = element_blank(), 
          axis.line.x = element_line(), axis.line.y = element_line())


Figure1 = ggdraw() + 
    draw_plot(Figure1A, 0,0.7,0.45,0.3) + 
    draw_plot_label(c("A","B","C"), c(0,0.45,0), c(1,1,0.75), size=18)

pdf("Figures/Figure1_part1.pdf", useDingbats = F, height=4.5, width=6.69)
Figure1
dev.off()

rm(Mercer_branchpoints_hc,Mercer_branchpoints_hc_chrom,
   Mercer_branchpoints_lc,Mercer_branchpoints_lc_chrom)

load("data/Figure_files.Rdata")

###### Figure 2A&B - Branchpoint motifs ######

Figure2A=ggplot(fivemer_summary, aes(x=percent_BP, y=median_branchpoint_prob, size=num_BP, col=BP_nt)) + 
  geom_point() +scale_color_manual(values = nt_cols, name="BP nucleotide")+theme_bw()  +
  scale_x_continuous(name="Relative motif frequency\n(BPs/Negatives)") + scale_y_continuous(name="Median branchpointer\nprobability score") + 
  theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches")) + scale_size_continuous(name="Number of BPs\nin training dataset", range=c(0.5,3))

Figure2B=ggplot(fivemer_summary, aes(y=median_branchpoint_prob, x=mean_U2,size=num_BP, col=BP_nt)) + geom_point() +
  scale_color_manual(values = nt_cols, name="BP nucleotide")+theme_bw()+
  scale_y_continuous(name="Median branchpointer\nprobability score") + 
  scale_x_continuous(name="Branchpoint strength\n(Mean U2 binding energy)") +
  theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches"))+ scale_size_continuous(name="Number of BPs\nin training dataset", range=c(0.5,3))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend_Figure2=g_legend(Figure2A)

###### Figure 2C&D - Prediction of splicing branchpoints in GENCODE introns ######

#want number of introns annotated by mercer capture for the trest/train data
branchpoints_introns$bps <- branchpoints_introns$annotated_BPs
branchpoints_introns$bps[branchpoints_introns$status == 'test/train'] <- 
  branchpoints_introns$known_BPs[branchpoints_introns$status == 'test/train']
branchpoints_introns$bps[branchpoints_introns$bps >= 4] ="4+"

Figure2C=ggplot(branchpoints_introns, aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black") + 
  scale_fill_manual(values = BP_multi_colors, name="Branchpoint\ncount") + 
  theme_bw() + 
  theme(legend.position = c(0.8,0.7),text=element_text(size=10),legend.key.size=unit(0.2, "inches")) + 
  scale_y_continuous(name="introns") +
  scale_x_discrete(labels=c("predicted","known","unknown"), name="Branchpoint status\nx")

Figure2D=ggplot(branchpoints_introns[branchpoints_introns$status!="unknown",], aes(x=dex_mean_log, fill=status)) + 
  geom_density(alpha=0.5) +
  scale_fill_brewer(palette = "Paired", name="Branchpoint\nstatus") + 
  theme_bw() + 
  theme(legend.position = c(0.8,0.8),text=element_text(size=10),legend.key.size=unit(0.2, "inches")) + 
  scale_x_continuous(name="Exon expression\nlog10(exon normalised counts + 0.1)")

Figure2=
    ggdraw() + draw_plot(Figure2A + theme(legend.position="none"), 0,0.5,0.4,.5) + 
    draw_plot(Figure2B+ theme(legend.position="none"), 0.6,0.5,0.4,.5) + 
    draw_grob(legend_Figure2,0.4,0.5,0.2,.5)+
    draw_plot(Figure2C, 0,0,0.5,0.5) + 
    draw_plot(Figure2D, 0.5,0,0.5,0.5) + 
    draw_plot_label(c("A","B","C","D"), c(0,0.6,0,0.5), c(1,1,.5,.5), size=18)

pdf("Figures/Figure2.pdf", useDingbats = F, height=5.5, width=6.69)
Figure2
dev.off()

####### Figure 2 Supplement #######
#Figure2C by gene biotype
FigureS4=ggplot(branchpoints_introns, aes(fill=factor(bps), x=status)) + 
  geom_bar(color="black") + 
  facet_wrap(~gene_biotype_broad, scales = "free")+ 
  scale_fill_manual(values = BP_multi_colors, name="Branchpoint\ncount") +
  theme_bw() +
  theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches")) +
  scale_x_discrete(drop=FALSE,labels=c("predicted","known","unknown"))

pdf("Figures/FigureS4.pdf", 6.69,6,useDingbats = F)
FigureS4
dev.off()


###### Figure 3 - Features of introns with branchpoint multiplicity ######

quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "intron_size")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted" & branchpoints_introns$gene_biotype == "protein_coding",], "dex_mean_log")
quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "cons_intron_mean")

max_U2 <- aggregate(U2_binding_energy ~ id, G24_all, max)
branchpoints_introns$max_U2 <- max_U2$U2_binding_energy[match(branchpoints_introns$introns, max_U2$id)]

quick_wilcox(branchpoints_introns[branchpoints_introns$status=="predicted",], "max_U2")

ggplot(branchpoints_introns[branchpoints_introns$status=="predicted",], 
       aes(gene_biotype_broad,max_U2, fill=factor(gene_biotype_broad))) + 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) +
  scale_y_log10(name="Intron size (nt)")+ 
  scale_fill_manual(values = BP_multi_colors, name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none",text=element_text(size=10))+ scale_x_discrete(name="Annotated branchpoints") +
  facet_wrap(~predicted_BPs_factor)

aggregate(max_U2 ~ gene_biotype_broad, branchpoints_introns[branchpoints_introns$predicted_BPs==2,], median)

wilcox.test(branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs!=0 & 
                                          branchpoints_introns$gene_biotype_broad == "lncRNA"],
            branchpoints_introns$max_U2[branchpoints_introns$predicted_BPs!=0 & 
                                          branchpoints_introns$gene_biotype_broad == "protein_coding"])

Figure3A=ggplot(branchpoints_introns[branchpoints_introns$status=="predicted",], 
                       aes(predicted_BPs_factor,intron_size, fill=factor(predicted_BPs_factor))) + 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) +
  scale_y_log10(name="Intron size (nt)")+ 
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_bw() + 
  theme(legend.position = "none",text=element_text(size=10))+ scale_x_discrete(name="Annotated branchpoints")

Figure3B=ggplot(branchpoints_introns[branchpoints_introns$gene_biotype =="protein_coding" & 
                                       branchpoints_introns$status=="predicted",], 
                       aes(y=dex_mean_log, fill=factor(predicted_BPs_factor), x=predicted_BPs_factor)) + 
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_bw()+ 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) + 
  scale_x_discrete(name="Annotated branchpoints") + 
  scale_y_continuous(name="log10\n(exon normalised counts + 0.1)") + 
  theme(legend.position = "none",text=element_text(size=10)) 

Figure3C=ggplot(branchpoints_introns[which(branchpoints_introns$status=="predicted" & !is.na(branchpoints_introns$cons_intron_mean)),], 
  aes(y=cons_intron_mean, fill=factor(predicted_BPs_factor), x=predicted_BPs_factor)) + 
  geom_violin(scale="width", lwd=0.25) + 
  geom_boxplot(alpha=0, width=0.5, outlier.size = 0.25, lwd=0.25) +
  scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
  theme_bw() +
  scale_x_discrete(name="Annotated branchpoints") + 
  scale_y_continuous(name="mean intron conservation") +
  theme(legend.position = "none",text=element_text(size=10))

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

alt_5exons <- data.frame(branchpoint_multiplicity=rep(c("1","2+"),each=7), 
                         alternative_5exons=rep(rownames(tbl),2),
                         raw_count=c(tbl[,1],tbl[,2]),
                         total = rep(colSums(tbl),each=7))
alt_5exons$percent_count <- alt_5exons$raw_count/alt_5exons$total

Figure3D=ggplot(alt_5exons, aes(x=alternative_5exons, fill=branchpoint_multiplicity, y=percent_count)) + 
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
    theme_bw() +
    scale_x_discrete(name="number of 5' exons associated with branchpoint") + 
    scale_y_continuous(name="percent of introns") +
    theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches"))


#and SJ reads
branchpoints_introns$`alternative5'annotated_factor`=branchpoints_introns$`alternative5'annotated`
branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$`alternative5'annotated` > 4] ="5+"
tbl <- table(branchpoints_introns$`alternative5'annotated_factor`[branchpoints_introns$status == "predicted"], 
             branchpoints_introns$multiplicity[branchpoints_introns$status == "predicted"])
chisq.test(tbl)

wilcox.test((branchpoints_introns$var_maxentscan[branchpoints_introns$predicted_BPs > 1 & branchpoints_introns$gencode_alternative5 > 1]),
            (branchpoints_introns$var_maxentscan[branchpoints_introns$predicted_BPs == 1 & branchpoints_introns$gencode_alternative5 > 1]))

Figure3E=ggplot(branchpoints_introns[branchpoints_introns$gencode_alternative5 > 1 & branchpoints_introns$annotated_BPs_factor!="0",], 
       aes(x=annotated_BPs_factor, y=var_maxentscan, fill=annotated_BPs_factor)) + 
    geom_boxplot( outlier.size = 0.25, lwd=0.25) + 
    scale_y_log10(name="5' splice site strength variation") + 
    scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") + 
    theme_bw() +
    theme(text=element_text(size=10),legend.position="none") + 
    scale_x_discrete(name="Number of branchpoints in intron") 

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

head(snp_locs)
tbl <- table(snp_locs_all$snp_in[snp_locs_all$set=="GTEx" & snp_locs_all$multi_BP !="0" & !is.na(snp_locs_all$snp_in)], 
             snp_locs_all$multi_BP[snp_locs_all$set=="GTEx" & snp_locs_all$multi_BP !="0" & !is.na(snp_locs_all$snp_in)])

c <- chisq.test(tbl)
c$p.value


Figure3F=ggplot(snp_locs, aes(x=Variant_location , fill=Var2, y=percent_Freq)) + 
    geom_bar(stat="identity",position="dodge") +  
    geom_text(aes(x=Variant_location , y=percent_Freq + 0.05 * sign(percent_Freq), label=Freq), hjust=0.5, size=2) +
    theme_bw() +
    theme(text=element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position="none") + 
    scale_fill_manual(values = BP_multi_colors[c(2,5,4,3)],name="Annotated branchpoints",labels=c("1","2+")) + 
    scale_x_discrete(name="variant location", labels=rep(c("3'SS", "3'SS to BP", "BP","5'SS to BP"), each=2)) + 
    scale_y_continuous(name="relative common SNP frequency")

legend_Figure3=g_legend(Figure3D)

Figure3 <- ggdraw() + 
    draw_plot(Figure3A, 0.00,0.70, 0.33,0.3) + 
    draw_plot(Figure3B, 0.33,0.70, 0.33,0.3) + 
    draw_plot(Figure3C, 0.66,0.70, 0.33,0.3) + 
    draw_plot(Figure3D + theme(legend.position="none"), 0.00,0.00, 0.66,0.3)+ 
    draw_plot(Figure3E, 0.66,0.00, 0.33,0.3) + 
    draw_plot(Figure3F, 0.33,0.30, 0.66,0.4) +
    draw_plot(legend_Figure3, 0.00,0.35, 0.33,0.4) +
    draw_plot_label(c("A","B","C", "D","E","F"), c(0,0.33,0.66,0,0,0.66), c(1,1,1,0.7,0.3,0.3), size=18)
Figure3
pdf("Figures/Figure3.pdf", height=7,width=6.69, useDingbats = FALSE)
Figure3
dev.off()

####### Figure 3 Supplement #######

#Intron size trend present in Mercer annotation
#non-annotated introns are longer

Mercer_branchpoints$intron_size <- Mercer_branchpoints$dist.1 + Mercer_branchpoints$dist.2
#all
Mercer_introns <- as.data.frame(table(Mercer_branchpoints$exon_id.2[Mercer_branchpoints$set=="HC"]))
Mercer_introns$size <- Mercer_branchpoints$intron_size[match(Mercer_introns$Var1, Mercer_branchpoints$exon_id.2)]
Mercer_introns$Freq[Mercer_introns$Freq >= 4] <- "4+"

FigureS5=ggplot(Mercer_introns, aes(x=Freq, y=size, group=Freq, fill=Freq)) + 
    geom_violin(scale="width") + 
    geom_boxplot(alpha=0.1, width=0.5) +
    scale_y_log10(name="Intron size (nt)") + 
    scale_fill_manual(values = BP_multi_colors[-1], name="Annotated branchpoints") + 
    theme_bw() + 
    theme(legend.position = "none") + 
    scale_x_discrete(name="Annotated branchpoints") 

wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 1], 
            Mercer_introns$size[Mercer_introns$Freq == 2])

wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 2], 
            Mercer_introns$size[Mercer_introns$Freq == 3])

wilcox.test(Mercer_introns$size[Mercer_introns$Freq == 3], 
            Mercer_introns$size[Mercer_introns$Freq == "4+"])


pdf("Figures/FigureS5.pdf", 3.35,3,useDingbats = F)
FigureS5
dev.off()


gene_types=as.data.frame(table(branchpoints_introns$gene_biotype))
keep_gene_types=gene_types$Var1[gene_types$Freq > 1000]
FigureS6=ggplot(branchpoints_introns[(branchpoints_introns$gene_biotype %in% keep_gene_types & branchpoints_introns$status=="predicted"),], 
                aes(y=intron_size, fill=factor(annotated_BPs_factor), x=annotated_BPs_factor)) + 
    geom_violin(scale="width") + 
    geom_boxplot(alpha=0.1, width=0.5) +
    scale_y_log10(name="Intron size (nt)") + 
    scale_fill_manual(values = BP_multi_colors[c(2,5)], name="Annotated branchpoints") +
    theme_bw() + 
    theme(legend.position = "none") + 
    scale_x_discrete(name="Annotated branchpoints") + 
    facet_wrap(~gene_biotype)

pdf("Figures/FigureS6.pdf", 6.69,6,useDingbats = F)
FigureS6
dev.off()

for(type in keep_gene_types){
  message(type)
  print(wilcox.test(branchpoints_introns$intron_size[branchpoints_introns$gene_biotype == type & branchpoints_introns$annotated_BPs_factor ==1],
                    branchpoints_introns$intron_size[branchpoints_introns$gene_biotype == type & branchpoints_introns$annotated_BPs_factor =="2+"]))
}


####### Supplementary Figure 7 -- exon skipping #######
load("data/exon_skipping.Rdata")

FigureS7 <- ggplot(splicing_strength, aes(x=variable, lower=lower, upper=upper, 
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

pdf("Figures/FigureS7.pdf",useDingbats = F, height=3.5, width=6.69)
FigureS7
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


###### Figure 4 - Disease causing variants affecting BPs #######
#locations of clinVar intronic variants and GTEX variants

load("data/diseaseVariants.Rdata")

clinVar_processed$snp_in <- NA
clinVar_processed$snp_in[clinVar_processed$dist_to_BP <= -1] <- "3'SS to BP"
clinVar_processed$snp_in[clinVar_processed$dist_to_BP >= 3] <- "BP to 5'SS"
clinVar_processed$snp_in[clinVar_processed$dist_to_BP == 0 | clinVar_processed$dist_to_BP == 0 | clinVar_processed$dist_to_BP == 0] <- "BP"
clinVar_processed$snp_in[clinVar_processed$dist_to_exon < 3] <- "3'SS"

clinVar_processed$set <- "ClinVar"
processed_GTEX$set <- "GTEx"

snp_locs_all <- rbind(processed_GTEX, clinVar_processed[,match(colnames(processed_GTEX), colnames(clinVar_processed))])

FigureS8=ggplot(snp_locs_all[snp_locs_all$multi_BP!="0" & !is.na(snp_locs_all$snp_in),], 
       aes(x=snp_in, fill=multi_BP)) + 
    geom_bar() +  
    theme_bw() +
    scale_fill_manual(name="annotated\nbranchpoints",values = BP_multi_colors[c(2,5)]) +
    theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches")) +
    facet_wrap(~set, scales="free") +scale_x_discrete(name="variant location", labels=c("3'SS", "3'SS to BP", "BP","5'SS to BP"))

pdf("Figures/FigureS8.pdf", width=6.69, height=3)
FigureS8
dev.off()

chisq.test(table(clinVar_processed$snp_in[clinVar_processed$multi_BP!=0], clinVar_processed$multi_BP[clinVar_processed$multi_BP!=0]))

chisq.test(table(processed_GTEX$snp_in[processed_GTEX$multi_BP!=0], processed_GTEX$multi_BP[processed_GTEX$multi_BP!=0]))

source("scripts/analysis/variation/plotBranchpointWindow_Figure.R")

#FECH branchpoint mutation
pdf("Figures/Figure4.pdf", width=6.69, height=6)
plotBranchpointWindow_Figure(clinVars_filtered$id[clinVars_filtered$GeneSymbol=="FECH"], clinvar_predictions,clinvar_attributes,exons)
dev.off()