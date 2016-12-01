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

load("data/Figure_files.Rdata")

####### Conservation ######
#see: get_cons_locs.R

chroms <- unique(G24_all$chromosome)
for(c in seq_along(chroms)){
  conservation_df_chrom <- as.data.frame(fread(paste0("data/conservation/conservation_df_",chroms[c],".txt")))
  
  if(exists("conservation_df")){
    conservation_df <- rbind(conservation_df,conservation_df_chrom)
  }else{
    conservation_df<- conservation_df_chrom
  }
}

m <- match(G24_all$id, conservation_df$id)
G24_all <- cbind(G24_all, conservation_df[m,-1])

m <- match(branchpoints_introns$introns, G24_all$id[which(G24_all$branchpoint_prob > 0.5)])
branchpoints_introns <-
  cbind(branchpoints_introns, G24_all[which(G24_all$branchpoint_prob > 0.5)[m],c("cons_SS_p0","cons_SS_n1",
                                          "cons_exon_mean","cons_exon_med",
                                          "cons_intron_mean","cons_intron_med")])

rm(conservation_df, conservation_df_chrom)

###### multi-branchpoint distances ######

w <- which(branchpoints_introns$predicted_BPs > 1)

m <- match(G24_all$id[G24_all$branchpoint_prob > 0.5], branchpoints_introns$introns[w])
mt <- as.data.frame(table(m))
mt$name <- branchpoints_introns$introns[w[mt$m]]

small_G24 <- G24_all[G24_all$branchpoint_prob > 0.5,c("id","to_3prime")]

for(i in 1:max(mt$Freq)){
    n <- match(mt$name, small_G24$id)
    mt <- cbind(mt, small_G24$to_3prime[n])
    small_G24 <- small_G24[-n[which(!is.na(n))],]
}

todim <- dim(mt)[2]

mt$min_dist_diff <- apply(mt[,c(4:todim)],1,function(x){min(diff(sort(as.numeric(x))))})
mt$max_dist_diff <- apply(mt[,c(4:todim)],1,function(x){max(diff(sort(as.numeric(x))))})

mt$max_dist <- apply(mt[,c(4:todim)],1,max, na.rm=T)
mt$min_dist <- apply(mt[,c(4:todim)],1,min, na.rm=T)
# mt$max_dist_diff <- mt$max_dist - mt$min_dist

m <- match(branchpoints_introns$introns, mt$name) 
branchpoints_introns <- cbind(branchpoints_introns, mt[m, c("max_dist","min_dist","max_dist_diff","min_dist_diff")])    
    
save(G24_all, branchpoints_introns, fivemer_summary,
     file = "data/Figure_files.Rdata")

