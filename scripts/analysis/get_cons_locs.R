################################################################################
#         Annotate conservation of branchpoints using phypoP                   #
################################################################################

options(stringsAsFactors = F)
library(data.table)
library(parallel)
library(plyr)

###### make intervals for conservation ######
load("data/gencode_v26.RData")

G26_all <- gencode_v26
rm(gencode_v26)

chroms <- unique(G26_all$chromosome)

if(exists("cons_tab")){rm(cons_tab)}
for(c in seq_along(chroms)){
  #BP +/- 100nt
  m <- which(G26_all$chromosome == chroms[c])
  sites <- sort(unique(G26_all$end[m]))
  
  w <- which(sites[-length(sites)] != (sites[-1]-1))
  starts <- sites[c(1, (w+1))]
  ends <- sites[c(w, length(sites))]
  
  starts <- starts-101
  ends <- ends+100
 
  if(exists("cons_tab")){
    cons_tab <- rbind(cons_tab,data.frame(chrom=chroms[c], start=starts,end=ends))
  }else{
    cons_tab <- data.frame(chrom= chroms[c], start=starts,end=ends)
  }
}

write.table(cons_tab, file=paste0("data/conservation/G26_cons_tab.csv"), row.names=F,quote=F,col.names=F,sep=",")

#make table covering intron regions from phylop bigwig
#python python scripts/analysis/genome_wide_predictions/get_bw_entries.py -bw data/HG38/hg38.phyloP100way.bw -csv data/conservation/G24_cons_tab.csv -o data/conservation/phyloP_introns100.csv
cmd <- "python scripts/analysis/genome_wide_predictions/get_bw_entries.py -bw data/conservation/hg38.phyloP100way.bw -csv data/conservation/G26_cons_tab.csv -o data/conservation/phyloP_introns100.csv"
system(cmd)

cons_scores <- as.data.frame(fread("data/conservation/phyloP_introns100.csv"))
cons_tab <- as.data.frame(fread("data/conservation/G26_cons_tab.csv"))
colnames(cons_tab) <- c("chrom","start","end")

#convert from wide to long format
chrom <- rep(cons_tab$chrom, 255)
p <- 1:255
score <- unlist(lapply(p, function(x) cons_scores[,x]))
position <- unlist(lapply(p, function(x) cons_tab$start + (x-1)))
    
conservation_long <- data.frame(chrom, score, position)

rm <- which(is.na(conservation_long$score))
conservation_long <- conservation_long[-rm,]

chroms <- unique(conservation_long$chrom)
for(c in seq_along(chroms)){
  ind <- which(conservation_long$chrom == chroms[c])
  dups <- which(duplicated(conservation_long$position[ind]))
  
  if(length(ind[dups]) > 0){
    conservation_long <- conservation_long[-ind[dups]]
  }
  
  message(chroms[c])
}


write.table(conservation_long, "data/conservation/conservation_in_introns.txt", row.names = F,quote=F, sep="\t")
conservation <- conservation_long
rm(conservation_long)
colnames(conservation)[1] <- "chromosome"

###### annotate branchpoint predictions with conservation ######

conservation <- as.data.frame(fread("data/conservation/conservation_in_introns.txt"))

chroms <- unique(G26_all$chromosome)
prob_score_cutoff <- 0.52

#only for predicted bps to save time
G26_BP <- G26_all[G26_all$branchpoint_prob >= prob_score_cutoff | G26_all$in_testtrain == "HC",]
G26_BP$id <- with(G26_BP, paste0(exon_id, "_", to_3prime, "_",branchpoint_nt))


for(c in seq_along(chroms)){
  
  c1 <- which(G26_BP$chromosome==chroms[c])
  c2 <- which(conservation$chromosome==chroms[c])
  
  #G24_df is id,(end) position, strand, dist2
  #cons_df is conservation
  cons_df <- conservation[c2,]
  G26_df <- G26_BP[c1,c("id","end","strand","to_3prime")]
  
  strand <- G26_df[,3]
  to_exon <- as.numeric(G26_df[,4])
  
  cons_SS_p0 <- cons_df$score[match((as.numeric(G26_df[,2])+to_exon), cons_df$pos)]
  cons_SS_n1 <- cons_df$score[match((as.numeric(G26_df[,2])+to_exon-1), cons_df$pos)]
  
  cons_SS_p0[strand=="-"] <- (cons_df$score[match((as.numeric(G26_df[,2])-(to_exon)+1), cons_df$pos)])[strand=="-"]
  cons_SS_n1[strand=="-"] <- (cons_df$score[match((as.numeric(G26_df[,2])-(to_exon)+2), cons_df$pos)])[strand=="-"]
  
  #conservation in exon + 50
  locs <- matrix(nrow=length(to_exon), ncol=50)
  for(i in 1:50){
    locs[,i] <- G26_df[,2] + to_exon + i
  }
  locsn <- matrix(nrow=length(to_exon), ncol=50)
  for(i in 1:50){
    locsn[,i] <- (G26_df[,2] - to_exon + 1) - i
  }
  locs[strand=="-",] <- locsn[strand=="-",]
  m <- match(locs, cons_df$pos)
  cons_locs <- matrix(cons_df$score[m], ncol=50)
  
  cons_exon_mean <- rowMeans(cons_locs, na.rm = T)
  cons_exon_med <- apply(cons_locs,1,median, na.rm=T)
  #conservation in intron
  
  locs <- matrix(nrow=length(to_exon), ncol=50)
  for(i in 1:50){
    locs[,i] <- G26_df[,2] + (to_exon -1) - (i-1)
  }
  locsn <- matrix(nrow=length(to_exon), ncol=50)
  for(i in 1:50){
    locsn[,i] <- (G26_df[,2] - to_exon) + i
  }
  locs[strand=="-",] <- locsn[strand=="-",]

  m <- match(locs, cons_df$pos)
  cons_locs <- matrix(cons_df$score[m], ncol=50)
  cons_intron_mean <- rowMeans(cons_locs, na.rm = T)
  cons_intron_med <- apply(cons_locs,1,median, na.rm=T)
  
  #cons_pos0
  m <- match(G26_BP$end[c1], conservation$pos[c2])
  cons_pos0 <- conservation$score[c2][m]
  
  #split into pos/negative
  pos <- which(G26_BP$strand[c1] == "+")
  neg <- which(G26_BP$strand[c1] == "-")
  
  #for pos
  cons_pos1 <- conservation$score[c2][match((G26_BP$end[c1] + 1), conservation$pos[c2])]
  cons_pos2 <- conservation$score[c2][match((G26_BP$end[c1] + 2), conservation$pos[c2])]
  cons_pos3 <- conservation$score[c2][match((G26_BP$end[c1] + 3), conservation$pos[c2])]
  cons_pos4 <- conservation$score[c2][match((G26_BP$end[c1] + 4), conservation$pos[c2])]
  cons_pos5 <- conservation$score[c2][match((G26_BP$end[c1] + 5), conservation$pos[c2])]
  cons_neg1 <- conservation$score[c2][match((G26_BP$end[c1] - 1), conservation$pos[c2])]
  cons_neg2 <- conservation$score[c2][match((G26_BP$end[c1] - 2), conservation$pos[c2])]
  cons_neg3 <- conservation$score[c2][match((G26_BP$end[c1] - 3), conservation$pos[c2])]
  cons_neg4 <- conservation$score[c2][match((G26_BP$end[c1] - 4), conservation$pos[c2])]
  cons_neg5 <- conservation$score[c2][match((G26_BP$end[c1] - 5), conservation$pos[c2])]
  
  conservation_df <- data.frame(id = G26_BP$id[c1],
                                cons_neg5,cons_neg4,cons_neg3,
                                cons_neg2, cons_neg1,cons_pos0,
                                cons_pos1,cons_pos2,cons_pos3,
                                cons_pos4,cons_pos5,
                                cons_SS_p0,cons_SS_n1,
                                cons_exon_mean,cons_exon_med,
                                cons_intron_mean,cons_intron_med)
  
  conservation_df$cons_neg5[neg] <- cons_pos5[neg]
  conservation_df$cons_neg4[neg] <- cons_pos4[neg]
  conservation_df$cons_neg3[neg] <- cons_pos3[neg]
  conservation_df$cons_neg2[neg] <- cons_pos2[neg]
  conservation_df$cons_neg1[neg] <- cons_pos1[neg]
  conservation_df$cons_pos5[neg] <- cons_neg5[neg]
  conservation_df$cons_pos4[neg] <- cons_neg4[neg]
  conservation_df$cons_pos3[neg] <- cons_neg3[neg]
  conservation_df$cons_pos2[neg] <- cons_neg2[neg]
  conservation_df$cons_pos1[neg] <- cons_neg1[neg]
  
  conservation_df <- cbind(conservation_df, index = c1)
  write.table(conservation_df, 
              file = paste0("data/conservation/conservation_df_",chroms[c],".txt"), 
              row.names = FALSE, quote = FALSE, sep = "\t")
}