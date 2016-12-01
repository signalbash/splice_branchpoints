args <- commandArgs(trailingOnly = TRUE)

k1=as.character(args[1])
k2=as.character(args[2])


#set directory
setwd(k2)

options(stringsAsFactors = F)
chrom=k1

#branchpoint_df_lc=read.csv(file=paste0("data/outputs/branchpoint_df_gene_lc_", chrom,".csv"))
#branchpoint_df_lc$X=NULL
branchpoint_df_hc=read.csv(file=paste0("data/outputs/branchpoint_df_gene_hc_", chrom,".csv"))
branchpoint_df_hc$X=NULL
branchpoint_df_N=read.csv(file=paste0("data/outputs/branchpoint_df_gene_N", chrom,".csv"))
branchpoint_df_N$X=NULL
branchpoint_df_N$X.1=NULL

branchpoint_df_hc=cbind(branchpoint_df_hc, set="HC")
branchpoint_df_hc=branchpoint_df_hc[(branchpoint_df_hc$dist.2 >=18 & branchpoint_df_hc$dist.2 <=44),]

#branchpoint_df_lc=cbind(branchpoint_df_lc, set="LC")
branchpoint_df=rbind(branchpoint_df_hc,branchpoint_df_N)
branchpoint_df=branchpoint_df[(branchpoint_df$dist.2 >=18 & branchpoint_df$dist.2 <=44),]

#make new id as chrZ_start_strand_stop_SET
new_ID=(paste0(branchpoint_df$Chromosome, "_", branchpoint_df$Start, branchpoint_df$Strand, branchpoint_df$Stop,"_",branchpoint_df$set))
branchpoint_df=cbind(new_ID, branchpoint_df)
#gets rid of any duplicated branchpoints (in case of multiple HC near same exon)
branchpoint_df=branchpoint_df[!(duplicated(new_ID)),]

write.csv(branchpoint_df, file=paste0("data/outputs/branchpoint_df_", chrom,".csv"))

branchpoint_df_x=read.csv(file=paste0("data/outputs/branchpoint_df_", chrom,".csv"))

#make a bed file to get fasta seqs
branchpoint_df_bed=branchpoint_df_x[,c(3,4,5,2,7,8)]
branchpoint_df_bed$Chromosome=gsub("chr","", branchpoint_df_bed$Chromosome)
#within 501nt window
branchpoint_df_bed$Start=branchpoint_df_bed$Start-250
branchpoint_df_bed$Stop=branchpoint_df_bed$Stop+250

write.table(branchpoint_df_bed, sep="\t", file=paste0("data/outputs/branchpoint_df_501_", chrom,".bed"), row.names = F,col.names = F,quote = F)

