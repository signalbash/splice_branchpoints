args <- commandArgs(trailingOnly = TRUE)

k1=as.character(args[1])
k2=as.character(args[2])

message(paste0("Finding nearest exons for  branchpoints on ",k1))

#set directory
setwd(k2)
options(stringsAsFactors = F)

chrom=k1
for(c in 1:length(chrom)){
#mismatch at bp(59359)
branchpoints_S1b=read.csv("data/inputs/branchpoints_S1b.csv", header = T)
#lower conf
branchpoints_S1e=read.csv("data/inputs/branchpoints_S1e.csv", header = T)

#load gencode annotations
gencode.v12.condensed <- read.table("data/inputs/altered_inputs/gencode.v12.condensed.txt", quote="\"", comment.char="")
colnames(gencode.v12.condensed) <- c("chromosome", "type","start","end",
                                     "strand","gene_id","gene_type",
                                     "transcript_id", "transcript_type")


#match to chromosome
gencode.v12.condensed=gencode.v12.condensed[gencode.v12.condensed$chromosome==chrom[c],]

#gencode genes/exons/transcripts
gencode.v12.condensed.genes=gencode.v12.condensed[gencode.v12.condensed$type=="gene",]
gencode.v12.condensed.exons=gencode.v12.condensed[gencode.v12.condensed$type=="exon",]
gencode.v12.condensed.trans=gencode.v12.condensed[gencode.v12.condensed$type=="transcript",]
unique_transcripts=unique(gencode.v12.condensed.exons$transcript_id)
gencode.v12.condensed.exons=cbind(gencode.v12.condensed.exons, exon_id="")

#make exon names
for(i in 1:length(unique_transcripts)){
  x=match(gencode.v12.condensed.exons$transcript_id, unique_transcripts[i])
  y=which(!(is.na(x)))
  gencode.v12.condensed.exons$exon_id[y] <- paste0(unique_transcripts[i],".",1:length(y))
  if(i %% 1000==0){
    message(paste(i, "of", length(unique_transcripts), " exons renamed."))
  }
}

rm=vector()
for(i in 1:length(gencode.v12.condensed.exons$chromosome)){
  e_start=gencode.v12.condensed.exons$start[i]
  e_end=gencode.v12.condensed.exons$end[i]
  e_strand=gencode.v12.condensed.exons$strand[i]
  e_chr=gencode.v12.condensed.exons$chromosome[i]
  
  x=which(gencode.v12.condensed.exons$start==e_start &
            gencode.v12.condensed.exons$end==e_end &
            gencode.v12.condensed.exons$strand==e_strand &
            gencode.v12.condensed.exons$chromosome==e_chr)
  if(length(x) > 1){
    rm=append(rm, x[-1])
    #used_name=gencode.v12.condensed.exons$exon_id[i]
    #replaced_names=paste(gencode.v12.condensed.exons$exon_id[x([-1])], collapse = ";")
    #line=data.frame(used_name=used_name, replaced_names=replaced_names)
    #duplicated_exons=rbind(duplicated_exons,line)
  }
  if(i %% 1000==0){
    message(paste(i, "of", length(gencode.v12.condensed.exons$chromosome), " exons processed for uniqueness."))
  }
}
rm=unique(rm)

#remove all duplcated exons (start&end&strand&chr are same)
gencode.v12.condensed.exons.u=gencode.v12.condensed.exons[-rm,]
write.csv(gencode.v12.condensed.exons.u, file=paste0("data/outputs/gencode.v12.condensed.exons_", chrom[c],".csv"))


#find alternative names for removed duplicates
duplicated_exons=data.frame(used_name=NA, replaced_names=NA)
for(i in 1:length(gencode.v12.condensed.exons.u$chromosome)){
  e_start=gencode.v12.condensed.exons.u$start[i]
  e_end=gencode.v12.condensed.exons.u$end[i]
  e_strand=gencode.v12.condensed.exons.u$strand[i]
  e_chr=gencode.v12.condensed.exons.u$chromosome[i]
  
  x=which(gencode.v12.condensed.exons$start==e_start &
            gencode.v12.condensed.exons$end==e_end &
            gencode.v12.condensed.exons$strand==e_strand &
            gencode.v12.condensed.exons$chromosome==e_chr)
  if(length(x) > 1){
    used_name=gencode.v12.condensed.exons.u$exon_id[i]
    replaced_names=paste(gencode.v12.condensed.exons$exon_id[(x[-1])], collapse = ";")
    line=data.frame(used_name=used_name, replaced_names=replaced_names)
    duplicated_exons=rbind(duplicated_exons,line)
  }
  if(i %% 1000==0){
    message(paste(i, "of", length(gencode.v12.condensed.exons.u$chromosome), "unique exon names processed."))
  }
}
duplicated_exons=duplicated_exons[-1,]
write.csv(duplicated_exons, file=paste0("data/outputs/duplicated_exons_", chrom[c],".csv"))


#find nearest exons for each branchpoint
branchpoint_df=branchpoints_S1b[branchpoints_S1b$Chromosome==chrom[c],]

for(i in 1:length(branchpoint_df$Chromosome)){
  branch_start=branchpoint_df$Start[i]
  branch_end=branchpoint_df$Stop[i]
  branch_chr=branchpoint_df$Chromosome[i]
  branch_strand=branchpoint_df$Strand[i]
  
  gencode.v12.condensed.exons.subset=gencode.v12.condensed.exons.u[which(gencode.v12.condensed.exons.u$chromosome==branch_chr & 
                                                                         gencode.v12.condensed.exons.u$strand==branch_strand ),]
  
  to_start=(gencode.v12.condensed.exons.subset$start - branch_end)
  to_start=min(to_start[to_start > 0])
  
  from_end=(branch_end - gencode.v12.condensed.exons.subset$end)
  from_end=min(from_end[from_end > 0])
  
  right_exon=which(gencode.v12.condensed.exons.subset$start == branch_end + to_start)
  left_exon=which(gencode.v12.condensed.exons.subset$end == branch_end - from_end)
  
  right_line=(cbind(gencode.v12.condensed.exons.subset[right_exon,], dist=to_start))[1,]
  left_line=(cbind(gencode.v12.condensed.exons.subset[left_exon,], dist=from_end))[1,]
  
  

  #reverse start/stop bed fields  if strand is negative
  if(branch_strand== "-"){
    whole_line=cbind(branchpoint_df[i,], right_line, left_line)
  }else{
    whole_line=cbind(branchpoint_df[i,], left_line, right_line)
  }
  if(exists("branchpoint_df2")){
    branchpoint_df2=rbind(branchpoint_df2, whole_line)
  }else{
    branchpoint_df2=whole_line
  }
}

#.1 is towards 5'SS
#.2 is towards 3'SS
colnames(branchpoint_df2)[7:17]=paste0(colnames(branchpoint_df2)[7:17],".1")
colnames(branchpoint_df2)[18:28]=paste0(colnames(branchpoint_df2)[18:28],".2")

write.csv(branchpoint_df2, file=paste0("data/outputs/branchpoint_df_all_hc_", chrom,".csv"))

#remove those where nearest 3' gene is not the same as the nearest 5' gene (i.e. not in an intron)
keep=which(branchpoint_df2$gene_id.1==branchpoint_df2$gene_id.2)
branchpoint_df3=branchpoint_df2[keep,]
write.csv(branchpoint_df3, file=paste0("data/outputs/branchpoint_df_gene_hc_", chrom,".csv"))

#For low confidence branchpoints - redo the same process
#find nearest exons for each branchpoint

branchpoint_df=branchpoints_S1e[branchpoints_S1e$Chromosome==chrom,]

rm(branchpoint_df2)
for(i in 1:length(branchpoint_df$Chromosome)){
  branch_start=branchpoint_df$Start[i]
  branch_end=branchpoint_df$Stop[i]
  branch_chr=branchpoint_df$Chromosome[i]
  branch_strand=branchpoint_df$Strand[i]
  
  gencode.v12.condensed.exons.subset=gencode.v12.condensed.exons.u[which(gencode.v12.condensed.exons.u$chromosome==branch_chr & 
                                                                           gencode.v12.condensed.exons.u$strand==branch_strand ),]
  
  to_start=(gencode.v12.condensed.exons.subset$start - branch_end)
  to_start=min(to_start[to_start > 0])
  
  from_end=(branch_end - gencode.v12.condensed.exons.subset$end)
  from_end=min(from_end[from_end > 0])
  
  right_exon=which(gencode.v12.condensed.exons.subset$start == branch_end + to_start)
  left_exon=which(gencode.v12.condensed.exons.subset$end == branch_end - from_end)
  
  right_line=(cbind(gencode.v12.condensed.exons.subset[right_exon,], dist=to_start))[1,]
  left_line=(cbind(gencode.v12.condensed.exons.subset[left_exon,], dist=from_end))[1,]
  
  
  #reverse start/stop bed fields  if strand is negative

  if(branch_strand== "-"){
    whole_line=cbind(branchpoint_df[i,], right_line, left_line)
  }else{
    whole_line=cbind(branchpoint_df[i,], left_line, right_line)
  }
  if(exists("branchpoint_df2")){
    branchpoint_df2=rbind(branchpoint_df2, whole_line)
  }else{
    branchpoint_df2=whole_line
  }
}

colnames(branchpoint_df2)[7:17]=paste0(colnames(branchpoint_df2)[7:17],".1")
colnames(branchpoint_df2)[18:28]=paste0(colnames(branchpoint_df2)[18:28],".2")
write.csv(branchpoint_df2, file=paste0("data/outputs/branchpoint_df_all_lc_", chrom,".csv"))

keep=which(branchpoint_df2$gene_id.1==branchpoint_df2$gene_id.2)
branchpoint_df3=branchpoint_df2[keep,]
write.csv(branchpoint_df3, file=paste0("data/outputs/branchpoint_df_gene_lc_", chrom,".csv"))
}

