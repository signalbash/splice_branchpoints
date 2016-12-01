args <- commandArgs(trailingOnly = TRUE)

k1=as.character(args[1])
k2=as.character(args[2])
setwd(k2)

#contains funcitons for ppt calculations
source("scripts/Functions.R")
options(stringsAsFactors = F)
chrom=k1
branchpoint_df_x=read.csv(file=paste0("data/outputs/branchpoint_df_", chrom,".csv"))

library(Biostrings)
library(plyr)
library(stringr)

s <- Biostrings::readDNAStringSet(paste0("data/outputs/branchpoint_df_501_", chrom,".fa"))
transcript_id <- names(s)
seq=vector()
for (i in seq(along=s)){
  seq[i] <- toString(s[i])
  if(i %%1000==0){
    message(paste(i, "sequences processed"))
  }
}
df <- data.frame(transcript_id,seq)

#get nucleotide at position -5 to +5 (relative to BP)
#already stranded correctly

seq_pos0=substr(df[,2],251,251)
seq_pos1=substr(df[,2],252,252)
seq_pos2=substr(df[,2],253,253)
seq_pos3=substr(df[,2],254,254)
seq_pos4=substr(df[,2],255,255)
seq_pos5=substr(df[,2],256,256)
seq_neg1=substr(df[,2],250,250)
seq_neg2=substr(df[,2],249,249)
seq_neg3=substr(df[,2],248,248)
seq_neg4=substr(df[,2],247,247)
seq_neg5=substr(df[,2],246,246)

branchpoint_df_x$X=NULL

#find canonical AG splice dinucleotides 
f=(gregexpr("AG",substr(df[,2], 252,501),perl=TRUE))
le=lapply(f, length)

#get position of 1st-5th canonical SS
canon_hit1=vector()
canon_hit2=vector()
canon_hit3=vector()
canon_hit4=vector()
canon_hit5=vector()
num_canon=vector()
for(i in 1:length(f)){
  canon=sort(f[[i]])
  num_canon[i]=length(canon)
  canon_hits=vector()
  for(j in 1:5){
    if(j <= num_canon[i]){
      canon_hits=append(canon_hits, canon[j])
    }else{
      #if no hit distance=300
      canon_hits=append(canon_hits, 300)
    }
  }
  canon_hit1[i]=canon_hits[1]
  canon_hit2[i]=canon_hits[2]
  canon_hit3[i]=canon_hits[3]
  canon_hit4[i]=canon_hits[4]
  canon_hit5[i]=canon_hits[5]
}

#get dist to and length of largest ppt
if(exists("pyra_df")){rm(pyra_df)}
for(i in 1:length(f)){
  line=get_ppt(i, df, branchpoint_df_x)
  if(exists("pyra_df")){
    pyra_df=rbind(pyra_df, line)
  }else{
    pyra_df=line
  }
}

pyra_df=as.data.frame(pyra_df)
colnames(pyra_df)=c("ppt_start","ppt_run_length","ppt_ppPercent")
branchpoint_df_x=cbind(branchpoint_df_x, pyra_df, canon_hit1, canon_hit2,
                       canon_hit3,canon_hit4, canon_hit5, seq_neg5,seq_neg4,seq_neg3,
                       seq_neg2,seq_neg1,seq_pos0,seq_pos1,seq_pos2,seq_pos3,seq_pos4,seq_pos5)

write.csv(branchpoint_df_x, file=paste0("data/outputs/branchpoint_df_with_seq_r2_", chrom,".csv"))

