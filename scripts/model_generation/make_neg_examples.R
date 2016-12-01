args <- commandArgs(trailingOnly = TRUE)

k1=as.character(args[1])
k2=as.character(args[2])

options(stringsAsFactors = F)

chrom=k1

message(paste0("Making negative examples for branchpoints on ",k1))
setwd(k2)

gencode.v12.condensed.exons.u=read.csv(file=paste0("data/outputs/gencode.v12.condensed.exons_", chrom,".csv"))

branchpoint_df_lc=read.csv(file=paste0("data/outputs/branchpoint_df_gene_lc_", chrom,".csv"))
branchpoint_df_hc=read.csv(file=paste0("data/outputs/branchpoint_df_gene_hc_", chrom,".csv"))

branchpoint_df_hc=cbind(branchpoint_df_hc, set="HC")
branchpoint_df_lc=cbind(branchpoint_df_lc, set="LC")

branchpoint_df_hc=branchpoint_df_hc[(branchpoint_df_hc$dist.2 >=18 & branchpoint_df_hc$dist.2 <=44),]

#data.frame with high confidence and low confidence branchpoints
branchpoint_df=rbind(branchpoint_df_hc,branchpoint_df_lc)

##Generate random set of negatives
if(exists("branchpoint_df_N")){rm(branchpoint_df_N)}
if(exists("branchpoint_df_N_large")){rm(branchpoint_df_N_large)}

#for each HC branchpoint
for(i in 1:length(branchpoint_df_hc$Chromosome)){
  #branchpoint_df[i,]
  pos_dist=branchpoint_df_hc$dist.2[i]
  neg_dist=branchpoint_df_hc$dist.1[i]

  #distance to 3'SS
  values=18:44
  
  for(j in seq(along=values)){
    line=branchpoint_df_hc[i,]
    if(line$Strand=="-"){
      line$Start=line$Start -(line$dist.2 - values[j])
      line$Stop=line$Stop - (line$dist.2 - values[j])
      line$dist.1=line$dist.1 + (line$dist.2 - values[j])
      line$dist.2=values[j]
    }else{
      line$Start=line$Start +(line$dist.2 - values[j])
      line$Stop=line$Stop + (line$dist.2 - values[j])
      line$dist.1=line$dist.1 + (line$dist.2 - values[j])
      line$dist.2=values[j]
    }
    line$set="NEG"
    
    #check if branchpoint made is in HC or LC branchpoint sets
    x=which(branchpoint_df$Start==line$Start & 
            branchpoint_df$Stop==line$Stop & 
            branchpoint_df$Chromosome==line$Chromosome &
            branchpoint_df$Strand==line$Strand)
    
    #if it is, remove it
    if(length(x)==0){
      if(exists("branchpoint_df_N")){
        branchpoint_df_N=rbind(branchpoint_df_N,line)
      }else{
      branchpoint_df_N=line
      }
    }
  }

  if(i %% 50==0){
    message(paste(i, "of", length(branchpoint_df_hc$Chromosome)))
    message(paste(length(branchpoint_df_N_large[,1]), " negative examples made."))
  }

  #every 10 exons, move the simulated branchpoint_df_N to a larger df (speeds up processing)
  if(i %% 10==0){
    if(exists("branchpoint_df_N_large")){
      branchpoint_df_N_large=rbind(branchpoint_df_N_large, branchpoint_df_N)
    }else{
      branchpoint_df_N_large=branchpoint_df_N
    }
    rm(branchpoint_df_N)
  }
}

#make sure all are put into branchpoint_df_N_large
if(exists("branchpoint_df_N_large") & exists("branchpoint_df_N")){
  branchpoint_df_N_large=rbind(branchpoint_df_N_large, branchpoint_df_N)
}else if(!exists("branchpoint_df_N_large") & exists("branchpoint_df_N")){
  branchpoint_df_N_large=branchpoint_df_N
}

write.csv(branchpoint_df_N_large, file=paste0("data/outputs/branchpoint_df_gene_N", chrom,".csv"))
