################################################################################
#                Annotate skipped and control exon triplets                    #
################################################################################

options(stringsAsFactors = FALSE)
library(data.table)
library(branchpointer)
library(stringr)
library(plyr)

#read in junction read counts
flux_junction <- as.data.frame(fread("data/GTeX/GTEx_Analysis_V4_RNA-seq_Flux1.6_junction_reads.txt"))
flux_junction <- flux_junction[-which(flux_junction$TargetID=="TargetID"),]
exons <- readExonAnnotation("data/genome_annotations/gencode.v19.annotation.gtf")
exons$newid <- with(exons, paste(chromosome, start,end,strand, sep="_"))

#remove any exons which are first/last
internal_exons <- arrange(exons, chromosome, gene_id, transcript_id,plyr::desc(exon_number))
internal_exons <- internal_exons[which(duplicated(internal_exons$transcript_id)),]
internal_exons <- internal_exons[-which(internal_exons$exon_number == 1),]
#internal_exons <- internal_exons[-which(duplicated(internal_exons$newid)),]

flux_junction$flux_start <- as.numeric(sapply(str_split(flux_junction$TargetID, "_"), "[[", 2))
flux_junction$flux_end <- as.numeric(sapply(str_split(flux_junction$TargetID, "_"), "[[", 3))

m <- match(str_sub(flux_junction$Gene_Symbol,1,15), str_sub(exons$gene_id,1,15))
flux_junction$strand <- exons$strand[m]

#convert to numeric value
for(i in 5:2920){
  flux_junction[,i] <- as.numeric(flux_junction[,i])
}

###### skipped exons ######

#junctions covering an exon

chroms <- unique(flux_junction$Chr)
#remove blank
chroms <- chroms[-2]

for(c in chroms){
  
  w <- which(flux_junction$Chr==c)
  message(c)
  #pos strand
  #keep junctions with duplicated starts and ends
  # ::|-------------------|::  - original
  # ::|----------|::          - same start
  #                 ::|---|:: - same end
  #              |::::|       - skipped exon
  
  dups_e <- which(duplicated(flux_junction$flux_end[w]) & (flux_junction$strand == "+")[w])
  dups_s <- which(duplicated(flux_junction$flux_start[w]) & (flux_junction$strand == "+")[w])
  keep <- which(dups_e %in% dups_s)
  
  #get read counts
  flux_junction_w <- flux_junction[w,]
  
  start <- flux_junction_w$flux_start[dups_e[keep]]
  end <- flux_junction_w$flux_end[dups_e[keep]]
    
  flux_junction_w$match_start <- match(flux_junction_w$flux_start, start)
  flux_junction_w$match_end <- match(flux_junction_w$flux_end, end)
  flux_junction_w$max_reads <- apply(flux_junction_w[,5:2920], 1, max)
  flux_junction_w$match <- flux_junction_w$match_start
  flux_junction_w$match[which((flux_junction_w$match_start != flux_junction_w$match_end) | is.na(flux_junction_w$match_start) |
                          is.na(flux_junction_w$match_end))] <- NA
  
  skipped_reads <- aggregate(max_reads ~ match, flux_junction_w, max)
  
  flux_junction_w$match <- NA
  flux_junction_w$match[which(flux_junction_w$flux_end < end[flux_junction_w$match_start])] <- 
    flux_junction_w$match_start[which(flux_junction_w$flux_end < end[flux_junction_w$match_start])] 
  
  reads1 <- aggregate(max_reads ~ match, flux_junction_w, max)
  
  flux_junction_w$match <- NA
  flux_junction_w$match[which(flux_junction_w$flux_start > start[flux_junction_w$match_end])] <- 
    flux_junction_w$match_end[which(flux_junction_w$flux_start > start[flux_junction_w$match_end])] 
  
  reads2 <- aggregate(max_reads ~ match, flux_junction_w, max)
  
  keep_df <- data.frame(n=1:length(keep))
  keep_df$skipped_reads <- skipped_reads$max_reads[match(keep_df$n, skipped_reads$match)]
  keep_df$reads1 <- reads1$max_reads[match(keep_df$n, reads1$match)]
  keep_df$reads2 <- reads2$max_reads[match(keep_df$n, reads2$match)]
  
  y <- which(keep_df$skipped_reads > 10 & keep_df$reads1 >10 & keep_df$reads2 >10)
  skipped_exons <- rep(NA, length(keep_df$n))
  
  #make sure 1 skipped exon and not multiple
  for(n in y){
    start <- flux_junction_w$flux_start[dups_e[keep[n]]]
    end <- flux_junction_w$flux_end[dups_e[keep[n]]]
    skipped_exons[n] <- length(which(internal_exons[-which(duplicated(internal_exons$newid)),]$start > (start+1) & 
                                       internal_exons[-which(duplicated(internal_exons$newid)),]$end < (end-1) & 
                                       internal_exons[-which(duplicated(internal_exons$newid)),]$chromosome == paste0("chr",c)))
  }
  
  y <- which(keep_df$skipped_reads > 10 & keep_df$reads1 >10 & keep_df$reads2 >10 & skipped_exons == 1)
  
  if(exists("junctions_with_se")){
    junctions_with_se <- rbind(junctions_with_se,flux_junction_w[dups_e[keep[y]],c(1:4,2921:2923)])
  }else{
    junctions_with_se <- flux_junction_w[dups_e[keep[y]],c(1:4,2921:2923)]
  }
  
  #uniques
  ids <- internal_exons[internal_exons$chromosome == paste0("chr",c) & internal_exons$strand == "+",]$newid
  
  ids <- ids[which(!(ids %in% ids[duplicated(ids)]))]
  
  x <- which(internal_exons$newid %in% ids)
  #no alternative annotated exons
  internal_exons_notskipped <- internal_exons[x,][!(internal_exons$start[x] %in% internal_exons$start[x][duplicated(internal_exons$start[x])]) &
                       !(internal_exons$end[x] %in% internal_exons$end[x][duplicated(internal_exons$end[x])]) ,]
  
  #find constitutive exons with at least 10 reads supporting junctions and none supporting skipping
  constit_exons <- vector()
  for(n in seq_along(internal_exons_notskipped$start)){
    l1 <- length(flux_junction_w[which(flux_junction_w$flux_end == internal_exons_notskipped$start[n]),]$flux_start)
    l2 <- length(flux_junction_w[which(flux_junction_w$flux_start == internal_exons_notskipped$end[n]),]$flux_end)
      if(l1 == 1 & l2 == 1){
      l1 <- length(which(flux_junction_w$flux_end == flux_junction_w[which(flux_junction_w$flux_start == internal_exons_notskipped$end[n]),]$flux_end))
      l2 <- length(which(flux_junction_w$flux_start == flux_junction_w[which(flux_junction_w$flux_end == internal_exons_notskipped$start[n]),]$flux_start))
      
      if(l1 == 1 & l2 == 1){
        m1 <- max(flux_junction_w[(which(flux_junction_w$flux_end == flux_junction_w[which(flux_junction_w$flux_start == internal_exons_notskipped$end[n]),]$flux_end)),5:2920])
        m2 <- max(flux_junction_w[(which(flux_junction_w$flux_start == flux_junction_w[which(flux_junction_w$flux_end == internal_exons_notskipped$start[n]),]$flux_start)),5:2920])
        if(m1 > 10 & m2 > 10){
          constit_exons <- append(constit_exons, internal_exons_notskipped$newid[n])
        }
      }
      #message(n)
    }
  }
  m <- match(exons$newid, constit_exons)
  constit_exons <- constit_exons[!(constit_exons %in% exons$newid[which(duplicated(m) & !is.na(m))])]  
  m <- match(constit_exons,exons$newid)
  
  if(exists("junctions_without_se")){
    junctions_without_se <- rbind(junctions_without_se,exons[c(m-1,m, m+1),])
  }else{
    junctions_without_se <- exons[c(m-1,m, m+1),]
  }
  
  #same on negative strand
  
  dups_e <- which(duplicated(flux_junction_w$flux_end) & (flux_junction_w$strand == "-"))
  dups_s <- which(duplicated(flux_junction_w$flux_start) & (flux_junction_w$strand == "-"))
  keep <- which(dups_e %in% dups_s)
  
  start <- flux_junction_w$flux_start[dups_e[keep]]
  end <- flux_junction_w$flux_end[dups_e[keep]]
  
  flux_junction_w$match_start <- match(flux_junction_w$flux_start, start)
  flux_junction_w$match_end <- match(flux_junction_w$flux_end, end)
  flux_junction_w$max_reads <- apply(flux_junction_w[,5:2920], 1, max)
  flux_junction_w$match <- flux_junction_w$match_start
  flux_junction_w$match[which((flux_junction_w$match_start != flux_junction_w$match_end) | is.na(flux_junction_w$match_start) |
                                is.na(flux_junction_w$match_end))] <- NA
  
  skipped_reads <- aggregate(max_reads ~ match, flux_junction_w, max)
  
  flux_junction_w$match <- NA
  flux_junction_w$match[which(flux_junction_w$flux_end < end[flux_junction_w$match_start])] <- 
    flux_junction_w$match_start[which(flux_junction_w$flux_end < end[flux_junction_w$match_start])] 
  
  reads1 <- aggregate(max_reads ~ match, flux_junction_w, max)
  
  flux_junction_w$match <- NA
  flux_junction_w$match[which(flux_junction_w$flux_start > start[flux_junction_w$match_end])] <- 
    flux_junction_w$match_end[which(flux_junction_w$flux_start > start[flux_junction_w$match_end])] 
  
  reads2 <- aggregate(max_reads ~ match, flux_junction_w, max)
  
  keep_df <- data.frame(n=1:length(keep))
  keep_df$skipped_reads <- skipped_reads$max_reads[match(keep_df$n, skipped_reads$match)]
  keep_df$reads1 <- reads1$max_reads[match(keep_df$n, reads1$match)]
  keep_df$reads2 <- reads2$max_reads[match(keep_df$n, reads2$match)]
  
  y <- which(keep_df$skipped_reads > 10 & keep_df$reads1 >10 & keep_df$reads2 >10)
  skipped_exons <- rep(NA, length(keep_df$n))

  for(n in y){
    start <- flux_junction_w$flux_start[dups_e[keep[n]]]
    end <- flux_junction_w$flux_end[dups_e[keep[n]]]
    skipped_exons[n] <- length(which(internal_exons$start > (start+1) & internal_exons$end < (end-1) & internal_exons$chromosome == paste0("chr",c)))
  }
  
  y <- which(keep_df$skipped_reads > 10 & keep_df$reads1 >10 & keep_df$reads2 >10 & skipped_exons == 1)
  
  if(exists("junctions_with_se")){
    junctions_with_se <- rbind(junctions_with_se,flux_junction_w[dups_e[y],c(1:4,2921:2923)])
  }else{
    junctions_with_se <- flux_junction_w[dups_e[y],c(1:4,2921:2923)]
  }

  #constitutive exons (negative strand)
  ids <- internal_exons[internal_exons$chromosome == paste0("chr",c) & internal_exons$strand == "-",]$newid
  ids <- ids[which(!(ids %in% ids[duplicated(ids)]))]
  x <- which(internal_exons$newid %in% ids)
  internal_exons_notskipped <- internal_exons[x,][!(internal_exons$start[x] %in% internal_exons$start[x][duplicated(internal_exons$start[x])]) &
                                                    !(internal_exons$end[x] %in% internal_exons$end[x][duplicated(internal_exons$end[x])]) ,]
  
  constit_exons <- vector()
  for(n in seq_along(internal_exons_notskipped$start)){
    l1 <- length(flux_junction_w[which(flux_junction_w$flux_end == internal_exons_notskipped$start[n]),]$flux_start)
    l2 <- length(flux_junction_w[which(flux_junction_w$flux_start == internal_exons_notskipped$end[n]),]$flux_end)
    if(l1 == 1 & l2 == 1){
      l1 <- length(which(flux_junction_w$flux_end == flux_junction_w[which(flux_junction_w$flux_start == internal_exons_notskipped$end[n]),]$flux_end))
      l2 <- length(which(flux_junction_w$flux_start == flux_junction_w[which(flux_junction_w$flux_end == internal_exons_notskipped$start[n]),]$flux_start))
      
      if(l1 == 1 & l2 == 1){
        m1 <- max(flux_junction_w[(which(flux_junction_w$flux_end == flux_junction_w[which(flux_junction_w$flux_start == internal_exons_notskipped$end[n]),]$flux_end)),5:2920])
        m2 <- max(flux_junction_w[(which(flux_junction_w$flux_start == flux_junction_w[which(flux_junction_w$flux_end == internal_exons_notskipped$start[n]),]$flux_start)),5:2920])
        if(m1 > 10 & m2 > 10){
          constit_exons <- append(constit_exons, internal_exons_notskipped$newid[n])
        }
      }
      #message(n)
    }
  }
  m <- match(exons$newid, constit_exons)
  constit_exons <- constit_exons[!(constit_exons %in% exons$newid[which(duplicated(m) & !is.na(m))])]  
  m <- match(constit_exons,exons$newid)
  
  if(exists("junctions_without_se")){
    junctions_without_se <- rbind(junctions_without_se,exons[c(m-1,m, m+1),])
  }else{
    junctions_without_se <- exons[c(m-1,m, m+1),]
  }
}

exons$exon_number <- as.numeric(exons$exon_number)

# match junctions to exon annotations
for(n in seq_along(junctions_with_se$TargetID)){
  
  exons1 <- exons[exons$start == junctions_with_se$flux_end[n] & 
                    exons$chromosome == paste0("chr", junctions_with_se$Chr[n]) & 
                    exons$strand==junctions_with_se$strand[n],]
  exons2 <- exons[exons$end == junctions_with_se$flux_start[n]& 
                    exons$chromosome == paste0("chr", junctions_with_se$Chr[n]) & 
                    exons$strand==junctions_with_se$strand[n],]
  
  exons1 <- exons1[exons1$transcript_id %in% exons2$transcript_id,]
  exons2 <- exons2[exons2$transcript_id %in% exons1$transcript_id,]
 
  tid <- exons2$transcript_id[which(abs(exons1$exon_number[match(exons1$transcript_id, exons2$transcript_id)] - exons2$exon_number) == 2)][1]
  if(!is.na(tid)){
    maxn <- (max(c(exons2$exon_number[exons2$transcript_id==tid],exons1$exon_number[exons1$transcript_id==tid])))
    
    subset <- exons[exons$transcript_id==tid & (exons$exon_number %in% (maxn-2):maxn),]
    
    subset$annotation <- "skipped_exon"
    
    if(exists("junctions_with_se_e")){
      junctions_with_se_e <- rbind(junctions_with_se_e, subset)
    }else{
      junctions_with_se_e <- subset
    }
  }
  message(n)
}

unique_transcripts <- unique(junctions_with_se_e$transcript_id)

for(t in seq_along(unique_transcripts)){
  #unique_transcripts[t]
  junctions_with_se_e_part <- arrange(junctions_with_se_e[junctions_with_se_e$transcript_id== unique_transcripts[t],], exon_number)
  junctions_with_se_e_part$annotation_2 <- NA
    
  #sets of 3
  junctions_with_se_e_part$annotation_2[1] <- "E1"
  junctions_with_se_e_part$annotation_2[match((junctions_with_se_e_part$exon_number[1] + 1), junctions_with_se_e_part$exon_number)] <- "SE"
  junctions_with_se_e_part$annotation_2[match((junctions_with_se_e_part$exon_number[1] + 2), junctions_with_se_e_part$exon_number)] <- "E2"
  
  while(any(is.na(junctions_with_se_e_part$annotation_2))){
    ind <- which(is.na(junctions_with_se_e_part$annotation_2))
    junctions_with_se_e_part$annotation_2[ind][1] <- "E1"
    junctions_with_se_e_part$annotation_2[ind][match((junctions_with_se_e_part$exon_number[ind][1] + 1), junctions_with_se_e_part$exon_number[ind])] <- "SE"
    junctions_with_se_e_part$annotation_2[ind][match((junctions_with_se_e_part$exon_number[ind][1] + 2), junctions_with_se_e_part$exon_number[ind])] <- "E2"
  }
  if(exists("junctions_with_se_e_parts")){
    junctions_with_se_e_parts <- rbind(junctions_with_se_e_parts, junctions_with_se_e_part)
  }else{
    junctions_with_se_e_parts <- junctions_with_se_e_part
  }
}
switch1 <- which(junctions_with_se_e_parts$strand =="-" & junctions_with_se_e_parts$annotation_2=="E1")
switch2 <- which(junctions_with_se_e_parts$strand =="-" & junctions_with_se_e_parts$annotation_2=="E2")

junctions_with_se_e_parts$annotation_2[switch1] <- "E2"
junctions_with_se_e_parts$annotation_2[switch2] <- "E1"

unique_transcripts <- unique(junctions_without_se$transcript_id)

junctions_without_se$exon_number <- as.numeric(junctions_without_se$exon_number)

for(t in seq_along(unique_transcripts)){
  #unique_transcripts[t]
  junctions_part <- arrange(junctions_without_se[junctions_without_se$transcript_id== unique_transcripts[t],], exon_number)
  junctions_part$annotation_2 <- NA
  
  #sets of 3
  junctions_part$annotation_2[1] <- "E1"
  junctions_part$annotation_2[match((junctions_part$exon_number[1] + 1), junctions_part$exon_number)] <- "SE"
  junctions_part$annotation_2[match((junctions_part$exon_number[1] + 2), junctions_part$exon_number)] <- "E2"
  
  while(any(is.na(junctions_part$annotation_2))){
    ind <- which(is.na(junctions_part$annotation_2))
    junctions_part$annotation_2[ind][1] <- "E1"
    junctions_part$annotation_2[ind][match((junctions_part$exon_number[ind][1] + 1), junctions_part$exon_number[ind])] <- "SE"
    junctions_part$annotation_2[ind][match((junctions_part$exon_number[ind][1] + 2), junctions_part$exon_number[ind])] <- "E2"
  }
  if(exists("junctions_without_se_parts")){
    junctions_without_se_parts <- rbind(junctions_without_se_parts, junctions_part)
  }else{
    junctions_without_se_parts <- junctions_part
  }
}

switch1 <- which(junctions_without_se_parts$strand =="-" & junctions_without_se_parts$annotation_2=="E1")
switch2 <- which(junctions_without_se_parts$strand =="-" & junctions_without_se_parts$annotation_2=="E2")

junctions_without_se_parts$annotation_2[switch1] <- "E2"
junctions_without_se_parts$annotation_2[switch2] <- "E1"
junctions_without_se_parts$annotation = "control"

junctions_without_se_parts <- junctions_without_se_parts[,c(1:12,14,13)]

skipped_exons_annotations <- rbind(junctions_with_se_e_parts,junctions_without_se_parts)
write.table(skipped_exons_annotations, file="data/skipped_exons_and_control_annotations_gencode19.txt", sep="\t", quote=FALSE, row.names=FALSE)
