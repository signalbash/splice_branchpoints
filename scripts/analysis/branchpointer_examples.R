library(branchpointer)
library(GenomicRanges)
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

exons <- gtfToExons("gencode.v26.annotation.gtf")

# functions and data from branchpointer
load("branchpointer/R/sysdata.rda")

getCanonical3SS <- function(ag){
  ag <- sort(ag)
  #get first five matches
  ag <- ag[1:5]
  #if less than 5 matches, replace remainder with 300
  ag[is.na(ag)] <- 300
  
  return(ag)
}
getPPT <- function(attributes){
  
  dist3prime <- as.numeric(attributes$to_3prime_point)
  
  #get sequence between the tested site and 3'exon
  seq <- substr(attributes$seq, 251, (250 + dist3prime))
  seq <- unlist(strsplit(seq, ""))
  
  pyramidines <- which(seq == "T" | seq == "C")
  
  if(length(pyramidines) >1){
    pyraDist <- rep(NA, length(pyramidines)-1)
    
    #get distance to the next pyramidine
    for(p in 2:length(pyramidines)){
      pyraDist[p-1] <- pyramidines[p] - pyramidines[p-1]
    }
    
    longestRun <- 0
    bestRun <- 0
    bestStart <- 0
    percentPyra <- 0
    start <- 1
    sets <- length(pyraDist) + 1
    startVec <- rep(NA, sets)
    runVec <- rep(NA, sets)
    percentPyraVec <- rep(NA, sets)
    
    #find best PPT
    for(p in seq_along(pyraDist)){
      
      #For the first instance of at least 2 sequential pyramidines
      if(pyraDist[p] == 1 & longestRun == 0){
        longestRun <- 1
        start <- p
        seq.p <- seq[(pyramidines[start]):(pyramidines[start+longestRun])]
        percentPyra <- length(which(seq.p=="C"|seq.p=="T"))/length(seq.p)
        
        if(longestRun > bestRun){
          bestRun <- longestRun
          bestStart <- start
          bestPercent <- percentPyra
        }
        
        
        #allowing gaps with one purine, grow the tract
      }else if(pyraDist[p] == 1 | pyraDist[p] == 2){
        longestRun <- longestRun + 1
        seq.p <- seq[(pyramidines[start]):(pyramidines[start+longestRun])]
        percentPyra <- length(which(seq.p == "C"|seq.p == "T"))/length(seq.p)
        
        if(longestRun > bestRun){
          bestRun <- longestRun
          bestStart <- start
          bestPercent <- percentPyra
        }
        
        #once too many purines are encontered, save tract
      }else{
        startVec[p] <- bestStart
        runVec[p] <- longestRun
        percentPyraVec[p] <- percentPyra
        longestRun <- 0
        start <- p + 1
        bestStart <- start
      }
    }
    
    #add the last PPT
    startVec[sets] <- start
    runVec[sets] <- longestRun
    percentPyraVec[sets] <- percentPyra
    
    #make into data.frame
    pyra.df <- data.frame(set=c(1:length(runVec)), runVec,startVec,percentPyraVec)
    pyra.df <- pyra.df[which(!is.na(pyra.df$startVec)),]
    pyra.df <- plyr::arrange(pyra.df, plyr::desc(runVec))
    
    bestStart <- pyra.df$startVec[which.max(pyra.df$runVec)]
    bestRun <- pyra.df$runVec[which.max(pyra.df$runVec)]
    
    #If PPTs are long enough, shorten to get pyramidine content above 80%
    pyra.df <- pyra.df[pyra.df$runVec >=10,]
    if(dim(pyra.df)[1] != 0){
      run.extend <- rep(NA, length(pyra.df[,1]))
      start.extend <- rep(NA, length(pyra.df[,1]))
      polyper.extend <- rep(NA, length(pyra.df[,1]))
      
      for(py in seq_along(length(pyra.df[,1]))){
        seq.p <- seq[(pyramidines[pyra.df$startVec[py]]):(pyramidines[pyra.df$startVec[py]+pyra.df$runVec[py]])]
        startp <- 1
        
        run <- length(seq.p)
        
        polyper.extend[py] <- pyra.df$percentPyraVec[py]
        while(polyper.extend[py] < 0.8 & run >= 10){
          purine <- which(seq.p == "A"| seq.p == "G")
          purine.neg <- (length(seq.p) +1 -purine)*-1
          purine.both <- c(purine,purine.neg)
          minPurine <- which.min(abs(purine.both))
          
          point <- purine.both[minPurine]
          
          if(point < 0){
            rm <- length(seq.p)-point*-1
            newSeq <- seq.p[1:rm]
            run <- length(newSeq)
          }else{
            newSeq <- seq.p[(point+1):length(seq.p)]
            run <- length(newSeq)
            startp <- startp+point
          }
          polyper.extend[py] <-
            length(which(newSeq == "C" | newSeq == "T")) / length(newSeq)
          seq.p <- newSeq
        }
        
        start.extend[py] <- (pyramidines[pyra.df$startVec[py]])+startp-1
        run.extend[py] <- run
      }
      
      bestStart <- start.extend[which.max(run.extend)]
      bestRun <- run.extend[which.max(run.extend)]
    }
    
    
    line <- c(bestStart,bestRun)
    
    #if only 1 pyramidine
  }else if(length(pyramidines) == 1){
    line <- c(pyramidines, 1)
    #if no pyramidines
  }else{
    line <- c(0,0)
  }
  return(line)
}

getBranchpointSequence <- function(query, uniqueId = "test",
                                   queryType,
                                   workingDirectory = ".",
                                   genome = NA,
                                   bedtoolsLocation = NA,
                                   BSgenome = NULL,
                                   useParallel = FALSE,
                                   cores = 1,
                                   rmChr = FALSE) {
  
  if(is.na(genome) & is.null(BSgenome)){
    stop("please specify a genome .fa file for sequence extraction or specify a BSgenome object")
  }
  
  if((is.na(genome) | is.na(bedtoolsLocation))& is.null(BSgenome)){
    stop("please specify a genome .fa file for sequence extraction and a bedtools binary location")
  }
  
  if(useParallel){
    maxCores <- parallel::detectCores()
    
    if(maxCores < cores){
      message(paste0("specified cores (", cores,") is greater than available cores(", maxCores,")"))
      message(paste0("using all available cores"))
      cores <- maxCores
    }
  }
  
  #make bed format file
  if(missing(queryType) | !(queryType %in% c("SNP", "region"))){
    stop("please specify queryType as \"region\" or \"SNP\"")
  }else if (queryType == "SNP") {
    bed <- query
    start(ranges(bed)) <- as.integer(start(ranges(bed)) - 1)
  }else if (queryType == "region") {
    bed <- query
    start(ranges(bed))[as.logical(strand(bed)== "+")] <-
      end(ranges(bed))[as.logical(strand(bed)== "+")]
    start(ranges(bed)) <- start(ranges(bed)) - 1
    width(ranges(bed)) <- 2
  }
  
  #extend bed file to cover +/- 250 nt from each query point
  start(ranges(bed))[which(as.logical(strand(bed) == "+"))] <-
    as.integer(start(ranges(bed)) - 250 -
                 (44 - query$to_3prime))[which(as.logical(strand(bed) == "+"))]
  start(ranges(bed))[which(as.logical(strand(bed) == "-"))] <-
    as.integer(start(ranges(bed)) - 277 +
                 (44 - query$to_3prime))[which(as.logical(strand(bed) == "-"))]
  width(ranges(bed))  <- 529
  
  
  if(!is.na(genome) & !missing(bedtoolsLocation)){
    #convert to .fasta using bedtools
    bedTable <- data.frame(seqnames(bed),
                           start(ranges(bed)),
                           end(ranges(bed)),
                           bed$id,
                           score=0,
                           strand(bed))
    
    if (rmChr == TRUE) {
      bedTable[,1] <- gsub("chr","", bedTable[,1])
    }
    
    utils::write.table(
      bedTable, sep = "\t", file = paste0(workingDirectory,"/mutation_",
                                          uniqueId,".bed"),
      row.names = FALSE,col.names = FALSE,quote = FALSE)
    cmd <- paste0(
      bedtoolsLocation," getfasta -fi ", genome,
      " -bed ",workingDirectory,"/mutation_",uniqueId,".bed -fo ",
      workingDirectory,"/mutation_",uniqueId,".fa -name -s")
    system(cmd)
    
    info <- file.info(paste0(workingDirectory,"/mutation_",uniqueId,".fa"))
    if(info$size > 0){
      fasta <-
        data.table::fread(paste0(workingDirectory,"/mutation_",uniqueId,".fa"),
                          header = FALSE, stringsAsFactors = FALSE)
      fasta <- as.data.frame(fasta)
      system(paste0("rm -f ",workingDirectory,"/mutation_",uniqueId,"*"))
      
      s <- fasta[seq(2,dim(fasta)[1],by = 2),1]
      mcols(query)$seq <- s
    }else{
      stop("cannot retrieve fa sequence, please check your exon annotation and genome .fa")
    }
    
  }else{
    # need to +1 for BSgenomes sequence retreval
    # ranges given are bedtools (legacy)
    
    start(ranges(bed)) <- start(ranges(bed)) +1
    width(ranges(bed))  <- 528
    
    # check if chromosomes need to renamed
    seqlevels.genome <- GenomeInfoDb::seqlevels(BSgenome)
    seqlevels.bed <- GenomeInfoDb::seqlevels(bed)
    
    # used chromosomes 
    seqnames.bed <- as.character(GenomeInfoDb::seqnames(bed))
    
    if(!(all(seqnames.bed %in% seqlevels.genome))){
      
      # try adding/removing chr from bed
      if(!all(stringr::str_sub(seqlevels.bed, 1, 3) == "chr")){
        GenomeInfoDb::seqlevels(bed) <- paste0("chr", seqlevels.bed)
      }else if(all(stringr::str_sub(seqlevels.bed, 1, 3) == "chr")){
        GenomeInfoDb::seqlevels(bed) <- gsub("chr", "", seqlevels.bed)
      }
      
      seqlevels.bed <- GenomeInfoDb::seqlevels(bed)
      seqnames.bed <- as.character(GenomeInfoDb::seqnames(bed))
      
      # if that doesn't fix the issue, break
      if(!(all(seqnames.bed %in% seqlevels.genome))){
        stop("Chromosome names of query and genome do not match")
      }
    }
    
    bed.seq <- suppressWarnings(Biostrings::getSeq(BSgenome, bed))
    mcols(query)$seq <- as.character(bed.seq)
  }
  
  ##mutate at SNP location
  if (queryType == "SNP") {
    #location of SNP
    loc <- 44 - query$to_3prime
    nt.ref <- as.character(query$ref_allele)
    nt.alt <- as.character(query$alt_allele)
    
    #change to compliment if on negative strand
    nt.ref[which(as.logical(strand(query) == "-"))] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(nt.ref[which(as.logical(strand(query) == "-"))])))
    nt.alt[which(as.logical(strand(query) == "-"))] <-
      as.character(Biostrings::complement(Biostrings::DNAStringSet(nt.alt[which(as.logical(strand(query) == "-"))])))
    
    #check ref allele
    refAlleleCorrect <- substr(query$seq, 251 + (loc),251 + (loc)) == nt.ref
    
    if (any(!refAlleleCorrect)) {
      rm <- which(refAlleleCorrect == FALSE)
      if (all(refAlleleCorrect[which(as.logical(strand(query) == "-"))] == FALSE) &
          all(refAlleleCorrect[which(as.logical(strand(query) == "+"))])) {
        message("reference alleles are incorrect for all negative strand introns")
        message("please input alleles as positive strand sequences")
      }else{
        message("reference alleles do not match sequence for:")
        message(paste(query$id[rm], collapse = "\n"))
      }
      message("removing from analysis")
      query <- query[-rm]
      nt.ref <- nt.ref[-rm]
      nt.alt <- nt.alt[-rm]
    }
  }
  
  #create 501nt sequences with each query point centered at 251
  seqs <- data.frame(seq = rep(query$seq, 27))
  seqs$i <- rep(18:44, each = length(query$seq))
  seqs <- apply(seqs, 1, function(x) substr(x[1],
                                            (as.numeric(x[2])-17),
                                            (as.numeric(x[2])-17) + 500))
  
  if (queryType == "SNP") {
    #create mutated sequence
    s.mut <-
      paste0(substr(query$seq, 1,250 + (loc)), nt.alt,
             substr(query$seq, 252 + (loc),528))
    
    
    seqs.mut <- data.frame(seq = rep(s.mut, 27))
    seqs.mut$i <- rep(18:44, each = length(s.mut))
    seqs.mut <- apply(seqs.mut, 1, function(x) substr(x[1],
                                                      (as.numeric(x[2])-17),
                                                      (as.numeric(x[2])-17) + 500))
    
    queryAllPoints <- rep(query, 27)
    queryAllPoints <- rep(queryAllPoints, 2)
    
    mcols(queryAllPoints)$status <- c(rep("REF", length(seqs)),
                                      rep("ALT", length(seqs.mut)))
    
    queryAllPoints$seq <- c(seqs, seqs.mut)
    mcols(queryAllPoints)$to_3prime_point <-
      rep(rep(44:18,each = length(query$id)),2)
  }else{
    queryAllPoints <- do.call("c",as.list(rep(query, 27)))
    mcols(queryAllPoints)$status <- rep("REF", length(seqs))
    queryAllPoints$seq <- seqs
    mcols(queryAllPoints)$to_3prime_point <-
      rep(44:18,each = length(query$seq))
  }
  
  mcols(queryAllPoints)$to_5prime_point <-
    (queryAllPoints$to_3prime + queryAllPoints$to_5prime) -
    queryAllPoints$to_3prime_point
  
  testSite <- start(ranges(queryAllPoints))
  posStrand <- which(as.logical(strand(queryAllPoints) == "+"))
  negStrand <- which(as.logical(strand(queryAllPoints) == "-"))
  testSite[posStrand] <- end(ranges(queryAllPoints))[posStrand]
  testSite[posStrand] <- (testSite + queryAllPoints$to_3prime -
                            queryAllPoints$to_3prime_point)[posStrand]
  testSite[negStrand] <- (testSite - queryAllPoints$to_3prime +
                            queryAllPoints$to_3prime_point)[negStrand]
  mcols(queryAllPoints)$test_site <- testSite
  
  #get sequence identity at position -5 to +5 relative to testing point
  mcols(queryAllPoints)$seq_pos0 <-
    factor(substr(queryAllPoints$seq,251,251), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_pos1 <-
    factor(substr(queryAllPoints$seq,252,252), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_pos2 <-
    factor(substr(queryAllPoints$seq,253,253), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_pos3 <-
    factor(substr(queryAllPoints$seq,254,254), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_pos4 <-
    factor(substr(queryAllPoints$seq,255,255), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$ seq_pos5 <-
    factor(substr(queryAllPoints$seq,256,256), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_neg1 <-
    factor(substr(queryAllPoints$seq,250,250), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_neg2 <-
    factor(substr(queryAllPoints$seq,249,249), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_neg3 <-
    factor(substr(queryAllPoints$seq,248,248), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_neg4 <-
    factor(substr(queryAllPoints$seq,247,247), levels = c("A","C","G","T"))
  mcols(queryAllPoints)$seq_neg5 <-
    factor(substr(queryAllPoints$seq,246,246), levels = c("A","C","G","T"))
  
  #find canonical AG splice dinucleotides
  f <- gregexpr("AG",substr(queryAllPoints$seq, 252,501),perl = TRUE)
  
  if (useParallel) {
    cluster <- parallel::makeCluster(cores)
    canonHits <- parallel::parLapply(cluster,f, getCanonical3SS)
    pyra <-
      parallel::parLapply(cluster,queryAllPoints, getPPT)
    parallel::stopCluster(cluster)
  }else{
    canonHits <- lapply(f, getCanonical3SS)
    pyra <- lapply(queryAllPoints, getPPT)
  }
  
  canon <- matrix(unlist(canonHits), ncol = 5, byrow = TRUE)
  canon <- as.data.frame(canon, stringsAsFactors=FALSE)
  colnames(canon) <-
    c("canon_hit1", "canon_hit2", "canon_hit3", "canon_hit4", "canon_hit5")
  
  mcols(queryAllPoints) <- cbind(mcols(queryAllPoints), canon)
  mcols(queryAllPoints)$ppt_start <- unlist(lapply(pyra, "[[", 1))
  mcols(queryAllPoints)$ppt_run_length <- unlist(lapply(pyra, "[[", 2))
  
  mcols(queryAllPoints)$seq <-
    stringr::str_sub(queryAllPoints$seq, (251 +
                                            queryAllPoints$to_3prime_point - 50),(250 +
                                                                                    queryAllPoints$to_3prime_point))
  
  return(queryAllPoints)
}

brca2 <- makeBranchpointWindowForExons(id="ENSE00000939171.1",idType = "exon_id", exons=exons)
brca2_attributes <- getBranchpointSequence(brca2,queryType = "region",BSgenome = g)


brca2_attributes.forModel <- as.data.frame(mcols(brca2_attributes))

colNamesForDataFrame <- HCN_names[-1]
colNamesForDataFrame <- gsub("new_ID", "id",
                             colNamesForDataFrame)
colNamesForDataFrame <- gsub("dist.1", "to_5prime_point",
                             colNamesForDataFrame)
colNamesForDataFrame <- gsub("dist.2", "to_3prime_point",
                             colNamesForDataFrame)

brca2_attributes.forModel <- brca2_attributes.forModel[,colNamesForDataFrame]
brca2_attributes.forModel <- cbind(set = "UNKNOWN", brca2_attributes.forModel)
colnames(brca2_attributes.forModel) <-  HCN_names

#remove any rows with NA values
isNA <- apply(brca2_attributes.forModel[,-c(1,2)], 2, is.na)
isNAv <- apply(isNA,1, any)
rm <- which(isNAv == TRUE)
if(length(rm) > 0){
  brca2_attributes.forModel <- brca2_attributes.forModel[-rm,]
  brca2_attributes <- brca2_attributes[-rm,]
}

#convert sequences to dummy vars
brca2_attributes.dummies <-
  predict(dummies, newdata = brca2_attributes.forModel[,-2])
brca2_attributes.forModel <-
  cbind(set = brca2_attributes.forModel$set, brca2_attributes.dummies)
brca2_attributes.forModel <- apply(brca2_attributes.forModel[,-1],2,as.numeric)

#make sure all dummy vars are accounted for
if(any(is.na(match(names(preProcValues$mean), colnames(brca2_attributes.forModel))))){
  newColNames <- names(preProcValues$mean)[which(is.na(match(names(preProcValues$mean),
                                                             colnames(brca2_attributes.forModel))))]
  dfAdd <- as.data.frame(matrix(data=0,nrow=nrow(brca2_attributes.forModel),
                                ncol=length(newColNames)))
  colnames(dfAdd) <- newColNames
  brca2_attributes.forModel <- cbind(brca2_attributes.forModel, dfAdd)
  brca2_attributes.forModel <-
    brca2_attributes.forModel[,match(names(preProcValues$mean),
                                    colnames(brca2_attributes.forModel))]
}

#pre-process values
brca2_attributes.forModel <- predict(preProcValues, brca2_attributes.forModel)
brca2_attributes.forModel <- cbind(brca2_attributes.forModel, Class="UNKNOWN")
brca2_attributes.forModel <-
  as.data.frame(brca2_attributes.forModel, stringsAsFactors=FALSE)
brca2_attributes.forModel$Class <- as.factor(brca2_attributes.forModel$Class)
for(n in 1:(length(colnames(brca2_attributes.forModel)) -1)){
  brca2_attributes.forModel[,n] <- as.numeric(brca2_attributes.forModel[,n])
}

#gbm prediction
p <- predict(object = branchpointer2.gbm,
             brca2_attributes.forModel,"prob")

#reconfigure
mcols(brca2_attributes)$branchpoint_prob <- p[,1]

#### U2 binding energy###
m <- match(colnames(U2_binding_df)[-c(1:3)],colnames(mcols(brca2_attributes)))
U2Eightmers <- do.call(paste0, as.data.frame(mcols(brca2_attributes)[,m]))

x <- match(U2Eightmers, U2_binding_df$eightmers)
mcols(brca2_attributes)$U2_binding_energy <- U2_binding_df$energy[x]
brca2_attributes.forModel
pheat_attributes_brca2 <- as.data.frame(brca2_attributes.forModel)
pheat_attributes_brca2 <- pheat_attributes_brca2[,c(1,2,3,4,5,6,7,30,25,19)]

###

mart.snp <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
query_rs587776767 <- makeBranchpointWindowForSNP("rs587776767", mart.snp = mart.snp,exons = exons, filter=F)

rs587776767_attributes <- getBranchpointSequence(query_rs587776767,queryType = "SNP",BSgenome = g)

rs587776767_attributes.forModel <- as.data.frame(mcols(rs587776767_attributes))

colNamesForDataFrame <- HCN_names[-1]
colNamesForDataFrame <- gsub("new_ID", "id",
                             colNamesForDataFrame)
colNamesForDataFrame <- gsub("dist.1", "to_5prime_point",
                             colNamesForDataFrame)
colNamesForDataFrame <- gsub("dist.2", "to_3prime_point",
                             colNamesForDataFrame)

rs587776767_attributes.forModel <- rs587776767_attributes.forModel[,colNamesForDataFrame]
rs587776767_attributes.forModel <- cbind(set = "UNKNOWN", rs587776767_attributes.forModel)
colnames(rs587776767_attributes.forModel) <-  HCN_names

#remove any rows with NA values
isNA <- apply(rs587776767_attributes.forModel[,-c(1,2)], 2, is.na)
isNAv <- apply(isNA,1, any)
rm <- which(isNAv == TRUE)
if(length(rm) > 0){
  rs587776767_attributes.forModel <- rs587776767_attributes.forModel[-rm,]
  rs587776767_attributes <- rs587776767_attributes[-rm,]
}

#convert sequences to dummy vars
rs587776767_attributes.dummies <-
  predict(dummies, newdata = rs587776767_attributes.forModel[,-2])
rs587776767_attributes.forModel <-
  cbind(set = rs587776767_attributes.forModel$set, rs587776767_attributes.dummies)
rs587776767_attributes.forModel <- apply(rs587776767_attributes.forModel[,-1],2,as.numeric)

#make sure all dummy vars are accounted for
if(any(is.na(match(names(preProcValues$mean), colnames(rs587776767_attributes.forModel))))){
  newColNames <- names(preProcValues$mean)[which(is.na(match(names(preProcValues$mean),
                                                             colnames(rs587776767_attributes.forModel))))]
  dfAdd <- as.data.frame(matrix(data=0,nrow=nrow(rs587776767_attributes.forModel),
                                ncol=length(newColNames)))
  colnames(dfAdd) <- newColNames
  rs587776767_attributes.forModel <- cbind(rs587776767_attributes.forModel, dfAdd)
  rs587776767_attributes.forModel <-
    rs587776767_attributes.forModel[,match(names(preProcValues$mean),
                                     colnames(rs587776767_attributes.forModel))]
}

#pre-process values
rs587776767_attributes.forModel <- predict(preProcValues, rs587776767_attributes.forModel)
rs587776767_attributes.forModel <- cbind(rs587776767_attributes.forModel, Class="UNKNOWN")
rs587776767_attributes.forModel <-
  as.data.frame(rs587776767_attributes.forModel, stringsAsFactors=FALSE)
rs587776767_attributes.forModel$Class <- as.factor(rs587776767_attributes.forModel$Class)
for(n in 1:(length(colnames(rs587776767_attributes.forModel)) -1)){
  rs587776767_attributes.forModel[,n] <- as.numeric(rs587776767_attributes.forModel[,n])
}

#gbm prediction
p <- predict(object = branchpointer2.gbm,
             rs587776767_attributes.forModel,"prob")

#reconfigure
mcols(rs587776767_attributes)$branchpoint_prob <- p[,1]

#### U2 binding energy###
m <- match(colnames(U2_binding_df)[-c(1:3)],colnames(mcols(rs587776767_attributes)))
U2Eightmers <- do.call(paste0, as.data.frame(mcols(rs587776767_attributes)[,m]))

x <- match(U2Eightmers, U2_binding_df$eightmers)
mcols(rs587776767_attributes)$U2_binding_energy <- U2_binding_df$energy[x]
rs587776767_attributes.forModel

pheat_attributes_rs587776767 <- as.data.frame(rs587776767_attributes.forModel)
pheat_attributes_rs587776767 <- pheat_attributes_rs587776767[,c(1,2,3,4,5,6,7,30,25,19)]

###

pheat_attributes_rs587776767 <- t(pheat_attributes_rs587776767)
pheat_attributes_rs587776767 <- cbind(pheat_attributes_rs587776767, c(rep(c(-2,2),5)))

###
pheat_attributes_brca2 <- as.data.frame(brca2_attributes.forModel)
pheat_attributes_brca2 <- pheat_attributes_brca2[,c(1,2,3,4,5,6,7,30,25,19)]

pheat_attributes_brca2 <- t(pheat_attributes_brca2)
pheat_attributes_brca2 <- cbind(pheat_attributes_brca2, c(rep(c(-2,2),5)))

rs587776767_attributes_df <- as.data.frame(rs587776767_attributes)
brca2_attributes_df <- as.data.frame(brca2_attributes)

save(pheat_attributes_brca2, pheat_attributes_rs587776767,rs587776767_attributes_df, brca2_attributes_df, 
     file="data/branchpointer_example_figures.Rdata")

  