################################################################################
#                branchpointer SNP predicitons to site summary                 #
################################################################################


predictions_to_stats <- function(predictions, snp_ids, snp_info){
    BP_prob_cutoff <- 0.5
    U2_cutoff_weak <- 0
    BP_prob_change <- 0.2
    
    BP_REF_num <- vector()
    BP_ALT_num <- vector()
    dist_to_exon <- vector()
    dist_to_BP_REF <- vector()
    dist_to_BP_ALT <- vector()
    max_prob_REF <- vector()
    max_prob_ALT <- vector()
    max_U2_REF <- vector()
    max_U2_ALT <- vector()
    deleted_n <- vector()
    created_n <- vector()
    
    for(z in seq(along = snp_ids)){
      branchpoint_mutated <- predictions[which(!is.na(match(predictions$id, snp_ids[z]))),]
    
      info <- snp_info[which(!is.na(match(snp_info$id, snp_ids[z]))),]
      
      if(min(abs(branchpoint_mutated$end - info$chrom_start)) == 0){
        dist_to_exon[z] <- branchpoint_mutated$distance[which.min(abs(branchpoint_mutated$end - info$chrom_start))]
      }else{
        if(info$strand == "+"){
        exon_start <- branchpoint_mutated$end[1] + branchpoint_mutated$distance[1]
        dist_to_exon[z] <- exon_start - info$chrom_start
        }else{
          exon_start <- branchpoint_mutated$end[1] - branchpoint_mutated$distance[1]
          dist_to_exon[z] <- info$chrom_start - exon_start
        }
      }
    
      #get distances from SNP to BPs
      diffs <- dist_to_exon[z]-branchpoint_mutated$distance[branchpoint_mutated$allele_status=="REF"][which(
        branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="REF"] >= BP_prob_cutoff)]
      
      BP_REF_num[z] <- length(diffs)
      
      if(BP_REF_num[z] > 0){
      dist_to_BP_REF[z] <- diffs[which.min(abs(diffs))]
      }else{
        dist_to_BP_REF[z] <- NA
      }
      
      diffsC <- dist_to_exon[z]-branchpoint_mutated$distance[branchpoint_mutated$allele_status=="ALT"][which(
        branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="ALT"] >= BP_prob_cutoff)]
      
      BP_ALT_num[z] <- length(diffsC)
      
      if(BP_ALT_num[z] >0){
        dist_to_BP_ALT[z] <- diffsC[which.min(abs(diffsC))]
      }else{
        dist_to_BP_ALT[z] <- NA
      }
      
      
      BP_Norm <- which(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="REF"] > BP_prob_cutoff)
      max_prob_REF[z] <- max(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="REF"])
      
      if(BP_REF_num[z] > 0){
        max_U2_REF[z] <- max(branchpoint_mutated$U2_binding_energy[branchpoint_mutated$allele_status=="REF"][BP_Norm])
      }else{
        max_U2_REF[z] <- NA
      }
      
      BP_Mut <- which(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="ALT"] > BP_prob_cutoff)
      max_prob_ALT[z] <- max(branchpoint_mutated$branchpoint_prob[branchpoint_mutated$allele_status=="ALT"])
      
      if(BP_ALT_num[z] > 0){
        max_U2_ALT[z] <- max(branchpoint_mutated$U2_binding_energy[branchpoint_mutated$allele_status=="ALT"][BP_Mut] )
      }else{
        max_U2_ALT[z] <- NA
      }
    
      dists <- unique(branchpoint_mutated$distance[branchpoint_mutated$branchpoint_prob > BP_prob_cutoff])
      
      ref_p <- branchpoint_mutated[branchpoint_mutated$allele_status=="REF" & branchpoint_mutated$distance %in% dists,]$branchpoint_prob
      alt_p <- branchpoint_mutated[branchpoint_mutated$allele_status=="ALT" & branchpoint_mutated$distance %in% dists,]$branchpoint_prob
      
      deleted_n[z] = length(which((ref_p - alt_p) > 0.2))
      created_n[z] = length(which((ref_p - alt_p) < -0.2))
      
      if(z%%100==0){
        message(z)
      }
    }

  m <- match(snp_ids, snp_info$id)
  n <- match(snp_ids, predictions$id)

  snp_info_matches <- cbind(snp_info[m,],
                            BP_REF_num,BP_ALT_num,
                            created_n,deleted_n,
                            dist_to_exon,dist_to_BP_REF,dist_to_BP_ALT,
                            max_prob_REF,max_prob_ALT,max_U2_REF,max_U2_ALT)
  
  return(snp_info_matches)
}
