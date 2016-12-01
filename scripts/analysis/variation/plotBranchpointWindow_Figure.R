transcript_overlap = function(query_name,branchpoint_predictions,
                                 exons,
                                 tid1=NA,tid2=NA) {
  #find the query
  ind = which(!is.na(match(
    branchpoint_predictions$id, query_name
  )))
  
  branchpoint_mutated = branchpoint_predictions[ind,]

  branchpoint_mutated$distance = branchpoint_mutated$distance * -1
  
  #plot gene/transcript structure

    x = match(exons$gene_id, exons$gene_id[grep(branchpoint_mutated$exon_3prime[1], exons$exon_id)][1])
    transcripts = exons[which(!is.na(x)),]

    transcripts_gene = transcripts
    transcripts_gene$transcript_id = transcripts_gene$gene_id
    transcripts = rbind(transcripts,transcripts_gene)
    transcripts$length = transcripts$end - transcripts$start

    if (transcripts$strand[1] == "+") {
      bp_exons = which(transcripts$start ==
                         transcripts$start[which(transcripts$exon_id ==
                                                   branchpoint_mutated$exon_3prime[1])][1])

    }else{
      bp_exons = which(transcripts$end ==
                         transcripts$end[which(transcripts$exon_id ==
                                                 branchpoint_mutated$exon_3prime[1])][1])
    }

    transcripts$SNP_exon = "isoform exon"
    transcripts$SNP_exon[bp_exons] = "isoform exon\nassociated with BP"
    transcripts$SNP_exon[transcripts$transcript_id == transcripts$gene_id[1]] = "gene"

    #make smaller data.frame for plotting connecting lines
    transcripts_lines = transcripts[!duplicated(transcripts$transcript_id),]

    maxs = vector()
    mins = vector()
    maxe = vector()
    mine = vector()
    for (t in seq_along(transcripts_lines$chromosome)) {
      ind = grep(transcripts_lines$transcript_id[t], transcripts$transcript_id)
      maxs[t] = max(transcripts$start[ind])
      mine[t] = min(transcripts$end[ind])
      maxe[t] = max(transcripts$end[ind])
      mins[t] = min(transcripts$start[ind])

    }

    #only plot transcripts overlapping the query
    keep = which(str_sub(transcripts_lines$transcript_id,1,15) %in% str_sub(c(tid1,tid2),1,15))
  
    if(length(keep) >0 ){
    transcripts_lines$max = maxs
    transcripts_lines$min = mine
    transcripts_lines = transcripts_lines[keep,]
    BP_transcripts = unique(transcripts_lines$transcript_id)
    transcripts = transcripts[transcripts$transcript_id %in% BP_transcripts,]
    transcripts$transcript_id_num = as.numeric(as.factor(transcripts$transcript_id)) * -1
    transcripts_lines$transcript_id_num = as.numeric(as.factor(transcripts_lines$transcript_id)) * -1

    start_value <- transcripts[transcripts$SNP_exon=="isoform exon\nassociated with BP",]$start
    end_value <- transcripts[transcripts$SNP_exon=="isoform exon\nassociated with BP",]$end
    
    transcripts_lines$transcript_id <-grep("ENS", unlist(str_split(transcripts_lines$transcript_id,"[.]")), value=TRUE)

    transcripts$SNP_exon <- factor(transcripts$SNP_exon, levels =c("gene","isoform exon","isoform exon\nassociated with BP" ))

    structure_plot = ggplot(transcripts, aes(xmin = start, xmax = start + length,
                                             ymin = transcript_id_num - 0.4,
                                             ymax = transcript_id_num + 0.4, fill = factor(SNP_exon))) +
      #geom_rect(xmin = start_value - width*2,xmax = end_value + width*2,ymin=min(transcripts$transcript_id_num)-0.5,
      #          ymax=max(transcripts$transcript_id_num)+0.5,col = NA, fill = "grey90") +
      geom_segment(data = transcripts_lines, aes(x = min,xend = max, y = transcript_id_num,
                                                 yend = transcript_id_num)) +
      geom_rect() +
      scale_fill_manual(values = c("black","grey60",nt_cols[1]), name="annotation") +
      scale_y_continuous(breaks = seq(-1, min(transcripts_lines$transcript_id_num),-1),
                         labels = transcripts_lines$transcript_id[
                           match(seq(-1, min(transcripts_lines$transcript_id_num),-1),
                                 transcripts_lines$transcript_id_num)]) +
      scale_x_continuous(name=transcripts_lines$chromosome[1]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = "white"), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),text=element_text(size=8),legend.key.size=unit(0.15, "inches"))

    structure_plot
    }
}


plotBranchpointWindow_Figure=function(query_name,
         branchpoint_predictions,
         query_attributes,
         exons,
         probability_cutoff = 0.5,structure = "plot_overlap",
         tid1=NA,tid2=NA) {
  #find the query
  ind = which(!is.na(match(
    branchpoint_predictions$id, query_name
  )))
  
  branchpoint_mutated = branchpoint_predictions[ind,]
  
  branchpoint_mutated$distance = branchpoint_mutated$distance * -1
  branchpoint_mutated$allele_status=factor(branchpoint_mutated$allele_status, levels=c("REF","ALT"))
  branchpoint_mutated$U2_binding_energy[branchpoint_mutated$branchpoint_prob < probability_cutoff] <- 0
  
  #make plots for reference sequence
  #U2 binding energy by sequence position
  plot_U2 = ggplot(branchpoint_mutated,
                       aes(x = distance, y = U2_binding_energy, fill=allele_status)) +
    geom_bar(stat = "identity",width = 1, position="dodge") +
    scale_x_continuous(limits = c(-27,-20), labels = seq(26,20,-2),
                       breaks = seq(-26,-20, 2), name ="Distance to 3' exon (nt)") +
    theme_bw() +
    scale_fill_manual(values=nt_cols[c(2,4)]) +
    theme(legend.position = "none",text=element_text(size=10)) +
    scale_y_continuous(name = "U2 binding energy",limits = c(0,5),
                       breaks = c(1,3,5), labels = c("1.0","3.0","5.0"))
  
  
  #Sequence identity
  plot_seq_ref = ggplot(branchpoint_mutated_ref, aes(x = distance, y = 1, col = seq_pos0,label = seq_pos0)) +
    geom_text(size = 6, family = "Courier") +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = levels(branchpoint_mutated_ref$seq_pos0)) +
    theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
          axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank(),text=element_text(size=8)) +
    scale_x_continuous(limits = c(-45,-17))
  
  #Sequence identity
  plot_seq_alt = ggplot(branchpoint_mutated_alt, aes(x = distance, y = 1, col = seq_pos0,label = seq_pos0)) +
    geom_text(size = 6, family = "Courier") +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = levels(branchpoint_mutated_ref$seq_pos0)) +
    theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
          axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank(),text=element_text(size=8)) +
    scale_x_continuous(limits = c(-45,-17))
  
  ######Get reference vs. alt sequence for 0-50 window#####
  
  seq = query_attributes[which(grepl(query_name, query_attributes$id) &
                                 grepl("REF",query_attributes$id))[1],]$seq
  pos = query_attributes[which(grepl(query_name, query_attributes$id) &
                                 grepl("REF",query_attributes$id))[1],]$to_3prime
  seq_1 = unlist(str_split(str_sub(seq,251 + (pos - 50),250 + pos),""))
  
  seq = query_attributes[which(grepl(query_name, query_attributes$id) &
                                 grepl("ALT",query_attributes$id))[1],]$seq
  pos = query_attributes[which(grepl(query_name, query_attributes$id) &
                                 grepl("ALT",query_attributes$id))[1],]$to_3prime
  seq_2 = unlist(str_split(str_sub(seq,251 + (pos - 50),250 + pos),""))
  
  whole_seq = data.frame(nt = c(seq_1, seq_2),
                         pos = rep(seq(-50,-1,1),2),
                         set = rep(c("REF","ALT"), each =50))
  mut_pos = whole_seq$pos[which(whole_seq$nt[whole_seq$set == "REF"] !=
                                  whole_seq$nt[whole_seq$set == "ALT"])]
  
  plot_seq_comparison = ggplot(whole_seq, aes(x = pos, y = set, col = nt,label = nt)) +
    geom_rect(xmin = mut_pos - 0.5,xmax = mut_pos + 0.5, ymin = 0.5, ymax = 0.5 + 2,
              col = NA, fill = "grey90") +
    geom_text(family = "Courier") +
    scale_color_manual(values =nt_cols) +
    theme_bw() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    scale_x_continuous(limits = c(-50,-1)) 
  
  #plot gene/transcript structure
  x = match(exons$gene_id, exons$gene_id[grep(branchpoint_mutated$exon_3prime[1], exons$exon_id)][1])
  transcripts = exons[which(!is.na(x)),]
  
  transcripts_gene = transcripts
  transcripts_gene$transcript_id = transcripts_gene$gene_id
  transcripts = rbind(transcripts,transcripts_gene)
  transcripts$length = transcripts$end - transcripts$start
  
  if (transcripts$strand[1] == "+") {
    bp_exons = which(transcripts$start ==
                       transcripts$start[which(transcripts$exon_id ==
                                                 branchpoint_mutated$exon_3prime[1])][1])
    
  }else{
    bp_exons = which(transcripts$end ==
                       transcripts$end[which(transcripts$exon_id ==
                                               branchpoint_mutated$exon_3prime[1])][1])
  }
  
  transcripts$SNP_exon = "isoform exon"
  transcripts$SNP_exon[bp_exons] = "isoform exon\nassociated with BP"
  transcripts$SNP_exon[transcripts$transcript_id == transcripts$gene_id[1]] = "gene"
  
  #make smaller data.frame for plotting connecting lines
  transcripts_lines = transcripts[!duplicated(transcripts$transcript_id),]
  
  maxs = vector()
  mins = vector()
  maxe = vector()
  mine = vector()
  for (t in seq_along(transcripts_lines$chromosome)) {
    ind = grep(transcripts_lines$transcript_id[t], transcripts$transcript_id)
    maxs[t] = max(transcripts$start[ind])
    mine[t] = min(transcripts$end[ind])
    maxe[t] = max(transcripts$end[ind])
    mins[t] = min(transcripts$start[ind])
    
  }
  
  #only plot transcripts overlapping the query
  if(structure == "plot_overlap"){
    keep = which((maxe >= max(transcripts$end[transcripts$SNP_exon == "isoform exon\nassociated with BP"])) &
                   (mins <= min(transcripts$start[transcripts$SNP_exon == "isoform exon\nassociated with BP"])))
  }else if(structure == "plot_tids"){
    keep = which(str_sub(transcripts_lines$transcript_id,1,15) %in% str_sub(c(tid1,tid2),1,15))
  }else{
    keep= 1:length(transcripts_lines$chromosome)
  }
  transcripts_lines$max = maxs
  transcripts_lines$min = mine
  transcripts_lines = transcripts_lines[keep,]
  BP_transcripts = unique(transcripts_lines$transcript_id)
  transcripts = transcripts[transcripts$transcript_id %in% BP_transcripts,]
  transcripts$transcript_id_num = as.numeric(as.factor(transcripts$transcript_id)) * -1
  transcripts_lines$transcript_id_num = as.numeric(as.factor(transcripts_lines$transcript_id)) * -1
  
  start_value <- transcripts[transcripts$SNP_exon=="isoform exon\nassociated with BP",]$start
  end_value <- transcripts[transcripts$SNP_exon=="isoform exon\nassociated with BP",]$end
  
  #width <- (max(c(transcripts$start, transcripts$end)) - min(c(transcripts$start, transcripts$end))) / 100
  
  
  #     if (transcripts$strand[1] == "+") {
  #       transcripts_lines$transcript_id = paste0(transcripts_lines$transcript_id, " (+)")
  #     }else{
  #       transcripts_lines$transcript_id = paste0(transcripts_lines$transcript_id, " (-)")
  #     }
  transcripts_lines$transcript_id <-grep("ENS", unlist(str_split(transcripts_lines$transcript_id,"[.]")), value=TRUE)
  #transcripts_lines$transcript_id <- with(transcripts_lines, paste(transcript_id, transcript_type, sep="\n"))
  
  transcripts$SNP_exon <- factor(transcripts$SNP_exon, levels =c("gene","isoform exon","isoform exon\nassociated with BP" ))
  
  structure_plot = ggplot(transcripts, aes(xmin = start, xmax = start + length,
                                           ymin = transcript_id_num - 0.4,
                                           ymax = transcript_id_num + 0.4, fill = factor(SNP_exon))) +
    #geom_rect(xmin = start_value - width*2,xmax = end_value + width*2,ymin=min(transcripts$transcript_id_num)-0.5,
    #          ymax=max(transcripts$transcript_id_num)+0.5,col = NA, fill = "grey90") +
    geom_segment(data = transcripts_lines, aes(x = min,xend = max, y = transcript_id_num,
                                               yend = transcript_id_num)) +
    geom_rect() +
    scale_fill_manual(values = c("black","grey60",nt_cols[1]), name="annotation") +
    theme_bw() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())  
  #structure_plot
  
  
  ggdraw() +
    draw_plot(structure_plot,0,0.7,1,0.3) + draw_plot(plot_seq_comparison,0,0.6,1,0.1) +
    draw_plot(plot_U2,0,0.3,0.5,.3) +
    draw_plot_label(c("A","B","C"), c(0,0,0), c(1,0.7,0.6), size=18)
  
}
