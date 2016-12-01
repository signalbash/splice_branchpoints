################################################################################
#                      Process Genome-wide predictions                         #
################################################################################

options(scipen = 999)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(plyr)

#process files in parts
species <- c("gencode.vM10", "Danio_rerio","Xenopus_tropicalis",
	"Drosophila_melanogaster","Gallus_gallus")

pattern <- paste0(species, "_split")


for(sp in seq_along(pattern)){

	if(exists("p")){rm(p)}

	files <- list.files("genome_predictions/", pattern = pattern[sp])
	  
	files_pred <- files[grep("pred", files)]
	files_attr <- files[grep("attr", files)]
	  
	for(f in seq_along(files_attr)){
	  p_part <- as.data.frame(fread(paste0("genome_predictions/", files_pred[f])))

	  if(exists("p")){
	    p <- rbind(p, p_part)
	  }else{
	    p <- p_part
	  }
	}

	#unique introns
	new_id <- with(p, paste0(chromosome,strand,end))

	message(species[sp])

	p <- p[!duplicated(new_id),]
	message("introns tested:")
	message(length(unique(p$id)))
	P_BP <- p[p$branchpoint_prob >= 0.5,]

	introns <- unique(P_BP$id)

	message("introns with bp:")
	message(length(introns))
	message("bps:")
	message(length(P_BP$id))
}
