################################################################################
#                      GTEx eQTL and sQTL file processing                      #
################################################################################

options(stringsAsFactors = F)
library(data.table)

#downloaded from the GTEx website
file=list.files("data/GTEX/GTEx_Analysis_V6_eQTLs/")

for(f in file){
    eQTLs <- as.data.frame(fread(paste0("data/GTEX/GTEx_Analysis_V6_eQTLs/",f)))
    w <- which(eQTLs$snp_pos > eQTLs$gene_start & eQTLs$snp_pos < eQTLs$gene_stop)
    
    keep_eQTLs <- eQTLs[w, c(1,2,6,22,23)]
    keep_eQTLs$tissue <- gsub("_Analysis.snpgenes","", f)
    if(exists("eQTLs_all")){
        eQTLs_all <- rbind(eQTLs_all, keep_eQTLs)
    }else{
        eQTLs_all <- keep_eQTLs
    }
   
    message(f) 
}

write.csv(eQTLs_all, "data/GTEX/eQTL_ALL.csv")

#downloaded from the GTEx website
file=list.files("data/GTEX/sQTLs-sQTLseeker-merged/", full.names=TRUE)
sQTLs.AdiposeTissue <- read.delim(file[1])
sQTLs.Blood <- read.delim(file[2])
sQTLs.BloodVessel <- read.delim(file[3])
sQTLs.Heart <- read.delim(file[4])
sQTLs.Lung <- read.delim(file[5])
sQTLs.Muscle <- read.delim(file[6])
sQTLs.Nerve <- read.delim(file[7])
sQTLs.Skin <- read.delim(file[8])
sQTLs.Thyroid <- read.delim(file[9])

sQTLs.AdiposeTissue$Tissue <- "AdiposeTissue"
sQTLs.Blood$Tissue <- "Blood"
sQTLs.BloodVessel$Tissue <- "BloodVessel"
sQTLs.Heart$Tissue <- "Heart"
sQTLs.Lung$Tissue <- "Lung"
sQTLs.Muscle$Tissue <- "Muscle"
sQTLs.Nerve$Tissue <- "Nerve"
sQTLs.Skin$Tissue <- "Skin"
sQTLs.Thyroid$Tissue <- "Thyroid"

sQTLs.ALL=rbind(sQTLs.AdiposeTissue, sQTLs.Blood,sQTLs.BloodVessel,
                sQTLs.Heart,sQTLs.Lung,sQTLs.Muscle,sQTLs.Nerve,
                sQTLs.Skin,sQTLs.Thyroid)

write.csv(sQTLs.ALL, "data/GTEX/sQTL_ALL.csv")




