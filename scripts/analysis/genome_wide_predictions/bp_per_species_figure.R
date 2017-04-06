options(scipen = 999)
options(stringsAsFactors = F)
library(data.table)
library(ggplot2)

files <- list.files("data/genome_predictions/", pattern="predictions.txt")
species <- c("gencode.vM10", "gencode.v24","Danio_rerio","Xenopus_tropicalis",
             "Drosophila_melanogaster","Gallus_gallus")
keep <- which(gsub("_predictions.txt","",files) %in% species)

files <- files[keep]

for(f in seq_along(files)){
  preds <- fread(paste0("data/genome_predictions/",files[f]), data.table=FALSE)
  
  branchpoints_by_intron_summary <- aggregate(branchpoint_prob ~ exon_3prime, preds, function(x) length(which(x >=0.5)))
  branchpoints_by_intron_summary$species <- gsub("_predictions.txt","",files[f])
  
  if(exists("all_bp_summary")){
    all_bp_summary <- rbind(all_bp_summary, branchpoints_by_intron_summary)
  }else{
    all_bp_summary <- branchpoints_by_intron_summary
  }
}

all_bp_summary$number_bp_factor <- all_bp_summary$branchpoint_prob
all_bp_summary$number_bp_factor[all_bp_summary$number_bp_factor > 1] <- "2+"

all_bp_table <- as.data.frame(table(all_bp_summary$number_bp_factor, all_bp_summary$species))
sums <- aggregate(Freq ~ Var2, all_bp_table, sum)
all_bp_table$Freq <- all_bp_table$Freq/sums$Freq[match(all_bp_table$Var2, sums$Var2)]

colnames(all_bp_summary)[4] <- "BP_per_intron"
colnames(all_bp_table)[1] <- "BP_per_intron"

pdf("Figures/Figure6.pdf", 6.69,6, useDingbats = FALSE)
ggplot(all_bp_summary, aes(x=species, fill=BP_per_intron)) + geom_bar()+ theme_bw()+
  theme(text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1),
        legend.key.size=unit(0.2, "inches")) + scale_y_continuous(name="Introns") + scale_x_discrete(name="Species")


ggplot(all_bp_table, aes(x=Var2, fill=BP_per_intron, y=Freq)) + geom_bar(stat="identity")+ theme_bw()+
  theme(text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1),
        legend.key.size=unit(0.2, "inches")) + scale_y_continuous(name="Introns") + scale_x_discrete(name="Species")
dev.off()
