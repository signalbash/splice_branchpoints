################################################################################
#                   Convert ensembl gtf to exon annotation                     #
################################################################################

options(scipen=999)
suppressMessages(library(argparse))
library(data.table)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("-a", "--annot", type="character", help="gtf file")
args <- parser$parse_args()

file <- args$annot

gtf <- fread(file)

gtf <- as.data.frame(gtf)

get_attribute <-  function(gtf, attribute){
  values <- unlist(str_split(gtf$V9, ";"))
  
  values <- grep(attribute, values, value = TRUE)
  
  values <- gsub(paste0(attribute," "), "", values)
  
  values <- gsub(" ", "", values)
  
  values <- gsub('"', "", values)
  
  if(length(values) == length(gtf[,9])){
    return(values)
  }else{
    return(NA)
  }
  
}

gtf <- gtf[gtf$V3 == "exon",]

gtf$gene_id <- get_attribute(gtf, "gene_id")
gtf$gene_biotype <- get_attribute(gtf, "gene_biotype")
gtf$transcript_id <- get_attribute(gtf, "transcript_id")
gtf$transcript_biotype <- get_attribute(gtf, "transcript_biotype")
gtf$exon_number <- get_attribute(gtf, "exon_number")
gtf$exon_id <- get_attribute(gtf, "exon_id")

write.table(gtf[,c(1,3,4,5,7,10:15)], file=gsub(".gtf", ".exons.txt", file), 
            quote=FALSE, col.names =FALSE, row.names = FALSE, sep ="\t")
