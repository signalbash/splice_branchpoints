################################################################################
#             Predict branchpoints from exon annotation + fasta                #
################################################################################

options(scipen=999)
suppressMessages(library(argparse))
suppressMessages(library(branchpointer))
suppressMessages(library(stringr))

parser <- ArgumentParser()
parser$add_argument("-s", "--split", type="character", help="index number for split")
parser$add_argument("-ss", "--splitsize", type="character", help="split size")
parser$add_argument("-fa", "--genomefasta", type="character", help="genome .fa file")
#from condense_annotation.sh x.gtf x.exons.txt
parser$add_argument("-e", "--genomeexons", type="character", help="genome exon annotation file")
parser$add_argument("-c", "--cores", type="character", help="cores")
parser$add_argument("-b", "--bedtools", type="character", help="location of the bedtools binary")
parser$add_argument("-p", "--prefix", type="character", help="prefix for output files")

args <- parser$parse_args()

core_n <- as.numeric(args$cores)
if(core_n < 2){
  use_P=FALSE
}else{
  use_P=TRUE
}

exons <- readExonAnnotation(args$genomeexons)

#by definition first exons shouldn't have branchpoints
keep <- which(exons$exon_number > 1)
exons_subset <-  exons[keep,]

#make windows
window_starts <- (exons_subset$start - 50)
window_ends <- (exons_subset$start - 10)
  
window_starts[exons_subset$strand=="-"] <- (exons_subset$end + 10)[exons_subset$strand=="-"]
window_ends[exons_subset$strand=="-"] <- (exons_subset$end + 50)[exons_subset$strand=="-"]
  
window_df <- data.frame(
  id = exons_subset$exon_id,
  chromosome = exons_subset$chromosome,
  chrom_start = window_starts,
  chrom_end = window_ends,
  strand = exons_subset$strand
)

query <- window_df[!duplicated(with(window_df, paste0(chromosome,chrom_start,strand))),]

message(length(query[,1]))

splitn <- as.numeric(args$splitsize)
ind <-as.numeric(args$split)
index <- (splitn*(ind-1) +1):min((splitn*ind),length(query[,1]))
query=query[index,]
chromnum <- paste0("_split_", ind)


query <- getQueryLoc(query,query_type="region",exons = exons, use_parallel = use_P, cores=core_n)

message("gotloc")

query_attributes <- getBranchpointSequence(query,
                                        query_type = "region",
                                        genome = args$genomefasta,
                                        bedtools_location=args$bedtools,
                                        unique_id = chromnum,
                                        use_parallel = use_P, 
                                        cores=core_n)

message("gotatt")

file_prefix <- args$prefix

write.table(query_attributes, file=paste0(file_prefix,chromnum,"_attributes.txt"), sep="\t",
            row.names=FALSE, quote=FALSE)

branchpoint_predictions <- predictBranchpoints(query_attributes)

write.table(branchpoint_predictions, file=paste0(file_prefix,chromnum,"_predictions.txt"), sep="\t",
            row.names=FALSE, quote=FALSE)

