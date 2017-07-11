################################################################################
#             Predict branchpoints from exon annotation + fasta                #
################################################################################

options(scipen=999)
suppressMessages(library(argparse))
suppressMessages(library(branchpointer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))

parser <- ArgumentParser()
parser$add_argument("-s", "--split", type="character", help="index number for split")
parser$add_argument("-ss", "--splitsize", type="character", help="split size")
#from condense_annotation.sh x.gtf x.exons.txt
parser$add_argument("-f", "--fa", type="character", help="genome .fasta annotation file")
parser$add_argument("-g", "--gtf", type="character", help="genome gtf annotation file")
parser$add_argument("-b", "--bedtools", type="character", help="bedtools binary location")
parser$add_argument("-c", "--cores", type="character", help="cores")
parser$add_argument("-p", "--prefix", type="character", help="prefix for output files")

args <- parser$parse_args()

core_n <- as.numeric(args$cores)
if(core_n < 2){
  use_P=FALSE
}else{
  use_P=TRUE
}

exons <- gtfToExons(args$gtf)

gene_ids <- unique(exons@elementMetadata$gene_id)

query <- makeBranchpointWindowForExons(gene_ids, "gene_id", exons)

message(length(query))

splitn <- as.numeric(args$splitsize)
ind <-as.numeric(args$split)
index <- (splitn*(ind-1) +1):min((splitn*ind),length(query[,1]))
query=query[index]
chromnum <- paste0("_split_", ind)

branchpointPredictions <- predictBranchpoints(query, queryType = "region", useParallel = use_P, cores = core_n,
                                                    genome = args$fa, bedtoolsLocation = args$bedtools)

table <- (as.data.frame(branchpointPredictions))
table$seq <- with(mcols(branchpointPredictions), 
                  paste0(seq_neg5,seq_neg4,seq_neg3,seq_neg2,seq_neg1,
                         seq_pos0,seq_pos1,seq_pos2,seq_pos3,seq_pos4,seq_pos5))

table <- table[,c('seqnames','test_site','test_site','strand' ,'exon_id','exon_number','exon_5prime','to_5prime_point','to_3prime_point',
                  'seq','seq_pos0', 'branchpoint_prob',"U2_binding_energy")]

colnames(table) <- c("chromosome","start","end","strand",
                     "exon_id", "exon_number", "exon_5prime","to_5prime",
                     "to_3prime","seq_motif","branchpoint_nt",
                     "branchpoint_prob","U2_binding_energy")

file_prefix <- args$prefix

write.table(table, file=paste0(file_prefix,chromnum,"_predictions.txt"), sep="\t",
            row.names=FALSE, quote=FALSE)


