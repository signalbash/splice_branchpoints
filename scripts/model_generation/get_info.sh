#!/bin/sh
#run from project level folder

Rscript scripts/get_branchpoint_loc.R "$1" "$2"

Rscript scripts/make_neg_examples.R "$1" "$2"

Rscript scripts/make_branchpoint_bed.R "$1" "$2"

bedtools getfasta -fi data/inputs/genome.fa -bed data/outputs/branchpoint_df_501_"$1".bed -fo data/outputs/branchpoint_df_501_"$1".fa -name -s

Rscript scripts/get_canon_ppt_seq_info.R "$1" "$2"

