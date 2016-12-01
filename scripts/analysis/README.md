### model_performance

`evaluate_performance.R`
Evaluate performance metrics of branchpointer on the testing dataset
Compares performance with SVM-BPfinder and HSF

*Figures* 
Figure 1D-F. Development and erformance of the branchpoint detection model. 
Figure S1. branchpointer classification performance metrics for probability cutoffs 0.01-0.99. 
Figure S2. Precision recall curve for branchpointer.
Figure S3. Frequency of each two nucleotide branchpoint motif at each probability score above 0.5.  

### genome_wide_predictions

```format_cdsutr_annotations.sh data/genome_predictions/gencode.v24.annotation.gtf data/genome_predictions/gencode.v24.condensed.cdsutr.txt```
formats GTF as a table with nly CDS and UTR annotations

`annotate_gencode_part1.R`
Get branchpoint and intronic features for gencodev24 annoations
Saves Rdata as `"data/Figure_files.Rdata"`

`get_cons_loc.R`
Get phyloP conservation scores for branchpoint windows
Calls `get_bw_entries.py` to extract intronic conservation scores from a bigwig file

`annotate_gencode_part2.R`
Get branchpoint and intronic features for gencodev24 annotations
Saves Rdata as `"data/Figure_files.Rdata"`

### variation

`combine_GTEx_QTL.R`
Combines the GTEx eQTL and sQTL results files into single tables

`predictions_to_stats.R`
Function for converting branchpointer SNP predictions to site stats

`common_variants.R`
Formats GTEx variants.
Predicts branchpoints in reference and alternative sequences using branchpointer.
Saves RData as `"data/commonVariants.Rdata"`

`disease_variants.R`
Formats ClinVar variants.
Predicts branchpoints in reference and alternative sequences using branchpointer.
Saves RData as `"data/diseaseVariants.Rdata"`

### Alternative splicing annotations 

`alternative_splicing_annotate.R`
Generates an annotation of exon triplets with evidence for both canonical splicing and exon skipping (skipped) and canonical use only (control)
Uses Gencodev19 annotations due to splice jjunction reads being aligned to this transcriptome.

`gencodev19_alternative_splicing.R`


`plot_figures.R`
Makes figures and computes statistics for genome wide predictions.

*Figures*
Figure 2. Prediction of splicing branchpoints in GENCODE introns.
Figure 3. Features of introns with branchpoint multiplicity.
Figure 4. Effect of the SNP rs2269219 on branchpoints in Fech.

Figure S4. Introns with annotated branchpoints from the Mercer annotation (known), and the branchpoint detection model (predicted) for the gene biotypes long noncoding RNA, protein coding, pseudogenes, and all other biotypes. The default cut-off score of 0.5 was used for predictions. 
Figure S5. Size of introns with single or multiple annotated branchpoints from the Mercer branchpoint annotation. 
Figure S6. Size of introns with single or multiple annotated branchpoints in parent gene biotypes. Only biotypes with at least 1000 branchpoint annotated introns are shown. 
Figure S7. Splicing element strength at constitutively spliced and skipped exon triplets. 
Figure S8. Locations of all intronic (1-50nt from the 3’ exon) ClinVar and GTEx SNPs.