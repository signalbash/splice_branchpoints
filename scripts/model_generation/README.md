## Model Generation

**1. Reconfigure gencode annotations to table format**

```
scripts/model_generation/condenseGTF.sh
```
condenses gencode 12 annotation (GTF) into columns, readable by R
Needs to be run in the same directory as the GTF
output in inputs/altered_inputs

**2. Generate Testing/Training dataset**
Run by chromosome
```
while read c; do scripts/model_generation/get_info.sh $c $PWD; done < scripts/chroms.txt
```

`get_info.sh` calls the following scripts:

* `get_branchpoint_loc.R` makes: 
  * `gencode.v12.condensed.exons_chrZ.csv` Unique exons from the GENCODE12 annotation, renamed as transcriptid_exonnum
  * `duplicated_exons_chrZ.csv` Corresponding exon names when same exon used by multiple transcripts
  * `branchpoint_df_all_lc_chrZ.csv` All low confidence (LC) branchpoints with distances to GENOCDE12 exons
  * `branchpoint_df_all_hc_chrZ.csv` All high confidence (HC) branchpoints with distances to GENOCDE12 exons
  * `branchpoint_df_gene_lc_chrZ.csv` Low confidence (LC) branchpoints with distances to GENOCDE12 exons - 3' and 5' exons share same parent gene 
  * `branchpoint_df_gene_hc_chrZ.csv` High confidence (HC) branchpoints with distances to GENOCDE12 exons - 3' and 5' exons share same parent gene 

* `make_neg_examples_1-50.R` Generates negative examples for each exon with a known branchpoint
	Negatives are all positions 18-44nt from 3' exon which are not contained in HC or LC dataset

Makes `branchpoint_df_gene_N_chrZ.csv`

* `make_branchpoint_bed.R` Combines HC branchpoints and generated negatives to generate bed file covering 501 surrounding nucleotides for each point.

Makes `branchpoint_df_chrZ.csv` Combines file for HC and negatives and renames points.
	`branchpoint_df_501_chrZ.bed` BED file (+/- 250nt)

* `bedtools getfasta -fi data/inputs/genome.fa -bed data/outputs/branchpoint_df_501_"$1".bed -fo data/outputs/branchpoint_df_501_"$1".fa -name -s` To convert to fasta sequence

* `get_canon_ppt_seq_info.R` Generates additional features for each point
	Makes `branchpoint_df_with_seq_r2_chrZ.csv` 

**3. Model Optimisation**

**Preprocess data from all chromosomes**
```
Rscript scripts/model_generation/make_model_svm_prepro.R
```
Makes: `Preprocessed_data.RData` - All Rdata, `Preprocessed_data_small.RData` Smaller file with only R objects required for making models.

**Optimise training positive:negative ratios**
```
Rscript scripts/model_generation/make_model_svm_trainsizes.R 1234 4
```
Generates multiple svm models with training true positives=1000 and varying true negative vales (1x - 20x)
makes: `training_size_sel_Cval4_1234.RData`
```
Rscript scripts/model_generation/process_size_selection.R data/training_size_sel_Cval4_1234.RData 
```
Calculates F1 values for each model, plots against number of training true negatives size
Creates `best_ratio.txt` = 8

**Optimise C value for SVM**
```
for i in 0.5 1 2 3 4 8; do for j in 1 2 3; do Rscript scripts/model_generation/make_model_svm_small.R 5000 10000 $i $j;done;done;
```
```
process_cval_bp.R
```
Calculates F1 values for each model, plots against c value
Creates `best_c.txt` = 4

**Recursive feature selection for SVM**

```
for i in 1 2 3 4 5 6 7 8 9 10; do Rscript scripts/model_generation/make_model_svm_rfe_combined.R 2000 16000 $i; done;
```

`make_model_svm_rfe_combined.R`
Performs recursive feature elimination on subset of data.
Makes: `bp44_2000.16000-5000.1e+05_Cval2_repX_splicing_model_rfe.RData` which contains the rfe profile
Then systematically removes each ordered variable and evaluates performance of new model
Makes: `bp44_2000.16000-10000.2e+05_Cval2_repX_rfe16000_splicing_model_ALL_opt_elim_optdf.RData` &
`bp44_2000.16000-10000.2e+05_Cval2_repX_rfe16000_splicing_model_ALL_opt_elim.RData`

```
Rscript scripts/model_generation/process_rfe_bp.R
```
**4. Final model generation**
```
Rscript scripts/model_generation/Make_final_model.R
```

**5. Feature removal
```
#single variables
Rscript scripts/model_generation/Make_model_removeVars1.R
#grouped variables
Rscript scripts/model_generation/Make_model_removeVarsSet.R
```

**6. Naive Bayes Model
```
Rscript scripts/model_generation/Make_model_nb.R
```

