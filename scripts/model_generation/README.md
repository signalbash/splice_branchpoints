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
Rscript scripts/model_generation/preProcessData.R
```
Makes: `Preprocessed_data.RData` - All Rdata, `Preprocessed_data_small.RData` Smaller file with only R objects required for making models.
Also contain indexes for testing and train subsets

**Optimise training positive:negative ratios**
```
Rscript scripts/model_generation/trainingRatioTesting.R 2;
```
Generates multiple svm models with training true positives=500 and varying true negative vales (1x - 20x)
makes: 
`gbm_training_size_sel_2.RData`

**Optimise parameters for gbm training**

First grid search
```
for i in 1 2 3 4 5 6; 
do Rscript scripts/model_generation/trainingGridSearch.R $i;
done;
```
Makes:
`gbm_gridSearch_1.RData`
`gbm_gridSearch_2.RData`
`gbm_gridSearch_3.RData`
`gbm_gridSearch_4.RData`
`gbm_gridSearch_5.RData`
`gbm_gridSearch_6.RData`

Second grid search
Using results from first grid search to expand parameter values
```
for i in 1 2 3 4 5 6; 
do Rscript scripts/model_generation/trainingGridSearch2.R $i;
done;
```
Makes:
`gbm_gridSearch2_1.RData`
`gbm_gridSearch2_2.RData`
`gbm_gridSearch2_3.RData`
`gbm_gridSearch2_4.RData`
`gbm_gridSearch2_5.RData`
`gbm_gridSearch2_6.RData`

**4. Final model generation**
```
Rscript scripts/model_generation/trainFinalModel.R 42
```
Makes:
`gbm_final_model_42.RData`

**5. Feature removal
Removal of all variables used in the final model one at a time, and in related sets.
Note that Sequence identity features at the same site (i.e. seq_pos0A/seq_pos0C/seq_pos0G/seq_pos0T) must be removed in groups, as it is possible to still infer sequence identity when only one is removed.

```
#single variables
for i in `seq 0 53`;
do Rscript scripts/model_generation/removeVariable.R $i;
done;

#grouped variables
for i in `seq 0 15`;
do Rscript scripts/model_generation/removeVariableSets.R $i;
done;
```
Makes:
`gbm_models_removeVar_i.RData` 
(i = 0:53)
and
`gbm_models_removeVarSets_i.RData`
(i = 0:15)

**6. Naive Bayes Model
Generates a naive bayes model with the same training dataset as the final model
```
Rscript scripts/model_generation/trainNBModel.R 42
```
Makes:
`nb_final_model_42.RData`
and
`nb_performance.RData`

