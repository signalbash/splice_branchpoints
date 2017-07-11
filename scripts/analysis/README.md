**1. Evaluate the performance ofthe branchpointer model**

`evaluate_performance.R`

**2. Get conservation for branchpoints and surrounding regions**
`get_cons_locs.R`
calls - `get_bw_entries.py`

**3. Annotate the Gencode v26 predictions with attributes for figures/stats**
`annotate_gencode_v26.R`

**4. Combine GTEX QTL calls into single files**
`combine_GTEx_QTL.R`

**5. Run branchpointer on GTEx variants and do some filtering**
`common_variants.R`

**6. Run branchpointer on ClinVar variants and do some filtering**
`disease_variants.R`

**7. Plot figures for manuscript**
`plot_figures.R`
 - calls `quick_wilcox.R`