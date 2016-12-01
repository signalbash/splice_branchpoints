# splice_branchpoints

Scripts used in:
Signal, B., Gloss, B.S., Dinger, M.E., & Mercer, T.E. Machine-Learning annotation of human splicing branchpoints. 2016. 

##branchpointer
The machine learning model for predicting intronic branchpoints is available as an R package at githib.com/betsig/branchpointer


##Abstract
The branchpoint element is required for the first lariat-forming reaction in splicing. We have developed a machine-learning algorithm trained with empirical human branchpoint annotations to identify branchpoint elements from primary genome sequence alone. Using this approach, we can accurately locate branchpoints elements in 85% of introns in current gene annotations. Consistent with branchpoints as basal genetic elements, we find our annotation is unbiased towards gene type and expression levels. A fraction of introns was found to encode multiple branchpoints raising the prospect that mutational redundancy is encoded in key genes. We also identify cases of deleterious branchpoint mutations in clinical variant databases that may explain disease pathogenicity. We propose the broad annotation of branchpoints constitutes a valuable resource for interpreting the impact of common- and disease-causing human genetic variation on gene splicing.


##Model Generation (model_generation)
Scripts for training the branchpoint machine learning model

##Prediction of branchpoints in genome annotations (genome_predictions)
Scripts for predicting branchpoints in Human (Gencodev24, Gencodev12, Genocdv19), Mouse (GencodevM10), Zebrafish, Xenopus, Drosophila, and Chicken.

##Analysis and Figures (analysis)
Scripts for analysing model performance, genomic attributes of branchpoints, and impact of variation on branchpoints.



