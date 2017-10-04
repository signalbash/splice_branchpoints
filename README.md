# splice_branchpoints

Scripts used in:
Signal, B., Gloss, B.S., Dinger, M.E., & Mercer, T.R. Machine-Learning annotation of human splicing branchpoints. 2016. 
https://www.biorxiv.org/content/early/2016/12/14/094003

##branchpointer
The machine learning model for predicting intronic branchpoints is available as an R package at githib.com/betsig/branchpointer

##Abstract
Motivation: The branchpoint element is required for the first lariat-forming reaction in splicing. How-ever current catalogues of human branchpoints remain incomplete due to the difficulty in experimen-tally identifying these splicing elements. To address this limitation, we have developed a machine-learning algorithm - branchpointer - to identify branchpoint elements solely from gene annotations and genomic sequence.
Results: Using branchpointer, we annotate branchpoint elements in 85% of human gene introns with sensitivity (61.8%) and specificity (97.8%). In addition to annotation, branchpointer can evaluate the impact of SNPs on branchpoint architecture to inform functional interpretation of genetic variants. Branchpointer identifies all published deleterious branchpoint mutations annotated in clinical variant databases, and finds thousands of additional clinical and common genetic variants with similar pre-dicted effects. This genome-wide annotation of branchpoints provides a reference for the genetic analysis of splicing, and the interpretation of noncoding variation.
Availability: Branchpointer is written and implemented in the statistical programming language R and is freely available under a BSD license as a package through Bioconductor.


##Model Generation (model_generation)
Scripts for training the branchpoint machine learning model

##Prediction of branchpoints in genome annotations (genome_predictions)
Scripts for predicting branchpoints in Human (Gencodev24, Gencodev12, Gencodev19), Mouse (GencodevM10), Zebrafish, Xenopus, Drosophila, and Chicken.

##Analysis and Figures (analysis)
Scripts for analysing model performance, genomic attributes of branchpoints, and impact of variation on branchpoints.



