# splice_branchpoints

Scripts used in:
Signal, B., Gloss, B.S., Dinger, M.E., & Mercer, T.E. Machine-Learning annotation of human splicing branchpoints. 2016. 

## branchpointer
The machine learning model for predicting intronic branchpoints is available as an R package at github.com/betsig/branchpointer and github.com/betsig/branchpointer_dev (development version)

## Abstract
**Motivation:** The branchpoint element is required for the first lariat-forming reaction in splicing. However current catalogues of human branchpoints remain incomplete due to the difficulty in experimentally identifying these splicing elements. To address this limitation, we have developed a machine-learning algorithm - branchpointer - to identify branchpoint elements solely from gene annotations.

**Results:** Using branchpointer, we annotate branchpoint elements in 85% of human gene introns with best-in-class accuracy (96.26%). In addition to annotation, branchpointer can evaluate the impact of SNPs on branchpoint architecture to inform functional interpretation of genetic variants. Branchpointer identifies all published deleterious branchpoint mutations annotated in clinical variant databases, and finds thousands of additional clinical and common genetic variants with similar predicted effects. This genome-wide annotation of branchpoints provides a reference for the genetic analysis of splicing, and the interpretation of noncoding variation. 

## Model Generation (model_generation)
Scripts for training the branchpoint machine learning model

## Prediction of branchpoints in genome annotations (genome_predictions)
Scripts for predicting branchpoints in Human (Gencodev24, Gencodev19)

## Analysis and Figures (analysis)
Scripts for analysing model performance, genomic attributes of branchpoints, and impact of variation on branchpoints.
