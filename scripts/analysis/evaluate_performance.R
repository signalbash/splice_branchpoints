################################################################################
#     Evaluate performance of the branchpointer model on the testing data      #
################################################################################

options(stringsAsFactors = F)

library(stringr)
library(data.table)
library(caret)
library(ggplot2)
library(cowplot)
library(reshape2)
library(PRROC)
library(branchpointer)
library(GenomicRanges)

# redo fastas?
rerunFasta <- FALSE


###### load in testing data and model ######

load("data/gbm_final_model_42.RData")
load("data/nb_model_only.RData")

n = predict(model, testing, "prob")
testing$branchpoint_prob=n$HC
n = predict(model, testingOptim, "prob")
testingOptim$branchpoint_prob=n$HC

n = predict(model.nb, testing, "prob")
testing$branchpoint_prob_NB=n$HC
n = predict(model.nb, testingOptim, "prob")
testingOptim$branchpoint_prob_NB=n$HC

#write training / testing indexes to file
training_ids <- branchpoint_df_HCN$new_ID[c(trainIndex.pos, trainIndex.neg)]
testing_ids <- branchpoint_df_HCN$new_ID[c(testIndex.pos, testIndex.neg)]
test_train_sites <- data.frame(id=c(training_ids, testing_ids),
                               testtrain=c(rep("train", length(training_ids)), rep("test", length(testing_ids))))
test_train_sites$chromosome <- unlist(lapply(str_split(test_train_sites$id, "_"), "[[",1))
test_train_sites$strand <- "+"
test_train_sites$strand[grep("-", test_train_sites$id)] <- "-"
test_train_sites$site <- unlist(lapply(str_split(test_train_sites$id, "_"), "[[",2))
test_train_sites$site <- gsub("[+]","_", gsub("[-]","_", test_train_sites$site))
test_train_sites$site <- unlist(lapply(str_split(test_train_sites$site, "_"), "[[",2))
test_train_sites$set <- unlist(lapply(str_split(test_train_sites$id, "_"), "[[",3))

write.csv(test_train_sites,"data/test_train_sites.csv", row.names = F, quote=F)

###### Find optimal probability cutoff ######

cutoff_performance <- data.frame(vals=seq(from=0.01, to=0.99, by=0.01),
                                 Accuracy=NA,
                                 Balanced_Accuracy=NA,
                                 Sensitivity=NA,
                                 Specificity=NA,
                                 PPV=NA,
                                 NPV=NA,
                                 F1=NA)

for(v in seq(along=cutoff_performance$vals)){
  testingOptim$Pred_Class="NEG"
  testingOptim$Pred_Class[testingOptim$branchpoint_prob > cutoff_performance$vals[v]] <-  "HC"
  keep=which(testingOptim$branchpoint_prob > cutoff_performance$vals[v] & 
               testingOptim$branchpoint_prob < cutoff_performance$vals[v]+0.01)
  c=confusionMatrix(testingOptim$Pred_Class, testingOptim$Class)
  
  cutoff_performance$Accuracy[v] <- as.numeric(c$overall['Accuracy'])
  cutoff_performance$Balanced_Accuracy[v] <- as.numeric(c$byClass['Balanced Accuracy'])
  cutoff_performance$Sensitivity[v] <- as.numeric(c$byClass['Sensitivity'])
  cutoff_performance$Specificity[v] <- as.numeric(c$byClass['Specificity'])
  cutoff_performance$PPV[v] <- as.numeric(c$byClass['Pos Pred Value'])
  cutoff_performance$NPV[v] <- as.numeric(c$byClass['Neg Pred Value'])
  cutoff_performance$F1[v] <- as.numeric(c$byClass['F1'])
}

write.csv(cutoff_performance, "data/cutoff_performance.csv")
BP_prob_cutoff <- cutoff_performance$vals[which.max(cutoff_performance$F1)]

###### Find optimal probability cutoff ######

cutoff_performance.NB <- data.frame(vals=seq(from=0.01, to=0.99, by=0.01),
                                 Accuracy=NA,
                                 Balanced_Accuracy=NA,
                                 Sensitivity=NA,
                                 Specificity=NA,
                                 PPV=NA,
                                 NPV=NA,
                                 F1=NA)

for(v in seq(along=cutoff_performance.NB$vals)){
  testingOptim$Pred_Class="NEG"
  testingOptim$Pred_Class[testingOptim$branchpoint_prob_NB > cutoff_performance$vals[v]] <-  "HC"
  keep=which(testingOptim$branchpoint_prob_NB > cutoff_performance.NB$vals[v] & 
               testingOptim$branchpoint_prob_NB < cutoff_performance.NB$vals[v]+0.01)
  c=confusionMatrix(testingOptim$Pred_Class, testingOptim$Class)
  
  cutoff_performance.NB$Accuracy[v] <- as.numeric(c$overall['Accuracy'])
  cutoff_performance.NB$Balanced_Accuracy[v] <- as.numeric(c$byClass['Balanced Accuracy'])
  cutoff_performance.NB$Sensitivity[v] <- as.numeric(c$byClass['Sensitivity'])
  cutoff_performance.NB$Specificity[v] <- as.numeric(c$byClass['Specificity'])
  cutoff_performance.NB$PPV[v] <- as.numeric(c$byClass['Pos Pred Value'])
  cutoff_performance.NB$NPV[v] <- as.numeric(c$byClass['Neg Pred Value'])
  cutoff_performance.NB$F1[v] <- as.numeric(c$byClass['F1'])
}

BP_prob_cutoff.NB <- cutoff_performance.NB$vals[which.max(cutoff_performance.NB$F1)]

testing$Pred_Class <- "NEG"
testing$Pred_Class[which(testing$branchpoint_prob >= BP_prob_cutoff)] <- "HC"

testingOptim$Pred_Class <- "NEG"
testingOptim$Pred_Class[which(testingOptim$branchpoint_prob >= BP_prob_cutoff)] <- "HC"

testing$Pred_Class_NB <- "NEG"
testing$Pred_Class_NB[which(testing$branchpoint_prob_NB >= BP_prob_cutoff.NB)] <- "HC"

testingOptim$Pred_Class_NB <- "NEG"
testingOptim$Pred_Class_NB[which(testingOptim$branchpoint_prob_NB >= BP_prob_cutoff.NB)] <- "HC"

#set branchpointer classes
training <- cbind(as.data.frame(branchpoint_df_HCN)[c(trainIndex.pos, trainIndex.neg), ],
                  Class=training$Class)

# combine two testing sets
testing <- cbind(as.data.frame(branchpoint_df_HCN)[c(testIndex.pos, testIndex.neg), ],
                 Class=testing$Class,
                 branchpoint_prob=testing$branchpoint_prob,
                 Pred_Class=testing$Pred_Class,
                 branchpoint_prob_NB=testing$branchpoint_prob_NB,
                 Pred_Class_NB=testing$Pred_Class_NB)
testing$Class <-
  c(rep("HC", testingPositives), rep("NEG", testingNegatives))
testing$new_ID <- branchpoint_df_HCN$new_ID[c(testIndex.pos, testIndex.neg)]
testing$set_name_long <- 
  gsub("NEG","Negatives", gsub("HC", "Mercer branchpoints",testing$Class))

testingOptim <- cbind(as.data.frame(branchpoint_df_HCN)[c(testIndexOptim.pos, testIndexOptim.neg), ],
                 Class=testingOptim$Class,
                 branchpoint_prob=testingOptim$branchpoint_prob,
                 Pred_Class=testingOptim$Pred_Class,
                 branchpoint_prob_NB=testingOptim$branchpoint_prob_NB,
                 Pred_Class_NB=testingOptim$Pred_Class_NB)
testingOptim$Class <-
  c(rep("HC", testingPositives), rep("NEG", testingNegatives))
testingOptim$new_ID <- branchpoint_df_HCN$new_ID[c(testIndexOptim.pos, testIndexOptim.neg)]
testingOptim$set_name_long <- 
  gsub("NEG","Negatives", gsub("HC", "Mercer branchpoints",testingOptim$Class))

testing <- rbind(testing, testingOptim[,match(colnames(testing), colnames(testingOptim))])
testing$test_set <- c(rep("T", length(c(testIndex.pos, testIndex.neg))),
                      rep("O", length(c(testIndexOptim.pos, testIndexOptim.neg))))


###### variable importance ######

#get variable importance from gbm
gbm_imp <- varImp(model,scale = F)
gbm_importance <- data.frame(gbm_imp$importance)
gbm_importance$variable <- rownames(gbm_importance)

#rename features to be more readable
gbm_importance$variableName <- gsub("seq_pos0", "BP: ",gbm_importance$variable)
gbm_importance$variableName <- gsub("seq_pos", "nt +",gbm_importance$variableName)
gbm_importance$variableName <- gsub("seq_neg", "nt -",gbm_importance$variableName)
gbm_importance$variableName <- gsub("dist.1", 
                                "5'exon distance",gbm_importance$variableName)
gbm_importance$variableName <- gsub("dist.2", 
                                "3'exon distance",gbm_importance$variableName)
gbm_importance$variableName <- gsub("canon_hit", 
                                "AG distance ",gbm_importance$variableName)
gbm_importance$variableName <- gsub("newFeat", 
                                "SVM model probability score",
                                gbm_importance$variableName)
gbm_importance$variableName <- gsub("ppt_start", 
                                "polypyrimidine tract distance",
                                gbm_importance$variableName)
gbm_importance$variableName <- gsub("ppt_run_length", 
                                "polypyrimidine tract length",
                                gbm_importance$variableName)

gbm_importance$color <- NA
gbm_importance$color[grepl("nt", gbm_importance$variableName) | grepl("BP", gbm_importance$variableName)] <- "Nucleotide identity"
gbm_importance$color[gbm_importance$variableName %in% 
                       c("SVM model probability score", "polypyrimidine tract distance","polypyrimidine tract length")] <- "Polypyrimidine tract"
gbm_importance$color[gbm_importance$variableName %in% c("AG distance 1", "AG distance 2","AG distance 5")] <- "Splice site locations"
gbm_importance$color[gbm_importance$variableName %in% c("3'exon distance", "5'exon distance")] <- "Exon locations"

#### Variable importance by removal
# performance metrics for each variable (/set) left out
# files <- list.files("data/removeVars/", full.names = TRUE)
# removed_vars <- list()
# meanAccuracy <- vector()
# meanF1 <- vector()
# medianAccuracy <- vector()
# medianF1 <- vector()
# for(f in 64:length(files)){
#   
#   load(files[f])
#   meanAccuracy[f] <- mean(unlist(lapply(confusionMatrixList, function(x) x$overall['Accuracy'])))
#   medianAccuracy[f] <- median(unlist(lapply(confusionMatrixList, function(x) x$overall['Accuracy'])))
#   meanF1[f] <- mean(unlist(lapply(confusionMatrixList, function(x) x$byClass['F1'])))
#   medianF1[f] <- median(unlist(lapply(confusionMatrixList, function(x) x$byClass['F1'])))
#   
#   removed_vars[[f]] <- removeVars
#   message(f)
# }
# removed_vars[[1]] <- NA
# 
# names(removed_vars) <- c("base",lapply(removed_vars,"[[",1)[2:54], "base2",
#                          "distance","seq_neg2","seq_neg3","seq_neg4", "seq_neg5",
#                          "canon" ,"dist+ppt","ppt","seq_pos0" ,"seq_pos1",
#                          "seq_pos2", "seq_pos3" ,"seq_pos4",
#                          "seq_pos5" ,"seq_neg1")
# 
# removedVarPerf <- data.frame(variable=names(removed_vars), 
#                              meanAccuracy, medianAccuracy,
#                              meanF1,medianF1)
# 
# #difference from model with all variables
# removedVarPerf$meanAccuracyChange <- 
#   removedVarPerf$meanAccuracy - removedVarPerf$meanAccuracy[removedVarPerf$variable=="base"] 
# removedVarPerf$meanF1Change <- 
#   removedVarPerf$meanF1 - removedVarPerf$meanF1[removedVarPerf$variable=="base"]
# 
# #remove single nt features -- shhould have no sig effect as can be imputed by values of other 3 nts
# rm <- (str_sub(removedVarPerf$variable, -1,-1) %in% c("A","T","C","G"))
# removedVarPerf <- removedVarPerf[!rm,]
# removedVarPerf <- removedVarPerf[-which(removedVarPerf$variable =="base2"),]
# 
# removedVarPerf$grouped <- "yes"
# removedVarPerf$grouped[removedVarPerf$variable %in% c("dist.1","dist.2",
#                                                       "ppt_start","ppt_run_length",
#                                                       "canon_hit1","canon_hit2","canon_hit3",
#                                                       "canon_hit4","canon_hit5")] <- "no"
# 
# #rename vars
# removedVarPerf$variableName <- gsub("seq_pos0", "BP nt",removedVarPerf$variable)
# removedVarPerf$variableName <- gsub("seq_pos", "nt +",removedVarPerf$variableName)
# removedVarPerf$variableName <- gsub("seq_neg", "nt -",removedVarPerf$variableName)
# removedVarPerf$variableName <- gsub("dist.1",
#                                  "5'exon distance",removedVarPerf$variableName)
# removedVarPerf$variableName <- gsub("dist.2",
#                                  "3'exon distance",removedVarPerf$variableName)
# removedVarPerf$variableName <- gsub("canon_hit",
#                                  "AG distance ",removedVarPerf$variableName)
# removedVarPerf$variableName <- gsub("newFeat",
#                                  "SVM model probability score",
#                                  removedVarPerf$variableName)
# removedVarPerf$variableName <- gsub("ppt_start",
#                                  "PPT distance",
#                                  removedVarPerf$variableName)
# removedVarPerf$variableName <- gsub("ppt_run_length",
#                                  "PPT length",
#                                  removedVarPerf$variableName)
# removedVarPerf$variableName[removedVarPerf$variableName == "dist+ppt"] <-
#   "PPT + exon distance variables"
# removedVarPerf$variableName[removedVarPerf$variableName == "distance"] <-
#   "exon distance variables"
# removedVarPerf$variableName[removedVarPerf$variableName == "distance"] <-
#   "exon distance variables"
# removedVarPerf$variableName[removedVarPerf$variableName == "canon"] <-
#   "AG distance variables"
# removedVarPerf$variableName[removedVarPerf$variableName == "ppt"] <-
#   "PPT variables"
# 
# # Order by accuracy for plotting
# removedVarPerf$order[order(removedVarPerf$meanAccuracy)] <- 1:nrow(removedVarPerf)
# FigureS3A <- ggplot(removedVarPerf, aes(x=factor(order), y=meanAccuracyChange, fill=grouped)) + geom_bar(stat="identity") +
#   scale_x_discrete(name="Variable removed",
#                    labels=removedVarPerf$variableName[order(removedVarPerf$meanAccuracy)]) +
#   scale_y_continuous(name="Change in Accuracy") +
#   scale_fill_manual(name="Variables\ngrouped", values=c("gray40", "gray60")) +
#   theme_figure + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# # Order by F1 for plotting
# removedVarPerf$order[order(removedVarPerf$meanF1Change)] <- 1:nrow(removedVarPerf)
# FigureS3B <- ggplot(removedVarPerf, aes(x=factor(order), y=meanF1Change, fill=grouped)) + geom_bar(stat="identity") +
#   scale_x_discrete(name="Variable removed",
#                    labels=removedVarPerf$variable[order(removedVarPerf$meanF1Change)]) +
#   scale_y_continuous(name="Change in F1") +
#   scale_fill_manual(name="Variables\ngrouped", values=c("gray40", "gray60")) +
#   theme_figure + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# save(removedVarPerf, removed_vars, file="data/removedVarsSummary.Rdata")
load("data/removedVarsSummary.Rdata")

# match leaveOut performance change and gbmImp
m1 <- match(gbm_importance$variable, removedVarPerf$variable)
m2 <- match(stringr::str_sub(gbm_importance$variable, 1,-2), removedVarPerf$variable)
m <- c(m1[which(!is.na(m1))], m2[which(!is.na(m2))])
gbm_importance$meanAccuracyChange <- removedVarPerf$meanAccuracyChange[m]
gbm_importance$meanF1Change <- removedVarPerf$meanF1Change[m]
gbm_importance$grouped <- removedVarPerf$grouped[m]
gbm_importance$variableLO <- removedVarPerf$variable[m]

removedVarPerf[which(!(1:nrow(removedVarPerf) %in% m)),]
gbm_importance <- plyr::arrange(gbm_importance, plyr::desc(Overall))

gbm_importance$variable_order <- factor(gbm_importance$variableName, levels=arrange(gbm_importance, plyr::desc(Overall))$variableName)


#### motif nucleotide importance

keep <- grepl("BP:", gbm_importance$variableName) | grepl("nt ", gbm_importance$variableName)
U2_df <- gbm_importance[keep,]
U2_df$pos <- stringr::str_sub(U2_df$variableName, 4,5)
U2_df$pos[grep("BP", U2_df$variableName)] <- 0
U2_df$nt <- stringr::str_sub(U2_df$variableName, -1,-1)
U2_df$Overall[is.na(U2_df$Overall)] <- 0
U2_df$pos <- as.numeric(U2_df$pos)

U2_df$importanceScaled <- rowMeans(cbind(scale((U2_df$Overall)), scale(U2_df$meanF1Change *-1)))
U2_df$importanceScaled <- U2_df$importanceScaled - min(U2_df$importanceScaled)

###### get sequences to send to SVM-BP ######

#read in gencode v19 exon annotation
gencode_v19.gtf <- gtfToExons("data/gencode.v19.annotation.gtf")
gencode_v19.gtf$chr_loc_start <- paste(seqnames(gencode_v19.gtf), start(ranges(gencode_v19.gtf)), sep="_")
gencode_v19.gtf$chr_loc_end <- paste(seqnames(gencode_v19.gtf), end(ranges(gencode_v19.gtf)), sep="_")

#find the nearest 3'/5' exons
testing$chromosome <- unlist(lapply(str_split(testing$new_ID, "_"), "[[", 1))
testing$strand <- "-"
testing$strand[grep("[+]", testing$new_ID)] <- "+"
testing$loc <- unlist(lapply(str_split(testing$new_ID, "_"), "[[", 2))
testing$start[testing$strand == "-"] <- as.numeric(unlist(lapply(str_split(testing$loc[testing$strand == "-"] , "-"), "[[", 1)))
testing$start[testing$strand == "+"] <- as.numeric(unlist(lapply(str_split(testing$loc[testing$strand == "+"] , "[+]"), "[[", 1)))

testing$exon_starts <- testing$start + 1 + as.numeric(testing$dist.2)
testing$exon_ends <- testing$start + 1 - as.numeric(testing$dist.2)

chr_loc <- with(testing, paste(chromosome, exon_starts, sep="_"))
m <- match(chr_loc, gencode_v19.gtf$chr_loc_start)

chr_loc <- with(testing, paste(chromosome, exon_ends, sep="_"))
m[testing$strand=="-"] <- match(chr_loc, gencode_v19.gtf$chr_loc_end)[testing$strand=="-"]
testing$gene_id <- gencode_v19.gtf$gene_id[m]
testing$gene_type <- gencode_v19.gtf$gene_type[m]
testing$transcript_id <- gencode_v19.gtf$transcript_id[m]
testing$transcript_type <- gencode_v19.gtf$transcript_type[m]
testing$exon_id <- gencode_v19.gtf$exon_id[m]
testing$exon_number <- gencode_v19.gtf$exon_number[m]

#make bed files for testing SVM-BP
#positive strand
testing.pos <- testing[testing$strand == "+",]
testing.pos <- testing.pos[!duplicated(testing.pos$exon_id),]

#make bed format file
bed <- testing.pos[,c('chromosome','exon_starts','exon_starts','exon_id','dist.2','strand')]
bed[,5] <-  0
bed[,2] <- bed[,3] - 101
bed[,3] <- bed[,3] - 1
uid <- "pos_test"

if(rerunFasta == TRUE){
  #write bed file
  write.table(
    bed, sep = "\t", file = paste0("data/intron_",uid,".bed"),
    row.names = F,col.names = F,quote = F
  )
  #convert to fasta using bedtools
  cmd <- paste0(
    "/Applications/apps/bedtools2/bin/bedtools getfasta -fi ",
    "data/GRCh37.p13.genome.fa",
    " -bed data/intron_",uid,".bed -fo data/intron_",uid,".fa -name -s"
  )
  system(cmd)
}
#read .fa
fasta <- fread(paste0("data/intron_",uid,".fa"), header = F)
fasta_pos <- as.data.frame(fasta)
fasta_pos <- as.data.frame(matrix(fasta_pos$V1, ncol = 2,byrow = T))
colnames(fasta_pos) <- c("id", "seq")
fasta_pos$id <- gsub(">","",fasta_pos$id)

#same again for negative strand
testing.neg <- testing[testing$strand == "-",]
testing.neg <- testing.neg[!duplicated(testing.neg$exon_id),]

#make bed format file
bed <- testing.neg[,c('chromosome','exon_ends','exon_ends','exon_id','dist.2','strand')]
bed[,5] <- 0
bed[,3] <- bed[,2] + 100
uid <- "neg_test"

if(rerunFasta == TRUE){
  
  #write bed file
  write.table(
    bed, sep = "\t", file = paste0("data/intron_",uid,".bed"),
    row.names = F,col.names = F,quote = F
  )
  
  #convert to fasta using bedtools
  cmd <- paste0(
    "/Applications/apps/bedtools2/bin/bedtools getfasta -fi ",
    "data/GRCh37.p13.genome.fa",
    " -bed data/intron_",uid,".bed -fo data/intron_",uid,".fa -name -s"
  )
  system(cmd)
}
#read .fa
fasta <- fread(paste0("data/intron_",uid,".fa"), header = F)
fasta_neg <- as.data.frame(fasta)
fasta_neg <- as.data.frame(matrix(fasta_neg$V1, ncol = 2,byrow = T))
colnames(fasta_neg) <- c("id", "seq")
fasta_neg$id <- gsub(">","",fasta_neg$id)

fasta <- rbind(fasta_pos,fasta_neg)

rm(fasta_pos, fasta_neg)

m <- match(testing$exon_id, fasta$id)
fasta_seqs = fasta[,c('id','seq')]
fasta_seqs$id = paste0(">",fasta_seqs$id)
fasta_seqs = fasta_seqs[!duplicated(fasta_seqs$id),]

if(rerunFasta == TRUE){
  #fa seqs for svm-bp
  write.table(
    fasta_seqs, file = paste0("data/testing_SVM_seqs.fa"),
    row.names = F,col.names = F,quote =
      F, sep = "\n"
  )
  
  cmd <- "python ~/Applications/svm-bpfinder/svm_bpfinder.py -i data/testing_SVM_seqs.fa -s Hsap > data/svm-BP_output.txt"
  system(cmd)
}

###### Read in other predictions ######
#csv of all heptamer scores run through HSF branchpoints
HSF_hept <- read.csv("data/HSF_heptamers.csv", header = TRUE)
colnames(HSF_hept) <- c("heptamer","HSF_score")
HSF_hept$HSF_score <- as.numeric(HSF_hept$HSF_score)
write.csv(HSF_hept, file="Tables/TableS1.csv", row.names = F)

#SVM BP output
SVM_BP_finder_output <-
  read.delim("data/svm-BP_output.txt")

testing$svmBPscore <- NA
testing$HSFBPscore <- NA

id_dist <- with(testing, paste(exon_id, dist.2, sep="_"))
id_dist_SVM <- apply(SVM_BP_finder_output[,c(1,3)],1,paste,collapse = "_")

m <- match(id_dist,id_dist_SVM)
testing$svmBPscore <- SVM_BP_finder_output$svm_scr[m]

#SVM BP doesn't provide values for non UNA motifs.
#we'll assign a mimimum score to allow us to use scores for performance metrics
testing$svmBPscore[is.na(testing$svmBPscore)] <-
  min(testing$svmBPscore, na.rm = T) - 0.1

#no appropriate cutoff given in text, pick out best using same method as branchpointer
SVM_BP_cutoff <- data.frame(score_range=seq(min(testing$svmBPscore),
                                            max(testing$svmBPscore),
                                            length.out = 100),
                            Accuracy = NA,
                            F1 = NA)

for (v in seq_along(SVM_BP_cutoff$score_range)) {

  testing[testing$test_set == "O",]$svmBP_class <- "NEG"
  testing$svmBP_class[testing[testing$test_set == "O",]$svmBPscore >= SVM_BP_cutoff$score_range[v]] <- "HC"
  c <- confusionMatrix(testing[testing$test_set == "O",]$svmBP_class, 
                       testing[testing$test_set == "O",]$Class)
  SVM_BP_cutoff$Accuracy[v] <- c$overall['Accuracy']
  SVM_BP_cutoff$F1[v] <- c$byClass['F1']
  
}

SVM_BP_BestCutoff <- SVM_BP_cutoff$score_range[which.max(SVM_BP_cutoff$F1)]

#use that cutoff for classification
testing$svmBP_class <- "NEG"
testing$svmBP_class[testing$svmBPscore >= SVM_BP_BestCutoff] <- "HC"

####heptmer scores
testing$heptamers <- with(testing, paste0(seq_neg5,seq_neg4,seq_neg3,seq_neg2,seq_neg1, seq_pos0, seq_pos1))
m <- match(testing$heptamers, HSF_hept$heptamer)
testing$HSFBPscore <- HSF_hept$HSF_score[m]
testing$HSFBP_class <- NA

#no appropriate cutoff given in text, pick out best using same method as branchpointer
HSF_BP_cutoff <- data.frame(score_range=seq(min(testing$HSFBPscore),
                                            max(testing$HSFBPscore),
                                            length.out = 100),
                            Accuracy = NA,
                            F1 = NA)

for (v in seq_along(HSF_BP_cutoff$score_range)) {
  
  testing[testing$test_set == "O",]$HSFBP_class <- "NEG"
  testing[testing$test_set == "O",]$HSFBP_class[testing[testing$test_set == "O",]$HSFBPscore >= 
                                                  HSF_BP_cutoff$score_range[v]] <- "HC"
  c <- confusionMatrix(testing[testing$test_set == "O",]$HSFBP_class, 
                       testing[testing$test_set == "O",]$Class)
  HSF_BP_cutoff$Accuracy[v] <- c$overall['Accuracy']
  HSF_BP_cutoff$F1[v] <- c$byClass['F1']
  
}

HSF_BP_bestCutoff <- HSF_BP_cutoff$score_range[which.max(HSF_BP_cutoff$F1)]

#use that cutoff for classification
testing$HSFBP_class <- "NEG"
testing$HSFBP_class[testing$HSFBPscore >= HSF_BP_bestCutoff] <- "HC"

#use basic UNA motif for classification
testing$UNA_class <- "NEG"
testing$UNA_class[which(testing$seq_pos0 == "A" &
                               testing$seq_neg2 == "T")] <- "HC"

#use branchpointer for classification
testing$branchpointer_class <- "NEG"
testing$branchpointer_class[which(testing$branchpoint_prob >= BP_prob_cutoff)] <- "HC"

#use branchpointer for classification
testing$branchpointer_class_NB <- "NEG"
testing$branchpointer_class_NB[which(testing$branchpoint_prob_NB >= BP_prob_cutoff.NB)] <- "HC"


###### compare methods for classification ######

#create confusion matrices
c1 <- confusionMatrix(testing$HSFBP_class[testing$test_set=="T"], testing$set[testing$test_set=="T"])
c2 <- confusionMatrix(testing$svmBP_class[testing$test_set=="T"], testing$set[testing$test_set=="T"])
c3 <- confusionMatrix(testing$branchpointer_class[testing$test_set=="T"], testing$set[testing$test_set=="T"])
c4 <- confusionMatrix(testing$UNA_class[testing$test_set=="T"], testing$set[testing$test_set=="T"])
c5 <- confusionMatrix(testing$branchpointer_class_NB[testing$test_set=="T"], testing$set[testing$test_set=="T"])

df <- data.frame(
  HSF = c1$byClass, svmBP = c2$byClass,
  branchpointer = c3$byClass, UNA = c4$byClass, branchpointer_NB = c5$byClass
)
accuracy_row <- data.frame(
  HSF = c1$overall, svmBP = c2$overall,
  branchpointer = c3$overall,UNA = c4$overall, branchpointer_NB = c5$overall
)
df <- rbind(df, accuracy_row)
df <- as.data.frame(t(df))

classification_performance <- df[,c(1:4,7,11,12)]

for(i in 1:ncol(classification_performance)){
  classification_performance[,i] <- round(classification_performance[,i], 3)
}

classification_performance <- cbind(Method=rownames(classification_performance), classification_performance)
classification_performance <- classification_performance[,c(1,6,8,2,3,4,5)]

write.table(
  classification_performance, file = "Tables/Table1.csv", row.names = F,quote = F,sep = ","
)

###### Compare methods using scores for pr/roc curves ######

#computing curves takes a long time
#saved curve objects in branchpointer_curves.Rdata

roc.1 <- roc.curve(scores.class0 =
                   testing$svmBPscore[testing$set=="HC" & testing$test_set=="T"],
                 scores.class1 =
                   testing$svmBPscore[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)
pr.1 <- pr.curve(scores.class0 =
                 testing$svmBPscore[testing$set=="HC"& testing$test_set=="T"],
               scores.class1 =
                 testing$svmBPscore[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)
svm_BP.AUC <- roc.1$auc
svm_BP.pr <- pr.1$auc.integral

roc.2 <- roc.curve(scores.class0 =
                   testing$HSFBPscore[testing$set=="HC"& testing$test_set=="T"],
                 scores.class1 =
                   testing$HSFBPscore[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)
pr.2 <- pr.curve(scores.class0 =
                 testing$HSFBPscore[testing$set=="HC"& testing$test_set=="T"],
               scores.class1 =
                 testing$HSFBPscore[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)

HSFBP.AUC <- roc.2$auc
HSFBP.pr <- pr.2$auc.integral

roc.3 <- roc.curve(scores.class0 =
                 testing$branchpoint_prob[testing$set=="HC"& testing$test_set=="T"],
               scores.class1 =
                 testing$branchpoint_prob[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)
pr.3 <- pr.curve(scores.class0 =
               testing$branchpoint_prob[testing$set=="HC"& testing$test_set=="T"],
             scores.class1 =
               testing$branchpoint_prob[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)
branchpointer.AUC <- roc.3$auc
branchpointer.pr <- pr.3$auc.integral

# ADD NB ROC
roc.4 <- roc.curve(scores.class0 =
                     testing$branchpoint_prob_NB[testing$set=="HC"& testing$test_set=="T"],
                   scores.class1 =
                     testing$branchpoint_prob_NB[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)
pr.4 <- pr.curve(scores.class0 =
                   testing$branchpoint_prob_NB[testing$set=="HC"& testing$test_set=="T"],
                 scores.class1 =
                   testing$branchpoint_prob_NB[testing$set=="NEG"& testing$test_set=="T"], curve=TRUE)
branchpointer_NB.AUC <- roc.4$auc
branchpointer_NB.pr <- pr.4$auc.integral

save(roc.1,roc.2,roc.3,roc.4,
     pr.1,pr.2,pr.3,pr.4, file="data/branchpointer_curves.Rdata")

load("data/branchpointer_curves.Rdata")

roc_curve_SVMBP <- data.frame(roc.1$curve)
colnames(roc_curve_SVMBP) <- c("FPR","TPR","cut")
roc_curve_SVMBP$method <- "SVM_BP"
roc_curve_HSFBP <- data.frame(roc.2$curve)
colnames(roc_curve_HSFBP) <- c("FPR","TPR","cut")
roc_curve_HSFBP$method <- "HSF_BP"
roc_curve_BP <- data.frame(roc.3$curve)
colnames(roc_curve_BP) <- c("FPR","TPR","cut")
roc_curve_BP$method <- "BP"
roc_curve_BP_NB <- data.frame(roc.4$curve)
colnames(roc_curve_BP_NB) <- c("FPR","TPR","cut")
roc_curve_BP_NB$method <- "BP_NB"
roc_curves=rbind(roc_curve_SVMBP,roc_curve_HSFBP,roc_curve_BP,roc_curve_BP_NB)

pr_curve_SVMBP <- data.frame(pr.1$curve)
pr_curve_SVMBP$method <- "SVM_BP"
pr_curve_HSFBP <- data.frame(pr.2$curve)
pr_curve_HSFBP$method <- "HSF_BP"
pr_curve_BP <- data.frame(pr.3$curve)
pr_curve_BP$method <- "BP"
pr_curve_BP_NB <- data.frame(pr.4$curve)
pr_curve_BP_NB$method <- "BP_NB"
pr_curves=rbind(pr_curve_SVMBP,pr_curve_HSFBP,pr_curve_BP,pr_curve_BP_NB)

rm(roc_curve_BP,roc_curve_HSFBP,roc_curve_SVMBP,
   pr_curve_BP,pr_curve_HSFBP,pr_curve_SVMBP)

roc_curves$method <- gsub("SVM_BP", "SVM-BPFinder",roc_curves$method )
roc_curves$method <- gsub("HSF_BP", "HSF",roc_curves$method )
roc_curves$method[roc_curves$method == "BP"] <- "branchpointer"
roc_curves$method[roc_curves$method == "BP_NB"] <- "Naive Bayes"

pr_curves$method <- gsub("SVM_BP", "SVM-BPFinder",pr_curves$method )
pr_curves$method <- gsub("HSF_BP", "HSF",pr_curves$method )
pr_curves$method[pr_curves$method == "BP"] <- "branchpointer"
pr_curves$method[pr_curves$method == "BP_NB"] <- "Naive Bayes"

save(roc_curves, pr_curves, testing, U2_df, cutoff_performance,gbm_importance, file="data/performance_objects.Rdata")
