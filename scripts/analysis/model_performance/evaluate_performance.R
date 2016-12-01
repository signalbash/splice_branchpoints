################################################################################
#     Evaluate performance of the branchpointer model on the testing data      #
################################################################################

#Creates Figure 1, S1, 3A

options(stringsAsFactors = F)

library(stringr)
library(data.table)
library(caret)
library(ggplot2)
library(cowplot)
library(reshape2)
library(PRROC)

nt_cols <- c("#359646","#4D7ABE","#FAA859","#CB3634")

###### load in testing data and model and reformat ######
load("data/Preprocessed_data_44_small.RData")
load("data/branchpoint44_ALL.Rdata")

#format data
testing_remade <- rbind(testing[keepAt,], testing[-keepAt,])
testing_remade$Class <- test_slim_remade$Class
testing_remade$newFeat <- test_slim_remade$newFeat

#predict probabilities 
predgbm <- predict(object = objGBMCplus, 
                     test_slim_remade[,predictorNames],"prob")
testing_remade$HC <- predgbm$HC

#recreate testing set
#positives (BPs)
xtesting1 <- (branchpoint_df_HCN[branchpoint_df_HCN$set == "HC",])[-inTrain1,]
xboost1 <- xtesting1[inboost1[c(1:boost_size_pos)],]
xtesting1x <- xtesting1[-inboost1[1:boost_size_pos],]
test_size_pos <- dim(testing1x)[1]
#negatives (non-BPs)
xtesting2 <- (branchpoint_df_HCN[branchpoint_df_HCN$set != "HC",])[-inTrain2,]
xboost2 <- xtesting2[inboost2[c(1:boost_size_neg)],]
xtesting2x <- xtesting2[-inboost2[c(1:boost_size_neg)],]
t2 <- as.numeric(rownames(testing2))
xtesting2 <- xtesting2x[match(t2, as.numeric(rownames(xtesting2x))),]

xtesting <- rbind(xtesting1x,xtesting2)

test_dataset <- rbind(xtesting[keepAt,], xtesting[-keepAt,])
test_dataset <- cbind(test_dataset, testing_remade[,c(54:56)])

rm(branchpoint_df_HCN, boost,boost_N,boost_slim,boost_slimA,boost_slimrmA,
   boost1,boost2,boosting,boostTransformed,filteredDescr,
   test,test_slim,test_slimA,test_slimrmA,testing_remade,
   testing,testing1,testing2,testing1x,testing2x,
   testTransformed,train,train_slim,training, training1,training2,
   trainTransformed,df, predgbm, test_slim_remade)
rm(xtesting2,xtesting1,xtesting2x,xtesting1x,xboost2,xboost1,xtesting)

#grep out the genomic locations from ids
chrom <- grep("chr",unlist(str_split(test_dataset$new_ID,"_")), value = T)
pos <- matrix(unlist(str_split(test_dataset$new_ID,"_")), 
              byrow = T, ncol = 3)[,2]
strand <- rep("-", length(pos))
strand[grepl("[+]", pos)] <- "+"
pos <- gsub("[+]","-", pos)
start <- as.numeric(matrix(unlist(str_split(pos,"-")), 
                           byrow = T, ncol = 2)[,1])
end <- as.numeric(matrix(unlist(str_split(pos,"-")), byrow = T, ncol = 2)[,2])
class <- matrix(unlist(str_split(test_dataset$new_ID,"_")), 
                byrow = T, ncol = 3)[,3]
test_dataset <- cbind(chrom, start,end,strand,test_dataset)

###### Find optimal probability cutoff ######

sens <- vector()
accuracy <- vector()
accuracy_b <- vector()
ppv <- vector()
vals <- seq(from=0.01, to=0.99, by=0.01)

for(v in seq(along=vals)){
  
  test_dataset$Class <- "NEG"
  test_dataset$Class[test_dataset$HC > vals[v]] <- "HC"
  keep <- which(test_dataset$HC > vals[v] & test_dataset$HC < vals[v]+0.01)
  c <- confusionMatrix(test_dataset$Class, test_dataset$set)
  sens[v] <- as.numeric(c$byClass[1])
  ppv[v] <- as.numeric(c$byClass[3])
  accuracy[v] <- as.numeric(c$overall[1])
  accuracy_b[v] <- as.numeric(c$byClass[11])
  
}
F1 <- 2*((ppv*sens)/(ppv+sens))
cutoff_performance <- data.frame(vals,accuracy,accuracy_b,sens,ppv,F1)

write.csv(cutoff_performance, "data/cutoff_performance.csv")
BP_prob_cutoff <- cutoff_performance$vals[which.max(cutoff_performance$F1)]

Figure1F = ggplot(cutoff_performance, aes(x=vals, y=F1))  +
  geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=1, shape=3)+ 
  theme_bw() +
  scale_x_continuous(name="branchpointer probability score cut-off") +
  scale_y_continuous(name="F1") + 
  theme(text=element_text(size=10))

SuppFigure_1A = ggplot(cutoff_performance, aes(x=vals, y=ppv))  +
  geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=1, shape=3)+ 
  theme_bw() +
  scale_x_continuous(name="branchpointer\nprobability score cut-off") +
  scale_y_continuous(name="positive predictive value") + 
  theme(text=element_text(size=10))

SuppFigure_1B = ggplot(cutoff_performance, aes(x=vals, y=sens))  +
  geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=1, shape=3)+ 
  theme_bw() +
  scale_x_continuous(name="branchpointer\nprobability score cut-off") +
  scale_y_continuous(name="sensitivity") + 
  theme(text=element_text(size=10))

SuppFigure_1C=ggplot(cutoff_performance, aes(x=vals, y=accuracy))  +
  geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=1, shape=3)+ 
  theme_bw() +
  scale_x_continuous(name="branchpointer\nprobability score cut-off") +
  scale_y_continuous(name="accuracy") + 
  theme(text=element_text(size=10))

SuppFigure_1D = ggplot(cutoff_performance, aes(x=vals, y=accuracy_b))  +
  geom_smooth(col=nt_cols[4], se=FALSE)+ 
  geom_point(size=1, shape=3)+ 
  theme_bw() +
  scale_x_continuous(name="branchpointer\nprobability score cut-off") +
  scale_y_continuous(name="balanced accuracy") + 
  theme(text=element_text(size=10))

#set branchpointer classes
test_dataset$Class <- "NEG"
test_dataset$Class[test_dataset$HC > BP_prob_cutoff] <-  "HC"
test_dataset$set_name_long <- 
  gsub("NEG","Negatives", gsub("HC", "Mercer branchpoints",test_dataset$set))

###### variable importance ######

#get variable importance from the two models
svm_importance <- data.frame(varImp(rfeProfile1))
svm_importance$variable <- rownames(svm_importance)
svm_importance$used <- 0
svm_importance$used[
  match(rfeProfile1$optVariables, svm_importance$variable)
  ] <- 1

gbm_imp <- varImp(objGBMCplus,scale = F)
gbm_importance <- data.frame(gbm_imp$importance)
gbm_importance$variable <- rownames(gbm_importance)

#combine
svm_importance <- svm_importance[svm_importance$used==1,]
m <- match(svm_importance$variable ,gbm_importance$variable)
svm_importance$gbm_importance <- NA
svm_importance$gbm_importance <- gbm_importance$Overall[m]

#add svm probability score as a feature
svm_importance <- rbind(svm_importance, 
                        c(0,"newFeat",1,
                          gbm_importance$Overall[gbm_importance$variable==
                                                   "newFeat"]))
svm_importance$Overall <- as.numeric(svm_importance$Overall)
svm_importance$gbm_importance <- as.numeric(svm_importance$gbm_importance)

#rename features to be more readable
svm_importance$variable <- gsub("seq_pos0", "BP: ",svm_importance$variable)
svm_importance$variable <- gsub("seq_pos", "nt +",svm_importance$variable)
svm_importance$variable <- gsub("seq_neg", "nt -",svm_importance$variable)
svm_importance$variable <- gsub("dist.1", 
                                "5'exon distance",svm_importance$variable)
svm_importance$variable <- gsub("dist.2", 
                                "3'exon distance",svm_importance$variable)
svm_importance$variable <- gsub("canon_hit", 
                                "AG distance ",svm_importance$variable)
svm_importance$variable <- gsub("newFeat", 
                                "SVM model probability score",
                                svm_importance$variable)
svm_importance$variable <- gsub("ppt_start", 
                                "polypyrimidine tract distance",
                                svm_importance$variable)
svm_importance$variable <- gsub("ppt_run_length", 
                                "polypyrimidine tract length",
                                svm_importance$variable)

svm_importance$color <- "5"
svm_importance$color[svm_importance$variable %in% c("BP: A", "nt -2T")] <- "1"
svm_importance$color[svm_importance$variable %in% c("SVM model probability score", "polypyrimidine tract distance","polypyrimidine tract length")] <- "2"
svm_importance$color[svm_importance$variable %in% c("AG distance 1", "AG distance 2","AG distance 5")] <- "3"
svm_importance$color[svm_importance$variable %in% c("3'exon distance", "5'exon distance")] <- "4"

Figure1D = ggplot(svm_importance, aes(x=Overall, y=gbm_importance, 
                                      label=variable,color=color)) + 
  geom_point() + 
  geom_text(hjust=0,vjust=1, size=2) +
  scale_y_log10(name="GBM variable importance") + 
  theme_bw() + 
  scale_color_manual(values=c(nt_cols, "grey60")) +
  scale_x_continuous(name="SVM variable importance") + 
  theme(text=element_text(size=10), legend.position="none")

rm(gbm_importance)

###### get sequences to send to SVM-BP ######

#read in gencode v12 exon annotation
G12_exons <- as.data.frame(fread("data/genome_annotations/gencode.v12.annotation.exons.txt"))
colnames(G12_exons)[1:9] <- c("chromosome","entry_type","exon_start",
                              "exon_end","strand","gene_id","gene_type",
                              "transcript_id","transcript_type")

#find the nearest 3'/5' exons
chroms <- sort(unique(test_dataset$chrom))
for (c in seq_along(chroms)) {
  
  exon_subset <- G12_exons[which(G12_exons$chromosome == chroms[c] &
                                  G12_exons$strand == "+"),]
  testing_subset <- test_dataset[(test_dataset$chrom == chroms[c] &
                               test_dataset$strand == "+"),]
  exon_starts <- testing_subset$end + testing_subset$dist.2
  m <- match(exon_starts, exon_subset$exon_start)
  testing_subset <- cbind(testing_subset, exon_subset[m,c(3,4,6:9)])
  
  if (exists("test_dataset_nearbygenes")) {
    
    test_dataset_nearbygenes <- rbind(test_dataset_nearbygenes,testing_subset)
    
  }else{
    
    test_dataset_nearbygenes <- testing_subset
    
  }
  
  exon_subset <- G12_exons[which(G12_exons$chromosome == chroms[c] &
                                  G12_exons$strand == "-"),]
  testing_subset <- test_dataset[(test_dataset$chrom == chroms[c] &
                               test_dataset$strand == "-"),]
  exon_ends <- testing_subset$end - testing_subset$dist.2
  m <- match(exon_ends, exon_subset$exon_end)
  testing_subset <- cbind(testing_subset, exon_subset[m,c(3,4,6:9)])
  
  if (exists("test_dataset_nearbygenes")) {
    
    test_dataset_nearbygenes <- rbind(test_dataset_nearbygenes,testing_subset)
    
  }else{
    
    test_dataset_nearbygenes <- testing_subset
    
  }
  
  message(chroms[c])
  
}

m <- match(rownames(test_dataset), rownames(test_dataset_nearbygenes))
test_dataset <- test_dataset_nearbygenes[m,]
rm(test_dataset_nearbygenes, testing_subset, exon_subset)

#make new ids for each intron (not site) to remove duplicates
ids <- apply(test_dataset[,c(35,31)], 1, function(x) {
  paste0(x, collapse = "_")
})
ids <- gsub(" ","",ids)
test_dataset$id <- ids

#make bed files for testing SVM-BP
#positive strand
test_dataset_pos <- test_dataset[test_dataset$strand == "+",]
test_dataset_pos <- test_dataset_pos[which(!duplicated(test_dataset_pos$id)),]

#make bed format file
bed <- test_dataset_pos[,c(1,31,31,6,7,4)]
bed[,5] <-  0
bed[,2] <- bed[,3] - 101
bed[,3] <- bed[,3] - 1
bed[,1] <- gsub("chr","", bed[,1])
uid <- "pos_test"

#write bed file
write.table(
  bed, sep = "\t", file = paste0("intron_",uid,".bed"),
  row.names = F,col.names = F,quote = F
)
#convert to fasta using bedtools
cmd <- paste0(
  "/Applications/apps/bedtools2/bin/bedtools getfasta -fi ", 
  "data/genome_predictions/GRCh37.p13.genome.fa",
  " -bed intron_",uid,".bed -fo intron_",uid,".fa -name -s"
)
system(cmd)

#read .fa
fasta <- fread(paste0("intron_",uid,".fa"), header = F)
fasta_pos <- as.data.frame(fasta)
fasta_pos <- as.data.frame(matrix(fasta_pos$V1, ncol = 2,byrow = T))
colnames(fasta_pos) <- c("new_ID", "seq")
fasta_pos$new_ID <- gsub(">","",fasta_pos$new_ID)
fasta_pos$id <- test_dataset_pos$id

#same again for negative strand
test_dataset_neg <- test_dataset[test_dataset$strand == "-",]
test_dataset_neg <- test_dataset_neg[!duplicated(test_dataset_neg$id),]

#make bed format file
bed <- test_dataset_neg[,c(1,32,32,6,7,4)]
bed[,5] <- 0
bed[,3] <- bed[,2] + 100
bed[,1] <- gsub("chr","", bed[,1])
uid <- "neg_test"
#write bed file
write.table(
  bed, sep = "\t", file = paste0("intron_",uid,".bed"), 
  row.names = F,col.names = F,quote = F
)
#convert to fasta using bedtools
cmd <- paste0(
  "/Applications/apps/bedtools2/bin/bedtools getfasta -fi ", 
  "data/genome_annotations/GRCh37.p13.genome.fa",
  " -bed intron_",uid,".bed -fo intron_",uid,".fa -name -s"
)
system(cmd)

#read .fa
fasta <- fread(paste0("intron_",uid,".fa"), header = F)
fasta_neg <- as.data.frame(fasta)
fasta_neg <- as.data.frame(matrix(fasta_neg$V1, ncol = 2,byrow = T))
colnames(fasta_neg) <- c("new_ID", "seq")
fasta_neg$new_ID <- gsub(">","",fasta_neg$new_ID)
fasta_neg$id <- test_dataset_neg$id

fasta <- rbind(fasta_pos,fasta_neg)

rm(fasta_pos, fasta_neg)

m <- match(test_dataset$id, fasta$id)
test_dataset <- cbind(test_dataset, fasta=fasta$seq[m])
fasta_seqs = test_dataset[,c(37,38)]
fasta_seqs$id = paste0(">",fasta_seqs$id)
fasta_seqs = fasta_seqs[!duplicated(fasta_seqs$id),]

rm(test_dataset_pos,test_dataset_neg, fasta, bed)

#fa seqs for svm-bp
write.table(
  fasta_seqs, file = "data/testing_SVM_seqs.fa", 
  row.names = F,col.names = F,quote =
    F, sep = "\n"
)

#need to add \n between id and seq
#run 1000 sequences at a time on webserver:
#http://regulatorygenomics.upf.edu/Software/SVM_BP/

###### Read in other predictions ######
#csv of all heptamer scores run through HSF branchpoints
HSF_hept <- read.csv("data/HSF_scores.csv", header = FALSE)
colnames(HSF_hept) <- c("heptamer","HSF_score")

#SVM BP output
SVM_BP_finder_output <-
  read.delim("data/svm-BP_output.txt")

test_dataset$svmBPscore <- NA
test_dataset$HSFBPscore <- NA

id_dist <- apply(test_dataset[,c(37,8)],1,paste,collapse = "_")
id_dist_SVM <- apply(SVM_BP_finder_output[,c(1,3)],1,paste,collapse = "_")

m <- match(id_dist,id_dist_SVM)
test_dataset$svmBPscore <- SVM_BP_finder_output$svm_scr[m]

#SVM BP doesn't provide values for non UNA motifs.
#we'll assign a mimimum score to allow us to use scores for performance metrics
test_dataset$svmBPscore[is.na(test_dataset$svmBPscore)] <- 
  min(test_dataset$svmBPscore, na.rm = T) - 0.1

#no appropriate cutoff given in text, pick out best using same method as branchpointer
sens <- vector()
accuracy <- vector()
accuracy_b <- vector()
ppv <- vector()
score_range <- seq(min(test_dataset$svmBPscore), 
                   max(test_dataset$svmBPscore), 
                   length.out = 100)

for (v in seq_along(score_range)) {
  
  test_dataset$svmBP_class <- "NEG"
  test_dataset$svmBP_class[test_dataset$svmBPscore >= score_range[v]] <- "HC"
  c <- confusionMatrix(test_dataset$svmBP_class, test_dataset$set)
  sens[v] <- as.numeric(c$byClass[1])
  ppv[v] <- as.numeric(c$byClass[3])
  accuracy[v] <- as.numeric(c$overall[1])
  accuracy_b[v] <- as.numeric(c$byClass[8])
  
}
F1 = 2 * ((ppv * sens) / (ppv + sens))
cutoff_performance <- data.frame(score_range,accuracy,
                                 accuracy_b,sens,ppv,F1)
SVM_BP_cutoff <- cutoff_performance$score_range[
  which.max(cutoff_performance$F1)]

#use that cutoff for classification
test_dataset$svmBP_class <- "NEG"
test_dataset$svmBP_class[test_dataset$svmBPscore >= SVM_BP_cutoff] <- "HC"

####heptmer scores
heptamers <- apply(test_dataset[,c(16:22)], 1, paste0,collapse = "")
test_dataset$heptamers <- heptamers
m <- match(test_dataset$heptamers, HSF_hept$heptamer)
test_dataset$HSFBPscore <- HSF_hept$HSF_score[m]

#no appropriate cutoff given in text, pick out best using same method as branchpointer
sens <- vector()
accuracy <- vector()
accuracy_b <- vector()
ppv <- vector()
score_range <- seq(min(test_dataset$HSFBPscore), 
                   max(test_dataset$HSFBPscore), 
                   length.out = 100)

for (v in seq_along(score_range)) {
  
  test_dataset$HSFBP_class <- "NEG"
  test_dataset$HSFBP_class[test_dataset$HSFBPscore >= score_range[v]] <- "HC"
  c <- confusionMatrix(test_dataset$HSFBP_class, test_dataset$set)
  sens[v] <- as.numeric(c$byClass[1])
  ppv[v] <- as.numeric(c$byClass[3])
  accuracy[v] <- as.numeric(c$overall[1])
  accuracy_b[v] <- as.numeric(c$byClass[8])
  
  }
F1 <- 2 * ((ppv * sens) / (ppv + sens))
cutoff_performance <- data.frame(score_range,accuracy,accuracy_b,sens,ppv,F1)
HSF_BP_cutoff <- cutoff_performance$score_range[
  which.max(cutoff_performance$F1)]

#use that cutoff for classification
test_dataset$HSFBP_class <- "NEG"
test_dataset$HSFBP_class[test_dataset$HSFBPscore >= HSF_BP_cutoff] <- "HC"

#use basic UNA motif for classification
test_dataset$UNA_class <- "NEG"
test_dataset$UNA_class[which(test_dataset$seq_pos0 == "A" &
                               test_dataset$seq_neg2 == "T")] <- "HC"

#use branchpointer for classification
test_dataset$branchpointer_class <- "NEG"
test_dataset$branchpointer_class[which(test_dataset$HC >= 0.5)] <- "HC"

###### compare methods for classification ######

#create confusion matrices
c1 <- confusionMatrix(test_dataset$HSFBP_class, test_dataset$set)
c2 <- confusionMatrix(test_dataset$svmBP_class, test_dataset$set)
c3 <- confusionMatrix(test_dataset$branchpointer_class, test_dataset$set)
c4 <- confusionMatrix(test_dataset$UNA_class, test_dataset$set)

df <- data.frame(
  HSF = c1$byClass, svmBP = c2$byClass, 
  branchpointer = c3$byClass, UNA = c4$byClass
)
accuracy_row <- data.frame(
  HSF = c1$overall, svmBP = c2$overall, 
  branchpointer = c3$overall,UNA = c4$overall
)
df <- rbind(df, accuracy_row)
df <- as.data.frame(t(df[c(9,1,2,3,4,8),]))
F1 <- 2 * ((df$`Pos Pred Value` * df$Sensitivity) / (df$`Pos Pred Value` +
                                                      df$Sensitivity))
classification_performance <- cbind(df, F1)

rm(df,accuracy_row, cutoff_performance)

write.table(
  classification_performance, file = "data/comparison_table.txt", row.names = T,quote = F,sep = "\t"
)

###### Compare methods using scores for pr/roc curves ######

#computing curves takes a long time
#saved curve objects in branchpointer_curves.Rdata

roc.1 <- roc.curve(scores.class0 = 
                   test_dataset$svmBPscore[test_dataset$set=="HC"], 
                 scores.class1 = 
                   test_dataset$svmBPscore[test_dataset$set=="NEG"], curve=TRUE)
pr.1 <- pr.curve(scores.class0 = 
                 test_dataset$svmBPscore[test_dataset$set=="HC"], 
               scores.class1 = 
                 test_dataset$svmBPscore[test_dataset$set=="NEG"], curve=TRUE)
svm_BP.AUC <- roc.1$auc
svm_BP.pr <- pr.1$auc.integral

roc.2 <- roc.curve(scores.class0 = 
                   test_dataset$HSFBPscore[test_dataset$set=="HC"], 
                 scores.class1 = 
                   test_dataset$HSFBPscore[test_dataset$set=="NEG"], curve=TRUE)
pr.2 <- pr.curve(scores.class0 = 
                 test_dataset$HSFBPscore[test_dataset$set=="HC"], 
               scores.class1 =
                 test_dataset$HSFBPscore[test_dataset$set=="NEG"], curve=TRUE)

HSFBP.AUC <- roc.2$auc
HSFBP.pr <- pr.2$auc.integral

roc.3 <- roc.curve(scores.class0 = 
                 test_dataset$HC[test_dataset$set=="HC"], 
               scores.class1 = 
                 test_dataset$HC[test_dataset$set=="NEG"], curve=TRUE)
pr.3 <- pr.curve(scores.class0 = 
               test_dataset$HC[test_dataset$set=="HC"], 
             scores.class1 = 
               test_dataset$HC[test_dataset$set=="NEG"], curve=TRUE)
branchpointer.AUC <- roc.3$auc
branchpointer.pr <- pr.3$auc.integral

save(roc.1,roc.2,roc.3,
     pr.1,pr.2,pr.3, file="data/branchpointer_curves.Rdata")

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
roc_curves=rbind(roc_curve_SVMBP,roc_curve_HSFBP,roc_curve_BP)

pr_curve_SVMBP <- data.frame(pr.1$curve)
pr_curve_SVMBP$method <- "SVM_BP"
pr_curve_HSFBP <- data.frame(pr.2$curve)
pr_curve_HSFBP$method <- "HSF_BP"
pr_curve_BP <- data.frame(pr.3$curve)
pr_curve_BP$method <- "BP"
pr_curves=rbind(pr_curve_SVMBP,pr_curve_HSFBP,pr_curve_BP)

rm(roc_curve_BP,roc_curve_HSFBP,roc_curve_SVMBP,
   pr_curve_BP,pr_curve_HSFBP,pr_curve_SVMBP)

roc_curves$method <- gsub("SVM_BP", "SVM-BPFinder",roc_curves$method )
roc_curves$method <- gsub("HSF_BP", "HSF",roc_curves$method )
roc_curves$method[roc_curves$method == "BP"] <- "branchpointer"

pr_curves$method <- gsub("SVM_BP", "SVM-BPFinder",pr_curves$method )
pr_curves$method <- gsub("HSF_BP", "HSF",pr_curves$method )
pr_curves$method[pr_curves$method == "BP"] <- "branchpointer"

#Figure1E and S2
#ROC curves
Figure1E=ggplot(roc_curves, aes(x = FPR,y = TPR, col = method)) + 
  geom_line() + 
  theme_bw() + 
  theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches"), 
        legend.position=c(0.6,0.3)) + 
  scale_color_manual(values=nt_cols)

FigureS2=ggplot(pr_curves, aes(x = X1,y = X2, col = method)) + 
  geom_line() + 
  theme_bw() + 
  theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches")) + 
  scale_color_manual(values=nt_cols) +  
  labs(x="Recall",y="Precision")


###### Make Figure pdfs ######

Figure1 = ggdraw() + 
  draw_plot(Figure1D, 0,0,0.45,0.5) + 
  draw_plot(Figure1E, 0.45,0,0.3,0.5) + 
  draw_plot(Figure1F,0.75,0,0.25,0.5) +
  draw_plot_label(c("D","E","F"), c(0,0.45,0.75), c(0.5,0.5,0.5), size=18)

FigureS1 = ggdraw() + draw_plot(SuppFigure_1A, 0,0.5,0.5,0.5) + 
  draw_plot(SuppFigure_1B, 0.5,0.5,0.5,0.5) + 
  draw_plot(SuppFigure_1C, 0,0,0.5,0.5) + 
  draw_plot(SuppFigure_1D, 0.5,0,0.5,0.5) + 
  draw_plot_label(c("A","B","C","D"), c(0,0.5,0,0.5), c(1,1,0.5,0.5), size=18)

pdf("Figures/Figure1.pdf", useDingbats = F, height=4.5, width=6.69)
Figure1
dev.off()

pdf("Figures/FigureS1.pdf", useDingbats = F, height=6.69, width=6.69)
FigureS1
dev.off()

pdf("Figures/FigureS2.pdf", 3.35,2.2,useDingbats = F)
FigureS2
dev.off()

test_dataset$short_motif <- factor(with(test_dataset, paste0(seq_neg2, "N", seq_pos0)), levels=c(
"ANA","CNA","GNA","TNA","ANC","CNC","GNC","TNC", "ANG","CNG","GNG","TNG", "ANT","CNT","GNT","TNT"))

Figure_S3 <- ggplot(test_dataset[test_dataset$HC >= 0.5,], aes(x=HC, fill = short_motif)) + 
  geom_histogram(bins=50) +
  theme_bw() + 
  theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches")) + 
  guides(fill = guide_legend(ncol = 4)) +
  scale_fill_manual(values = c("#00441b","#238b45","#74c476","#c7e9c0",
                               "#08306b","#2171b5",'#6baed6',"#c6dbef",
                               "#7f2704","#d94801",'#fd8d3c',"#fdd0a2",
                               "#67000d","#cb181d","#fb6a4a","#fcbba1"), drop=FALSE, name="motif")+
  scale_x_continuous(name="branchpointer probability score")
  

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend_FigureS3=g_legend(Figure_S3)

FigureS3=ggdraw() + draw_plot(Figure_S3 + theme(legend.position="none"), 0,0,0.7,1) + 
  draw_grob(legend_FigureS3,0.7,0,0.3,1)

pdf("Figures/FigureS3.pdf", 6.69,4,useDingbats = F)
FigureS3
dev.off()

