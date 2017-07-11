##ML testing

options(stringsAsFactors = F)

library(caret)
library(kernlab)
library(MLmetrics)


## See http://topepo.github.io/caret/training.html#metrics
f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
  c(F1 = f1_val)
}


load("data/Preprocessed_data_small.RData")
message("loaded preprocessed values")
args <- commandArgs(trailingOnly = TRUE)

r <- as.numeric(args[1])

li <- list(name=c("distance", 
                  "ppt", 
                  "seq_pos0", "seq_pos1", "seq_pos2", "seq_pos3", "seq_pos4", "seq_pos5",
                  "seq_neg1", "seq_neg2", "seq_neg3", "seq_neg4", "seq_neg5",
                  "canon",
                  "dist+ppt"),
           vars=list(c("dist.1", "dist.2"),
                     c("ppt_start", "ppt_run_length"),
                     c("seq_pos0A", "seq_pos0C","seq_pos0G", "seq_pos0T"),
                     c("seq_pos1A", "seq_pos1C","seq_pos1G", "seq_pos1T"),
                     c("seq_pos2A", "seq_pos2C","seq_pos2G", "seq_pos2T"),
                     c("seq_pos3A", "seq_pos3C","seq_pos3G", "seq_pos3T"),
                     c("seq_pos4A", "seq_pos4C","seq_pos4G", "seq_pos4T"),
                     c("seq_pos5A", "seq_pos5C","seq_pos5G", "seq_pos5T"),
                     c("seq_neg1A", "seq_neg1C","seq_neg1G", "seq_neg1T"),
                     c("seq_neg2A", "seq_neg2C","seq_neg2G", "seq_neg2T"),
                     c("seq_neg3A", "seq_neg3C","seq_neg3G", "seq_neg3T"),
                     c("seq_neg4A", "seq_neg4C","seq_neg4G", "seq_neg4T"),
                     c("seq_neg5A", "seq_neg5C","seq_neg5G", "seq_neg5T"),
                     c("canon_hit1","canon_hit2","canon_hit3","canon_hit4","canon_hit5"),
                     c("dist.1", "dist.2","ppt_start", "ppt_run_length")))

removeVarsName <- li$name
removeVars <- unlist(li$vars[r])

set.seed(1)

filteredDescr=filteredDescr[,-1]

trainingPositives=4000
trainingNegatives=trainingPositives * 8
testingPositives=5000
testingNegatives=testingPositives * 20

control <- trainControl(method="repeatedcv", number=10,summaryFunction = f1, 
                        classProbs = TRUE)

modelList <- list()
confusionMatrixList <- list()

for(reps in 1:5){
  trainIndex.pos <- sample(index.pos.train, trainingPositives)
  trainIndex.neg <- sample(index.neg.train, trainingNegatives)
  training <- filteredDescr[c(trainIndex.pos, trainIndex.neg), ]
  trainingClass <-
    c(rep("HC", trainingPositives), rep("NEG", trainingNegatives))
    
  testIndex.pos <- sample(index.pos.test, testingPositives)
  testIndex.neg <- sample(index.neg.test, testingNegatives)
  testing <- filteredDescr[c(testIndex.pos, testIndex.neg), ]
  testingClass <-
    c(rep("HC", testingPositives), rep("NEG", testingNegatives))
  
  training <- as.data.frame(training)
  testing <- as.data.frame(testing)
  for (n in 1:(length(colnames(training)))) {
    training[, n] <- as.numeric(training[, n])
    testing[, n] <- as.numeric(testing[, n])
  }
    
  nzv <- nearZeroVar(training)
  if (length(nzv > 0)) {
    training = training[, -nzv]
    testing = testing[, -nzv]
  }
    
  preProcValues <-
    preProcess(training, method = c("center", "scale"))
    
  training <- as.data.frame(predict(preProcValues, training))
  testing <- as.data.frame(predict(preProcValues, testing))
  
  training <- cbind(training , Class = as.factor(trainingClass))
  testing <- cbind(testing, Class = as.factor(testingClass))
    
  message("Building model ", reps)
  ptm <- proc.time()
    
  m <- which(colnames(training) %in% c("Class", removeVars))
  
  gbmGrid <-  expand.grid(interaction.depth = c(21), 
                          n.trees = 2000, 
                          shrinkage = 0.005,
                          n.minobsinnode = c(20))
    
  model <- train(x = as.data.frame(training[, -m]),y = training$Class,
                method="gbm",
                trControl = control,
                metric = "F1",tuneGrid=gbmGrid)
    
  message("Time to build: ",(proc.time() - ptm)[3])
    
  n = predict(model, testing)
  c = confusionMatrix(n, testing$Class)
  message("F1: ", round(c$byClass['F1'], 3))
  message("Acc: ", round(c$overall['Accuracy'], 3))
  modelList[[reps]] <- model
  confusionMatrixList[[reps]] <- c
  
  save(modelList,confusionMatrixList, file=paste0("data/gbm_models_removeVarSets_",as.numeric(args[1]),".RData"))
  
}
save.image(file=paste0("data/gbm_models_removeVarSets_",as.numeric(args[1]),".RData"))


