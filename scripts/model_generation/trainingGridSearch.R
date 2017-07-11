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


load("data/gbm_training_size_sel_2.RData")

trainingRatioSizes <- data.frame(negatives=trainingNegativesVector,
                                 F1 = unlist(lapply(confusionMatrixList, function(x) x$byClass['F1'])),
                                 Accuracy=unlist(lapply(confusionMatrixList, function(x) x$overall['Accuracy'])))

medianF1 <- aggregate(F1 ~ negatives, trainingRatioSizes, median)
bestTrainRatio <- medianF1$negatives[which.max(medianF1$F1)] / trainingPositives

args <- commandArgs(trailingOnly = TRUE)

set.seed(as.numeric(args[1]))

trainingPositives=1000
trainingNegatives=trainingPositives*bestTrainRatio
testingPositives=5000
testingNegatives=testingPositives*20

control <- trainControl(method="repeatedcv", number=10,summaryFunction = f1, 
                        classProbs = TRUE)

confusionMatrixList <- list()

  
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
  
  max_s <- max(0.01, 0.1*min(1, nrow(training)/10000))
  
  gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), 
                          n.trees = seq(from=1000, to=2500, by = 500), 
                          shrinkage = seq(0, max_s,.01)[-1],
                          n.minobsinnode = c(20))
  


  m <- which(colnames(training) == "Class")
  
  modelGridSearchList <- list()
  for(j in 1:nrow(gbmGrid)){
    ptm <- proc.time()
    
    message("Building model")
    
    modelGridSearchList[[j]] <- caret::train(x = as.data.frame(training[, -m]),y = training$Class,
                          trControl = control,
                          method = "gbm",
                          verbose = FALSE, 
                          tuneGrid=gbmGrid[j,])
    
    message("Time to build: ",(proc.time() - ptm)[3])
    
    n = predict(modelGridSearchList[[j]], testing)
    c = confusionMatrix(n, testing$Class)
    confusionMatrixList[[j]] <- confusionMatrix(n, testing$Class)
    
    message("Search ", j, " of ", nrow(gbmGrid))
    message("F1: ", round(c$byClass['F1'], 3))
    message("Acc: ", round(c$overall['Accuracy'], 3))
  }

save.image(file=paste0("data/gbm_gridSearch_",as.numeric(args[1]),".RData"))

