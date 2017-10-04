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

set.seed(as.numeric(args[1]))

filteredDescr=filteredDescr[,-1]

trainingPositives=500
trainingNegativesVector=c(500,1000,2000,4000,5000,6000,7000,8000,9000,10000,15000, 20000, 30000, 40000, 50000, 60000)
testingPositives=5000
testingNegatives=testingPositives*20

control <- trainControl(method="repeatedcv", number=10,summaryFunction = f1, 
                        classProbs = TRUE)

confusionMatrixList <- list()

# repeat 5x with subsampling  
trainingNegativesVector <- rep(trainingNegativesVector, each = 5)

for(i in seq_along(trainingNegativesVector)){
  
  trainingNegatives=trainingNegativesVector[i]
  
  
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
  
  message("Building model")
  ptm <- proc.time()
  
  m <- which(colnames(training) == "Class")
  
  model <-train(x = as.data.frame(training[, -m]),y = training$Class,
              method="nb",
              trControl = control,
              metric = "F1")
  
  message("Time to build: ",(proc.time() - ptm)[3])
  
  n = predict(model, testing)
  c = confusionMatrix(n, testing$Class)
  confusionMatrixList[[i]] <- confusionMatrix(n, testing$Class)
  
  message(trainingPositives, ":", trainingNegativesVector[i])
  message("F1: ", round(c$byClass['F1'], 3))
  message("Acc: ", round(c$overall['Accuracy'], 3))
  
}

save.image(file=paste0("data/nb_training_size_sel_",as.numeric(args[1]),".RData"))


#


load("~/clusterS/splice_branchpoints/data/nb_training_size_sel_1.RData")

df <- data.frame(size=trainingNegativesVector, f1=unlist(lapply(confusionMatrixList, function(x) x$byClass['F1'])),
                                                         acc=unlist(lapply(confusionMatrixList, function(x) x$byClass['Accuracy'])))


                                                         
                                                         plot(trainingNegativesVector, unlist(lapply(confusionMatrixList, function(x) x$byClass['F1']))



