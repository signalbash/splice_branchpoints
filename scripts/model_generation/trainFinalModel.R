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

trainingPositives=40000
trainingNegatives=trainingPositives * 8
testingPositives=5000
testingNegatives=testingPositives * 20

control <- trainControl(method="repeatedcv", number=10,summaryFunction = f1, 
                        classProbs = TRUE)


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

gbmGrid <-  expand.grid(interaction.depth = c(21), 
                        n.trees = 2000, 
                        shrinkage = 0.005,
                        n.minobsinnode = c(20))
  
model <-train(x = as.data.frame(training[, -m]),y = training$Class,
              method="gbm",
              trControl = control,
              metric = "F1",tuneGrid=gbmGrid)
  
message("Time to build: ",(proc.time() - ptm)[3])
  
n = predict(model, testing)
c = confusionMatrix(n, testing$Class)
message("F1: ", round(c$byClass['F1'], 3))
message("Acc: ", round(c$overall['Accuracy'], 3))

save.image(file=paste0("data/gbm_final_model_",as.numeric(args[1]),".RData"))

