##ML testing

options(stringsAsFactors = F)

library(caret)
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
trainingNegatives=trainingPositives * 100
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
  
model <-train(x = as.data.frame(training[, -m]),y = training$Class,
              method="nb",
              trControl = control,
              metric = "F1")
  
message("Time to build: ",(proc.time() - ptm)[3])
  
n = predict(model, testing)
c = confusionMatrix(n, testing$Class)
message("F1: ", round(c$byClass['F1'], 3))
message("Acc: ", round(c$overall['Accuracy'], 3))

save.image(file=paste0("data/nb_final_model_",as.numeric(args[1]),".RData"))

model.nb <- model
n = predict(model.nb, testing, "prob")
testing$branchpoint_prob_NB=n$HC

###### Second testing set for F1 optimisation ######

set.seed(as.numeric(args[1]) + 1)
testIndexOptim.pos <- sample(index.pos.test[which(!(index.pos.test %in% testIndex.pos))], 
                             testingPositives)
testIndexOptim.neg <- sample(index.neg.test[which(!(index.neg.test %in% testIndex.neg))], 
                             testingNegatives)
testingOptim <- filteredDescr[c(testIndexOptim.pos, testIndexOptim.neg), ]
testingClass <-
  c(rep("HC", testingPositives), rep("NEG", testingNegatives))

testingOptim <- as.data.frame(testingOptim)
for (n in 1:(length(colnames(testingOptim)))) {
  testingOptim[, n] <- as.numeric(as.character(testingOptim[, n]))
}

if (length(nzv > 0)) {
  testingOptim = testingOptim[, -nzv]
}

testingOptim <- as.data.frame(predict(preProcValues, testingOptim))
testingOptim <- cbind(testingOptim, Class = as.factor(testingClass))

n = predict(model.nb, testingOptim, "prob")
testingOptim$branchpoint_prob_NB=n$HC


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

nb_cutoff_performance <- cutoff_performance

testing$Pred_Class_NB <- "NEG"
testing$Pred_Class_NB[testing$branchpoint_prob_NB > cutoff_performance$vals[which.max(cutoff_performance$F1)]] <- "HC"
testing_nb <- testing[,c('Class','branchpoint_prob_NB','Pred_Class_NB')]

save(testing_nb, nb_cutoff_performance, model.nb, file="data/nb_performance.RData")
save(model.nb, file="data/nb_model_only.RData")

