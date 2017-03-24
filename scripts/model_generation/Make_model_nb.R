options(stringsAsFactors = F)

library(caret)
library()
library(pROC)
library(stringr)

load("data/branchpoint_ALL.Rdata")

set.seed(1)

model_nb1 <- train(train_slim[,-match("Class", colnames(train_slim))],
                         train_slim$Class,
                         method='nb',
                         trControl=myControl)

preds <- predict(object=model_nb1, test_slim)
c_nb1 <- confusionMatrix(preds, test_slim$Class)

preds_prob1 <- predict(object=model_nb1, test_slim, "prob")

# both training sets
train_and_boost <- rbind(train_slim, boost_slim[match(colnames(train_slim), colnames(boost_slim))])

model_nb2 <- train(train_and_boost[,-match("Class", colnames(train_and_boost))],
                   train_and_boost$Class,
                   method='nb',
                   trControl=myControl)
preds <- predict(object=model_nb2, test_slim)
c_nb2 <- confusionMatrix(preds, test_slim$Class)

preds_prob2 <- predict(object=model_nb2, test_slim, "prob")

test_slim$HC_nb1 <- preds_prob1$HC
test_slim$HC_nb2 <- preds_prob2$HC

test_slim$set <- test_slim$Class

###### Find optimal probability cutoff - nb ######

sens <- vector()
accuracy <- vector()
accuracy_b <- vector()
ppv <- vector()
vals <- seq(from=0.01, to=0.99, by=0.01)

for(v in seq(along=vals)){
  
  test_slim$Class <- "NEG"
  test_slim$Class[test_slim$HC_nb1 > vals[v]] <- "HC"
  keep <- which(test_slim$HC_nb1 > vals[v] & test_slim$HC_nb1 < vals[v]+0.01)
  c <- confusionMatrix(test_slim$Class, test_slim$set)
  sens[v] <- as.numeric(c$byClass[1])
  ppv[v] <- as.numeric(c$byClass[3])
  accuracy[v] <- as.numeric(c$overall[1])
  accuracy_b[v] <- as.numeric(c$byClass[11])
  
}
F1 <- 2*((ppv*sens)/(ppv+sens))
cutoff_performance_nb1 <- data.frame(vals,accuracy,accuracy_b,sens,ppv,F1)

sens <- vector()
accuracy <- vector()
accuracy_b <- vector()
ppv <- vector()
vals <- seq(from=0.01, to=0.99, by=0.01)

for(v in seq(along=vals)){
  
  test_slim$Class <- "NEG"
  test_slim$Class[test_slim$HC_nb2 > vals[v]] <- "HC"
  keep <- which(test_slim$HC_nb2 > vals[v] & test_slim$HC_nb2 < vals[v]+0.01)
  c <- confusionMatrix(test_slim$Class, test_slim$set)
  sens[v] <- as.numeric(c$byClass[1])
  ppv[v] <- as.numeric(c$byClass[3])
  accuracy[v] <- as.numeric(c$overall[1])
  accuracy_b[v] <- as.numeric(c$byClass[11])
  
}
F1 <- 2*((ppv*sens)/(ppv+sens))
cutoff_performance_nb2 <- data.frame(vals,accuracy,accuracy_b,sens,ppv,F1)

val <- cutoff_performance_nb2$vals[which.max(cutoff_performance_nb2$F1)]
test_slim$Class <- "NEG"
test_slim$Class[test_slim$HC_nb2 > val] <- "HC"
keep <- which(test_slim$HC_nb2 > val & test_slim$HC_nb2 < val+0.01)
c_nb <- confusionMatrix(test_slim$Class, test_slim$set)

save.image("data/nb_classifiers.Rdata")

