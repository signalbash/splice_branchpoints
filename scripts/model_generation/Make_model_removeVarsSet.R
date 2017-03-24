options(stringsAsFactors = F)

args <- commandArgs(trailingOnly = TRUE)

r <- as.numeric(args[1])

library(caret)
library(kernlab)
library(randomForest)
library(pROC)
library(stringr)

load("data/branchpoint_ALL.Rdata")

set.seed(1)

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


modelr1_list <- list()
modelr1c_list <- list()
modelr2_list <- list()
modelr2c_list <- list()

removeVarsName <- li$name

for(r in seq_along(removeVarsName)){
  set.seed(1)
  removeVars <- unlist(li$vars[r])

  train_slim.rm <-
      train_slim[, -match(removeVars, colnames(train_slim))]

  model1 <- ksvm(
      Class ~ .,
      data = train_slim.rm,
      kernal = "rbfdot",
      kpar = "automatic",
      C = 4,
      prob.model = T
  )
  message("SVM model made")

  preds <- predict(model1, test_slim)
  model1_c <- confusionMatrix(preds, test_slim$Class)

  #predict using SVM
  newFeat <- predict(model1, test_slim, "probabilities")
  test_slim$newFeat <-  newFeat[, 1]

  newFeat <- predict(model1, boost_slim, "probabilities")
  boost_slim$newFeat <- newFeat[, 1]

  boost_slim.rm <- boost_slim[, -match(c(removeVars), colnames(boost_slim))]

  model2 <-
      train(Class ~ ., boost_slim.rm,
            method = 'gbm',
            trControl = myControl)

  message("gbm made")

  preds <- predict(model2, test_slim)

  model2_c <- confusionMatrix(preds, test_slim$Class)

  modelr1_list[[r]] <- model1
  modelr1c_list[[r]] <- model1_c
  modelr2_list[[r]] <- model2
  modelr2c_list[[r]] <- model2_c

}

save(li, modelr1_list, modelr1c_list,modelr2_list, modelr2_list, file="data/removeSet_vars.Rdata")