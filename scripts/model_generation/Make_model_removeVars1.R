options(stringsAsFactors = F)

library(caret)
library(kernlab)
library(randomForest)
library(pROC)
library(stringr)

load("data/branchpoint_ALL.Rdata")

set.seed(1)

model1_allFeats <- ksvm(Class~., 
                        data=train_slim,
                        kernal="rbfdot",
                        kpar="automatic",
                        C=4,prob.model=T)
message("SVM model made")

preds <- predict(object=model1_allFeats, test_slim)
c_svm <- confusionMatrix(preds, test_slim$Class)

#predict using SVM
newFeat <- predict(model1_allFeats, test_slim, "probabilities")
test_slim$newFeat <- newFeat[,1]
newFeat <- predict(model1_allFeats, boost_slim, "probabilities")
boost_slim$newFeat <- newFeat[,1]

model2_allFeats <- train(Class~.,
                         boost_slim,
                         method='gbm',
                         trControl=myControl)
message("gbm made")

preds <- predict(object=model2_allFeats, test_slim)
c_boost <- confusionMatrix(preds, test_slim$Class)

######################################################################

model1_list <- list()
model1c_list <- list()
model2_list <- list()
model2c_list <- list()

removeVars <- colnames(train)
removeVars <- removeVars[which(!(removeVars == "Class"))]

removeVars <- removeVars[removeNum]

for(r in seq_along(removeVars)){

    set.seed(1)

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

    boost_slim.rm <- boost_slim[, -match(c(removeVars), 
        colnames(boost_slim))]

    model2 <-
        train(Class ~ ., boost_slim.rm,
              method = 'gbm',
              trControl = myControl)

    message("gbm made")

    preds <- predict(model2, test_slim)

    model2_c <- confusionMatrix(preds, test_slim$Class)

    model1_list[[r]] <- model1
    model1c_list[[r]] <- model1_c
    model2_list[[r]] <- model2
    model2c_list[[r]] <- model2_c

}

save(removeVars, model1_list,model1c_list, model2, model2c_list, 
     c_svm,c_boost,model1_allFeats,model2_allFeats, file="data/remove1_vars.Rdata")

