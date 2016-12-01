options(stringsAsFactors = F)

library(caret)
library(kernlab)
library(randomForest)
library(pROC)
library(stringr)

#load in shuffle optimised models
load("data/Preprocessed_data_small.RData")
HCN_names=colnames(branchpoint_df_HCN)


#cvalue=as.numeric(read.table("data/best_c.txt")[1,1])
#multiplier=as.numeric(read.table("data/best_ratio.txt")[1,1])
best_model=read.table("data/best_rfe_opt.txt")[1,1]

data1=paste0("data/models/",best_model)
data2=paste0("data/models/", gsub("_optdf","",best_model))
data3=paste0("data/bp_2000.16000-5000.1e+05_Cval4_rep1_splicing_model_rfe.RData") 

load(data1)
load(data2)
load(data3)

opt=c(rfeProfile1$optVariables, "Class")

ctrl <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10)

#get opt features
df=df[-2,]
F1=2*((df$pos_pred_value*df$sensitivity)/(df$pos_pred_value+df$sensitivity))
df=cbind(df, F1)
df$o_v[which.max(F1)]
opt=opt[1:df$o_v[which.max(F1)]]


rm=which(colnames(filteredDescr)=="set")
if(length(rm >0)){
  filteredDescr=filteredDescr[,-rm]
}


##center and scale and split

train_size_pos=20000
train_size_neg=train_size_pos*as.numeric(read.table("data/best_ratio.txt")[1,1])

#filteredDescr=filteredDescr[,-1]

inTrain1 <- sample(seq(along=branchpoint_df_HCN$set[branchpoint_df_HCN$set=="HC"]), train_size_pos)
inTrain2 <- sample(seq(along=branchpoint_df_HCN$set[branchpoint_df_HCN$set!="HC"]), train_size_neg)

boost_size_pos=20000
boost_size_neg=boost_size_pos*as.numeric(read.table("data/best_ratio.txt")[1,1])

training1=(filteredDescr[which(branchpoint_df_HCN$set=="HC"),])[inTrain1,]
training2=(filteredDescr[branchpoint_df_HCN$set!="HC",])[inTrain2,]
training=rbind(training1,training2)
trainingClass1 <- (branchpoint_df_HCN$set[branchpoint_df_HCN$set=="HC"])[inTrain1]
trainingClass2 <- (branchpoint_df_HCN$set[branchpoint_df_HCN$set!="HC"])[inTrain2]
trainingClass=c(trainingClass1,trainingClass2)

testing1=(filteredDescr[branchpoint_df_HCN$set=="HC",])[-inTrain1,]
inboost1 <- sample(seq(along=testing1[,1]), length(testing1[,1]))
boost1=testing1[inboost1[c(1:boost_size_pos)],]
testing1x=testing1[-inboost1[1:boost_size_pos],]
test_size_pos=dim(testing1x)[1]

testing2=(filteredDescr[branchpoint_df_HCN$set!="HC",])[-inTrain2,]
inboost2 <- sample(seq(along=testing2[,1]), length(testing2[,1]))
boost2=testing2[inboost2[c(1:boost_size_neg)],]
testing2x=testing2[-inboost2[c(1:boost_size_neg)],]

testing2=testing2x[sample(seq(along=testing2x[,1]), test_size_pos*20),]

testing=rbind(testing1x,testing2)

testingClass1 <- rep("HC", test_size_pos)
testingClass2 <- rep("NEG", test_size_pos*20)
testingClass=c(testingClass1,testingClass2)

boosting=rbind(boost1,boost2)
boostClass1 <- rep("HC", length(boost1[,1]))
boostClass2 <- rep("NEG",length(boost2[,1]))
boostClass=c(boostClass1,boostClass2)


training=as.data.frame(training)
testing=as.data.frame(testing)
boosting=as.data.frame(boosting)
for(n in 1:(length(colnames(training)))){
  training[,n] <- as.numeric(training[,n])
  testing[,n] <- as.numeric(testing[,n])
  boosting[,n] <- as.numeric(boosting[,n])
}
if(length(nzv >0)){
  training=training[,-nzv]
  testing=testing[,-nzv]
  boosting=boosting[,-nzv]
}

trainTransformed <- predict(preProcValues, training)
testTransformed <- predict(preProcValues, testing)
boostTransformed <- predict(preProcValues, boosting)

train=cbind(trainTransformed, Class=trainingClass)
test=cbind(testTransformed, Class=testingClass)
boost=cbind(boostTransformed, Class=boostClass)

train_slim=as.data.frame(train)
test_slim=as.data.frame(test)
boost_slim=as.data.frame(boost)

train_slim$Class=as.factor(train_slim$Class)
test_slim$Class=as.factor(test_slim$Class)
boost_slim$Class=as.factor(boost_slim$Class)

for(n in 1:(length(colnames(test_slim))-1)){
  train_slim[,n] <- as.numeric(train_slim[,n])
  test_slim[,n] <- as.numeric(test_slim[,n])
  boost_slim[,n] <- as.numeric(boost_slim[,n])
}


myControl <- trainControl(method='cv', number=10, returnResamp='none')

predictorNames=opt


outcomeName="Class"
opt2=append(opt, "Class")

gene=ksvm(Class~., data=train_slim[,match(opt2, colnames(train_slim))], 
  kernal="rbfdot",kpar="automatic",C=as.numeric(read.table("data/best_c.txt")[1,1]),
  prob.model=T)
message("SVM model made")

save(gene, file="data/branchpoint_SVM_gene.Rdata")

gene_p=predict(gene, test_slim)
c2=confusionMatrix(gene_p, test_slim$Class)

rm(branchpoint_df_HCN,filteredDescr)


#predict using SVM

newFeat=predict(gene, test_slim, "probabilities")
newFeat=newFeat[,1]
test_slim=cbind(test_slim, newFeat)
message("test new feet made")

predictorNames=c(predictorNames, "newFeat")
newFeat=predict(gene, boost_slim, "probabilities")
newFeat=newFeat[,1]
boost_slim=cbind(boost_slim, newFeat)


save.image(file="data/branchpoint_ALL.Rdata")

boost_N=boost_slim

message("boostn made")

objGBMCplus <- train(boost_N[,predictorNames],
                boost_N$Class,
                 method='gbm',
                 trControl=myControl)


message("gbm made")

predgbmCplus=predict(object=objGBMCplus, test_slim[,predictorNames])

c3=confusionMatrix(predgbmCplus, test_slim$Class)

keepAt=which(testing$seq_pos0A==1)
test_slimA=test_slim[keepAt,]
test_slimrmA=test_slim[-keepAt,]
test_slim_remade=rbind(test_slimA,test_slimrmA)

A_value=boost_slim$seq_pos0A[which(boosting$seq_pos0A==1)[1]]
boost_slimA=boost_slim[boost_slim$seq_pos0A==A_value,]
boost_slimrmA=boost_slim[boost_slim$seq_pos0A!=A_value,]


save.image(file="data/branchpoint_ALL.Rdata")

save(objGBMCplus,c3, file="data/branchpoint_gbms.Rdata")

save(HCN_names,c,dummies, gene,nzv,objGBMCplus,predictorNames,preProcValues,file="data/models/final_model.RData")

