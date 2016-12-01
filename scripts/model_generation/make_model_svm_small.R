##ML testing
options(stringsAsFactors = F)

library(caret)
library(kernlab)

load("data/Preprocessed_data_small.RData")
message("loaded preprocessed values")
args <- commandArgs(trailingOnly = TRUE)

trp=as.numeric(args[1])
tep=as.numeric(args[2])

cvalue=as.numeric(args[3])
repeatn=as.numeric(args[4])

##center and scale and split
multiplier=read.table("data/best_ratio.txt")[1,1]

train_size_pos=trp
train_size_neg=trp*multiplier
test_size_pos=tep
test_size_neg=tep*20

message(paste(train_size_pos,
train_size_neg,
test_size_pos,
test_size_neg))

table(branchpoint_df_HCN$set)
filteredDescr=filteredDescr[,-1]

inTrain1 <- sample(seq(along=branchpoint_df_HCN$set[branchpoint_df_HCN$set=="HC"]), train_size_pos)
inTrain2 <- sample(seq(along=branchpoint_df_HCN$set[branchpoint_df_HCN$set!="HC"]), train_size_neg)

training1=(filteredDescr[which(branchpoint_df_HCN$set=="HC"),])[inTrain1,]
training2=(filteredDescr[branchpoint_df_HCN$set!="HC",])[inTrain2,]
training=rbind(training1,training2)
trainingClass1 <- (branchpoint_df_HCN$set[branchpoint_df_HCN$set=="HC"])[inTrain1]
trainingClass2 <- (branchpoint_df_HCN$set[branchpoint_df_HCN$set!="HC"])[inTrain2]
trainingClass=c(trainingClass1,trainingClass2)

testing1=(filteredDescr[branchpoint_df_HCN$set=="HC",])[-inTrain1,]
inTest1 <- sample(seq(along=testing1[,1]), test_size_pos)
testing1=testing1[inTest1,]

testing2=(filteredDescr[branchpoint_df_HCN$set!="HC",])[-inTrain2,]
inTest2 <- sample(seq(along=testing2[,1]), test_size_neg)
testing2=testing2[inTest2,]

testing=rbind(testing1,testing2)
testingClass1 <- rep("HC", test_size_pos)
testingClass2 <- rep("NEG", test_size_neg)
testingClass=c(testingClass1,testingClass2)

training=as.data.frame(training)
testing=as.data.frame(testing)
for(n in 1:(length(colnames(training)))){
  training[,n] <- as.numeric(training[,n])
  testing[,n] <- as.numeric(testing[,n])
}


nzv <- nearZeroVar(training)
if(length(nzv >0)){
training=training[,-nzv]
testing=testing[,-nzv]
}


preProcValues <- preProcess(training, method = c("center", "scale"))


trainTransformed <- predict(preProcValues, training)
testTransformed <- predict(preProcValues, testing)

train=cbind(trainTransformed, Class=trainingClass)
test=cbind(testTransformed, Class=testingClass)

train_slim=as.data.frame(train)
test_slim=as.data.frame(test)

train_slim$Class=as.factor(train_slim$Class)
test_slim$Class=as.factor(test_slim$Class)

for(n in 1:(length(colnames(test_slim))-1)){
  train_slim[,n] <- as.numeric(train_slim[,n])
  test_slim[,n] <- as.numeric(test_slim[,n])
}


message("Building model")

gene=ksvm(Class~., data=train_slim, kernal="rbfdot",kpar="automatic",C=cvalue,prob.model=T)
n=predict(gene,test_slim)
c=confusionMatrix(n, test_slim$Class)


save.image(file=paste0("data/bp44" ,trp, "-",tep,"_Cval",cvalue,"_rep", repeatn,"_splicing_model_ALL.RData"))
save(c, file=paste0("data/bp44" ,trp, "-",tep,"_Cval",cvalue,"_rep", repeatn,"_splicing_model_C.RData"))
save(gene, c,branchpoint_df_HCN, nzv,preProcValues, dummies, file=paste0("data/bp44",trp, "-",tep,"_Cval",cvalue,"_rep", repeatn,"_splicing_model.RData"))
