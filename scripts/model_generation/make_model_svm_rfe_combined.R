##ML testing

options(stringsAsFactors = F)

library(caret)
library(kernlab)

load("data/Preprocessed_data_small.RData")
message("loaded preprocessed values")
args <- commandArgs(trailingOnly = TRUE)

trp=as.numeric(args[1])
trn=as.numeric(args[2])
tep=5000
ten=tep*20
cvalue=read.table("data/best_c.txt")[1,1]
repeatn=as.numeric(args[3])
multiplier=read.table("data/best_ratio.txt")[1,1]

message(paste(trp, trn,tep,ten))
##center and scale and split

train_size_pos=trp
train_size_neg=trn
test_size_pos=tep
test_size_neg=ten

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


#####################
#rfe
library(randomForest)
ctrl <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10)

rfeProfile1=rfe(x=train_slim[,-length(colnames(train_slim))], y=train_slim$Class, 
                 sizes=c(5,10,15,20,25,30,35,40,45,50), rfeControl=ctrl)

save(rfeProfile1, file=paste0("data/bp_" ,trp, ".", trn, "-",tep,".", ten,"_Cval",cvalue,"_rep", repeatn,"_splicing_model_rfe.RData"))

opt=c(rfeProfile1$optVariables, "Class")

trn2=trn
trp2=trp

trp=2000
trn=trp*multiplier
tep=10000
ten=tep*20

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
gene_optim=ksvm(Class~., data=train_slim[,match(opt, colnames(train_slim))], kernal="rbfdot",kpar="automatic",C=cvalue,prob.model=T)
#proc.time() - ptm
n_o=predict(gene_optim,test_slim[,match(opt, colnames(test_slim))])
c_o=confusionMatrix(n_o, test_slim$Class)


save.image(file=paste0("data/models/bp_",trp2, ".", trn2, "-",tep,".", ten,"_Cval",cvalue,"_rep", repeatn,"_splicing_model_ALL_optim.RData"))
save(c_o,file=paste0("data/models/bp_",trp2, ".", trn2, "-",tep,".", ten,"_Cval",cvalue,"_rep", repeatn,"_splicing_model_c_optim.RData"))

#opt rming
o_file=paste0("rfe",train_size_neg)

accuracy=vector()
accuracy_balanced=vector()
sensitivity=vector()
specificity=vector()
pos_pred_value=vector()
neg_pred_value=vector()
prevalence=vector()
detection_rate=vector()
detection_prevalence=vector()
o_v=vector()

optALL=colnames(train_slim)
gene_optim=ksvm(Class~., data=train_slim, kernal="rbfdot",kpar="automatic",C=cvalue,prob.model=T)
n_o=predict(gene_optim,test_slim)
c=confusionMatrix(n_o, test_slim$Class)

		accuracy=append(accuracy,as.numeric(c$overall[1]))
		accuracy_balanced=append(accuracy_balanced,as.numeric(c$byClass[8]))
		sensitivity=append(sensitivity,as.numeric(c$byClass[1]))
		specificity=append(specificity,as.numeric(c$byClass[2]))
		pos_pred_value=append(pos_pred_value,as.numeric(c$byClass[3]))
		neg_pred_value=append(neg_pred_value,as.numeric(c$byClass[4]))
		prevalence=append(prevalence,as.numeric(c$byClass[5]))
		detection_rate=append(detection_rate,as.numeric(c$byClass[6]))
		detection_prevalence=append(detection_prevalence,as.numeric(c$byClass[7]))
		o_v=append(o_v,length(colnames(train_slim))-1)


for(o in length(opt):15){
sub_opt=as.character(opt[1:o])
sub_opt=append(sub_opt, "Class")
message(paste0("Building model with ",o," variables"))
gene_optim=ksvm(Class~., data=train_slim[,match(sub_opt, colnames(train_slim))], kernal="rbfdot",kpar="automatic",C=cvalue,prob.model=T)
#proc.time() - ptm
n_o=predict(gene_optim,test_slim[,match(sub_opt, colnames(test_slim))])
c=confusionMatrix(n_o, test_slim$Class)

		accuracy=append(accuracy,as.numeric(c$overall[1]))
		accuracy_balanced=append(accuracy_balanced,as.numeric(c$byClass[8]))
		sensitivity=append(sensitivity,as.numeric(c$byClass[1]))
		specificity=append(specificity,as.numeric(c$byClass[2]))
		pos_pred_value=append(pos_pred_value,as.numeric(c$byClass[3]))
		neg_pred_value=append(neg_pred_value,as.numeric(c$byClass[4]))
		prevalence=append(prevalence,as.numeric(c$byClass[5]))
		detection_rate=append(detection_rate,as.numeric(c$byClass[6]))
		detection_prevalence=append(detection_prevalence,as.numeric(c$byClass[7]))
		o_v=append(o_v,o)

}
save.image(file=paste0("data/models/bp_",as.numeric(args[1]), ".", as.numeric(args[2]), "-",tep,".", ten,"_Cval",cvalue,"_rep", repeatn,"_",o_file,"_splicing_model_ALL_opt_elim.RData"))

df=data.frame(accuracy,accuracy_balanced,sensitivity,specificity,pos_pred_value,
	neg_pred_value,prevalence,detection_rate,detection_prevalence,o_v)

save(opt, df, file=paste0("data/models/bp_",as.numeric(args[1]), ".", as.numeric(args[2]), "-",tep,".", ten,"_Cval",cvalue,"_rep", repeatn,"_",o_file,"_splicing_model_ALL_opt_elim_optdf.RData"))



