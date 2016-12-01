##ML testing

options(stringsAsFactors = F)

library(caret)
library(kernlab)

load("data/Preprocessed_data_small.RData")
message("loaded preprocessed values")
args <- commandArgs(trailingOnly = TRUE)

set.seed(as.numeric(args[1]))
cvalue=(as.numeric(args[2]))

filteredDescr=filteredDescr[,-1]

trp=1000
trn_v=c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,20000)
tep=5000
#cvalue=2
#message(paste(trp, trn,tep,ten))
##center and scale and split

rep=vector()
cvalue_n=vector()
trp_n=vector()
trn_n=vector()
tep_n=vector()
ten_n=vector()
accuracy=vector()
accuracy_balanced=vector()
sensitivity=vector()
specificity=vector()
pos_pred_value=vector()
neg_pred_value=vector()
prevalence=vector()
detection_rate=vector()
detection_prevalence=vector()

for(j in 1:length(trn_v)){
	trn=trn_v[j]
	train_size_pos=trp
	train_size_neg=trn

	test_size_pos=tep
	test_size_neg=tep*20
	ten=test_size_neg

	for(r in 1:5){
	  	repeatn=r

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
		ptm <- proc.time()
		gene=ksvm(Class~., data=train_slim, kernal="rbfdot",kpar="automatic",C=cvalue,prob.model=T)
		#proc.time() - ptm
		n=predict(gene,test_slim)
		c=confusionMatrix(n, test_slim$Class)
		proc.time() - ptm
		c

		rep=append(rep,repeatn)
		cvalue_n=append(cvalue_n,cvalue)
		trp_n=append(trp_n,trp)
		trn_n=append(trn_n,trn)
		tep_n=append(tep_n,tep)
		ten_n=append(ten_n,ten)
		accuracy=append(accuracy,as.numeric(c$overall[1]))
		accuracy_balanced=append(accuracy_balanced,as.numeric(c$byClass[8]))
		sensitivity=append(sensitivity,as.numeric(c$byClass[1]))
		specificity=append(specificity,as.numeric(c$byClass[2]))
		pos_pred_value=append(pos_pred_value,as.numeric(c$byClass[3]))
		neg_pred_value=append(neg_pred_value,as.numeric(c$byClass[4]))
		prevalence=append(prevalence,as.numeric(c$byClass[5]))
		detection_rate=append(detection_rate,as.numeric(c$byClass[6]))
		detection_prevalence=append(detection_prevalence,as.numeric(c$byClass[7]))
		message(paste("TN:",trn, "sens:", as.numeric(c$byClass[1]), "ppv:", as.numeric(c$byClass[3])))

	}
}

save.image(file=paste0("data/training_size_sel_Cval",cvalue,"_",as.numeric(args[1]),".RData"))
