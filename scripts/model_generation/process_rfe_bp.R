files=list.files("data/models")
files=files[grep("bp_2000.16000", files)]
files=files[grep("optdf", files)]


best_rfe=data.frame(file=files, 
                    accuracy=NA,
                    sensitivity=NA,
                    specificity=NA,
                    pos_pred_value=NA,
                    neg_pred_value=NA,
                    o_v=NA,
                    F1=NA)

for(f in files){
  load(paste0("data/models/",f))
  df=df[-2,]
  df$F1=2*((df$pos_pred_value*df$sensitivity)/(df$pos_pred_value+df$sensitivity))
  best_rfe[match(f, best_rfe$file),-1] <- df[which.max(df$F1),c(1,3,4,5,6,10,11)]
}

best_rfe_opt=best_rfe[which.max(best_rfe$F1),c("file")]
write.table(best_rfe_opt, "data/best_rfe_opt.txt", quote=F,row.names=F,col.names=F)

######
library(caret)
files=list.files("data/")
files=files[grep("bp_2000.16000", files)]
files=files[grep("rfe", files)]

load(paste0("data/",files[1]))
allVars=data.frame(variable=unique(rfeProfile1$variables$var), importance=0)

for(f in files){
  load(paste0("data/",f))
  df=as.data.frame(varImp(rfeProfile1))
  x=match(rownames(df), allVars$variable)
  allVars$importance[x] <- allVars$importance[x] + df$Overall
}

load(paste0("data/models/",best_rfe_opt))

allVars$used=0
allVars$used[which(!is.na(match(allVars$variable, opt[1:best_rfe[which.max(best_rfe$F1),c("o_v")]])))] <- 1

pdf('plots/variable_importance.pdf',8,4)
ggplot(allVars, aes(x=reorder(variable, order(importance,decreasing=T)),y=importance, col=factor(used))) + geom_point()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
