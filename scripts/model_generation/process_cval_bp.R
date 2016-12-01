
file=list.files("data/")
file=file[grep("bp", file)]
file=file[grep("_C.R", file)]

library(ggplot2)
library(stringr)
accuracy=vector()
accuracy_balanced=vector()
sensitivity=vector()
specificity=vector()
pos_pred_value=vector()
neg_pred_value=vector()
prevalence=vector()
detection_rate=vector()
detection_prevalence=vector()

for( i in 1:length(file)){
  load(paste0("data/",file[i]))
  
  accuracy=append(accuracy,as.numeric(c$overall[1]))
  accuracy_balanced=append(accuracy_balanced,as.numeric(c$byClass[8]))
  sensitivity=append(sensitivity,as.numeric(c$byClass[1]))
  specificity=append(specificity,as.numeric(c$byClass[2]))
  pos_pred_value=append(pos_pred_value,as.numeric(c$byClass[3]))
  neg_pred_value=append(neg_pred_value,as.numeric(c$byClass[4]))
  prevalence=append(prevalence,as.numeric(c$byClass[5]))
  detection_rate=append(detection_rate,as.numeric(c$byClass[6]))
  detection_prevalence=append(detection_prevalence,as.numeric(c$byClass[7]))
  
}

df= data.frame(file, accuracy,accuracy_balanced, sensitivity, specificity, pos_pred_value, neg_pred_value, prevalence, detection_rate,detection_prevalence)

df$F1=2*((df$pos_pred_value*df$sensitivity)/(df$pos_pred_value+df$sensitivity))


filen=matrix(unlist(str_split(file, "Cval")), byrow=T, ncol=2)[,2]
df$cvalue=as.numeric(matrix(unlist(str_split(filen, "_rep")), byrow=T, ncol=2)[,1])
df$rep=as.numeric(str_sub(matrix(unlist(str_split(filen, "_rep")), byrow=T, ncol=2)[,2],1,1))

mean_df=aggregate(. ~ cvalue, df[,-1], mean)

best_c=mean_df[which.max(mean_df$F1),c("cvalue")]
write.table(best_c, "data/best_c.txt", quote=F,row.names=F,col.names=F)

pdf("plots/model_eval_cval_F1.pdf",8,8)
ggplot(df, aes(x=cvalue, y=F1, group=cvalue)) +geom_boxplot()+geom_point()
dev.off()
