args <- commandArgs(trailingOnly = TRUE)

k1=as.character(args[1]) #data/training_size_sel _***.Rdata

load(k1)
df= data.frame(rep, cvalue_n, trp_n, trn_n, tep_n,ten_n, accuracy,accuracy_balanced, sensitivity, specificity, pos_pred_value, neg_pred_value, prevalence, detection_rate,detection_prevalence)
df$F1=2*((df$pos_pred_value*df$sensitivity)/(df$pos_pred_value+df$sensitivity))
mean_df=aggregate(. ~ trn_n, df, mean)

best_ratio=mean_df[which.max(mean_df$F1),c("trn_n")]/mean_df[which.max(mean_df$F1),c("trp_n")]

library(ggplot2)

pdf("plots/model_eval_F1.pdf",8,8)
ggplot(df, aes(x=trn_n, y=F1, group=trn_n)) + geom_boxplot()
dev.off()
write.csv(df, "data/model_eval.csv")
