#variable_class == bps
quick_wilcox <- function(data_frame, variable_test,variable_class="bps"){
  class1 <- 1
  class2 <- 2
  class3 <- 3
  class4 <- "4+"
  
  #1 v 2
  c1 <- which(data_frame[,variable_class] %in% class1)
  c2 <- which(data_frame[,variable_class] %in% class2)
  
  w1 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  #1 v all
  c1 <- which(data_frame[,variable_class] %in% class1)
  c2 <- which(data_frame[,variable_class] %in% c(class2,class3,class4))
  
  w2 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  
  #2 v 3
  c1 <- which(data_frame[,variable_class] %in% class2)
  c2 <- which(data_frame[,variable_class] %in% class3)
  
  w3 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  #3 v 4
  c1 <- which(data_frame[,variable_class] %in% class3)
  c2 <- which(data_frame[,variable_class] %in% class4)
  
  w4 <- wilcox.test(data_frame[c1,variable_test],
                    data_frame[c2,variable_test])
  
  df <- data.frame(variable=variable_test,class1=c(1,1,2,3), class2=c(2,"2,3,4+", 3,"4+"), pval=c(w1$p.value, w2$p.value,w3$p.value,w4$p.value))
  
  return(df)
}
