####

start_time = Sys.time()


###################################################################
#### Read Tables for Traninig ####

sink("../results/result-1__rf_log.txt", append=FALSE, split=TRUE)

# Current work-space where this code presents in
#wd <- system("pwd",intern=T)
#setwd(wd)

args=commandArgs(trailingOnly=TRUE)
input <- args[1]
#data = read.table("Input-1__scrna_class.txt", header=T)
data = read.table(input, header=T)

### Class : Tumor (column)
data$Tumor = factor(data$Tumor, levels=c(0,1), labels=c(F,T))

head(data) #result
str(data) #result


###################################################################
#### Traning Split and Test Split ####

# library for data sampling
library(caret)

# script for table split
set.seed(1234)
in_train = createDataPartition(data$Tumor, p=0.8, list=F)
data_train = data[ in_train,]
data_test = data[-in_train,]

# result
prop.table(table(data_train$Tumor))
prop.table(table(data_test$Tumor))


###################################################################
#### Decision Tree (1.Learning) ####

# library for data sampling
library(rpart); library(rpart.plot)

# script for learning decision tree
DT_model = rpart(Tumor~. ,data= data_train)

# result
#png(filename="tto_Final__RF-Decision-Tree.surf.png", width=480, height=480)
#rpart.plot(DT_model, type=4, extra=1)
#dev.off()


###################################################################
#### Decision Tree (2.Testing) ####

DT_predict = predict(DT_model, data_test, type="class")
confusionMatrix(DT_predict, data_test$Tumor, positive="TRUE")
# predictive value VS practical value
# True Positive / True Negative / False Positive / False Negative
# Accuracy / Error rate / Sensitivity / Specificity


###################################################################
#### Decision Tree (3.ROC curve) ####

# ROC curve / AUC (Area under the ROC curve)
library(ROCR)
DT_predict_prob = predict(DT_model, data_test, type="prob")
DT_ROC_predict = prediction(predictions = DT_predict_prob[,"TRUE"], labels = data_test$Tumor)
DT_ROC_performance = performance(DT_ROC_predict, measure ="tpr", x.measure ="fpr")

DT_AUC = performance(DT_ROC_predict, measure="auc")

str(DT_AUC)
DT_AUC@y.values


###################################################################
#### Random Forest (1.Learning) ####

library(randomForest)

# script for learning random forest
set.seed(1234)
RF_model = randomForest(Tumor ~ ., data=data_train, ntree=100, maxnodes=9)

# result
RF_model


###################################################################
#### Random Forest (2.Testing) ####

# confusion matrix
RF_predict = predict(RF_model, data_test, type="class")
confusionMatrix(RF_predict, data_test$Tumor, positive="TRUE")


###################################################################
#### Random Forest (3.ROC curve) ####

# ROC curve
RF_predict_prob = predict(RF_model, data_test, type="prob")
RF_ROC_predict = prediction(predictions = RF_predict_prob[,"TRUE"], labels = data_test$Tumor)
RF_ROC_performance = performance(RF_ROC_predict, measure = "tpr", x.measure="fpr")

#png(filename="tto_Final__RF-ROC-Performance.surf.png", width=480, height=480)
#plot(DT_ROC_performance, main="ROC curve", col="blue")
#plot(RF_ROC_performance, col="red", add=T)
#abline(a=0, b=1, lwd=2, lty=2)

RF_AUC = performance(RF_ROC_predict, measure = "auc")
RF_AUC@y.values
value <- round(as.numeric(RF_AUC@y.values),3)
label_auc <- paste0('AUC = ',as.character(value))
#legend(0.75,0.1,label_auc,bty="n")

#dev.off()


###################################################################
#### Random Forest (4.Variable Importance) ####

# script for variable importance
RF_model = randomForest(Tumor ~ ., data=data, ntree=100, importance=TRUE)
RF_importance = importance(RF_model)

# result
RF_importance


###################################################################
#### Variable Importance (1.Sorting) ####

# script for variable importance
RF_importance_order = order(RF_importance[,"MeanDecreaseAccuracy"], decreasing=T)
RF_importance[RF_importance_order,]

# writing file of variable importance
write.table(RF_importance[RF_importance_order,], '../results/result-1__rf_var-imp.txt')

#png(filename="tto_Final__RF-Variable-Importance.surf.png", width=960, height=960)
#varImpPlot(RF_model)
#dev.off()


####

end_time = Sys.time()
print(end_time - start_time)

