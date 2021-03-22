### XGBoost

setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/2436 genes_Loess_Norm_Merge_FC/Constant DEG/Constant unique genes/')
# write.table(genes2,'ML data_174 genes.csv', sep = ';', col.names = T, row.names = T)

dataset= read.csv('ML data_174 genest.csv.', sep = ';', header = T, row.names = 1)

#Dataset with deceptive classes
dataset= dataset[,c(1:175)]
#Dataset with correct classes
dataset= dataset[,c(1:174,176)]
#Dataset with compounds
dataset= dataset[,c(1:174,177)]


# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class)
dataset$Class = as.numeric(dataset$Class)
library(xgboost)  # the main algorithm
library(archdata) # for the sample dataset
library(caret)    # for the confusionmatrix() function (also needs e1071 package)
library(dplyr)    # for some data preperation
library(Ckmeans.1d.dp) # for xgb.ggplot.importance
library(caTools)

set.seed(123)
#Compound
split = sample.split(dataset$Compound, SplitRatio = 0.80)
#Class
split = sample.split(dataset$Class, SplitRatio = 0.80)
# Full data set
data_variables <- as.matrix(dataset[,-175])
data_label <- dataset[,'Compound']
data_label <- dataset[,'Class']
data_matrix <- xgb.DMatrix(data = as.matrix(data_variables), label = data_label)
# split train data and make xgb.DMatrix
train_data   <- data_variables[split,]
train_label  <- data_label[split]
train_matrix <- xgb.DMatrix(data = train_data, label = train_label)
# split test data and make xgb.DMatrix
test_data  <- data_variables[-split,]
test_label <- data_label[-split]
test_matrix <- xgb.DMatrix(data = test_data, label = test_label)


numberOfClasses <- length(unique(train_label))
xgb_params <- list("objective" = "multi:softmax",
                   "eval_metric" = "mlogloss",
                   "num_class" = numberOfClasses+1)
nround    <- 50 # number of XGBoost rounds
cv.nfold  <- 5

# Fit cv.nfold * cv.nround XGB models and save OOF predictions
cv_model <- xgb.cv(params = xgb_params,
                   data = train_matrix, 
                   nrounds = nround,
                   nfold = cv.nfold,
                   verbose = FALSE,
                   prediction = TRUE)
OOF_prediction <- data.frame(cv_model$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = train_label + 1)
head(OOF_prediction)

# confusion matrix
confusionMatrix(factor(OOF_prediction$max_prob),
                factor(OOF_prediction$label),
                mode = "everything")

bst_model <- xgb.train(params = xgb_params,
                       data = train_matrix,
                       nrounds = nround)

# Predict hold-out test set
test_pred <- predict(bst_model, newdata = test_matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test_label+1,
         max_prob = max.col(., "last"))
# confusion matrix of test set
confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")