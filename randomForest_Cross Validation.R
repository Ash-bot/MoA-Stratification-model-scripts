################################175 genes
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/2436 genes_Loess_Norm_Merge_FC/Constant DEG/Constant unique genes/')
# write.table(genes2,'ML data_174 genes.csv', sep = ';', col.names = T, row.names = T)
dataset= read.csv('ML data_174 genest.csv.', sep = ';', header = T, row.names = 1)
#Dataset with deceptive classes
dataset= dataset[,c(1:175)]
#Dataset with correct classes
dataset= dataset[,c(1:174,176)]
#Dataset with compounds
dataset= dataset[,c(1:174,177)]
# Featurescaling
dataset[-175] = scale(dataset[-175])
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = as.factor(dataset$Class)
dataset$Class = as.numeric(dataset$Class)
# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
library(caTools)
library(e1071)
library(caret)
set.seed(123)

#Class
split = sample.split(dataset$Class, SplitRatio = 0.80)
training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)

# Fitting classifier to the Training set
# install.packages('randomForest')
library(randomForest)

#Create 10 equally size folds
folds <- cut(seq(1,nrow(training_set)),breaks=10,labels=FALSE)
y= data.frame()
#Perform 10 fold cross validation


control= trainControl(method = 'cv', savePredictions = 'all', summaryFunction = multiClassSummary)
ntree=800
mtry= 6
tune= expand.grid(mtry=mtry)
cv= train(training_set[-175],training_set$Class, method = 'rf',trControl = control, ntree=800,
          tuneGrid = tune)
# install.packages('caret')
cv$results$AccuracySD
cv$results$Accuracy
# > cv$results$AccuracySD
# [1] 0.1671147
# > cv$results$Accuracy
# [1] 0.8035173
##############################################Alll genes
#All genes ML
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/')
# all= read.csv('FC_LoessNorm_mergdata.csv', row.names = 1, sep = ';')
# all=t(all)
# write.table(all, 'allgeneFCNorm.csv', sep = ';')
dataset= as.data.frame(read.csv('allgeneFCNorm.csv', row.names = 1, sep = ';'))
#Dataset with deceptive classes
dataset= dataset[,c(1:2464)]
#Dataset with correct classes
dataset= dataset[,c(1:2463,2465)]
#Dataset with compounds
dataset= dataset[,c(1:2463,2466)]
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)
# Featurescaling
dataset[-2464] = scale(dataset[-2464])
library(caTools)
library(e1071)
library(caret)
set.seed(123)
#Compound
split = sample.split(dataset$Compound, SplitRatio = 0.80)
#Class
split = sample.split(dataset$Class, SplitRatio = 0.80)
training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)
library(randomForest)


control= trainControl(method = 'cv', savePredictions = 'all', summaryFunction = multiClassSummary)
ntree=800
mtry= 6
tune= expand.grid(mtry=mtry)
cv= train(training_set[-2464],training_set$Class, method = 'rf',trControl = control, ntree=800,
          tuneGrid = tune)
# install.packages('caret')
cv$results$AccuracySD
cv$results$Accuracy
# > cv$results$AccuracySD
# [1] 0.1734041
# > cv$results$Accuracy
# [1] 0.7792316