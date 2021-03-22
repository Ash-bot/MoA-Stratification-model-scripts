
################################175 genes
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/2436 genes_Loess_Norm_Merge_FC/Constant DEG/Constant unique genes/')
# write.table(genes2,'ML data_174 genes.csv', sep = ';', col.names = T, row.names = T)
dataset= read.csv('ML data_174 genest2.csv.', sep = ';', header = T, row.names = 1)
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
library(ggplot2) #Load package ggplot2 to use function qplot 
qplot(PF11_0527, PF13_0042, data=dataset, color = Compound)
#Plot data
plot( dataset[ , c(171:174)], col=dataset[ ,175 ])
# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
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

# Fitting classifier to the Training set
# install.packages('randomForest')
library(randomForest)
# we need x and y argument which is the matrix of features (predictors) 
#and the dependant variable (prediction = y)
#also the # of trees we want
classifier = randomForest( x= training_set[-175], #without the dependant variable
                           y= training_set$Class,
                           ntree = 500) #NB be cautious of the # of trees chosen
#you do not want to overfit your model so that it may fit the training set correctly
#but not be as accurate in predicting the test set or new dataset
pred =predict(classifier, test_set[,-175])
print(pred)
cm = table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
#    test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    6
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    9
# 20             14   14
# 21             15   15
Accuracy = 17/21*100 # 80.95% 
summary(classifier)
############################## Parameter tuning
RF_tune <- tune(randomForest, train.x=training_set[-175], train.y=training_set$Class,
                ranges=list(ntree=c(1,10,100,500,1000,5000)))
print(RF_tune)
# Parameter tuning of 'randomForest':
#   - sampling method: 10-fold cross validation 
# - best parameters:
#   ntree
# 500
# - best performance: 0.2819444 
RF_tune <- tune(randomForest, train.x=training_set[-175], train.y=training_set$Class,
                ranges=list(ntree=c(100,200,300,400,500, 600)))
print(RF_tune)
# Parameter tuning of 'randomForest':
#   - sampling method: 10-fold cross validation 
# - best parameters:
#   ntree
# 400
# - best performance: 0.3763889 

RF_tune <- tune(randomForest,train.x=training_set[-175], train.y=training_set$Class,
                ranges=list(ntree=c(380,400,420,440,460,480,500)))
print(RF_tune)
# Parameter tuning of 'randomForest':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   ntree
# 460
# 
# - best performance: 0.2569444

RF_tune= tuneRF(training_set[-175],training_set$Class, mtryStart= (sqrt(103)), ntreeTry=460, stepFactor=2, improve=0.05,
                trace=TRUE, plot=TRUE, doBest=FALSE)
print(RF_tune)

# mtry = 10.14889  OOB error = 26.83% 
# Searching left ...
# mtry = 6 	OOB error = 29.27% 
# -0.09090909 0.05 
# Searching right ...
# mtry = 20 	OOB error = 32.93% 
# -0.2272727 0.05 
# > print(RF_tune)
# mtry  OOBError
# 6.OOB                 6.00000 0.2926829
# 10.1488915650922.OOB 10.14889 0.2682927
# 20.OOB               20.00000 0.3292683

RF_model_after_tune <-randomForest( x= training_set[-175], #without the dependant variable
                                    y= training_set$Class,
                                    ntree = 460,
                                    mtry = 10.14889,
                                    keep.forest=TRUE,
                                    proximity=TRUE,
                                    importance = T)

Import=as.data.frame(importance(RF_model_after_tune)) 

pred =predict(RF_model_after_tune, test_set[,-175])
print(pred)
cm = table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
#    test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    6
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    9
# 20             14   14
# 21             15   15
Accuracy = 18/21*100 # 90% 
summary(RF_model_after_tune)

#####DO best RF tune
RF_tune= tuneRF(dataset[-175],dataset$Class, mtryStart= (sqrt(103)), ntreeTry=460, stepFactor=2, improve=0.05,
                trace=TRUE, plot=TRUE, doBest=T)
print(RF_tune)

# Call:
#   randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1]) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 6
# 
# OOB estimate of  error rate: 24.27%
# Confusion matrix:
#   1 2  3 4 5 6 7 8 9 10 11 12 13 14 15 class.error
# 1  17 0  1 0 0 0 0 0 0  0  0  0  0  0  0  0.05555556
# 2   2 3  0 0 0 0 0 0 0  0  0  0  0  0  0  0.40000000
# 3   1 0 13 0 0 0 0 0 0  0  0  0  0  0  0  0.07142857
# 4   0 0  1 6 0 0 0 0 0  0  0  0  0  0  0  0.14285714
# 5   0 0  0 0 8 0 0 0 0  0  0  0  0  0  0  0.00000000
# 6   0 0  0 0 0 5 0 0 0  0  0  0  0  0  0  0.00000000
# 7   1 0  0 0 0 0 4 0 0  0  0  0  0  0  0  0.20000000
# 8   5 0  0 0 0 0 0 0 0  0  0  0  0  0  0  1.00000000
# 9   0 0  0 0 0 1 0 0 4  0  0  1  0  0  0  0.33333333
# 10  1 0  0 0 0 0 0 0 0  5  0  0  0  0  0  0.16666667
# 11  1 0  0 0 0 0 0 0 0  0  2  3  0  0  0  0.66666667
# 12  0 0  0 0 0 0 0 0 1  0  3  2  0  0  0  0.66666667
# 13  0 0  0 0 0 1 0 0 0  0  0  0  2  1  0  0.50000000
# 14  0 0  0 0 0 0 0 0 0  0  0  1  0  2  0  0.33333333
# 15  0 0  0 0 0 0 0 0 0  0  0  0  0  0  5  0.00000000

RF_model_after_tune <-randomForest( x= training_set[-175], #without the dependant variable
                                    y= training_set$Class,
                                    ntree = 460,
                                    mtry = 6,
                                    keep.forest=TRUE,
                                    proximity=TRUE,
                                    importance = T)

Import=as.data.frame(importance(RF_model_after_tune)) 

pred =predict(RF_model_after_tune, test_set[,-175])
print(pred)
cm = table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 0 0  0  0  0  1  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 1 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
# > result= as.data.frame(test_set$Class)
# > result$pred = pred
# > print(result)
# test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    6
# 10              7    7
# 11              8    8
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    1
# 20             14   14
# 21             15   15
Accuracy = 19/21*100 # 90.48% 
summary(RF_model_after_tune)

###################################
##############################################Alll genes
#All genes ML
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/')
# all= read.csv('FC_LoessNorm_mergdata.csv', row.names = 1, sep = ';')
# all=t(all)
# write.table(all, 'allgeneFCNorm.csv', sep = ';')
dataset= as.data.frame(read.csv('allgeneFCNorm2.csv', row.names = 1, sep = ';'))
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
# we need x and y argument which is the matrix of features (predictors) 
#and the dependant variable (prediction = y)
#also the # of trees we want
classifier = randomForest( x= training_set[-2464], #without the dependant variable
                           y= training_set$Class,
                           ntree = 500) #NB be cautious of the # of trees chosen
#you do not want to overfit your model so that it may fit the training set correctly
#but not be as accurate in predicting the test set or new dataset
pred =predict(classifier, test_set[,-2464])
print(pred)
cm = table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
# print(result)
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 1 0  0  0  0  1  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 1  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
# > result= as.data.frame(test_set$Class)
# > result$pred = pred
# > print(result)
# test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    6
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9   11
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    1
# 20             14   14
# 21             15   15
Accuracy = 16/21*100 # 76.19% 
summary(classifier)
############################## Parameter tuning
RF_tune <- tune(randomForest, train.x=training_set[-2464], train.y=training_set$Class,
                ranges=list(ntree=c(1,10,100,500,1000)))

print(RF_tune)
# Parameter tuning of 'randomForest':
#  - sampling method: 10-fold cross validation 
# - best parameters:
#   ntree
# 100
# - best performance: 0.3027778 
RF_tune <- tune(randomForest, train.x=training_set[-2464], train.y=training_set$Class,
                ranges=list(ntree=c(500,1000,2000,3000,4000, 5000)))
print(RF_tune)
# Parameter tuning of 'randomForest':
#   
#   - sampling method: 10-fold cross validation 
# - best parameters:
#   ntree
# 4000
# - best performance: 0.2944444 


RF_tune= tuneRF(training_set[-2464],training_set$Class, mtryStart= (sqrt(103)), ntreeTry=4000, stepFactor=2, improve=0.05,
                trace=TRUE, plot=TRUE, doBest=T)
print(RF_tune)

# mtry = 10.14889  OOB error = 24.39% 
# Searching left ...
# mtry = 6 	OOB error = 23.17% 
# 0.05 0.05 
# Searching right ...
# mtry = 20 	OOB error = 26.83% 
# -0.1 0.05 
# > print(RF_tune)
# 
# Call:
#   randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1]) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 6
# 
# OOB estimate of  error rate: 20.73%
# Confusion matrix:
#   1 2  3 4 5 6 7 8 9 10 11 12 13 14 15 class.error
# 1  13 0  0 0 0 0 0 1 0  0  0  0  0  0  0  0.07142857
# 2   1 3  0 0 0 0 0 0 0  0  0  0  0  0  0  0.25000000
# 3   0 0 11 0 0 0 0 0 0  0  0  0  0  0  0  0.00000000
# 4   0 0  1 5 0 0 0 0 0  0  0  0  0  0  0  0.16666667
# 5   0 0  0 0 6 0 0 0 0  0  0  0  0  0  0  0.00000000
# 6   0 0  0 0 0 4 0 0 0  0  0  0  0  0  0  0.00000000
# 7   1 0  0 0 0 0 3 0 0  0  0  0  0  0  0  0.25000000
# 8   4 0  0 0 0 0 0 0 0  0  0  0  0  0  0  1.00000000
# 9   0 0  0 0 0 1 0 0 3  0  1  0  0  0  0  0.40000000
# 10  0 0  0 0 0 0 0 0 0  5  0  0  0  0  0  0.00000000
# 11  0 0  0 0 0 0 0 0 1  0  3  1  0  0  0  0.40000000
# 12  0 0  0 0 0 0 0 0 0  0  2  3  0  0  0  0.40000000
# 13  0 0  0 0 0 0 0 0 0  0  0  0  2  1  0  0.33333333
# 14  1 0  0 0 0 0 0 0 1  0  0  0  0  0  0  1.00000000
# 15  0 0  0 0 0 0 0 0 0  0  0  0  0  0  4  0.00000000

RF_model_after_tune <-randomForest( x= training_set[-2464], #without the dependant variable
                                    y= training_set$Class,
                                    ntree = 4000,
                                    mtry = 6,
                                    keep.forest=TRUE,
                                    proximity=TRUE,
                                    importance = T)

Import=as.data.frame(importance(RF_model_after_tune)) 

pred =predict(RF_model_after_tune, test_set[,-2464])
print(pred)
cm = table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 1 0  0  0  0  0  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  1  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
# test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    6
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    6
# 20             14   14
# 21             15   15
Accuracy = 17/21*100 # 80.95% 
summary(RF_model_after_tune)

# Length Class  Mode     
# call                8  -none- call     
# type                1  -none- character
# predicted          83  factor numeric  
# err.rate        14400  -none- numeric  
# confusion         306  -none- numeric  
# votes            1411  matrix numeric  
# oob.times          83  -none- numeric  
# classes            17  -none- character
# importance       3306  -none- numeric  
# importanceSD     3132  -none- numeric  
# localImportance     0  -none- NULL     
# proximity        6889  -none- numeric  
# ntree               1  -none- numeric  
# mtry                1  -none- numeric  
# forest             14  -none- list     
# y                  83  factor numeric  
# test                0  -none- NULL     
# inbag               0  -none- NULL   

#####DO best RF tune
RF_tune= tuneRF(dataset[-2464],dataset$Class, mtryStart= (sqrt(103)), ntreeTry=160, stepFactor=2, improve=0.05,
                trace=TRUE, plot=TRUE, doBest=T)
print(RF_tune)

# Call:
#   randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1]) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 6
# 
# OOB estimate of  error rate: 18.45%
# Confusion matrix:
#   1 2  3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 class.error
# 1  5 0  0 0 0 0 0 0 0  1  0  0  0  0  0  0  0   0.1666667
# 2  0 4  0 0 0 0 0 0 1  0  0  0  0  0  0  0  0   0.2000000
# 3  0 0 14 0 0 0 0 0 0  0  0  0  0  0  0  0  0   0.0000000
# 4  0 0  1 6 0 0 0 0 0  0  0  0  0  0  0  0  0   0.1428571
# 5  0 0  0 0 8 0 0 0 0  0  0  0  0  0  0  0  0   0.0000000
# 6  0 0  0 0 0 5 0 0 0  0  0  0  0  0  0  0  0   0.0000000
# 7  0 0  0 0 0 0 5 0 0  0  0  0  0  0  0  0  0   0.0000000
# 8  0 0  0 0 0 0 0 1 4  0  0  0  0  0  0  0  0   0.8000000
# 9  0 0  0 0 0 0 0 0 6  0  0  0  0  0  0  0  0   0.0000000
# 10 0 0  1 0 0 0 0 0 0  5  0  0  0  0  0  0  0   0.1666667
# 11 0 0  0 0 0 1 0 0 0  0  4  0  0  1  0  0  0   0.3333333
# 12 0 0  0 0 0 0 0 0 0  0  0  6  0  0  0  0  0   0.0000000
# 13 0 0  0 0 0 0 0 0 1  0  0  0  4  1  0  0  0   0.3333333
# 14 0 0  0 0 0 0 0 0 0  0  0  0  4  2  0  0  0   0.6666667
# 15 0 0  1 0 0 0 0 0 0  0  0  0  0  0  2  1  0   0.5000000
# 16 0 0  0 0 0 0 0 0 0  0  0  0  0  1  0  2  0   0.3333333
# 17 0 0  0 0 0 0 0 0 0  0  0  0  0  0  0  0  5   0.0000000

RF_model_after_tune <-randomForest( x= training_set[-2462], #without the dependant variable
                                    y= training_set$Class,
                                    ntree = 4000,
                                    mtry = 6,
                                    keep.forest=TRUE,
                                    proximity=TRUE,
                                    importance = T)

Import=as.data.frame(importance(RF_model_after_tune)) 

pred =predict(RF_model_after_tune, test_set[,-2464])
print(pred)
cm = table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
# test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    6
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    6
# 20             14   14
# 21             15   15
# > print(cm)
# 
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 1 0  0  0  0  0  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  1  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
Accuracy = 17/21*100 # 80.95% 
summary(RF_model_after_tune)
