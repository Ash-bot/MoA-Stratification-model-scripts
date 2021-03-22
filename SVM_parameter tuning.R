
################################################175 genes
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/2436 genes_Loess_Norm_Merge_FC/Constant DEG/Constant unique genes/')
# write.table(genes2,'ML data_174 genes.csv', sep = ';', col.names = T, row.names = T)
dataset= read.csv('ML data_174 genest2.csv', sep = ';', header = T, row.names = 1)
#Dataset with deceptive classes
# dataset= dataset[,c(1:175)]
#Dataset with correct classes
dataset= dataset[,c(1:174,176)]
#Dataset with compounds
# dataset= dataset[,c(1:174,177)]
# Featurescaling
dataset[-175] = scale(dataset[-175])
# rm(dataset)
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = as.factor(dataset$Class)
# dataset$Class = as.numeric(dataset$Class)
# library(ggplot2) #Load package ggplot2 to use function qplot 
# qplot(PF11_0527, PF13_0042, data=dataset, color = Compound)
#Plot data
# plot( dataset[ , c(171:174)], col=dataset[ ,175 ])
# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
library(caTools)
set.seed(123)
#Compound
# split = sample.split(dataset$Compound, SplitRatio = 0.80)
#Class
split = sample.split(dataset$Class, SplitRatio = 0.80)

training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)

library(e1071)

svm_model <- svm(Class ~ ., data=training_set, cross=10)
summary(svm_model)
# Call:
#   svm(formula = Class ~ ., data = training_set, cross = 10)
# 
# 
# Parameters:
#   SVM-Type:  C-classification 
# SVM-Kernel:  radial 
# cost:  1 
# 
# Number of Support Vectors:  82
# 
# ( 14 4 11 6 6 4 4 4 5 5 5 5 3 2 4 )
# 
# 
# Number of Classes:  15 
# 
# Levels: 
#   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 50 
# Single Accuracies:
#   25 75 37.5 37.5 55.55556 87.5 62.5 25 62.5 33.33333 
accuracy= mean(c(25, 75, 37.5, 37.5, 55.55556, 87.5, 62.5, 25, 62.5, 33.33333 ))
std= sd(c(25, 75, 37.5, 37.5, 55.55556, 87.5, 62.5, 25, 62.5, 33.33333))
# > print(accuracy)
# [1] 50.13889
# > print(std)
# [1] 21.64024
pred <- predict(svm_model,training_set[-175])
system.time(pred <- predict(svm_model,training_set[-175]))
# user  system elapsed 
# 0.01    0.00    0.01 
print(pred)

#testset results
pred <- predict(svm_model,test_set[-175])
system.time(predict(svm_model,test_set[-175]))

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 1 0 1 0  0  0  0  1  1  1
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0

# test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    1
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    1
# 20             14    1
# 21             15    1
Accuracy= 14/21*100 #66.67%

svm_tune <- tune(svm, train.x=dataset[-175], train.y=dataset$Class, 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))
print(svm_tune)
# Parameter tuning of 'svm':
#   - sampling method: 10-fold cross validation 
# - best parameters:
#   cost gamma
# 0.1   0.5
# - best performance: 0.8254545 

svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="radial",
                            cost=0.1, gamma=0.5,cross=10)
summary(svm_model_after_tune)
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 12.19512 
# Single Accuracies:
#   12.5 0 12.5 0 11.11111 12.5 12.5 25 12.5 22.22222
accuracy= mean(c(12.5, 0 ,12.5, 0, 11.11111, 12.5, 12.5, 25, 12.5, 22.22222))
std= sd(c(12.5, 0 ,12.5, 0, 11.11111, 12.5, 12.5, 25, 12.5, 22.22222))
# > print(accuracy)
# [1] 12.08333
# > print(std)
# [1] 7.912469
svm_tune <- tune(svm, train.x=dataset[-175], train.y=dataset$Class,
                 kernel= "radial", 
                 ranges=list(cost=10^(-3:10), gamma=c(0,0.1,0.3,.5,1,2,4,8,10)))
print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   cost gamma
# 0.001     0
# 
# - best performance: 0.8272727

svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="radial",
                            cost=0.001, gamma=0, cross=10)
summary(svm_model_after_tune)
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 17.07317 
# Single Accuracies:
#   25 25 25 12.5 11.11111 0 25 12.5 25 11.11111   
accuracy= mean(c(25, 25, 25, 12.5, 11.11111, 0, 25, 12.5, 25, 11.11111))
std= sd(c(25, 25, 25, 12.5, 11.11111, 0, 25, 12.5, 25, 11.11111))
# > print(accuracy)
# [1] 17.22222
# > print(std)
# [1] 8.934106

pred <- predict(svm_model_after_tune,training_set[-175])
#testset results
pred <- predict(svm_model_after_tune,test_set[-175])
system.time(predict(svm_model_after_tune,test_set[-175]))

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
Accuracy= 13/21*100 #61.9%
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 1 0 1 0  0  0  0  1  1  1
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
#    test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    1
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    1
# 20             14    1
# 21             15    1
#####################################################

svm_tune <- tune(svm, train.x=dataset[-175], train.y=dataset$Class,
                 kernel= "polynomial", 
                 ranges=list(cost=10^(-3:10), gamma=c(0,0.1,0.3,.5,1,2,4,8,10),
                             degree= c(1,2,3,4,5,6)))

print(svm_tune)
# Parameter tuning of 'svm':
#   - sampling method: 10-fold cross validation 
# - best parameters:
#   cost gamma degree
#     1   0.1      1
# - best performance: 0.1363636 
svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="polynomial",
                            degree=1, cost=1, gamma=0.1, cross=10)
summary(svm_model_after_tune)
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 77.63889
# Single Accuracies:
#   62.5 62.5 62.5 87.5 100 62.5 87.5 75 87.5 88.88889  
accuracy= mean(c(62.5, 62.5 ,62.5, 87.5, 100 ,62.5 ,87.5, 75 ,87.5, 88.88889 ))
std= sd(c(62.5, 62.5 ,62.5, 87.5, 100 ,62.5 ,87.5, 75 ,87.5, 88.88889))
# print(accuracy)
# [1] 77.63889
# print(std)
# [1] 14.30623
pred <- predict(svm_model_after_tune,training_set[-175])
system.time(predict(svm_model_after_tune,training_set[-175]))
#testset results
pred <- predict(svm_model_after_tune,test_set[-175])
system.time(predict(svm_model_after_tune,test_set[-175]))

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
print(result)
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 0 0  0  0  0  0  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  1  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 1 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  0  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  1  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
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
# 18             12   12
# 19             13    6
# 20             14   14
# 21             15   15
Accuracy = 19/21*100#90.48% too accurate

##################################################### Linear

svm_tune <- tune(svm, train.x=dataset[-175], train.y=dataset$Class,
                 kernel= "linear", 
                 ranges=list(cost=10^(-3:10), gamma=c(0,0.1,0.3,.5,1,2,4,8,10),
                             degree= c(1,2,3,4,5,6)))

print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   cost gamma degree
# 0.1     0      1
# 
# - best performance: 0.1463636
svm_model_after_tune <- svm(Class ~ ., data=training_set,
                            kernel="linear",degree=1, cost=0.1, gamma=0,cross = 10)
summary(svm_model_after_tune)
# # Call:
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 81.70732 
# Single Accuracies:
#   100 75 62.5 62.5 66.66667 100 75 100 87.5 88.88889 

accuracy= mean(c(100, 75, 62.5, 62.5, 66.66667, 100, 75, 100, 87.5, 88.88889 ))
std= sd(c(100, 75, 62.5, 62.5, 66.66667, 100, 75, 100, 87.5, 88.88889))
# > print(accuracy)
# [1] 81.80556
# > print(std)
# [1] 15.43063
pred <- predict(svm_model_after_tune,training_set[-175])
system.time(predict(svm_model_after_tune,training_set[-175]))
#testset results

pred <- predict(svm_model_after_tune,test_set[-175])
system.time(predict(svm_model_after_tune,test_set[-175]))

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 0 0  0  0  0  0  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  1  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 1 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  0  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  1  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
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
# 18             12   12
# 19             13    6
# 20             14   14
# 21             15   15
Accuracy = 19/21*100 #90.48%

##################################################### Sigmoid

svm_tune <- tune(svm, train.x=dataset[-175], train.y=dataset$Class,
                 kernel= "sigmoid", 
                 ranges=list(cost=10^(-3:10), gamma=c(0,0.1,0.3,.5,1,2,4,8,10),
                             coef0= c(0,1,2,3,4,5,6)))


print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   cost gamma coef0
# 10000   0.1     0
# 
# - best performance: 0.5036364
svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="sigmoid",
                            coef0=0, cost=1000, gamma=0.1,cross = 10)
summary(svm_model_after_tune)
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 50 
# Single Accuracies:
#   62.5 62.5 37.5 50 44.44444 37.5 50 50 37.5 66.66667 
accuracy= mean(c(62.5, 62.5, 37.5, 50, 44.44444, 37.5, 50, 50, 37.5, 66.66667))
std= sd(c(62.5, 62.5, 37.5, 50, 44.44444, 37.5, 50, 50, 37.5, 66.66667))
# > print(accuracy)
# [1] 49.86111
# > print(std)
# [1] 11.00353

pred <- predict(svm_model_after_tune,training_set[-175])
system.time(predict(svm_model_after_tune,training_set[-175]))
#testset results
pred <- predict(svm_model_after_tune,test_set[-175])
system.time(predict(svm_model_after_tune,test_set[-175]))

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 0 0 0 0 0 1 1 0  1  1  0  0  0  1
# 2  0 1 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  1  0  0
# 5  0 0 0 0 2 1 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  0  1  0  1  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# > print(result)
# test_set$Class pred
# 1               1    1
# 2               2    2
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    5
# 8               5    5
# 9               6    5
# 10              7    1
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10    1
# 17             11    1
# 18             12   11
# 19             13    4
# 20             14   11
# 21             15    1
Accuracy = 12/21*100 #57.14%

##############################################Alll genes
#All genes ML
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/')
# all= read.csv('FC_LoessNorm_mergdata.csv', row.names = 1, sep = ';')
# all=t(all)
# write.table(all, 'allgeneFCNorm.csv', sep = ';')
dataset= as.data.frame(read.csv('allgeneFCNorm2.csv', row.names = 1, sep = ';'))
#Dataset with deceptive classes
# dataset= dataset[,c(1:2464)]
#Dataset with correct classes
dataset= dataset[,c(1:2463,2465)]
#Dataset with compounds
# dataset= dataset[,c(1:2463,2466)]
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)
# Featurescaling
dataset[-2464] = scale(dataset[-2464])
# install.packages('caTools')
library(caTools)
set.seed(123)
#Compound
split = sample.split(dataset$Compound, SplitRatio = 0.80)
#Class
split = sample.split(dataset$Class, SplitRatio = 0.80)

training_set = subset(dataset, split == TRUE)
test_set = subset(dataset, split == FALSE)
library(e1071)

svm_model <- svm(Class~ ., data=training_set, cross=10)
summary(svm_model)
# Parameters:
#   SVM-Type:  C-classification 
# SVM-Kernel:  radial 
# cost:  1 
# 
# Number of Support Vectors:  82
# 
# ( 14 4 11 6 6 4 4 4 5 5 5 5 3 2 4 )
# 
# 
# Number of Classes:  15 
# 
# Levels: 
#   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 46.34146 
# Single Accuracies:
#   25 62.5 37.5 37.5 66.66667 87.5 50 25 37.5 33.33333 
accuracy= mean(c(25 ,62.5, 37.5 ,37.5, 66.66667, 87.5, 50, 25 ,37.5, 33.33333 ))
std= sd(c(25 ,62.5, 37.5 ,37.5, 66.66667, 87.5, 50, 25 ,37.5, 33.33333 ))
# > print(accuracy)
# [1] 46.25
# > print(std)
# [1] 20.26541
#testset results
pred <- predict(svm_model,test_set[-2464])

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
Accuracy= 14/21*100 #66.67%
# 
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 1 0 1 0  0  0  0  1  1  1
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 1  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
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
# 9               6    1
# 10              7    7
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    9
# 16             10   10
# 17             11   11
# 18             12   11
# 19             13    1
# 20             14    1
# 21             15    1
svm_tune <- tune(svm, train.x=training_set[-2464], train.y=training_set$Class, 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))
print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   cost gamma
# 0.1   0.5
# 
# - best performance: 0.8277778

svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="radial",
                            cost=0.1, gamma=0.5,cross=10)
summary(svm_model_after_tune)
svm_tune <- tune(svm,train.x=training_set[-2464], train.y=training_set$Class,
                 kernel= "radial", 
                 ranges=list(cost=10^(-3:10), gamma=c(0,0.1,0.3,.5,1,2,4,8,10)))
 print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
 # cost gamma
 # 0.001     0
 # 
 # - best performance: 0.8319444

svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="radial",
                            cost=0.001, gamma=0, cross=10)
summary(svm_model_after_tune)
# 10-fold cross-validation on training data:
#   
# Total Accuracy: 17.07317 
# Single Accuracies:
#   25 12.5 0 12.5 33.33333 0 12.5 25 12.5 33.33333 
accuracy= mean(c(25, 12.5, 0, 12.5, 33.33333, 0, 12.5, 25, 12.5, 33.33333))
std= sd(c(25, 12.5, 0, 12.5, 33.33333, 0, 12.5, 25, 12.5, 33.33333))
# > print(accuracy)
# [1] 16.66667
# > print(std)
# [1] 12.10805

pred <- predict(svm_model_after_tune,training_set[-2464])
#testset results
pred <- predict(svm_model_after_tune,test_set[-2464])
cm= table(pred, test_set$Class)
print(cm)
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 3 1 2 1 1 1 1  1  1  1  1  1  1
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# > print(result)
# test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    1
# 4               3    1
# 5               3    1
# 6               4    1
# 7               5    1
# 8               5    1
# 9               6    1
# 10              7    1
# 11              8    1
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    1
# 16             10    1
# 17             11    1
# 18             12    1
# 19             13    1
# 20             14    1
# 21             15    1
result= as.data.frame(test_set$Class)
result$pred = pred
Accuracy= 4/21*100 #19.05%
#####################################################

svm_tune <- tune(svm, train.x=training_set[-2464], train.y=training_set$Class,
                 kernel= "polynomial", 
                 ranges=list(cost=10^(-3:10), gamma=c(0,0.1,0.3,.5,1,2,4,8,10),
                             degree= c(1,2,3,4,5,6)))

print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
# cost gamma degree
# 0.1   0.1      1
# 
# - best performance: 0.2513889 

svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="polynomial",
                            degree=1, cost=0.1, gamma=0.1, cross=10)
summary(svm_model_after_tune)
# 10-fold cross-validation on training data:
#   
# Total Accuracy: 78.04878 
# Single Accuracies:
#   62.5 75 87.5 50 55.55556 100 75 87.5 100 88.88889


accuracy= mean(c(62.5, 75, 87.5, 50, 55.55556, 100, 75, 87.5, 100, 88.88889))
std= sd(c(62.5, 75, 87.5, 50, 55.55556, 100, 75, 87.5, 100, 88.88889))
# > print(accuracy)
# [1] 78.19445
# > print(std)
# [1] 17.67828
#testset results
pred <- predict(svm_model_after_tune,test_set[-2464])

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 1 0  0  0  0  0  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  1  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 1  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
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
# 19             13    6
# 20             14   14
# 21             15   15

Accuracy = 16/21*100#76.19% too accurate

##################################################### Linear

svm_tune <- tune(svm,  train.x=training_set[-2464], train.y=training_set$Class,
                 kernel= "linear", 
                 ranges=list(cost=10^(-3:1), gamma=c(0,0.1,0.3,.5,1),
                             degree= c(1,2,3)))

print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
#   cost gamma degree
# 0.01     0      1
# 
# - best performance: 0.1680556 

svm_model_after_tune <- svm(Class ~ ., data=training_set,
                            kernel="linear",degree=1, cost=0.01, gamma=0,cross = 10)
summary(svm_model_after_tune)
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 81.70732 
# Single Accuracies:
#   62.5 87.5 87.5 87.5 100 62.5 75 100 87.5 66.66667  

accuracy= mean(c(62.5, 87.5, 87.5, 87.5, 100, 62.5, 75, 100, 87.5, 66.66667))
std= sd(c(62.5, 87.5, 87.5, 87.5, 100, 62.5, 75, 100, 87.5, 66.66667))
# > print(accuracy)
# [1] 81.66667
# > print(std)
# [1] 14.19115
#testset results

pred <- predict(svm_model_after_tune,test_set[-2464])

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 0 1 0  0  0  0  0  0  0
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 0 0 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 2 0 0 0 0  0  0  0  0  0  0
# 6  0 0 0 0 0 1 0 0 0  0  0  0  1  0  0
# 7  0 0 0 0 0 0 1 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  1  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 1  0  1  1  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  1  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  1
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
# 19             13    6
# 20             14   14
# 21             15   15
Accuracy = 16/21*100 # 76.19%

##################################################### Sigmoid

svm_tune <- tune(svm, train.x=training_set[-2464], train.y=training_set$Class,
                 kernel= "sigmoid", 
                 ranges=list(cost=10^(-3:1), gamma=c(0,0.1,0.3,.5,1),
                             coef0= c(0,1,2)))


print(svm_tune)
# Parameter tuning of 'svm':
#   
#   - sampling method: 10-fold cross validation 
# 
# - best parameters:
# cost gamma coef0
# 0.1   0.1     2
# 
# - best performance: 0.7097222 
svm_model_after_tune <- svm(Class ~ ., data=training_set, kernel="sigmoid",
                            coef0=2, cost=0.1, gamma=0.1,cross = 10)
summary(svm_model_after_tune)

# 
# 10-fold cross-validation on training data:
#   
#   Total Accuracy: 29.26829 
# Single Accuracies:
#   25 37.5 12.5 25 33.33333 50 12.5 25 25 44.44444 

accuracy= mean(c(25, 37.5, 12.5, 25, 33.33333, 50, 12.5, 25, 25, 44.44444))
std= sd(c(25, 37.5, 12.5, 25, 33.33333, 50, 12.5, 25, 25, 44.44444))

# > print(accuracy)
# [1] 29.02778
# > print(std)
# [1] 12.37852
#testset results
pred <- predict(svm_model_after_tune,test_set[-2464])
system.time(predict(svm_model_after_tune,test_set[-2464]))

cm= table(pred, test_set$Class)
print(cm)
result= as.data.frame(test_set$Class)
result$pred = pred
# pred 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
# 1  4 1 0 0 0 0 1 0 1  1  1  1  1  1  1
# 2  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 3  0 0 3 0 2 1 0 0 0  0  0  0  0  0  0
# 4  0 0 0 1 0 0 0 0 0  0  0  0  0  0  0
# 5  0 0 0 0 0 0 0 1 0  0  0  0  0  0  0
# 6  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 7  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 9  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 10 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 11 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 12 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 13 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 14 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# 15 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0
# > print(result)
# test_set$Class pred
# 1               1    1
# 2               2    1
# 3               3    3
# 4               3    3
# 5               3    3
# 6               4    4
# 7               5    3
# 8               5    3
# 9               6    3
# 10              7    1
# 11              8    5
# 12              1    1
# 13              1    1
# 14              1    1
# 15              9    1
# 16             10    1
# 17             11    1
# 18             12    1
# 19             13    1
# 20             14    1
# 21             15    1
Accuracy = 8/21*100 #38.10%
