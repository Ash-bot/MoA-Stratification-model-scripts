#Logistic regression _Minimodel

library(h2o)
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/')
######Genes from total geneset narrowed down to 175 genes
dataset= as.data.frame(read.csv('VarImp_GLM_H_175SelectedGenes2.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
dataset= dataset[,c(1:175)]
#Dataset with correct classes
dataset= dataset[,c(1:174,176)]
#Dataset with compounds
dataset= dataset[,c(1:174,177)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
#   
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563658731799_2 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04051 )                       2975
# number_of_active_predictors number_of_iterations training_frame
# 1                         270                    6      train.hex
# # 
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid cv_5_valid
# accuracy                 0.7407265  0.12648924  0.8333333        1.0 0.84615386  0.6666667  0.7777778
# err                      0.2592735  0.12648924 0.16666667        0.0 0.15384616 0.33333334 0.22222222
# err_count                      2.0    1.264911        1.0        0.0        2.0        2.0        2.0
# logloss                  1.2288436   0.2726161  1.2841812 0.63086635  1.3257129  1.7382741  1.0195022
# max_per_class_error           0.85  0.22638462        1.0        0.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.90588236  0.05407287  0.9411765        1.0  0.9117647 0.88235295  0.9411765
# mean_per_class_error    0.09411765  0.05407287 0.05882353        0.0  0.0882353 0.11764706 0.05882353
# mse                     0.41074842  0.07331649  0.4630925 0.21546975 0.47992012 0.49881816 0.36788377
# null_deviance             48.08551   12.420051  34.881477    59.5407  79.605865  36.811638  53.431126
# r2                       0.9748187 0.010562429 0.97188646 0.98649937 0.96366197  0.9799582  0.9821565
# residual_deviance          19.1728    6.009468  15.410174  12.617327  34.468536   20.85929  18.351038
# rmse                       0.63532 0.059652798 0.68050903 0.46418718  0.6927627  0.7062706 0.60653424
# cv_6_valid  cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                       0.8   0.8333333       0.75        0.6         0.3
# err                            0.2  0.16666667       0.25        0.4         0.7
# err_count                      2.0         1.0        1.0        2.0         7.0
# logloss                  1.0204494   0.9210043  0.8604071  1.8083411   1.6796967
# max_per_class_error            1.0         0.5        1.0        1.0         1.0
# mean_per_class_accuracy 0.88235295   0.9705882  0.9411765 0.88235295   0.7058824
# mean_per_class_error    0.11764706 0.029411765 0.05882353 0.11764706  0.29411766
# mse                      0.3654512  0.33413604 0.30584976  0.5194155  0.55744725
# null_deviance            61.923798    33.93256  24.291977  29.882149    66.55382
# r2                       0.9879627   0.9809669 0.93475205 0.98221177  0.97813076
# residual_deviance        20.408989   11.052052   6.883257  18.083412   33.593933
# rmse                    0.60452557    0.578045 0.55303687 0.72070485  0.74662393

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(175)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  5  13 9  6  10 10 13 11 12 17
# Levels: 1 10 11 12 13 17 2 3 4 5 6 9
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
test_accHb4hp= 19/24*100#79.17

######## H-dataset-150
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 150/')
dataset= as.data.frame(read.csv('VarImp_GLM_H_150SelectedGenes2.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
dataset= dataset[,c(1:151)]
#Dataset with correct classes
dataset= dataset[,c(1:150,152)]
#Dataset with compounds
dataset= dataset[,c(1:150,153)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563658731799_3 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04051 )                       2567
# number_of_active_predictors number_of_iterations training_frame
# 1                         248                    6      train.hex
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid cv_5_valid
# accuracy                 0.7340598 0.109158054  0.6666667        1.0 0.84615386  0.6666667  0.7777778
# err                     0.26594016 0.109158054 0.33333334        0.0 0.15384616 0.33333334 0.22222222
# err_count                      2.0   1.0488088        2.0        0.0        2.0        2.0        2.0
# logloss                  1.2220039    0.275078  1.3121858 0.63916427  1.2087581  1.7950182    1.03164
# max_per_class_error           0.85  0.22638462        1.0        0.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.90588236  0.04342595 0.88235295        1.0  0.9117647 0.88235295  0.9411765
# mean_per_class_error    0.09411765  0.04342595 0.11764706        0.0  0.0882353 0.11764706 0.05882353
# mse                     0.41129297   0.0721348 0.47562188 0.21921268  0.4558588 0.51074314 0.37389055
# null_deviance             48.08551   12.420051  34.881477    59.5407  79.605865  36.811638  53.431126
# r2                       0.9747548 0.010746702 0.97112584  0.9862649  0.9654838  0.9794791 0.98186517
# residual_deviance        18.947453   5.6584244  15.746229  12.783286   31.42771  21.540218   18.56952
# rmse                    0.63593227 0.058664676 0.68965346 0.46820155 0.67517316 0.71466297 0.61146593
# cv_6_valid  cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                       0.8   0.8333333       0.75        0.6         0.4
# err                            0.2  0.16666667       0.25        0.4         0.6
# err_count                      2.0         1.0        1.0        2.0         6.0
# logloss                  1.0446688   0.8564996  0.8779725  1.7870237   1.6671085
# max_per_class_error            1.0         0.5        1.0        1.0         1.0
# mean_per_class_accuracy 0.88235295   0.9705882  0.9411765 0.88235295   0.7647059
# mean_per_class_error    0.11764706 0.029411765 0.05882353 0.11764706  0.23529412
# mse                     0.37707546  0.31950617 0.31286618  0.5213369   0.5468179
# null_deviance            61.923798    33.93256  24.291977  29.882149    66.55382
# r2                       0.9875799  0.98180026  0.9332552   0.982146  0.97854775
# residual_deviance        20.893375   10.277996    7.02378  17.870235    33.34217
# rmse                     0.6140647   0.5652487  0.5593444 0.72203666   0.7394714

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(151)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  5  13 9  6  10 10 13 14 12 17
# Levels: 1 10 12 13 14 17 2 3 4 5 6 9
test_accHb4hp= 20/24*100 #83.33

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

######## H-dataset-125
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 125/')
dataset= as.data.frame(read.csv('VarImp_GLM_H_125SelectedGenes2.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:126)]
#Dataset with correct classes
dataset= dataset[,c(1:125,127)]

#Dataset with compounds
# dataset= dataset[,c(1:125,128)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563658731799_4 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04051 )                       2142
# number_of_active_predictors number_of_iterations training_frame
# 1                         233                    6      train.hex

# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid cv_5_valid
# accuracy                 0.7697008   0.1232171  0.8333333        1.0  0.7692308  0.6666667  0.7777778
# err                     0.23029914   0.1232171 0.16666667        0.0 0.23076923 0.33333334 0.22222222
# err_count                      1.8    1.174734        1.0        0.0        3.0        2.0        2.0
# logloss                  1.2030711  0.26393202  1.1192789  0.6454218  1.2580942   1.785754  1.0710115
# max_per_class_error            0.8  0.28284273        1.0        0.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.91764706 0.046317693  0.9411765        1.0 0.88235295 0.88235295  0.9411765
# mean_per_class_error    0.08235294 0.046317693 0.05882353        0.0 0.11764706 0.11764706 0.05882353
# mse                     0.40893996  0.06787475  0.4347627 0.22399732 0.47233027 0.49616283   0.389867
# null_deviance             48.08551   12.420051  34.881477    59.5407  79.605865  36.811638  53.431126
# r2                      0.97463775 0.011270343  0.9736063  0.9859651  0.9642366 0.98006487  0.9810903
# residual_deviance         18.75747   5.7385974  13.431347  12.908436   32.71045  21.429047  19.278208
# rmse                     0.6346819  0.05531221 0.65936536 0.47328356  0.6872629 0.70438826 0.62439334
# cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                       0.9        1.0       0.75        0.6         0.4
# err                            0.1        0.0       0.25        0.4         0.6
# err_count                      1.0        0.0        1.0        2.0         6.0
# logloss                   1.027299 0.84811336 0.90828013   1.752156   1.6153028
# max_per_class_error            1.0        0.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765        1.0  0.9411765 0.88235295   0.7647059
# mean_per_class_error    0.05882353        0.0 0.05882353 0.11764706  0.23529412
# mse                     0.37317863 0.31805557 0.32491434 0.51585394   0.5402768
# null_deviance            61.923798   33.93256  24.291977  29.882149    66.55382
# r2                       0.9877082  0.9818829  0.9306849  0.9823338  0.97880435
# residual_deviance         20.54598  10.177361   7.266241   17.52156   32.306057
# rmse                     0.6108835  0.5639642 0.57001257  0.7182297  0.73503524

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(126)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  5  13 9  6  10 4  13 14 12 17
# Levels: 1 10 12 13 14 17 2 3 4 5 6 9

test_accHb4hp= 19/24*100 #79.17%

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

######## H-dataset-100
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 100/')
dataset= as.data.frame(read.csv('VarImp_GLM_H_100SelectedGenes2.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:101)]
#Dataset with correct classes
dataset= dataset[,c(1:100,102)]
#Dataset with compounds
# dataset= dataset[,c(1:100,103)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563658731799_5 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04051 )                       1717
# number_of_active_predictors number_of_iterations training_frame
# 1                         208                    6      train.hex
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid cv_5_valid
# accuracy                 0.7697008   0.1232171  0.8333333        1.0  0.7692308  0.6666667  0.7777778
# err                     0.23029914   0.1232171 0.16666667        0.0 0.23076923 0.33333334 0.22222222
# err_count                      1.8    1.174734        1.0        0.0        3.0        2.0        2.0
# logloss                   1.188375  0.24791656  1.0835927   0.718335  1.3162091  1.5977546  1.0361546
# max_per_class_error            0.8  0.28284273        1.0        0.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.91764706 0.046317693  0.9411765        1.0 0.88235295 0.88235295  0.9411765
# mean_per_class_error    0.08235294 0.046317693 0.05882353        0.0 0.11764706 0.11764706 0.05882353
# mse                     0.41038856 0.065967694 0.42478117  0.2572366  0.4943987 0.50535995 0.37754247
# null_deviance             48.08551   12.420051  34.881477    59.5407  79.605865  36.811638  53.431126
# r2                       0.9746904  0.01073582  0.9742123  0.9838824 0.96256566  0.9796954  0.9816881
# residual_deviance        18.717224   5.9374485  13.003112  14.366699  34.221436  19.173056  18.650784
# rmse                     0.6362925 0.052537594  0.6517524   0.507185 0.70313495  0.7108868 0.61444485
# cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                       0.9        1.0       0.75        0.6         0.4
# err                            0.1        0.0       0.25        0.4         0.6
# err_count                      1.0        0.0        1.0        2.0         6.0
# logloss                  1.0194277 0.82768077   0.880429  1.7690333   1.6351333
# max_per_class_error            1.0        0.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765        1.0  0.9411765 0.88235295   0.7647059
# mean_per_class_error    0.05882353        0.0 0.05882353 0.11764706  0.23529412
# mse                      0.3741656 0.31350362 0.31158796 0.50609875   0.5392106
# null_deviance            61.923798   33.93256  24.291977  29.882149    66.55382
# r2                       0.9876757  0.9821422  0.9335279 0.98266786   0.9788462
# residual_deviance        20.388554   9.932169   7.043432  17.690332   32.702667
# rmse                    0.61169076 0.55991393 0.55820066  0.7114062   0.7343096

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(101)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  17 3  3  3  3  3  4  4  5  5  5  6  5  13 9  6  10 4  13 14 12 17
# Levels: 1 10 12 13 14 17 2 3 4 5 6 9

test_accHb4hp= 18/24*100#75

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

######## H-dataset-75
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 75/')
dataset= as.data.frame(read.csv('VarImp_GLM_H_75SelectedGenes2.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:76)]
#Dataset with correct classes
dataset= dataset[,c(1:75,77)]
#Dataset with compounds
# dataset= dataset[,c(1:75,78)]
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563658731799_6 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04051 )                       1292
# number_of_active_predictors number_of_iterations training_frame
# 1                         174                    6      train.hex
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid cv_5_valid
# accuracy                 0.7857265  0.12461675  0.8333333        1.0 0.84615386  0.6666667  0.7777778
# err                      0.2142735  0.12461675 0.16666667        0.0 0.15384616 0.33333334 0.22222222
# err_count                      1.7   1.1423659        1.0        0.0        2.0        2.0        2.0
# logloss                  1.1582031  0.23133151  1.0705888  0.7258893  1.2537767  1.3390794  1.0563093
# max_per_class_error           0.75  0.28504387        1.0        0.0        0.5        1.0        1.0
# mean_per_class_accuracy  0.9235294 0.045753967  0.9411765        1.0  0.9411765 0.88235295  0.9411765
# mean_per_class_error    0.07647059 0.045753967 0.05882353        0.0 0.05882353 0.11764706 0.05882353
# mse                     0.40713435  0.06523331  0.4165941 0.26326346  0.4651536 0.48080158 0.38671565
# null_deviance             48.08551   12.420051  34.881477    59.5407  79.605865  36.811638  53.431126
# r2                       0.9750425 0.010282443  0.9747093  0.9835048 0.96478003  0.9806821 0.98124313
# residual_deviance        18.269777    5.720101  12.847066  14.517787   32.59819  16.068953  19.013567
# rmse                    0.63388574 0.051590756   0.645441 0.51309204 0.68202174 0.69339854  0.6218647
# cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                       0.9  0.8333333        1.0        0.6         0.4
# err                            0.1 0.16666667        0.0        0.4         0.6
# err_count                      1.0        1.0        0.0        2.0         6.0
# logloss                  1.0032411 0.85651684 0.85957193  1.7908553    1.626203
# max_per_class_error            1.0        1.0        0.0        1.0         1.0
# mean_per_class_accuracy  0.9411765  0.9411765        1.0 0.88235295   0.7647059
# mean_per_class_error    0.05882353 0.05882353        0.0 0.11764706  0.23529412
# mse                     0.36328396 0.31924948 0.30378303 0.50941795  0.56308085
# null_deviance            61.923798   33.93256  24.291977  29.882149    66.55382
# r2                       0.9880341  0.9818149 0.93519294  0.9825542  0.97790974
# residual_deviance        20.064823  10.278202  6.8765755  17.908554    32.52406
# rmse                    0.60273045 0.56502163 0.55116516  0.7137352  0.75038713

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(76)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  17 3  3  3  3  17 4  4  5  5  5  6  3  13 9  6  10 4  13 14 12 17
# Levels: 1 10 12 13 14 17 2 3 4 5 6 9

test_accHb4hp= 17/24*100# 70.83

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

######## H-dataset-50
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 50/')
dataset= as.data.frame(read.csv('VarImp_GLM_H_50SelectedGenes2.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:51)]
#Dataset with correct classes
dataset= dataset[,c(1:50,52)]
#Dataset with compounds
# dataset= dataset[,c(1:50,53)]
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563658731799_7 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.03448 )                        867
# number_of_active_predictors number_of_iterations training_frame
# 1                         160                    6      train.hex
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.75970083   0.11970574        1.0        1.0  0.7692308  0.6666667
# err                     0.24029915   0.11970574        0.0        0.0 0.23076923 0.33333334
# err_count                      1.9    1.1597414        0.0        0.0        3.0        2.0
# logloss                  1.2420244   0.22602962  1.2422305 0.90995693  1.2729523  1.4216607
# max_per_class_error           0.75   0.28504387        0.0        0.0        1.0        1.0
# mean_per_class_accuracy  0.9147059  0.045989692        1.0        1.0 0.88235295 0.88235295
# mean_per_class_error    0.08529412  0.045989692        0.0        0.0 0.11764706 0.11764706
# mse                      0.4373666   0.06683725  0.5002976 0.34366366 0.44978312 0.47117344
# null_deviance             48.08551    12.420051  34.881477    59.5407  79.605865  36.811638
# r2                       0.9743526 0.0080227405  0.9696278 0.97846717  0.9659438  0.9810689
# residual_deviance        19.670908    5.7208195  14.906766  18.199139   33.09676  17.059929
# rmse                     0.6572067  0.052182287  0.7073172  0.5862283  0.6706587 0.68642074
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.8  0.8333333       0.75        0.6         0.4
# err                     0.22222222        0.2 0.16666667       0.25        0.4         0.6
# err_count                      2.0        2.0        1.0        1.0        2.0         6.0
# logloss                  0.9788993  1.4490857 0.97066104    0.78544   1.887402   1.5019559
# max_per_class_error            0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765  0.9117647  0.9411765  0.9411765 0.88235295   0.7647059
# mean_per_class_error    0.05882353  0.0882353 0.05882353 0.05882353 0.11764706  0.23529412
# mse                     0.36166394  0.5227846  0.3676874  0.2616432  0.5813502  0.51361895
# null_deviance            53.431126  61.923798   33.93256  24.291977  29.882149    66.55382
# r2                      0.98245823 0.98278046 0.97905576  0.9441828 0.98009074   0.9798502
# residual_deviance        17.620188  28.981714  11.647933    6.28352   18.87402   30.039116
# rmse                      0.601385 0.72303843 0.60637236 0.51151067  0.7624632   0.7166721

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(51)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  8  10 3  6  3  17 4  4  5  5  5  6  3  13 9  1  10 6  13 14 12 17
# Levels: 1 10 12 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 15/24*100 #62.5

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)


###########################################SELECTED GENES MINIMODELS
library(h2o)
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 150/')
######Genes from total geneset narrowed down to 175 genes
dataset= as.data.frame(read.csv('VarImp_GLM_L_150SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:150,152)]
#Dataset with correct classes
dataset= dataset[,c(1:150,151)]
#Dataset with compounds
# dataset= dataset[,c(1:150,153)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563475663677_6 
# GLM Model: summary
# family        link                               regularization
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04629 )
# number_of_predictors_total number_of_active_predictors number_of_iterations
# 1                       2265                         172                    6
# training_frame
# 1      train.hex
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.70905983   0.07955801  0.8333333        0.8 0.84615386  0.6666667
# err                     0.29094017   0.07955801 0.16666667        0.2 0.15384616 0.33333334
# err_count                      2.1   0.49497473        1.0        2.0        2.0        2.0
# logloss                  1.2618718   0.20719978 0.95425695 0.83343905   1.236533  1.2808452
# max_per_class_error           0.95   0.10606602        1.0        0.5        1.0        1.0
# mean_per_class_accuracy 0.89111114  0.022443345 0.93333334 0.93333334        0.9  0.8666667
# mean_per_class_error    0.10888889  0.022443345 0.06666667 0.06666667        0.1 0.13333334
# mse                     0.43660754   0.06918597 0.33043218 0.30598992  0.4525112 0.44969916
# null_deviance             44.97518    11.931735  32.282913  57.517498  77.582664   33.40214
# r2                       0.9780123 0.0027329293 0.97729856 0.97898424  0.9721303  0.9785858
# residual_deviance        19.356289    5.1248193  11.451083  16.668781   32.14986  15.370142
# rmse                     0.6566302   0.05217428  0.5748323  0.5531636  0.6726895  0.6705961
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.8  0.6666667        0.5        0.6         0.6
# err                     0.22222222        0.2 0.33333334        0.5        0.4         0.4
# err_count                      2.0        2.0        2.0        2.0        2.0         4.0
# logloss                  1.2755147  0.8782212  1.2372599  1.6253194  1.6396896   1.6576385
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.93333334  0.8666667        0.9  0.8666667  0.8666667  0.84444445
# mean_per_class_error    0.06666667 0.13333334        0.1 0.13333334 0.13333334  0.15555556
# mse                      0.4215476 0.31437856 0.43351048   0.466715 0.60802436    0.583267
# null_deviance            53.431126  55.332127   33.93256   20.88248  27.858946    57.52935
# r2                       0.9759201  0.9871209  0.9753064   0.978965 0.98085564  0.97495633
# residual_deviance        22.959265  17.564425  14.847118  13.002555  16.396896    33.15277
# rmse                    0.64926696  0.5606947 0.65841514 0.68316543 0.77975917   0.7637192

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(151)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  1  1  1  11 12 12 15
# Levels: 1 11 12 15 2 3 4 5 6 8
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  1  1  6  6  1  1  11 12 12 10
# Levels: 1 10 11 12 2 3 4 5 6
test_accHb4hp= 20/24*100 #83.33

######## L-dataset-125
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 125/')
dataset= as.data.frame(read.csv('VarImp_GLM_L_125SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:126)]
#Dataset with correct classes
dataset= dataset[,c(1:125,127)]
#Dataset with compounds
# dataset= dataset[,c(1:125,128)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563475663677_7 
# GLM Model: summary
# family        link                               regularization
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04629 )
# number_of_predictors_total number_of_active_predictors number_of_iterations
# 1                       1890                         172                    5
# training_frame
# 1      train.hex
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.70905983   0.07955801  0.8333333        0.8 0.84615386  0.6666667
# err                     0.29094017   0.07955801 0.16666667        0.2 0.15384616 0.33333334
# err_count                      2.1   0.49497473        1.0        2.0        2.0        2.0
# logloss                  1.2408366   0.20992847 0.94075936 0.83092856  1.2230116  1.2543725
# max_per_class_error           0.95   0.10606602        1.0        0.5        1.0        1.0
# mean_per_class_accuracy 0.89111114  0.022443345 0.93333334 0.93333334        0.9  0.8666667
# mean_per_class_error    0.10888889  0.022443345 0.06666667 0.06666667        0.1 0.13333334
# mse                      0.4292535   0.07002736 0.32827163 0.30495167 0.44813004  0.4449057
# null_deviance             44.97518    11.931735  32.282913  57.517498  77.582664   33.40214
# r2                       0.9783868 0.0027408837   0.977447  0.9790555  0.9724001   0.978814
# residual_deviance        19.045702     5.136488  11.289112  16.618572  31.798302   15.05247
# rmse                      0.650848  0.053152677 0.57294995  0.5522243 0.66942513  0.6670125
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.8  0.6666667        0.5        0.6         0.6
# err                     0.22222222        0.2 0.33333334        0.5        0.4         0.4
# err_count                      2.0        2.0        2.0        2.0        2.0         4.0
# logloss                    1.27158  0.8428192  1.1534035  1.6250898  1.6215856   1.6448162
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.93333334  0.8666667        0.9  0.8666667  0.8666667  0.84444445
# mean_per_class_error    0.06666667 0.13333334        0.1 0.13333334 0.13333334  0.15555556
# mse                     0.41981354 0.29869035  0.3997058  0.4665502  0.6021517  0.57936436
# null_deviance            53.431126  55.332127   33.93256   20.88248  27.858946    57.52935
# r2                       0.9760191  0.9877636  0.9772319  0.9789724 0.98104054   0.9751239
# residual_deviance         22.88844  16.856384  13.840841  13.000718  16.215857   32.896324
# rmse                     0.6479302  0.5465257  0.6322229  0.6830448 0.77598435   0.7611599

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(126)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  1  1  1  11 12 12 15
# Levels: 1 11 12 15 2 3 4 5 6 8
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  1  1  6  6  1  1  11 12 12 10
# Levels: 1 10 11 12 2 3 4 5 6
test_accHb4hp= 20/24*100 #83.33

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

######## L-dataset-100
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 100/')
dataset= as.data.frame(read.csv('VarImp_GLM_L_100SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:101)]
#Dataset with correct classes
dataset= dataset[,c(1:100,102)]
#Dataset with compounds
# dataset= dataset[,c(1:100,103)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563475663677_8 
# GLM Model: summary
# family        link                               regularization
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04629 )
# number_of_predictors_total number_of_active_predictors number_of_iterations
# 1                       1515                         153                    6
# training_frame
# 1      train.hex
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.70905983   0.07955801  0.8333333        0.8 0.84615386  0.6666667
# err                     0.29094017   0.07955801 0.16666667        0.2 0.15384616 0.33333334
# err_count                      2.1   0.49497473        1.0        2.0        2.0        2.0
# logloss                  1.2174613   0.21081613  0.9398362   0.821252  1.2048882  1.2378756
# max_per_class_error           0.95   0.10606602        1.0        0.5        1.0        1.0
# mean_per_class_accuracy 0.89111114  0.022443345 0.93333334 0.93333334        0.9  0.8666667
# mean_per_class_error    0.10888889  0.022443345 0.06666667 0.06666667        0.1 0.13333334
# mse                      0.4210745    0.0700608 0.32778946 0.30015117 0.44033882  0.4405949
# null_deviance             44.97518    11.931735  32.282913  57.517498  77.582664   33.40214
# r2                       0.9788191 0.0027019642  0.9774801  0.9793852    0.97288  0.9790193
# residual_deviance        18.730867    5.1710663  11.278035   16.42504  31.327095 14.8545065
# rmse                     0.6444781  0.053490795   0.572529  0.5478605  0.6635803 0.66377324
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.8  0.6666667        0.5        0.6         0.6
# err                     0.22222222        0.2 0.33333334        0.5        0.4         0.4
# err_count                      2.0        2.0        2.0        2.0        2.0         4.0
# logloss                   1.268108  0.8397211  1.0237432   1.623902  1.5778258   1.6374618
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.93333334  0.8666667        0.9  0.8666667  0.8666667  0.84444445
# mean_per_class_error    0.06666667 0.13333334        0.1 0.13333334 0.13333334  0.15555556
# mse                      0.4165389 0.29746172 0.35538134 0.46687993  0.5885515  0.57705724
# null_deviance            53.431126  55.332127   33.93256   20.88248  27.858946    57.52935
# r2                       0.9762062 0.98781395  0.9797568 0.97895753  0.9814688  0.97522295
# residual_deviance        22.825945  16.794422  12.284918  12.991216  15.778258   32.749237
# rmse                    0.64539826  0.5454005  0.5961387 0.68328613 0.76717114  0.75964284

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(101)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  1  1  1  11 12 12 15
# Levels: 1 11 12 15 2 3 4 5 6 8
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  1  1  6  6  1  10 11 12 12 10
# Levels: 1 10 11 12 2 3 4 5 6

test_accHb4hp= 19/24*100 #79.17

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

######## L-dataset-75
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 75/')
dataset= as.data.frame(read.csv('VarImp_GLM_L_75SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:76)]
#Dataset with correct classes
dataset= dataset[,c(1:75,77)]
#Dataset with compounds
# dataset= dataset[,c(1:75,78)]
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563475663677_9 
# GLM Model: summary
# family        link                               regularization
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04629 )
# number_of_predictors_total number_of_active_predictors number_of_iterations
# 1                       1140                         128                    6
# training_frame
# 1      train.hex
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.7257265  0.082903095  0.8333333        0.8 0.84615386  0.6666667
# err                      0.2742735  0.082903095 0.16666667        0.2 0.15384616 0.33333334
# err_count                      2.0    0.5477226        1.0        2.0        2.0        2.0
# logloss                  1.1443878   0.18606153 0.93035424 0.80770665  1.1415792  1.1893346
# max_per_class_error           0.95   0.10606602        1.0        0.5        1.0        1.0
# mean_per_class_accuracy 0.89444447  0.024152294 0.93333334 0.93333334        0.9  0.8666667
# mean_per_class_error    0.10555556  0.024152294 0.06666667 0.06666667        0.1 0.13333334
# mse                      0.4059668   0.06938589 0.32236642 0.29457825 0.42242783   0.431166
# null_deviance             44.97518    11.931735  32.282913  57.517498  77.582664   33.40214
# r2                      0.97959185 0.0026294482  0.9778527   0.979768  0.9739831  0.9794683
# residual_deviance        17.900873      5.22654  11.164251  16.154133  29.681059  14.272016
# rmse                    0.63271385   0.05310356  0.5677732 0.54275066  0.6499445  0.6566323
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.8  0.8333333        0.5        0.6         0.6
# err                     0.22222222        0.2 0.16666667        0.5        0.4         0.4
# err_count                      2.0        2.0        1.0        2.0        2.0         4.0
# logloss                  1.2319121 0.84365654 0.95098466  1.1723863  1.5620469    1.613917
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.93333334  0.8666667 0.93333334  0.8666667  0.8666667  0.84444445
# mean_per_class_error    0.06666667 0.13333334 0.06666667 0.13333334 0.13333334  0.15555556
# mse                     0.40574265 0.29642338 0.33103752   0.402681 0.58417183  0.56907296
# null_deviance            53.431126  55.332127   33.93256   20.88248  27.858946    57.52935
# r2                      0.97682285  0.9878565  0.9811434   0.981851 0.98160666   0.9755658
# residual_deviance        22.174417   16.87313  11.411816    9.37909  15.620468    32.27834
# rmse                     0.6369793  0.5444478 0.57535857  0.6345715  0.7643113  0.75436926

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(76)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  1  1  1  11 12 12 15
# Levels: 1 11 12 15 2 3 4 5 6 8
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  1  1  6  6  1  3  11 12 12 10
# Levels: 1 10 11 12 2 3 4 5 6

test_accHb4hp= 19/24*100 #79.17

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

######## H-dataset-50
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 50/')
dataset= as.data.frame(read.csv('VarImp_GLM_L_50SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:51)]
#Dataset with correct classes
dataset= dataset[,c(1:50,52)]
#Dataset with compounds
# dataset= dataset[,c(1:50,53)]
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)

print(response)
ncol(test)
nrow(test)
###No HPT
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial',
             max_iterations = 1000)
print(lgr)
# Model Details:
#   ==============
# H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563475663677_10 
# GLM Model: summary
# family        link                               regularization
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04192 )
# number_of_predictors_total number_of_active_predictors number_of_iterations
# 1                        765                         114                    6
# training_frame
# 1      train.hex
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid  cv_2_valid  cv_3_valid cv_4_valid
# accuracy                 0.7434188  0.09581947  0.8333333         0.9   0.9230769  0.6666667
# err                      0.2565812  0.09581947 0.16666667         0.1  0.07692308 0.33333334
# err_count                      1.8   0.6164414        1.0         1.0         1.0        2.0
# logloss                  1.0737971  0.19540878  0.7824226  0.73422223   1.0366538  1.1564221
# max_per_class_error            0.9  0.14142136        1.0         0.5         0.5        1.0
# mean_per_class_accuracy 0.90444446  0.03126944 0.93333334  0.96666664  0.96666664  0.8666667
# mean_per_class_error    0.09555556  0.03126944 0.06666667 0.033333335 0.033333335 0.13333334
# mse                     0.38554305 0.072830915  0.2741087  0.26383147   0.3918536 0.42733043
# null_deviance             44.97518   11.931735  32.282913   57.517498   77.582664   33.40214
# r2                      0.98070574 0.002709503  0.9811681   0.9818797  0.97586614  0.9796509
# residual_deviance         16.86718   5.2230024   9.389071   14.684445   26.952997  13.877067
# rmse                    0.61551225    0.057826  0.5235539   0.5136453   0.6259821  0.6537052
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.8  0.8333333        0.5        0.6         0.6
# err                     0.22222222        0.2 0.16666667        0.5        0.4         0.4
# err_count                      2.0        2.0        1.0        2.0        2.0         4.0
# logloss                  1.2621362 0.80187696 0.84290653   1.103238  1.4290413   1.5890516
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.93333334  0.8666667 0.93333334  0.8666667  0.8666667  0.84444445
# mean_per_class_error    0.06666667 0.13333334 0.06666667 0.13333334 0.13333334  0.15555556
# mse                     0.41205555 0.27998284 0.30171221 0.39480767  0.5353255   0.5744224
# null_deviance            53.431126  55.332127   33.93256   20.88248  27.858946    57.52935
# r2                      0.97646224    0.98853 0.98281384 0.98220587 0.98314464   0.9753361
# residual_deviance        22.718452  16.037539  10.114878   8.825904  14.290413   31.781033
# rmse                    0.64191556 0.52913404  0.5492833  0.6283372  0.7316594   0.7579066

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(51)])[,1]
predict= as.data.frame(pred)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  1  1  1  11 12 12 15
# Levels: 1 11 12 15 2 3 4 5 6 8
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  1  1  8  8  1  4  11 12 7  10
# Levels: 1 10 11 12 2 3 4 5 6 7 8

test_accHb4hp= 18/24*100 #75

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)