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
# dataset$Class = as.factor(dataset$Class)
dataset$Class = as.numeric(dataset$Class)
# Splitting the dataset into the Training set and Test set
# install.packages('caTools')

# install.packages("h2o")
h2o.shutdown()
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.8), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 80%
test <- h2o.assign(splits[[2]], "test.hex") # 20%
response <- "Class"
# response <- "Compound"
predictors <- setdiff(names(dataset), response)
print(response)

###############################LR
lgr= h2o.glm(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234,
             family = 'multinomial')
print(lgr)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1562754524862_4 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.04629 )                       2625
# number_of_active_predictors number_of_iterations training_frame
# 1                         172                    6      train.hex
# 
# Coefficients: glm multinomial coefficients
# names coefs_class_0 coefs_class_1 coefs_class_2 coefs_class_3 coefs_class_4 coefs_class_5
# 1 Intercept     -1.535511     -2.334058     -2.134619     -2.443245     -2.764356     -3.709164
# 2 PF11_0527      0.083070      0.000000      0.411881      0.000000      0.000000      0.000000
# 3 PF13_0042      0.000000      0.000000      0.000000      0.000000      0.000000      0.000000
# 4  PFI0965w      0.000000      0.000000      0.000000      0.000000      0.000000      0.000000
# 5  PFA0180w      0.000000      0.000000      0.000000      0.000000      0.000000      0.000000
# coefs_class_6 coefs_class_7 coefs_class_8 coefs_class_9 coefs_class_10 coefs_class_11 coefs_class_12
# 1     -3.002509     -2.742743     -2.348712     -2.598422      -2.678757      -2.233380      -3.172991
# 2      0.000000     -0.060693      0.000000      0.000000       0.000000       0.000000       0.000000
# 3      0.000000      0.000000      0.000000      0.071606       0.000000       0.000000       0.000000
# 4      0.000000      0.000000      0.000000      0.000000      -0.238453       0.000000       0.000000
# 5      0.000000      0.000000      0.000000      0.000000       0.000000       0.000000       0.208463
# coefs_class_13 coefs_class_14 std_coefs_class_0 std_coefs_class_1 std_coefs_class_2 std_coefs_class_3
# 1      -2.160189      -2.417505         -1.530013         -2.388308         -2.203941         -2.415009
# 2       0.000000       0.000000          0.083141          0.000000          0.412235          0.000000
# 3       0.000000       0.000000          0.000000          0.000000          0.000000          0.000000
# 4       0.000000       0.000000          0.000000          0.000000          0.000000          0.000000
# 5       0.000000       0.000000          0.000000          0.000000          0.000000          0.000000
# std_coefs_class_4 std_coefs_class_5 std_coefs_class_6 std_coefs_class_7 std_coefs_class_8
# 1         -2.705010         -3.611473         -3.087589         -2.752050         -2.520614
# 2          0.000000          0.000000          0.000000         -0.060745          0.000000
# 3          0.000000          0.000000          0.000000          0.000000          0.000000
# 4          0.000000          0.000000          0.000000          0.000000          0.000000
# 5          0.000000          0.000000          0.000000          0.000000          0.000000
# std_coefs_class_9 std_coefs_class_10 std_coefs_class_11 std_coefs_class_12 std_coefs_class_13
# 1         -2.692790          -2.744672          -2.294233          -3.101243          -2.177072
# 2          0.000000           0.000000           0.000000           0.000000           0.000000
# 3          0.069085           0.000000           0.000000           0.000000           0.000000
# 4          0.000000          -0.191508           0.000000           0.000000           0.000000
# 5          0.000000           0.000000           0.000000           0.204200           0.000000
# std_coefs_class_14
# 1          -2.332533
# 2           0.000000
# 3           0.000000
# 4           0.000000
# 5           0.000000
# 
# ---
#   names coefs_class_0 coefs_class_1 coefs_class_2 coefs_class_3 coefs_class_4 coefs_class_5
# 170   PF10_0374      0.000000      0.000000      0.000000      0.000000      0.000000      0.000000
# 171    PFF0595c      0.000000      0.000000      0.000000      0.000000      0.000000      0.000000
# 172    PFF0200c      0.000000      0.155639      0.000000     -0.128463      0.000000      0.000000
# 173 MAL13P1.137      0.000000      0.016988      0.000000     -0.198539      0.000000      0.000000
# 174   MAL7P1.65      0.000000      0.000000      0.000000      0.000000      0.000000      0.000000
# 175   PF13_0090      0.000000      0.000000      0.000000      0.000000      0.000000      0.000000
# coefs_class_6 coefs_class_7 coefs_class_8 coefs_class_9 coefs_class_10 coefs_class_11 coefs_class_12
# 170     -0.658165      0.000000      0.000000      0.000000       0.000000       0.000000       0.000000
# 171      0.000000      0.000000      0.000000      0.000000       0.000000       0.000000       0.000000
# 172      0.221492      0.000000      0.000000      0.000000       0.000000       0.000000       0.000000
# 173      0.000000      0.000000      0.000000      0.000000       0.000000       0.000000       0.000000
# 174      0.000000      0.000000      0.000000      0.000000       0.000000       0.000000       0.000000
# 175      0.000000      0.000000      0.000000      0.000000       0.000000       0.000000       0.000000
# coefs_class_13 coefs_class_14 std_coefs_class_0 std_coefs_class_1 std_coefs_class_2
# 170       0.000000       0.000000          0.000000          0.000000          0.000000
# 171       0.000000       0.000000          0.000000          0.000000          0.000000
# 172       0.000000      -0.014208          0.000000          0.169112          0.000000
# 173       0.000000      -0.006794          0.000000          0.018780          0.000000
# 174       0.000000       0.000000          0.000000          0.000000          0.000000
# 175       0.000000       0.000000          0.000000          0.000000          0.000000
# std_coefs_class_3 std_coefs_class_4 std_coefs_class_5 std_coefs_class_6 std_coefs_class_7
# 170          0.000000          0.000000          0.000000         -0.659692          0.000000
# 171          0.000000          0.000000          0.000000          0.000000          0.000000
# 172         -0.139583          0.000000          0.000000          0.240666          0.000000
# 173         -0.219489          0.000000          0.000000          0.000000          0.000000
# 174          0.000000          0.000000          0.000000          0.000000          0.000000
# 175          0.000000          0.000000          0.000000          0.000000          0.000000
# std_coefs_class_8 std_coefs_class_9 std_coefs_class_10 std_coefs_class_11 std_coefs_class_12
# 170          0.000000          0.000000           0.000000           0.000000           0.000000
# 171          0.000000          0.000000           0.000000           0.000000           0.000000
# 172          0.000000          0.000000           0.000000           0.000000           0.000000
# 173          0.000000          0.000000           0.000000           0.000000           0.000000
# 174          0.000000          0.000000           0.000000           0.000000           0.000000
# 175          0.000000          0.000000           0.000000           0.000000           0.000000
# std_coefs_class_13 std_coefs_class_14
# 170           0.000000           0.000000
# 171           0.000000           0.000000
# 172           0.000000          -0.015438
# 173           0.000000          -0.007511
# 174           0.000000           0.000000
# 175           0.000000           0.000000
# 
# H2OMultinomialMetrics: glm
# ** Reported on training data. **
#   
#   Training Set Metrics: 
#   =====================
#   
#   Extract training frame with `h2o.getFrame("train.hex")`
# MSE: (Extract with `h2o.mse`) 0.1993412
# RMSE: (Extract with `h2o.rmse`) 0.4464764
# Logloss: (Extract with `h2o.logloss`) 0.5737886
# Mean Per-Class Error: 0.05
# Null Deviance: (Extract with `h2o.nulldeviance`) 416.9504
# Residual Deviance: (Extract with `h2o.residual_deviance`) 90.65859
# R^2: (Extract with `h2o.r2`) 0.9908527
# AIC: (Extract with `h2o.aic`) NaN
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>,train = TRUE)`)
# =========================================================================
#   Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
# 1 10 11 12 13 14 15 2 3 4 5 6 7 8 9  Error     Rate
# 1      12  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 = 0 / 12
# 10      0  6  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 6
# 11      0  0  5  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 5
# 12      0  0  0  3  0  0  0 0 0 0 0 0 0 0 1 0.2500 =  1 / 4
# 13      0  0  0  0  4  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 14      0  0  0  0  0  3  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 3
# 15      0  0  0  0  0  0  4 0 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 2       0  0  0  0  0  0  0 4 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 3       0  0  0  0  0  0  0 0 8 0 0 0 0 0 0 0.0000 =  0 / 8
# 4       0  0  0  0  0  0  0 0 0 5 0 0 0 0 0 0.0000 =  0 / 5
# 5       0  0  0  0  0  0  0 0 0 0 5 0 0 0 0 0.0000 =  0 / 5
# 6       0  0  0  0  0  0  0 0 0 0 0 4 0 0 0 0.0000 =  0 / 4
# 7       0  0  0  0  0  0  0 0 0 0 0 0 5 0 0 0.0000 =  0 / 5
# 8       2  0  0  0  0  0  0 0 0 0 0 0 0 2 0 0.5000 =  2 / 4
# 9       0  0  0  0  0  0  0 0 0 0 0 0 0 0 6 0.0000 =  0 / 6
# Totals 14  6  5  3  4  3  4 4 8 5 5 4 5 2 7 0.0380 = 3 / 79
# 
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,train = TRUE)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.962025
# 2   2  1.000000
# 3   3  1.000000
# 4   4  1.000000
# 5   5  1.000000
# 6   6  1.000000
# 7   7  1.000000
# 8   8  1.000000
# 9   9  1.000000
# 10 10  1.000000
# 
# 
# 
# H2OMultinomialMetrics: glm
# ** Reported on cross-validation data. **
#   ** 10-fold cross-validation on training data (Metrics computed for combined holdout predictions) **
#   
#   Cross-Validation Set Metrics: 
#   =====================
#   
#   Extract cross-validation frame with `h2o.getFrame("train.hex")`
# MSE: (Extract with `h2o.mse`) 0.4314468
# RMSE: (Extract with `h2o.rmse`) 0.6568461
# Logloss: (Extract with `h2o.logloss`) 1.232732
# Mean Per-Class Error: 0.3127778
# Null Deviance: (Extract with `h2o.nulldeviance`) 449.7518
# Residual Deviance: (Extract with `h2o.residual_deviance`) 194.7716
# R^2: (Extract with `h2o.r2`) 0.9802019
# AIC: (Extract with `h2o.aic`) NaN
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,xval = TRUE)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.734177
# 2   2  0.797468
# 3   3  0.835443
# 4   4  0.886076
# 5   5  0.924051
# 6   6  0.936709
# 7   7  0.949367
# 8   8  0.974684
# 9   9  0.987342
# 10 10  0.987342
# 
# 
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid cv_5_valid
# accuracy                0.70905983   0.07955801  0.8333333        0.8 0.84615386  0.6666667  0.7777778
# err                     0.29094017   0.07955801 0.16666667        0.2 0.15384616 0.33333334 0.22222222
# err_count                      2.1   0.49497473        1.0        2.0        2.0        2.0        2.0
# logloss                  1.2703218   0.20847753 0.95907193  0.8322017  1.2425241  1.2928679  1.2765213
# max_per_class_error           0.95   0.10606602        1.0        0.5        1.0        1.0        1.0
# mean_per_class_accuracy 0.89111114  0.022443345 0.93333334 0.93333334        0.9  0.8666667 0.93333334
# mean_per_class_error    0.10888889  0.022443345 0.06666667 0.06666667        0.1 0.13333334 0.06666667
# mse                     0.43896112   0.06958366 0.33096576 0.30550736 0.45473403 0.45368764 0.42167318
# null_deviance             44.97518    11.931735  32.282913  57.517498  77.582664   33.40214  53.431126
# r2                       0.9779014 0.0027344343  0.9772619  0.9790174 0.97199345  0.9783958 0.97591287
# residual_deviance        19.477158    5.1383867  11.508863  16.644033  32.305626  15.514416  22.977385
# rmse                     0.6583922   0.05234908  0.5752962  0.5527272  0.6743397 0.67356336  0.6493637
# cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                       0.8  0.6666667        0.5        0.6         0.6
# err                            0.2 0.33333334        0.5        0.4         0.4
# err_count                      2.0        2.0        2.0        2.0         4.0
# logloss                  0.8898781  1.2586261  1.6321157  1.6524968   1.6669143
# max_per_class_error            1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.8666667        0.9  0.8666667  0.8666667  0.84444445
# mean_per_class_error    0.13333334        0.1 0.13333334 0.13333334  0.15555556
# mse                     0.31899592  0.4390319 0.46740445 0.61103565   0.5865753
# null_deviance            55.332127   33.93256   20.88248  27.858946    57.52935
# r2                      0.98693174 0.97499186  0.9789339  0.9807608   0.9748143
# residual_deviance        17.797562  15.103514  13.056926  16.524967   33.338287
# rmse                     0.5647972  0.6625948 0.68366987  0.7816877   0.7658821
lgr_perf <- h2o.performance(model = lgr, newdata = test)
print(lgr_perf)
# H2OMultinomialMetrics: glm
# 
# Test Set Metrics: 
#   =====================
#   
#   MSE: (Extract with `h2o.mse`) 0.3348994
# RMSE: (Extract with `h2o.rmse`) 0.578705
# Logloss: (Extract with `h2o.logloss`) 0.9019858
# Mean Per-Class Error: 0.1555556
# Null Deviance: (Extract with `h2o.nulldeviance`) 119.0125
# Residual Deviance: (Extract with `h2o.residual_deviance`) 43.29532
# R^2: (Extract with `h2o.r2`) 0.9806148
# AIC: (Extract with `h2o.aic`) NaN
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>, <data>)`)
# =========================================================================
#   Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
# 1 10 11 12 13 14 15 2 3 4 5 6 7 8 9  Error     Rate
# 1      4  0  0  0  0  0  0 0 0 0 0 2 0 0 0 0.3333 =  2 / 6
# 10     0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 11     0  0  1  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 1
# 12     0  0  0  2  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 2
# 13     0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 14     0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 15     0  1  0  0  0  0  0 0 0 0 0 0 0 0 0 1.0000 =  1 / 1
# 2      0  0  0  0  0  0  0 1 0 0 0 0 0 0 0 0.0000 =  0 / 1
# 3      0  0  0  0  0  0  0 0 6 0 0 0 0 0 0 0.0000 =  0 / 6
# 4      0  0  0  0  0  0  0 0 0 2 0 0 0 0 0 0.0000 =  0 / 2
# 5      0  0  0  0  0  0  0 0 0 0 3 0 0 0 0 0.0000 =  0 / 3
# 6      0  0  0  0  0  0  0 0 0 0 0 1 0 0 0 0.0000 =  0 / 1
# 7      0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 8      1  0  0  0  0  0  0 0 0 0 0 0 0 0 0 1.0000 =  1 / 1
# 9      0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# Totals 5  1  1  2  0  0  0 1 6 2 3 3 0 0 0 0.1667 = 4 / 24
# 
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>, <data>)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.833333
# 2   2  0.916667
# 3   3  1.000000
# 4   4  1.000000
# 5   5  1.000000
# 6   6  1.000000
# 7   7  1.000000
# 8   8  1.000000
# 9   9  1.000000
# 10 10  1.000000

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(175)])[,1]
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
test_accLb4hp= 20/24*100 #83.33%
#save model as MLRb4hpt_H
h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/" )

# Logistic regression has no hyperparamters to tune

GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)
Variables= h2o.varimp(lgr)
# Variables= h2o.varimp(model)

setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/')
# write.csv(Variables,file = 'GLM_L_Model_VI.csv')
# Variables= read.csv('GLM_L_Model_VI.csv')
top174=head(Variables[order(Variables$percentage,decreasing = TRUE), ], 174)
top= top174$variable
# write.csv(top175,file = 'GLM_L_Model_175VI.csv')
dataset2= t(dataset)
dataset2= as.data.frame(dataset2)
names <- rownames(dataset2)
rownames(dataset2) = NULL
data <- cbind(names,dataset2)
# rm(data)
# data= as.data.frame(read.csv('GLM_H_Model_175VI.csv'))
topg= subset(data, data$names %in% top)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound

# write.csv(HmodelVI_genes, file = 'VarImp_GLM_L_175SelectedGenes.csv')

##################Top 150
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/')
variables= as.data.frame(read.csv('GLM_L_Model_VI.csv', sep = ','))
top150=head(variables[order(variables$percentage,decreasing = TRUE), ], 150)
top150s= as.list(top150$names)
top125=head(variables[order(variables$percentage,decreasing = TRUE), ], 125)
top125s= as.list(top125$names)
top100=head(variables[order(variables$percentage,decreasing = TRUE), ], 100)
top100s= as.list(top100$names)
top75=head(variables[order(variables$percentage,decreasing = TRUE), ], 75)
top75s= as.list(top75$names)
top50=head(variables[order(variables$percentage,decreasing = TRUE), ], 50)
top50s= as.list(top50$names)

#All genes ML
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/')
dataset= as.data.frame(read.csv('allgeneFCNorm2.csv', row.names = 1, sep = ';'))
dataset2= t(dataset)
dataset2= as.data.frame(dataset2)
names <- rownames(dataset2)
rownames(dataset2) = NULL
data <- cbind(names,dataset2)
# rm(data)
# data= as.data.frame(read.csv('GLM_H_Model_175VI.csv'))

#Top150
topg= subset(data, data$names %in% top150$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 150/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_L_150SelectedGenes.csv')

#Top100
topg= subset(data, data$names %in% top100$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 100/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_L_100SelectedGenes.csv')

#Top125
topg= subset(data, data$names %in% top125$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 125/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_L_125SelectedGenes.csv')

#Top75
topg= subset(data, data$names %in% top75$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 75/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_L_75SelectedGenes.csv')

#Top50
topg= subset(data, data$names %in% top50$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 50/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_L_50SelectedGenes.csv')
#########################All genes
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
# Featurescaling
dataset[-2464] = scale(dataset[-2464])
# rm(dataset)
#Compound
dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)

# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
h2o.shutdown()
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
#   H2OMultinomialModel: glm
# Model ID:  GLM_model_R_1563658731799_1 
# GLM Model: summary
# family        link                               regularization number_of_predictors_total
# 1 multinomial multinomial Elastic Net (alpha = 0.5, lambda = 0.05229 )                      36960
# number_of_active_predictors number_of_iterations training_frame
# 1                       36945                 1001      train.hex
# 
# Coefficients: glm multinomial coefficients
# names coefs_class_0 coefs_class_1 coefs_class_2 coefs_class_3 coefs_class_4 coefs_class_5
# 1   Intercept     -1.866623     -2.609706     -2.587832     -2.844414     -2.972156     -3.507493
# 2 MAL13P1.100      0.000347      0.000363      0.000049      0.000007     -0.000201     -0.000795
# 3 MAL13P1.105     -0.000326      0.000717      0.000200      0.000104     -0.000460      0.000247
# 4 MAL13P1.112     -0.000001     -0.000380      0.002233      0.000011     -0.000095     -0.000283
# 5 MAL13P1.116     -0.002427      0.000239     -0.000031     -0.000002      0.000037      0.000228
# coefs_class_6 coefs_class_7 coefs_class_8 coefs_class_9 coefs_class_10 coefs_class_11
# 1     -2.876749     -3.033927     -2.336962     -2.850016      -2.837190      -2.744752
# 2     -0.000256      0.000040      0.001769      0.000993       0.000131       0.000100
# 3     -0.000057      0.000742      0.000641     -0.000063      -0.000296      -0.000277
# 4     -0.001263      0.000547      0.000641     -0.000018      -0.000743       0.000294
# 5      0.000046     -0.000322      0.003034     -0.007819       0.001433      -0.000084
# coefs_class_12 coefs_class_13 coefs_class_14 std_coefs_class_0 std_coefs_class_1 std_coefs_class_2
# 1      -3.033344      -2.758252      -2.608685         -1.841842         -2.657229         -2.570232
# 2      -0.003787      -0.000145       0.000973          0.000315          0.000329          0.000044
# 3      -0.000331      -0.000130      -0.000250         -0.000346          0.000761          0.000212
# 4       0.001380       0.000061      -0.000217         -0.000001         -0.000398          0.002344
# 5       0.000246       0.000216       0.000238         -0.002301          0.000227         -0.000029
# std_coefs_class_3 std_coefs_class_4 std_coefs_class_5 std_coefs_class_6 std_coefs_class_7
# 1         -2.817743         -2.918985         -3.399906         -2.922232         -3.012054
# 2          0.000007         -0.000182         -0.000721         -0.000232          0.000036
# 3          0.000110         -0.000489          0.000262         -0.000060          0.000787
# 4          0.000011         -0.000099         -0.000297         -0.001326          0.000574
# 5         -0.000002          0.000036          0.000216          0.000044         -0.000306
# std_coefs_class_8 std_coefs_class_9 std_coefs_class_10 std_coefs_class_11 std_coefs_class_12
# 1         -2.501354         -2.929158          -2.935079          -2.783800          -2.932913
# 2          0.001605          0.000902           0.000119           0.000091          -0.003437
# 3          0.000680         -0.000067          -0.000314          -0.000294          -0.000351
# 4          0.000673         -0.000019          -0.000780           0.000308           0.001448
# 5          0.002875         -0.007411           0.001358          -0.000080           0.000233
# std_coefs_class_13 std_coefs_class_14
# 1          -2.776732          -2.540463
# 2          -0.000132           0.000883
# 3          -0.000138          -0.000265
# 4           0.000064          -0.000228
# 5           0.000205           0.000225
# 
# ---
#   names coefs_class_0 coefs_class_1 coefs_class_2 coefs_class_3 coefs_class_4 coefs_class_5
# 2459 PFL2490c      0.000483     -0.000171     -0.000002     -0.000189      0.000047      0.000974
# 2460 PFL2505c      0.000333     -0.000220     -0.000061      0.000164     -0.000022      0.000003
# 2461 PFL2520w     -0.002675      0.000419     -0.000030      0.000045      0.000826      0.000043
# 2462 PFL2535w      0.000361      0.000011      0.000200      0.000140     -0.000027     -0.004626
# 2463 PFL2540w     -0.000367     -0.000683      0.001245     -0.000090      0.000194      0.000281
# 2464 PFL2550w      0.000090      0.000237     -0.000009     -0.000007     -0.000190     -0.003671
# coefs_class_6 coefs_class_7 coefs_class_8 coefs_class_9 coefs_class_10 coefs_class_11
# 2459     -0.000024     -0.000016     -0.006991      0.000283       0.001008       0.000024
# 2460      0.000056     -0.000182      0.000617     -0.000291       0.001646      -0.000317
# 2461      0.000084     -0.000228     -0.000006     -0.000569       0.002531      -0.000028
# 2462     -0.000004      0.000131      0.005010     -0.000167       0.000133       0.000001
# 2463      0.001262      0.000074     -0.002342      0.000013      -0.000007       0.000241
# 2464     -0.000001      0.000180      0.000514      0.000177       0.002951      -0.000058
# coefs_class_12 coefs_class_13 coefs_class_14 std_coefs_class_0 std_coefs_class_1
# 2459      -0.000063       0.000009      -0.000080          0.000470         -0.000167
# 2460      -0.004192      -0.000005       0.000072          0.000363         -0.000240
# 2461       0.000316      -0.000118       0.000488         -0.002697          0.000422
# 2462      -0.001499      -0.000003      -0.000200          0.000397          0.000012
# 2463       0.002916      -0.000217      -0.000066         -0.000366         -0.000680
# 2464       0.000205      -0.000651      -0.000239          0.000086          0.000227
# std_coefs_class_2 std_coefs_class_3 std_coefs_class_4 std_coefs_class_5 std_coefs_class_6
# 2459         -0.000002         -0.000184          0.000046          0.000949         -0.000023
# 2460         -0.000066          0.000178         -0.000024          0.000003          0.000061
# 2461         -0.000030          0.000045          0.000833          0.000043          0.000085
# 2462          0.000220          0.000154         -0.000030         -0.005091         -0.000005
# 2463          0.001240         -0.000090          0.000193          0.000280          0.001256
# 2464         -0.000008         -0.000007         -0.000182         -0.003515         -0.000001
# std_coefs_class_7 std_coefs_class_8 std_coefs_class_9 std_coefs_class_10 std_coefs_class_11
# 2459         -0.000016         -0.006810          0.000275           0.000982           0.000023
# 2460         -0.000198          0.000672         -0.000317           0.001794          -0.000345
# 2461         -0.000230         -0.000007         -0.000574           0.002552          -0.000028
# 2462          0.000144          0.005513         -0.000183           0.000146           0.000001
# 2463          0.000073         -0.002332          0.000013          -0.000007           0.000240
# 2464          0.000173          0.000492          0.000169           0.002826          -0.000056
# std_coefs_class_12 std_coefs_class_13 std_coefs_class_14
# 2459          -0.000061           0.000009          -0.000078
# 2460          -0.004569          -0.000006           0.000079
# 2461           0.000319          -0.000119           0.000492
# 2462          -0.001650          -0.000003          -0.000220
# 2463           0.002903          -0.000216          -0.000065
# 2464           0.000196          -0.000623          -0.000229
# 
# H2OMultinomialMetrics: glm
# ** Reported on training data. **
#   
#   Training Set Metrics: 
#   =====================
#   
#   Extract training frame with `h2o.getFrame("train.hex")`
# MSE: (Extract with `h2o.mse`) 0.1074175
# RMSE: (Extract with `h2o.rmse`) 0.3277461
# Logloss: (Extract with `h2o.logloss`) 0.3581745
# Mean Per-Class Error: 0
# Null Deviance: (Extract with `h2o.nulldeviance`) 416.9504
# Residual Deviance: (Extract with `h2o.residual_deviance`) 56.59158
# R^2: (Extract with `h2o.r2`) 0.9950709
# AIC: (Extract with `h2o.aic`) NaN
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>,train = TRUE)`)
# =========================================================================
#   Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
# 1 10 11 12 13 14 15 2 3 4 5 6 7 8 9  Error     Rate
# 1      12  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 = 0 / 12
# 10      0  6  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 6
# 11      0  0  5  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 5
# 12      0  0  0  4  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 13      0  0  0  0  4  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 14      0  0  0  0  0  3  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 3
# 15      0  0  0  0  0  0  4 0 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 2       0  0  0  0  0  0  0 4 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 3       0  0  0  0  0  0  0 0 8 0 0 0 0 0 0 0.0000 =  0 / 8
# 4       0  0  0  0  0  0  0 0 0 5 0 0 0 0 0 0.0000 =  0 / 5
# 5       0  0  0  0  0  0  0 0 0 0 5 0 0 0 0 0.0000 =  0 / 5
# 6       0  0  0  0  0  0  0 0 0 0 0 4 0 0 0 0.0000 =  0 / 4
# 7       0  0  0  0  0  0  0 0 0 0 0 0 5 0 0 0.0000 =  0 / 5
# 8       0  0  0  0  0  0  0 0 0 0 0 0 0 4 0 0.0000 =  0 / 4
# 9       0  0  0  0  0  0  0 0 0 0 0 0 0 0 6 0.0000 =  0 / 6
# Totals 12  6  5  4  4  3  4 4 8 5 5 4 5 4 6 0.0000 = 0 / 79
# 
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,train = TRUE)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  1.000000
# 2   2  1.000000
# 3   3  1.000000
# 4   4  1.000000
# 5   5  1.000000
# 6   6  1.000000
# 7   7  1.000000
# 8   8  1.000000
# 9   9  1.000000
# 10 10  1.000000
# 
# 
# 
# H2OMultinomialMetrics: glm
# ** Reported on cross-validation data. **
#   ** 10-fold cross-validation on training data (Metrics computed for combined holdout predictions) **
#   
#   Cross-Validation Set Metrics: 
#   =====================
#   
#   Extract cross-validation frame with `h2o.getFrame("train.hex")`
# MSE: (Extract with `h2o.mse`) 0.3964632
# RMSE: (Extract with `h2o.rmse`) 0.6296532
# Logloss: (Extract with `h2o.logloss`) 1.147209
# Mean Per-Class Error: 0.2733333
# Null Deviance: (Extract with `h2o.nulldeviance`) 449.7518
# Residual Deviance: (Extract with `h2o.residual_deviance`) 181.259
# R^2: (Extract with `h2o.r2`) 0.9818072
# AIC: (Extract with `h2o.aic`) NaN
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,xval = TRUE)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.784810
# 2   2  0.835443
# 3   3  0.873418
# 4   4  0.898734
# 5   5  0.898734
# 6   6  0.911392
# 7   7  0.936709
# 8   8  0.936709
# 9   9  0.962025
# 10 10  0.974684
# 
# 
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                  0.7780342  0.115458995  0.8333333        1.0  0.7692308  0.6666667
# err                       0.2219658  0.115458995 0.16666667        0.0 0.23076923 0.33333334
# err_count                       1.7    0.7778175        1.0        0.0        3.0        2.0
# logloss                   1.1495508   0.26369718  0.9159138  0.6060301  1.4420564   1.207748
# max_per_class_error             0.8   0.28284273        1.0        0.0        1.0        1.0
# mean_per_class_accuracy   0.9066667  0.043204937 0.93333334        1.0  0.8666667  0.8666667
# mean_per_class_error    0.093333334  0.043204937 0.06666667        0.0 0.13333334 0.13333334
# mse                      0.39645723  0.083731264  0.2904171 0.21688579  0.5286338 0.44653195
# null_deviance              44.97518    11.931735  32.282913  57.517498  77.582664   33.40214
# r2                       0.98019767 0.0035661815  0.9800477   0.985104 0.96744204  0.9787366
# residual_deviance         18.125904    6.1423955  10.990966  12.120603  37.493465  14.492976
# rmse                      0.6227696   0.06563249  0.5389036 0.46570998 0.72707206  0.6682305
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.8  0.8333333        1.0        0.4         0.7
# err                     0.22222222        0.2 0.16666667        0.0        0.6         0.3
# err_count                      2.0        2.0        1.0        0.0        3.0         3.0
# logloss                   1.046625  1.0226486 0.97300744 0.86036605  2.0112462   1.4098668
# max_per_class_error            1.0        1.0        1.0        0.0        1.0         1.0
# mean_per_class_accuracy 0.93333334  0.8666667 0.93333334        1.0        0.8   0.8666667
# mean_per_class_error    0.06666667 0.13333334 0.06666667        0.0        0.2  0.13333334
# mse                      0.3488966  0.3297317  0.3624291   0.329112  0.6379902  0.47394395
# null_deviance            53.431126  55.332127   33.93256   20.88248  27.858946    57.52935
# r2                      0.98007005  0.9864919  0.9793553  0.9851668 0.97991216   0.9796503
# residual_deviance         18.83925   20.45297  11.676089  6.8829284  20.112461   28.197336
# rmse                     0.5906747  0.5742227 0.60202086 0.57368284  0.7987429  0.68843585

pred= predict(lgr, test)            
labels <- as.data.frame(test[,c(2464)])[,1]
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
# [1] 1  2  3  4  3  3  3  3  4  4  5  5  5  6  1  1  1  6  1  1  11 9  12 15
# Levels: 1 11 12 15 2 3 4 5 6 9
test_accHb4hp= 20/24*100 #83.33%

#save model as MLRb4hpt_H
GLM_H_model=h2o.saveModel(lgr, path = "/Users/ashle/Desktop/LAbbook/R files/models")
model= h2o.loadModel(GLM_H_model)

Variables= h2o.varimp(lgr)
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/')
# write.csv(Variables,file = 'GLM_H_Model_VI2.csv')

top174=head(Variables[order(Variables$percentage,decreasing = TRUE), ], 174)
top=top174$variable
# write.csv(top174,file = 'GLM_H_Model_175VI2.csv')
# top174= read.csv('GLM_H_Model_175VI2.csv')
dataset2= t(dataset)
dataset2= as.data.frame(dataset2)
names <- rownames(dataset2)
rownames(dataset2) = NULL
data <- cbind(names,dataset2)
# rm(data)
# data= as.data.frame(read.csv('GLM_H_Model_175VI.csv'))
topg= subset(data, data$names %in% top)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound

# write.csv(HmodelVI_genes, file = 'VarImp_GLM_H_175SelectedGenes2.csv')

###################Top 150
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/')
# variables= as.data.frame(read.csv('GLM_H_Model_175VI2.csv', sep = ','))
top150=head(Variables[order(Variables$percentage,decreasing = TRUE), ], 150)
top150s= as.list(top150$names)
top125=head(Variables[order(Variables$percentage,decreasing = TRUE), ], 125)
top125s= as.list(top125$names)
top100=head(Variables[order(Variables$percentage,decreasing = TRUE), ], 100)
top100s= as.list(top100$names)
top75=head(Variables[order(Variables$percentage,decreasing = TRUE), ], 75)
top75s= as.list(top75$names)
top50=head(Variables[order(Variables$percentage,decreasing = TRUE), ], 50)
top50s= as.list(top50$names)

#All genes ML
setwd('/Users/ashle/Desktop/LAbbook/Additional datasets/Old gene ID/1356 genes/')
dataset= as.data.frame(read.csv('allgeneFCNorm.csv', row.names = 1, sep = ';'))
dataset2= t(dataset)
dataset2= as.data.frame(dataset2)
names <- rownames(dataset2)
rownames(dataset2) = NULL
data <- cbind(names,dataset2)
# rm(data)
# data= as.data.frame(read.csv('GLM_H_Model_175VI.csv'))

#Top150
topg= subset(data, data$names %in% top150$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 150/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_H_150SelectedGenes2.csv')

#Top100
topg= subset(data, data$names %in% top100$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 100/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_H_100SelectedGenes2.csv')

#Top125
topg= subset(data, data$names %in% top125$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 125/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_H_125SelectedGenes2.csv')

#Top75
topg= subset(data, data$names %in% top75$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 75/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_H_75SelectedGenes2.csv')

#Top50
topg= subset(data, data$names %in% top50$variable)
topg= data.frame(topg, row.names = 1)
HmodelVI_genes= as.data.frame(t(topg))
HmodelVI_genes$D_class= dataset$Classd
HmodelVI_genes$Class= dataset$Class
HmodelVI_genes$Compound= dataset$Compound
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 50/')
write.csv(HmodelVI_genes, file = 'VarImp_GLM_H_50SelectedGenes2.csv')
