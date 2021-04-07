#Random forest _Minimodel
######All genes filtered minimodels VAr IMP
library(h2o)
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/')
######Genes from total geneset narrowed down to 175 genes
dataset= as.data.frame(read.csv('VarImp_RF_H_175SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:174)]
#Dataset with correct classes
dataset= dataset[,c(1:174,176)]
#Dataset with compounds
# dataset= dataset[,c(1:174,177)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)
h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size =  '5G')#
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
rfh= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
                  0.08861
print(rfh)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564049536577_7 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               90989         1         6    2.55647
# min_leaves max_leaves mean_leaves
# 1          2          8     3.83176
# 
# 
# H2OMultinomialMetrics: drf
# ** Reported on training data. **
#   ** Metrics reported on Out-Of-Bag training samples **
#   
#   Training Set Metrics: 
#   =====================
#   
#   Extract training frame with `h2o.getFrame("train.hex")`
# MSE: (Extract with `h2o.mse`) 0.3391148
# RMSE: (Extract with `h2o.rmse`) 0.5823357
# Logloss: (Extract with `h2o.logloss`) 0.8926575
# Mean Per-Class Error: 0.1
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>,train = TRUE)`)
# =========================================================================
#   Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
# 1 10 11 12 13 14 15 16 17 2 3 4 5 6 7 8 9  Error     Rate
# 1      4  0  0  1  0  0  0  0  0 0 0 0 0 0 0 0 0 0.2000 =  1 / 5
# 10     0  4  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 11     0  0  6  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 6
# 12     0  0  0  5  0  0  0  0  0 0 0 0 0 0 0 1 0 0.1667 =  1 / 6
# 13     0  0  0  0  5  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 5
# 14     0  0  0  0  1  3  0  0  0 0 0 0 0 0 0 0 0 0.2500 =  1 / 4
# 15     0  0  0  0  0  0  4  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 16     0  0  0  0  0  0  0  3  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 3
# 17     0  0  0  1  0  0  0  0  3 0 0 0 0 0 0 0 0 0.2500 =  1 / 4
# 2      0  0  0  0  0  0  0  0  0 4 0 0 0 0 0 0 0 0.0000 =  0 / 4
# 3      0  0  0  0  0  0  0  0  0 0 8 0 0 0 0 0 0 0.0000 =  0 / 8
# 4      0  0  0  0  0  0  0  0  0 0 0 5 0 0 0 0 0 0.0000 =  0 / 5
# 5      0  0  0  0  0  0  0  0  0 0 0 0 5 0 0 0 0 0.0000 =  0 / 5
# 6      0  0  0  0  0  0  0  0  0 0 0 0 0 3 0 0 1 0.2500 =  1 / 4
# 7      0  0  0  0  0  0  0  0  0 0 0 0 0 0 5 0 0 0.0000 =  0 / 5
# 8      0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 3 1 0.2500 =  1 / 4
# 9      0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 1 2 0.3333 =  1 / 3
# Totals 4  4  6  7  6  3  4  3  3 4 8 5 5 3 5 5 4 0.0886 = 7 / 79
# 
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,train = TRUE)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.911392
# 2   2  0.974684
# 3   3  0.987342
# 4   4  0.987342
# 5   5  1.000000
# 6   6  1.000000
# 7   7  1.000000
# 8   8  1.000000
# 9   9  1.000000
# 10 10  1.000000
# 
# 
# 
# H2OMultinomialMetrics: drf
# ** Reported on cross-validation data. **
#   ** 10-fold cross-validation on training data (Metrics computed for combined holdout predictions) **
#   
#   Cross-Validation Set Metrics: 
#   =====================
#   
#   Extract cross-validation frame with `h2o.getFrame("train.hex")`
# MSE: (Extract with `h2o.mse`) 0.5562359
# RMSE: (Extract with `h2o.rmse`) 0.7458123
# Logloss: (Extract with `h2o.logloss`) 1.583815
# Mean Per-Class Error: 0.4102941
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>,xval = TRUE)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.607595
# 2   2  0.772152
# 3   3  0.848101
# 4   4  0.873418
# 5   5  0.898734
# 6   6  0.924051
# 7   7  0.949367
# 8   8  0.962025
# 9   9  0.974684
# 10 10  0.974684
# 
# 
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.5944017  0.12725464  0.6666667        0.7 0.53846157        0.5
# err                     0.40559828  0.12725464 0.33333334        0.3 0.46153846        0.5
# err_count                      3.1   1.0700468        2.0        3.0        6.0        3.0
# logloss                  1.5535641  0.23888822  1.3008568   1.368355  1.7594912  1.6163968
# max_per_class_error            0.9  0.14142136        1.0        0.5        1.0        1.0
# mean_per_class_accuracy 0.85294116  0.05094267 0.88235295  0.9117647  0.7352941  0.8235294
# mean_per_class_error    0.14705883  0.05094267 0.11764706  0.0882353  0.2647059  0.1764706
# mse                     0.55123436  0.06499044 0.49022743 0.48702657 0.60105807 0.61048096
# r2                       0.9660986 0.013460995 0.97023916 0.96948457 0.95448977 0.97547174
# rmse                    0.73982626 0.044110615  0.7001625  0.6978729 0.77527934 0.78133285
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                  0.8888889        0.7        0.5       0.75        0.2         0.5
# err                      0.11111111        0.3        0.5       0.25        0.8         0.5
# err_count                       1.0        3.0        3.0        1.0        4.0         5.0
# logloss                   1.1054196  1.4347843  1.5750908  1.2281473    1.81397    2.333129
# max_per_class_error             0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy   0.9705882 0.85294116 0.85294116  0.9411765  0.7647059   0.7941176
# mean_per_class_error    0.029411765 0.14705883 0.14705883 0.05882353 0.23529412  0.20588236
# mse                       0.4353043  0.5342865  0.5813437 0.40610778  0.6767628  0.68974566
# r2                        0.9788864  0.9824016  0.9668855  0.9133637  0.9768232  0.97294056
# rmse                       0.659776 0.73094904   0.762459 0.63726586   0.822656   0.8305093

pred= predict(rfh, test)            
labels <- as.data.frame(test[,c(175)])[,1]
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
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  17 6  9  10 10 13 11 11 17
# Levels: 1 10 11 13 17 2 3 4 5 6 8 9
test_accHb4hp= 20/24*100#83.33
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees     model_ids            logloss
# 1          5      6    500 rfh_model_121  0.992811866687745
# 2         20      6    500  rfh_model_83 0.9929583806062712
# 3         15      6    500  rfh_model_14 0.9929583806062712
# 4         10      6    500 rfh_model_123 0.9929583806062712
# 5          5      6    400 rfh_model_124  0.998265007027591
# 6         10      6    400  rfh_model_56 0.9984771185232193
# 7         15      6    400  rfh_model_84 0.9984771185232193
# 8         20      6    400 rfh_model_111 0.9984771185232193
# 9          5      6    100   rfh_model_9 1.0043912281031102
# 10         5      6    200  rfh_model_73 1.0050411552707113


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_121')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
rf= h2o.randomForest(predictors,response, train, model_id = 'best_rf',
                     nfolds = 10)
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
rf= h2o.randomForest(predictors,response, train,max_depth = 5, mtries = 6, ntrees = 500,
                     nfolds = 10)
print(best_rf)
print(rf)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564049536577_8 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1             500                     8500              964164         1         5    2.90047
# min_leaves max_leaves mean_leaves
# 1          2          9     4.26212

# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid  cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.72303754  0.10557732 0.54545456   0.8333333        0.8        0.5
# err                       0.2769625  0.10557732 0.45454547  0.16666667        0.2        0.5
# err_count                       2.4   1.1489125        5.0         2.0        1.0        4.0
# logloss                   1.5525213  0.12639315  1.6740603   1.4688768  1.5802083   1.697469
# max_per_class_error            0.85  0.22638462        1.0         1.0        1.0        1.0
# mean_per_class_accuracy   0.8892157  0.05104163  0.7647059  0.92156863  0.9411765  0.8235294
# mean_per_class_error    0.110784315  0.05104163 0.23529412 0.078431375 0.05882353  0.1764706
# mse                      0.58570546  0.04254124  0.6281036   0.5430168  0.5979363   0.619545
# r2                       0.96577984  0.01160435 0.95905143   0.9843235 0.97200674 0.98062974
# rmse                      0.7642686 0.028274976  0.7925299  0.73689675  0.7732634  0.7871118
# cv_5_valid  cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.6666667  0.85714287  0.6666667  0.7777778  0.5833333         1.0
# err                     0.33333334  0.14285715 0.33333334 0.22222222 0.41666666         0.0
# err_count                      2.0         1.0        2.0        2.0        5.0         0.0
# logloss                  1.4779177   1.2215717  1.6701679  1.2676749  1.7455761   1.7216908
# max_per_class_error            1.0         0.5        1.0        1.0        1.0         0.0
# mean_per_class_accuracy  0.9117647   0.9705882 0.88235295 0.88235295  0.7941176         1.0
# mean_per_class_error     0.0882353 0.029411765 0.11764706 0.11764706 0.20588236         0.0
# mse                      0.5638655   0.4704292 0.64351225 0.50669146 0.61371875   0.6702355
# r2                      0.94937867   0.9836749  0.9706382  0.9771481  0.9437815  0.93716544
# rmse                    0.75090975  0.68587846 0.80219215  0.7118226   0.783402   0.8186791
prob=predict(rf, test)
labels <- as.data.frame(test[,c(175)])[,1]
predict= as.data.frame(prob)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  6  8  9  10 3  13 11 13 17
# Levels: 1 10 11 13 17 2 3 4 5 6 8 9
nrow(test)
test_accL= 19/24*100#79.17
######## H-dataset-150
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 150/')
dataset= as.data.frame(read.csv('VarImp_RF_H_150SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:151)]
#Dataset with correct classes
dataset= dataset[,c(1:150,152)]
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
RF= h2o.randomForest(x = predictors, y =response,
                     training_frame = train,
                     nfolds = 10,
                     seed= 1234)
print(RF)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564049536577_9 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               90591         1         5    2.54706
# min_leaves max_leaves mean_leaves
# 1          2          7     3.79176
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.67311966  0.121314876  0.6666667        0.7  0.6923077  0.6666667
# err                     0.32688034  0.121314876 0.33333334        0.3 0.30769232 0.33333334
# err_count                      2.4    0.7211103        2.0        3.0        4.0        2.0
# logloss                  1.3917426   0.23774625  1.1647084  1.2371113  1.7447755  1.3963181
# max_per_class_error           0.95   0.10606602        1.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.88235295   0.03834825 0.88235295 0.88235295  0.8235294 0.88235295
# mean_per_class_error    0.11764706   0.03834825 0.11764706 0.11764706  0.1764706 0.11764706
# mse                      0.5066861  0.069524385  0.4455682 0.47490963 0.58943546   0.533166
# r2                      0.96891844 0.0120361885 0.97295034 0.97024375 0.95536983 0.97857815
# rmse                    0.70829856  0.049996182 0.66750896 0.68913686   0.767747  0.7301822
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                  0.8888889        0.8  0.6666667       0.75        0.2         0.7
# err                      0.11111111        0.2 0.33333334       0.25        0.8         0.3
# err_count                       1.0        2.0        2.0        1.0        4.0         3.0
# logloss                  0.90497464   1.276104  1.4704807 0.98298997  1.7129177   2.0270455
# max_per_class_error             0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy   0.9705882 0.88235295  0.9117647  0.9411765  0.7647059  0.88235295
# mean_per_class_error    0.029411765 0.11764706  0.0882353 0.05882353 0.23529412  0.11764706
# mse                      0.34836826 0.48172823 0.56616676 0.35757026  0.6338638  0.63608444
# r2                        0.9831031  0.9841328    0.96775 0.92371833 0.97829235  0.97504574
# rmse                      0.5902273  0.6940664  0.7524405  0.5979718 0.79615563    0.797549

pred= predict(RF, test)            
labels <- as.data.frame(test[,c(151)])[,1]
predict= as.data.frame(pred)
results= predict$predict
# results$true= list(labels)
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  8  12 6  10 3  13 13 14 17
# Levels: 1 10 12 13 14 17 2 3 4 5 6 8
test_accHb4hp= 19/24*100#79.17%
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
#    max_depth mtries ntrees     model_ids            logloss
# 1          5      5    100 rfh_model_137 0.9519759308236725
# 2         20      5    100 rfh_model_122 0.9530693010009622
# 3         10      5    100 rfh_model_131 0.9530693010009622
# 4         15      5    100 rfh_model_226 0.9530693010009622
# 5          5      6    500 rfh_model_123  0.956208604965141
# 6         15      6    500 rfh_model_135 0.9562135944785068
# 7         20      6    500 rfh_model_148 0.9562135944785068
# 8         20      6    400 rfh_model_141 0.9588620310770554
# 9          5      6    400 rfh_model_165 0.9588857631970644
# 10        15      6    300 rfh_model_147  0.964178215170615


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_193')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 5, mtries = 5,
                      ntrees = 100,
                     nfolds = 10)
# Error: water.exceptions.H2OModelBuilderIllegalArgumentException: Illegal argument(s) for DRF model: 
# rfh_model_193_cv_1.  Details: ERRR on field: _ntrees: 
# The tree model will not fit in the driver node's memory (1.7 KB per tree x 50 > Zero  ) - 
# try decreasing ntrees and/or max_depth or increasing min_rows!
print(best_rf)
print(drf)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564140389365_2 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1             100                     1700              193823         1         5    2.96294
# min_leaves max_leaves mean_leaves
# 1          2          9     4.37882
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.67300504   0.11312918 0.72727275  0.6666667        0.6      0.875
# err                     0.32699496   0.11312918 0.27272728 0.33333334        0.4      0.125
# err_count                      2.4    0.6480741        3.0        2.0        2.0        1.0
# logloss                  1.5032569   0.17147957  1.4125694  1.2466031  1.4314054  1.3303078
# max_per_class_error            1.0          0.0        1.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.88235295  0.030847318  0.8235294  0.9117647 0.88235295  0.9411765
# mean_per_class_error    0.11764706  0.030847318  0.1764706  0.0882353 0.11764706 0.05882353
# mse                      0.5616536   0.05264223  0.5241827 0.48091644 0.56883556 0.50570947
# r2                      0.97368264 0.0045259283  0.9785433  0.9793154 0.97189546  0.9716838
# rmse                     0.7479392  0.033470344  0.7240046  0.6934814 0.75421184  0.7111325
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                      0.25        0.7  0.6666667  0.6666667        0.8   0.7777778
# err                           0.75        0.3 0.33333334 0.33333334        0.2  0.22222222
# err_count                      3.0        3.0        4.0        3.0        1.0         2.0
# logloss                  2.1629944  1.3868624  1.6106602  1.5909778  1.4201378   1.4400512
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.8235294 0.88235295  0.8235294  0.9117647  0.9411765  0.88235295
# mean_per_class_error     0.1764706 0.11764706  0.1764706  0.0882353 0.05882353  0.11764706
# mse                      0.7605754 0.51772165  0.6097212 0.55597514 0.54920727   0.5436914
# r2                       0.9750631 0.96092665 0.97631514 0.96425873 0.97616285   0.9826618
# rmse                     0.8721098 0.71952873  0.7808465  0.7456374 0.74108523   0.7373543

prob=predict(drf, test)
labels <- as.data.frame(test[,c(151)])[,1]
predict= as.data.frame(prob)
results= predict$predict
results$true= list(labels)
print(labels)
print(predict$predict)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  13 9  6  10 10 13 13 11 17
# Levels: 1 10 11 13 17 2 3 4 5 6 8 9
# nrow(test)
test_accL= 20/24*100#83.33

######## H-dataset-125
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 125/')
dataset= as.data.frame(read.csv('VarImp_RF_H_125SelectedGenes.csv', row.names = 1, sep = ','))
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
RF= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(RF)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564140389365_3 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               90557         1         6    2.55176
# min_leaves max_leaves mean_leaves
# 1          2          8     3.79059

# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.67423075   0.1440152  0.6666667        0.8  0.6923077  0.8333333
# err                     0.32576925   0.1440152 0.33333334        0.2 0.30769232 0.16666667
# err_count                      2.4   0.9591663        2.0        2.0        4.0        1.0
# logloss                  1.3445415  0.18090361  1.2749567  1.2124014   1.404962  1.3588592
# max_per_class_error           0.85  0.22638462        1.0        0.5        1.0        1.0
# mean_per_class_accuracy 0.88529414 0.046920754 0.88235295  0.9411765  0.8235294  0.9411765
# mean_per_class_error    0.11470588 0.046920754 0.11764706 0.05882353  0.1764706 0.05882353
# mse                     0.49584073 0.065085016 0.47060215 0.44294536 0.53961366  0.5303416
# r2                       0.9696083 0.011719387 0.97143054  0.9722465 0.95914215 0.97869164
# rmse                    0.70092547 0.047666572 0.68600446  0.6655414 0.73458403 0.72824556
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                       1.0        0.7        0.5       0.75        0.2         0.6
# err                            0.0        0.3        0.5       0.25        0.8         0.4
# err_count                      0.0        3.0        3.0        1.0        4.0         4.0
# logloss                  0.8847052   1.326756  1.5736799  1.0099367  1.7034311   1.6957271
# max_per_class_error            0.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy        1.0 0.85294116 0.85294116  0.9411765  0.7647059  0.85294116
# mean_per_class_error           0.0 0.14705883 0.14705883 0.05882353 0.23529412  0.14705883
# mse                     0.33852485 0.48548594 0.58580095  0.3519418  0.6069346  0.60621643
# r2                       0.9835805   0.984009  0.9666316 0.92491907 0.97921455   0.9762175
# rmse                     0.5818289  0.6967682  0.7653763  0.5932468 0.77906007    0.778599

pred= predict(RF, test)            
labels <- as.data.frame(test[,c(126)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  8  11 9  10 3  13 11 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 20/24*100 #83.33
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1         20      6    100 rfh_model_10  0.887244391976268
# 2         20      6    400 rfh_model_12 0.9076333114555654
# 3          5      5    100 rfh_model_11 0.9108763511117426
# 4         15      5    100  rfh_model_3  0.911748184597554
# 5         15      6    200 rfh_model_20 0.9127660387548095
# 6         20      5   5000  rfh_model_8 0.9382812978075651
# 7          5      5   5000 rfh_model_22 0.9382844728808194
# 8          5      5   1000 rfh_model_14 0.9401272256560003
# 9         20      5    200 rfh_model_13 0.9497858793296996
# 10         5      4    400  rfh_model_6 0.9710234248664078


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_10')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 20, mtries = 6,
                      ntrees = 100,
                      nfolds = 10)

print(best_rf)
print(drf)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564140389365_4 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1             100                     1700              189969         1         6    2.86118
# min_leaves max_leaves mean_leaves
# 1          2         10     4.19765
# # Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.6820599  0.15705173        0.5        0.8  0.6666667 0.90909094
# err                     0.31794012  0.15705173        0.5        0.2 0.33333334 0.09090909
# err_count                      2.5   1.1067972        5.0        2.0        2.0        1.0
# logloss                  1.3848811  0.15211357  1.5993385  1.4641984  1.3533398  1.3469039
# max_per_class_error            0.9  0.21213204        1.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.87352943   0.0518265  0.7647059  0.9117647 0.88235295  0.9411765
# mean_per_class_error     0.1264706   0.0518265 0.23529412  0.0882353 0.11764706 0.05882353
# mse                      0.5245027  0.05482864  0.5947899  0.5656187  0.5018324 0.52763516
# r2                      0.97188747 0.009424315 0.96416926 0.97781014  0.9813561 0.94652945
# rmse                     0.7220096  0.04002996  0.7712262 0.75207627  0.7084013   0.726385
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.5555556 0.71428573        0.2      0.875        1.0         0.6
# err                     0.44444445  0.2857143        0.8      0.125        0.0         0.4
# err_count                      4.0        2.0        4.0        1.0        0.0         4.0
# logloss                  1.5046772  1.2827731  1.6828618  1.2246718  0.8813978   1.5086484
# max_per_class_error            1.0        1.0        1.0        1.0        0.0         1.0
# mean_per_class_accuracy  0.7941176 0.88235295  0.8235294  0.9411765        1.0   0.7941176
# mean_per_class_error    0.20588236 0.11764706  0.1764706 0.05882353        0.0  0.20588236
# mse                     0.57261956  0.4976872 0.62932456  0.4751081 0.33462322    0.545788
# r2                       0.9703059  0.9854495  0.9833688  0.9637582 0.98982567    0.956302
# rmse                     0.7567163  0.7054695  0.7932998 0.68928087  0.5784663  0.73877466
prob=predict(drf, test)
labels <- as.data.frame(test[,c(126)])[,1]
predict= as.data.frame(prob)
results= predict$predict
print(labels)
print(predict$predict)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  11 8  9  10 10 13 11 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8 9
nrow(test)
test_accL= 21/24*100#87.5


######## H-dataset-100
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 100/')
dataset= as.data.frame(read.csv('VarImp_RF_H_100SelectedGenes.csv', row.names = 1, sep = ','))
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
# h2o.shutdown()
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
RF= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(RF)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564140389365_5 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               90503         1         5    2.51529
# min_leaves max_leaves mean_leaves
# 1          2          7     3.78471

# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.65534186 0.083930716        0.5        0.8  0.6923077  0.6666667
# err                       0.3446581 0.083930716        0.5        0.2 0.30769232 0.33333334
# err_count                       2.6   0.6480741        3.0        2.0        4.0        2.0
# logloss                   1.4365638  0.26778883  1.3028479  1.2064646   1.656355  1.4092295
# max_per_class_error             0.9  0.14142136        1.0        0.5        1.0        1.0
# mean_per_class_accuracy  0.87941176  0.03411004  0.8235294  0.9411765 0.85294116 0.88235295
# mean_per_class_error    0.120588236  0.03411004  0.1764706 0.05882353 0.14705883 0.11764706
# mse                      0.50594753 0.074264534  0.4883278 0.45048508   0.584793  0.5495541
# r2                         0.969193 0.011404801 0.97035444  0.9717741  0.9557213  0.9779197
# rmse                     0.70722556  0.05375676 0.69880456 0.67118186  0.7647176  0.7413192
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.7  0.6666667       0.75        0.4         0.6
# err                     0.22222222        0.3 0.33333334       0.25        0.6         0.4
# err_count                      2.0        3.0        2.0        1.0        3.0         4.0
# logloss                 0.88326764   1.200982  1.5315603  1.1310881    1.73463   2.3092134
# max_per_class_error            0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765 0.85294116  0.9117647  0.9411765  0.8235294   0.8235294
# mean_per_class_error    0.05882353 0.14705883  0.0882353 0.05882353  0.1764706   0.1764706
# mse                     0.33782107 0.45377183 0.57421064 0.34101734 0.63125163  0.64824295
# r2                       0.9836147  0.9850536  0.9672918  0.9272496  0.9783818   0.9745687
# rmse                     0.5812238  0.6736259  0.7577669  0.5839669 0.79451346  0.80513537

pred= predict(RF, test)            
labels <- as.data.frame(test[,c(101)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  4  4  4  5  5  5  6  6  9  9  8  10 3  13 11 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 20/24*100 #83.33
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1          5      6   5000 rfh_model_15 0.8801944160638108
# 2         15      6    400  rfh_model_8 0.8859208828563875
# 3         20      5    200 rfh_model_14 0.9242069083371385
# 4         15      5    200 rfh_model_20 0.9242069083371385
# 5         15      5    500 rfh_model_21 0.9287343715631812
# 6         10      5    500 rfh_model_13 0.9287343715631812
# 7         15      5    300 rfh_model_17 0.9380101145599292
# 8         10      5    100  rfh_model_5 0.9465894995852153
# 9         15      4   1000  rfh_model_9 0.9655025215384242
# 10        15      4    400 rfh_model_22 0.9681610239633858


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_15')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 5,
                      mtries = 5, ntrees = 5000,
                      nfolds = 10)

print(best_rf)
print(drf)
# Model Details:
#   ==============
# H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564142749921_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1            5000                    85000             9738498         1         5    2.90294
# min_leaves max_leaves mean_leaves
# 1          2         10     4.27342
# Cross-Validation Metrics Summary: 
#   mean          sd  cv_1_valid cv_2_valid  cv_3_valid cv_4_valid
# accuracy                 0.7006868 0.091127805   0.8333333        0.5       0.875 0.84615386
# err                      0.2993132 0.091127805  0.16666667        0.5       0.125 0.15384616
# err_count                      2.3  0.83964276         1.0        3.0         1.0        2.0
# logloss                  1.4141717  0.16206184   1.1546481  1.6528873   1.1726128  1.2106744
# max_per_class_error            0.9  0.14142136         0.5        1.0         0.5        1.0
# mean_per_class_accuracy        0.9  0.04180197   0.9705882 0.85294116   0.9705882  0.9117647
# mean_per_class_error           0.1  0.04180197 0.029411765 0.14705883 0.029411765  0.0882353
# mse                       0.524325 0.041007206  0.44373357 0.63155746  0.45584968 0.47130564
# r2                      0.95973766 0.021881208   0.9880877  0.9631506  0.98787934  0.9771119
# rmse                    0.72300416 0.028195564   0.6661333  0.7947059   0.6751664 0.68651706
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.6666667  0.5714286  0.5714286  0.6666667 0.64285713   0.8333333
# err                     0.33333334 0.42857143 0.42857143 0.33333334 0.35714287  0.16666667
# err_count                      2.0        3.0        3.0        2.0        5.0         1.0
# logloss                  1.5674185  1.5005212  1.8610421  1.2542943  1.5306811   1.2369374
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.88235295 0.88235295  0.9411765 0.88235295  0.7647059   0.9411765
# mean_per_class_error    0.11764706 0.11764706 0.05882353 0.11764706 0.23529412  0.05882353
# mse                     0.56975174  0.5378409  0.5557408 0.50036585  0.5828581   0.4942463
# r2                       0.9403748 0.96369946 0.87622136 0.97642255 0.95994383   0.9644853
# rmse                      0.754819  0.7333764 0.74548024 0.70736545 0.76345146  0.70302653
prob=predict(drf, test)
labels <- as.data.frame(test[,c(101)])[,1]
predict= as.data.frame(prob)
results= predict$predict
print(labels)
print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  13 9  9  10 3  13 13 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
nrow(test)
test_accL= 21/24*100# 87.5
######## H-dataset-75
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 75/')
dataset= as.data.frame(read.csv('VarImp_RF_H_75SelectedGenes.csv', row.names = 1, sep = ','))
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
# h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size = "5G")#
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
RF= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(RF)
# Model Details:
#   ==============
# H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564143424608_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               91715         1         5    2.61765
# min_leaves max_leaves mean_leaves
# 1          2          7     3.89529
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.68081194  0.10912813  0.8333333        0.8  0.7692308  0.6666667
# err                     0.31918803  0.10912813 0.16666667        0.2 0.23076923 0.33333334
# err_count                      2.4  0.84852815        1.0        2.0        3.0        2.0
# logloss                  1.4045041  0.20880207   1.190685   1.073402  1.4385175  1.3096356
# max_per_class_error            0.9  0.14142136        1.0        0.5        1.0        1.0
# mean_per_class_accuracy  0.8882353  0.03924019  0.9411765  0.9411765 0.88235295 0.88235295
# mean_per_class_error    0.11176471  0.03924019 0.05882353 0.05882353 0.11764706 0.11764706
# mse                      0.5009718  0.06869367 0.45275113 0.39891857 0.55508447  0.5092229
# r2                      0.96843547 0.014304555  0.9725143  0.9750051 0.95797074 0.97954017
# rmse                    0.70437574 0.049125247 0.67286783     0.6316 0.74503994 0.71359855
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                  0.8888889        0.7        0.5       0.75        0.4         0.5
# err                      0.11111111        0.3        0.5       0.25        0.6         0.5
# err_count                       1.0        3.0        3.0        1.0        3.0         5.0
# logloss                   0.9000228  1.2460148  1.6191347  1.8281164  1.6839746   1.7555373
# max_per_class_error             0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy   0.9705882 0.85294116 0.85294116  0.9411765  0.8235294   0.7941176
# mean_per_class_error    0.029411765 0.14705883 0.14705883 0.05882353  0.1764706  0.20588236
# mse                      0.33926272 0.47921497 0.59773314  0.4120017  0.6327376  0.63279057
# r2                       0.98354477 0.98421556  0.9659519 0.91210634  0.9783309  0.97517496
# rmse                      0.5824626 0.69225353   0.773132  0.6418736 0.79544806   0.7954814

pred= predict(RF, test)            
labels <- as.data.frame(test[,c(76)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  6  3  3  3  3  4  4  5  5  5  6  6  8  8  8  10 3  13 13 11 17
# Levels: 1 10 11 13 17 2 3 4 5 6 8

test_accHb4hp= 16/24*100#66.67
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1         15      6   1000  rfh_model_3 0.8525549279894086
# 2         10      6   1000 rfh_model_61 0.8525549279894086
# 3          5      6   1000 rfh_model_54  0.852714592423893
# 4         10      6    500 rfh_model_14 0.8554728437838492
# 5         15      6    500 rfh_model_46 0.8554728437838492
# 6         15      6    300 rfh_model_25 0.8565919027537608
# 7         20      6    300 rfh_model_58 0.8565919027537608
# 8         15      6    400 rfh_model_24 0.8586608468917437
# 9         10      6    400 rfh_model_15 0.8586608468917437
# 10        20      6    400 rfh_model_42 0.8586608468917437


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_15')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 15,
                      mtries = 6, ntrees = 1000,
                      nfolds = 10)

print(best_rf)
print(drf)
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                  0.6953355  0.10930698      0.625       0.75 0.85714287        1.0
# err                       0.3046645  0.10930698      0.375       0.25 0.14285715        0.0
# err_count                       2.6   1.1916375        3.0        3.0        1.0        0.0
# logloss                   1.4272196  0.14015761  1.6748301  1.7670211  1.3005669  1.0937778
# max_per_class_error             0.9  0.21213204        1.0        1.0        1.0        0.0
# mean_per_class_accuracy   0.8955882 0.050183486  0.9411765  0.9411765  0.9411765        1.0
# mean_per_class_error    0.104411766 0.050183486 0.05882353 0.05882353 0.05882353        0.0
# mse                       0.5184132  0.03822948   0.576275 0.53344744 0.52354085 0.41653574
# r2                        0.9723845 0.010384441 0.98152226  0.9722784  0.9802665 0.96358955
# rmse                     0.71900576 0.026869552  0.7591278  0.7303749  0.7235612  0.6453958
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                      0.75 0.45454547        0.6       0.75        0.5   0.6666667
# err                           0.25 0.54545456        0.4       0.25        0.5  0.33333334
# err_count                      1.0        6.0        4.0        2.0        4.0         2.0
# logloss                  1.3578796  1.4763372  1.3270632  1.4573846  1.5957944   1.2215414
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765  0.7352941  0.8382353 0.88235295 0.85294116  0.88235295
# mean_per_class_error    0.05882353  0.2647059 0.16176471 0.11764706 0.14705883  0.11764706
# mse                      0.5201697 0.52246386  0.4758696  0.5578409   0.605292  0.45269704
# r2                       0.9873129  0.9777401  0.9325008  0.9706158  0.9784666    0.979552
# rmse                     0.7212279  0.7228166 0.68983305 0.74688745  0.7780052  0.67282766

pred= predict(drf, test)            
labels <- as.data.frame(test[,c(76)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 3  13 13 13 17
# Levels: 1 10 13 17 2 3 4 5 6 8 9

test_accHb4hp= 21/24*100#87.5


######## H-dataset-50
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 50/')
dataset= as.data.frame(read.csv('VarImp_RF_H_50SelectedGenes.csv', row.names = 1, sep = ','))
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
# h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size = '3G')#
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
RF= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(RF)
# Model Details:
#   ==============
# H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564144989683_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               92056         1         6    2.66118
# min_leaves max_leaves mean_leaves
# 1          2          8     3.93059
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.6674786  0.10545035  0.6666667        0.7  0.7692308  0.8333333
# err                     0.33252138  0.10545035 0.33333334        0.3 0.23076923 0.16666667
# err_count                      2.6   1.0099505        2.0        3.0        3.0        1.0
# logloss                  1.3537263  0.18378358  1.3091987   1.247402   1.300117  1.5494792
# max_per_class_error           0.95  0.10606602        1.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.87058824 0.051112197 0.88235295 0.88235295 0.88235295  0.9411765
# mean_per_class_error    0.12941177 0.051112197 0.11764706 0.11764706 0.11764706 0.05882353
# mse                     0.49254444 0.056243736 0.50203395 0.42506135  0.5146818 0.57529604
# r2                       0.9701112 0.010984282  0.9695224  0.9733671 0.96102995 0.97688544
# rmse                     0.6993504 0.041553836 0.70854354  0.6519673 0.71741325 0.75848275
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                  0.8888889        0.5  0.6666667       0.75        0.4         0.5
# err                      0.11111111        0.5 0.33333334       0.25        0.6         0.5
# err_count                       1.0        5.0        2.0        1.0        3.0         5.0
# logloss                   0.9724292  1.4665736  1.3171619 0.96785706  1.5191598   1.8878846
# max_per_class_error             0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy   0.9705882  0.7352941 0.88235295  0.9411765  0.8235294   0.7647059
# mean_per_class_error    0.029411765  0.2647059 0.11764706 0.05882353  0.1764706  0.23529412
# mse                       0.3790751  0.5398498 0.50522906 0.34061214 0.55830956   0.5852957
# r2                        0.9816137  0.9822184 0.97122115  0.9273361  0.9808798   0.9770382
# rmse                     0.61569077  0.7347447  0.7107947 0.58361983 0.74720114   0.7650462

pred= predict(RF, test)            
labels <- as.data.frame(test[,c(51)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  17 3  3  3  3  3  3  4  4  5  5  5  6  6  8  12 8  10 3  13 11 2  17
# Levels: 1 10 11 12 13 17 2 3 4 5 6 8

test_accHb4hp= 14/24*100# 58.33
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
#    max_depth mtries ntrees    model_ids            logloss
# 1         10      6   1000  rfh_model_2 0.8475851562875105
# 2         15      6    100  rfh_model_5 0.8757082145266195
# 3          5      5    300 rfh_model_12 0.8855723201279415
# 4         10      5    100  rfh_model_4 0.8942397074135566
# 5         10      4    500  rfh_model_8 0.9263761404976976
# 6          5      3    300  rfh_model_7  0.980709134366816
# 7         15      3    400  rfh_model_6 0.9812840694175369
# 8         15      3   5000  rfh_model_1 0.9875478100292615
# 9         20      4    100  rfh_model_9 0.9963373321833088
# 10         5      2    400 rfh_model_11 1.0828228676873994


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_15')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 10,
                      mtries = 6, ntrees = 1000,
                      nfolds = 10)

print(best_rf)
print(drf)
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.64411616  0.14229976  0.8333333 0.72727275      0.875 0.33333334
# err                     0.35588384  0.14229976 0.16666667 0.27272728      0.125  0.6666667
# err_count                      2.6   0.9591663        1.0        3.0        1.0        4.0
# logloss                  1.3586637  0.19609265  1.0501792  1.5371517  1.1170111  1.8406888
# max_per_class_error           0.95  0.10606602        1.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.88235295  0.05008643  0.9411765  0.9411765  0.9411765  0.8235294
# mean_per_class_error    0.11764706  0.05008643 0.05882353 0.05882353 0.05882353  0.1764706
# mse                      0.5045191  0.06770588 0.40932557 0.49890095 0.43434817  0.6681387
# r2                       0.9299463  0.08351783 0.98439014 0.98509455  0.9583858  0.9548724
# rmse                      0.707161 0.047129862  0.6397856  0.7063292 0.65905094 0.81739753
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid  cv_9_valid cv_10_valid
# accuracy                       0.7       0.75 0.33333334        0.5   0.8888889         0.5
# err                            0.3       0.25  0.6666667        0.5  0.11111111         0.5
# err_count                      3.0        2.0        2.0        4.0         1.0         5.0
# logloss                  1.2528417 0.96277785  1.6815994   1.380063   1.1853014   1.5790223
# max_per_class_error            1.0        1.0        1.0        1.0         0.5         1.0
# mean_per_class_accuracy  0.9117647 0.88235295 0.88235295  0.7647059   0.9705882   0.7647059
# mean_per_class_error     0.0882353 0.11764706 0.11764706 0.23529412 0.029411765  0.23529412
# mse                     0.48070157  0.3619997 0.65595555  0.5076272  0.45541185   0.5727818
# r2                      0.98324496 0.98305196  0.5783143 0.95974207   0.9719693   0.9403973
# rmse                     0.6933265 0.60166407 0.80991083  0.7124796   0.6748421   0.7568235
pred= predict(drf, test)            
labels <- as.data.frame(test[,c(51)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  6  6  9  8  10 3  13 11 11 17
# Levels: 1 10 11 13 17 2 3 4 5 6 8 9

test_accHb4hp= 18/24*100#75


###########################################SELECTED GENES MINIMODELS
library(h2o)
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 150/')
######Genes from total geneset narrowed down to 175 genes
dataset= as.data.frame(read.csv('VarImp_RF_L_150SelectedGenes.csv', row.names = 1, sep = ','))
#Dataset with deceptive classes
# dataset= dataset[,c(1:151)]
#Dataset with correct classes
dataset= dataset[,c(1:150,152)]
#Dataset with compounds
# dataset= dataset[,c(1:150,153)]
# Featurescaling
# dataset[-2464] = scale(dataset[-2464]) #has already been feature scaled
# rm(dataset)
#Compound
# dataset$Compound = factor(dataset$Compound) # specify that the last classifier is a factor(not #)
#Classes
dataset$Class = factor(dataset$Class) # specify that the last classifier is a factor(not #)
# h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size = '5G')#
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
RF= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(RF)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564145715563_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               90643         1         5    2.54000
# min_leaves max_leaves mean_leaves
# 1          2          8     3.79765
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.6867094  0.13206667  0.8333333        0.8 0.46153846        0.5
# err                      0.3132906  0.13206667 0.16666667        0.2 0.53846157        0.5
# err_count                      2.5   1.3874437        1.0        2.0        7.0        3.0
# logloss                  1.3994712  0.28715375  1.0789123  1.2524871  1.6412867  1.4901443
# max_per_class_error           0.85  0.16201852        1.0        0.5        1.0        1.0
# mean_per_class_accuracy 0.88529414  0.06406771  0.9411765  0.9411765  0.6764706  0.8235294
# mean_per_class_error    0.11470588  0.06406771 0.05882353 0.05882353 0.32352942  0.1764706
# mse                     0.50465506  0.07529053 0.42336625 0.45588166  0.6042313 0.57192445
# r2                      0.96887505 0.012644292  0.9742982 0.97143596  0.9542495  0.9770209
# rmse                     0.7064217  0.05302579   0.650666  0.6751901  0.7773232  0.7562569
# cv_5_valid  cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                  0.8888889         0.9  0.8333333       0.75        0.4         0.5
# err                      0.11111111         0.1 0.16666667       0.25        0.6         0.5
# err_count                       1.0         1.0        1.0        1.0        3.0         5.0
# logloss                   0.9182518   1.1613648  1.3353152  1.0227739  1.7429589   2.3512166
# max_per_class_error             0.5         0.5        1.0        1.0        1.0         1.0
# mean_per_class_accuracy   0.9705882   0.9705882  0.9411765  0.9411765  0.8235294   0.8235294
# mean_per_class_error    0.029411765 0.029411765 0.05882353 0.05882353  0.1764706   0.1764706
# mse                      0.35331798  0.45452896 0.50752926 0.37079516  0.6144438  0.69053173
# r2                         0.982863   0.9850287  0.9710901   0.920897  0.9789574   0.9729097
# rmse                      0.5944056  0.67418766  0.7124109  0.6089295  0.7838646   0.8309824

pred= predict(RF, test)            
labels <- as.data.frame(test[,c(151)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  6  8  11 6  10 3  13 13 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8
test_accHb4hp= 18/24*100 #75%


##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1          5      6    200 rfh_model_34 0.9378731724947076
# 2         10      6    300 rfh_model_18 0.9403903894026979
# 3          5      6    300 rfh_model_30 0.9405185138840548
# 4         15      6    500 rfh_model_15 0.9457732666783094
# 5         20      5   1000 rfh_model_35 0.9727961826632047
# 6         15      5   1000 rfh_model_33 0.9727961826632047
# 7         10      5   1000 rfh_model_32 0.9727961826632047
# 8          5      5    500 rfh_model_31 0.9759366042227121
# 9         15      5    500  rfh_model_9 0.9759968504409777
# 10        15      5    300 rfh_model_38 0.9857279417370317


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_7')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 5,
                      mtries = 6, ntrees = 200,
                      nfolds = 10)

print(best_rf)
print(drf)
# H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564219599226_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1             200                     3400              381866         1         5    2.86647
# min_leaves max_leaves mean_leaves
# 1          2         10     4.19824
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.65641415   0.13493566  0.8181818        0.5  0.6666667  0.8181818
# err                     0.34358585   0.13493566 0.18181819        0.5 0.33333334 0.18181819
# err_count                      2.5    0.7905694        2.0        4.0        3.0        2.0
# logloss                  1.4485192    0.1587054  1.3964504  1.5984583  1.4158387  1.4358865
# max_per_class_error            1.0          0.0        1.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.87058824   0.04382254  0.9117647  0.8235294  0.9117647 0.88235295
# mean_per_class_error    0.12941177   0.04382254  0.0882353  0.1764706  0.0882353 0.11764706
# mse                     0.54325205  0.053439636  0.5336597  0.6024937  0.5322056  0.5628734
# r2                       0.9723532 0.0078551695  0.9697693 0.97624177 0.96626866  0.9771144
# rmse                     0.7354263  0.034642678  0.7305201  0.7762047  0.7295242  0.7502489
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.6666667        0.6        0.8       0.75  0.7777778  0.16666667
# err                     0.33333334        0.4        0.2       0.25 0.22222222   0.8333333
# err_count                      2.0        2.0        2.0        1.0        2.0         5.0
# logloss                  1.2703145  1.3246512  1.3497007  1.3737957  1.2566936   2.0634027
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.88235295 0.88235295 0.88235295  0.9411765 0.88235295   0.7058824
# mean_per_class_error    0.11764706 0.11764706 0.11764706 0.05882353 0.11764706  0.29411766
# mse                     0.47566617  0.4785064 0.52193147  0.5100331 0.47639957  0.73875165
# r2                      0.98627883  0.9667704  0.9833195  0.9832776 0.96800303   0.9464888
# rmse                     0.6896855  0.6917415  0.7224482   0.714166   0.690217  0.85950667
pred= predict(drf, test)            
labels <- as.data.frame(test[,c(151)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  6  8  8  10 3  13 13 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8

test_accHb4hp= 19/24*100#79.17
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
# h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size = "5G")#
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
rf= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(rf)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564145903623_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      750               82282         1         6    2.71067
# min_leaves max_leaves mean_leaves
# 1          2         10     4.05333
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.66876066   0.08300664  0.6666667        0.8 0.61538464        0.5
# err                      0.3312393   0.08300664 0.33333334        0.2  0.3846154        0.5
# err_count                      2.6    0.9591663        2.0        2.0        5.0        3.0
# logloss                  1.3641136   0.20166585  1.1779276   1.087362  1.7105938   1.621274
# max_per_class_error            0.9   0.14142136        1.0        0.5        1.0        1.0
# mean_per_class_accuracy  0.8611111   0.04820865  0.8666667 0.93333334 0.76666665        0.8
# mean_per_class_error     0.1388889   0.04820865 0.13333334 0.06666667 0.23333333        0.2
# mse                     0.49632764   0.06579136 0.42873746  0.4072897  0.5885232 0.60621345
# r2                      0.97450215 0.0043365443 0.97054476  0.9720268 0.96375346  0.9711327
# rmse                     0.7012844  0.047580738 0.65478045 0.63819253 0.76715267  0.7785971
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                  0.8888889        0.7  0.6666667       0.75        0.6         0.5
# err                      0.11111111        0.3 0.33333334       0.25        0.4         0.5
# err_count                       1.0        3.0        2.0        1.0        2.0         5.0
# logloss                   1.0860951  1.2878257  1.3721327 0.90705895  1.6489388   1.7419268
# max_per_class_error             0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.96666664        0.8        0.9 0.93333334  0.8666667   0.7777778
# mean_per_class_error    0.033333335        0.2        0.1 0.06666667 0.13333334  0.22222222
# mse                      0.42702854  0.4542144  0.5295148 0.33317986  0.5918149   0.5967601
# r2                         0.975607 0.98139226  0.9698378 0.98498344 0.98136604    0.974377
# rmse                      0.6534742  0.6739543 0.72767764 0.57721734 0.76929504   0.7725025

pred= predict(rf, test)            
labels <- as.data.frame(test[,c(126)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  1  1  1  11 12 12 15
# Levels: 1 11 12 15 2 3 4 5 6 8
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  6  15 1  8  1  3  11 9  12 15
# Levels: 1 11 12 15 2 3 4 5 6 8 9
test_accHb4hp= 19/24*100 #79,17
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1         10      6   1000  rfh_model_6 0.9204154824809099
# 2         15      6    500 rfh_model_11 0.9299376541100939
# 3         20      6    500 rfh_model_59 0.9299376541100939
# 4         20      6    400 rfh_model_30 0.9368587588319971
# 5         15      6    400 rfh_model_17 0.9368587588319971
# 6         20      6    300 rfh_model_60 0.9407269546938074
# 7         15      6    300 rfh_model_41 0.9407269546938074
# 8         20      5   1000 rfh_model_53 0.9429863743818299
# 9         15      5   1000 rfh_model_72 0.9429863743818299
# 10        10      5   1000 rfh_model_22 0.9429863743818299


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_7')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 10,
                      mtries = 6, ntrees = 1000,
                      nfolds = 10)

print(best_rf)
print(drf)
# H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564217439985_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1            1000                    15000             1737433         1         8    3.00093
# min_leaves max_leaves mean_leaves
# 1          2         13     4.45353
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.6800577  0.12317345        0.5        0.8 0.44444445        1.0
# err                     0.31994227  0.12317345        0.5        0.2  0.5555556        0.0
# err_count                      2.8   1.2961482        4.0        1.0        5.0        0.0
# logloss                  1.3597076  0.17145155  1.3977748  1.1421897   1.649085  1.1024573
# max_per_class_error            0.8   0.2345208        1.0        1.0        1.0        0.0
# mean_per_class_accuracy  0.8888889 0.052234042        0.8 0.93333334  0.9111111        1.0
# mean_per_class_error    0.11111111 0.052234042        0.2 0.06666667 0.08888889        0.0
# mse                     0.51100564  0.05015912 0.54436296 0.45593637 0.63147366 0.43787408
# r2                       0.9631938  0.02436086  0.9731179  0.9853113 0.98412496  0.9793334
# rmse                     0.7131813  0.03448221  0.7378096  0.6752306  0.7946532  0.6617205
# cv_5_valid cv_6_valid  cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                0.45454547  0.7777778  0.85714287  0.6666667        0.6         0.7
# err                     0.54545456 0.22222222  0.14285715 0.33333334        0.4         0.3
# err_count                      6.0        2.0         1.0        2.0        4.0         3.0
# logloss                  1.8725572  1.4660374    1.164291  1.1227777   1.439472   1.2404343
# max_per_class_error            1.0        0.5         0.5        1.0        1.0         1.0
# mean_per_class_accuracy 0.76666665 0.93333334  0.96666664  0.8666667        0.8   0.9111111
# mean_per_class_error    0.23333333 0.06666667 0.033333335 0.13333334        0.2  0.08888889
# mse                     0.63161725  0.5168674   0.4465423  0.4360003  0.5352692  0.47411293
# r2                       0.9482211  0.8640706  0.97854847 0.97241473  0.9719608  0.97483474
# rmse                    0.79474354  0.7189349  0.66823816  0.6603032  0.7316209   0.6885586
pred= predict(drf, test)            
labels <- as.data.frame(test[,c(126)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  1  1  1  11 12 12 15
# Levels: 1 11 12 15 2 3 4 5 6 8
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  1  1  8  1  1  11 9  12 15
# Levels: 1 11 12 15 2 3 4 5 6 8 9

test_accHb4hp= 22/24*100#91.67

######## L-dataset-100
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 100/')
dataset= as.data.frame(read.csv('VarImp_RF_L_100SelectedGenes.csv', row.names = 1, sep = ','))
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
# h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size = "5G")#
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
rf= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(rf)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564146047079_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               90565         1         5    2.55176
# min_leaves max_leaves mean_leaves
# 1          2          8     3.78824
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid  cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.6944017  0.12639537  0.8333333         0.9 0.53846157        0.5
# err                      0.3055983  0.12639537 0.16666667         0.1 0.46153846        0.5
# err_count                      2.4   1.2328829        1.0         1.0        6.0        3.0
# logloss                  1.4375714   0.2514449   1.163354   1.0974592  1.6107309  1.5281224
# max_per_class_error           0.95  0.10606602        1.0         0.5        1.0        1.0
# mean_per_class_accuracy 0.88529414 0.050473537  0.9411765   0.9705882  0.7647059  0.8235294
# mean_per_class_error    0.11470588 0.050473537 0.05882353 0.029411765 0.23529412  0.1764706
# mse                      0.4985258 0.057544004 0.44333962   0.4166715  0.5910395  0.5856707
# r2                      0.96826094 0.014818981 0.97308564  0.97389275 0.95524836 0.97646856
# rmse                     0.7037456 0.040422376  0.6658375    0.645501  0.7687909  0.7652913
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.8888889        0.8  0.8333333       0.75        0.4         0.5
# err                     0.11111111        0.2 0.16666667       0.25        0.6         0.5
# err_count                      1.0        2.0        1.0        1.0        3.0         5.0
# logloss                 0.98292613  1.2861464  1.2493305   2.029871  1.3486127   2.0791605
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765  0.9117647  0.9411765  0.9411765  0.8235294   0.7941176
# mean_per_class_error    0.05882353  0.0882353 0.05882353 0.05882353  0.1764706  0.20588236
# mse                     0.38484466  0.4849099  0.4950395  0.4232604 0.51626074   0.6442215
# r2                      0.98133385   0.984028  0.9718015 0.90970445 0.98231983   0.9747265
# rmse                     0.6203585 0.69635475  0.7035904 0.65058464  0.7185129   0.8026341

pred= predict(rf, test)            
labels <- as.data.frame(test[,c(101)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  13 9  6  10 4  13 11 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 20/24*100#83.44
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1          5      6   5000 rfh_model_14 0.8756171249989272
# 2         20      6    200 rfh_model_39 0.8836320603455067
# 3         15      6   1000 rfh_model_24 0.8838416803010948
# 4          5      6    200 rfh_model_26 0.8845212690701518
# 5         20      6    400 rfh_model_25 0.8965850456782981
# 6         15      5    200 rfh_model_34 0.9083078675626509
# 7          5      5    200  rfh_model_8 0.9084814898703458
# 8         15      5   1000 rfh_model_16 0.9107691777229058
# 9          5      5   5000 rfh_model_22 0.9126342876024472
# 10        15      5    400 rfh_model_23 0.9175875014457805


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_7')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 5,
                      mtries = 6, ntrees = 5000,
                      nfolds = 10)

print(best_rf)
print(drf)
# H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564216737587_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1            5000                    85000             9598945         1         5    2.80679
# min_leaves max_leaves mean_leaves
# 1          2         10     4.14304
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.75632674  0.07917146       0.75        0.6        1.0        0.8
# err                     0.24367326  0.07917146       0.25        0.4        0.0        0.2
# err_count                      1.8   0.7615773        1.0        2.0        0.0        1.0
# logloss                    1.34176  0.09477708  1.5429109  1.3241954  1.1761271  1.2437173
# max_per_class_error            0.9  0.21213204        1.0        1.0        0.0        1.0
# mean_per_class_accuracy 0.90588236  0.04342595  0.9411765 0.88235295        1.0  0.9411765
# mean_per_class_error    0.09411765  0.04342595 0.05882353 0.11764706        0.0 0.05882353
# mse                      0.5073022 0.037868317 0.60056376 0.50387716 0.46727064 0.49794096
# r2                      0.95895815 0.027148735  0.9164433 0.97510487  0.9833704 0.97972554
# rmse                    0.71125305 0.026658285  0.7749605 0.70984304   0.683572  0.7056493
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.8181818      0.625  0.6666667  0.6923077  0.7777778   0.8333333
# err                     0.18181819      0.375 0.33333334 0.30769232 0.22222222  0.16666667
# err_count                      2.0        3.0        1.0        4.0        2.0         2.0
# logloss                  1.2135553  1.4119643  1.2232085  1.2554017  1.5214819   1.5050384
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy 0.88235295 0.85294116  0.9411765  0.7647059  0.9117647   0.9411765
# mean_per_class_error    0.11764706 0.14705883 0.05882353 0.23529412  0.0882353  0.05882353
# mse                      0.4732866  0.5420971 0.40766612  0.4780887  0.5765403   0.5256908
# r2                       0.9856472  0.9750042  0.8588848  0.9746559  0.9791519   0.9615934
# rmse                     0.6879583  0.7362724 0.63848734 0.69143957  0.7593025   0.7250454
pred= predict(drf, test)            
labels <- as.data.frame(test[,c(101)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  8  10 3  13 13 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 21/24*100#87.5

######## L-dataset-75
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 75/')
dataset= as.data.frame(read.csv('VarImp_RF_L_75SelectedGenes.csv', row.names = 1, sep = ','))
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
# h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size = "5G")#
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
rf= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(rf)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564146184043_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               91985         1         6    2.63529
# min_leaves max_leaves mean_leaves
# 1          2          8     3.92824
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid  cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.72517097   0.1242497  0.8333333         0.9 0.84615386  0.8333333
# err                      0.27482906   0.1242497 0.16666667         0.1 0.15384616 0.16666667
# err_count                       2.0   0.8944272        1.0         1.0        2.0        1.0
# logloss                   1.3306059  0.23934646   1.059285   1.0277185  1.3835034  1.4815592
# max_per_class_error             0.9  0.14142136        1.0         0.5        1.0        1.0
# mean_per_class_accuracy  0.90294117 0.041646477  0.9411765   0.9705882  0.9117647  0.9411765
# mean_per_class_error    0.097058825 0.041646477 0.05882353 0.029411765  0.0882353 0.05882353
# mse                      0.48539022  0.06808908 0.41261777  0.39083934  0.5314021   0.538062
# r2                        0.9700893 0.012084997  0.9749507   0.9755113  0.9597639 0.97838145
# rmse                      0.6932425 0.049015645  0.6423533   0.6251715  0.7289733 0.73352706
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                  0.8888889        0.8        0.5       0.75        0.4         0.5
# err                      0.11111111        0.2        0.5       0.25        0.6         0.5
# err_count                       1.0        2.0        3.0        1.0        3.0         5.0
# logloss                  0.94475085  1.2991433  1.5750005 0.97541076  1.4516685    2.108019
# max_per_class_error             0.5        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy   0.9705882 0.88235295 0.85294116  0.9411765  0.8235294   0.7941176
# mean_per_class_error    0.029411765 0.11764706 0.14705883 0.05882353  0.1764706  0.20588236
# mse                       0.3656829  0.4615362 0.59947175  0.3584265 0.55193776  0.64392585
# r2                       0.98226327  0.9847979 0.96585286  0.9235357   0.981098   0.9747381
# rmse                      0.6047172 0.67936456 0.77425563  0.5986873  0.7429251   0.8024499

pred= predict(rf, test)            
labels <- as.data.frame(test[,c(76)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  8  9  10 3  13 11 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 21/24*100#87.5%
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1         20      6    100  rfh_model_8 0.8294647304209105
# 2         15      6    300  rfh_model_9  0.852851199763693
# 3         20      6   1000  rfh_model_4 0.8620860128921354
# 4         20      6    500  rfh_model_2 0.8636005704734541
# 5         15      5    100  rfh_model_5 0.9180293896392766
# 6         20      5    200 rfh_model_10 0.9211977061943984
# 7          5      4    100  rfh_model_6 0.9295864600950793
# 8          5      4    200  rfh_model_1 0.9440309568702447
# 9         10      4    400 rfh_model_11 0.9451424246183063
# 10        15      4   1000  rfh_model_7 0.9490967245765132


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_7')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 20,
                      mtries = 6, ntrees = 100,
                      nfolds = 10)

print(best_rf)
print(drf)
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.74348485  0.11995035       0.75 0.42857143        0.7  0.5714286
# err                     0.25651515  0.11995035       0.25  0.5714286        0.3 0.42857143
# err_count                      2.1   0.9192388        2.0        4.0        3.0        3.0
# logloss                  1.3626722  0.21778214  1.2079031  2.1442893  1.2652397  1.5486453
# max_per_class_error            0.8  0.28284273        1.0        1.0        1.0        1.0
# mean_per_class_accuracy        0.9 0.044798665  0.9117647  0.8235294 0.88235295  0.8235294
# mean_per_class_error           0.1 0.044798665  0.0882353  0.1764706 0.11764706  0.1764706
# mse                       0.503458 0.056636173 0.47229618   0.677778 0.49722052 0.57578313
# r2                       0.9745007 0.006663918 0.95958966   0.964518  0.9718129  0.9796146
# rmse                    0.70738083 0.039181434  0.6872381  0.8232727  0.7051386  0.7588037
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.8333333        1.0        1.0  0.8181818  0.6666667   0.6666667
# err                     0.16666667        0.0        0.0 0.18181819 0.33333334  0.33333334
# err_count                      1.0        0.0        0.0        2.0        3.0         3.0
# logloss                    1.14481  1.0242898  1.0645164  1.4368956  1.3203928   1.4697397
# max_per_class_error            1.0        0.0        0.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765        1.0        1.0  0.9117647  0.8235294  0.88235295
# mean_per_class_error    0.05882353        0.0        0.0  0.0882353  0.1764706  0.11764706
# mse                      0.4156886 0.40089506 0.42454198  0.5436421 0.48823527  0.53849965
# r2                       0.9688234 0.98839855  0.9814745  0.9645578 0.98068994  0.98552805
# rmse                     0.6447392 0.63316274  0.6515689  0.7373209 0.69873834   0.7338254
pred= predict(drf, test)            
labels <- as.data.frame(test[,c(76)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  6  12 8  10 3  13 13 14 17
# Levels: 1 10 12 13 14 17 2 3 4 5 6 8

test_accHb4hp= 19/24*100#79.17



######## H-dataset-50
setwd('/Users/ashle/Desktop/LAbbook/R files/models/Variable Importance/Top 50/')
dataset= as.data.frame(read.csv('VarImp_RF_L_50SelectedGenes.csv', row.names = 1, sep = ','))
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
h2o.shutdown()
# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1, max_mem_size = "5G")#
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
rf= h2o.randomForest(x = predictors, y =response,
             training_frame = train,
             nfolds = 10,
             seed= 1234)
print(rf)
# Model Details:
#   ==============
#   
#   H2OMultinomialModel: drf
# Model ID:  DRF_model_R_1564146360888_1 
# Model Summary: 
#   number_of_trees number_of_internal_trees model_size_in_bytes min_depth max_depth mean_depth
# 1              50                      850               91872         1         5    2.64235
# min_leaves max_leaves mean_leaves
# 1          2          9     3.92000
# Cross-Validation Metrics Summary: 
#   mean          sd cv_1_valid  cv_2_valid cv_3_valid cv_4_valid
# accuracy                 0.6630342  0.10293197        0.5         0.9  0.7692308  0.6666667
# err                      0.3369658  0.10293197        0.5         0.1 0.23076923 0.33333334
# err_count                      2.5   0.7905694        3.0         1.0        3.0        2.0
# logloss                  1.2965903  0.18691102  1.2674845   1.1015458   1.368708  1.4210575
# max_per_class_error           0.95  0.10606602        1.0         0.5        1.0        1.0
# mean_per_class_accuracy 0.87647057  0.04032738  0.8235294   0.9705882 0.88235295 0.88235295
# mean_per_class_error    0.12352941  0.04032738  0.1764706 0.029411765 0.11764706 0.11764706
# mse                      0.4751523 0.058746792 0.46953496  0.40406182  0.5100537  0.5164839
# r2                       0.9715951 0.009565656 0.97149533  0.97468287 0.96138036  0.9792484
# rmse                     0.6864897 0.044069555  0.6852262  0.63565856 0.71418047 0.71866816
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                 0.7777778        0.7  0.6666667       0.75        0.4         0.5
# err                     0.22222222        0.3 0.33333334       0.25        0.6         0.5
# err_count                      2.0        3.0        2.0        1.0        3.0         5.0
# logloss                  1.0119519  1.3115095  1.2227006  0.8404347  1.7427771   1.6777334
# max_per_class_error            1.0        1.0        1.0        1.0        1.0         1.0
# mean_per_class_accuracy  0.9411765  0.8235294 0.88235295  0.9411765  0.8235294   0.7941176
# mean_per_class_error    0.05882353  0.1764706 0.11764706 0.05882353  0.1764706  0.20588236
# mse                      0.3994835 0.49749193 0.47409207 0.30434558  0.5965687   0.5794071
# r2                      0.98062384 0.98361355 0.97299474 0.93507296 0.97956955  0.97726923
# rmse                    0.63204706  0.7053311 0.68854344 0.55167526  0.7723786   0.7611879

pred= predict(rf, test)            
labels <- as.data.frame(test[,c(51)])[,1]
predict= as.data.frame(pred)
results= predict$predict
print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  6  11 9  10 3  13 13 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 20/24*100# 83.33%
##################hyperparameter tuning
search_criteria <- list(strategy = 'RandomDiscrete')
grid_space <- list()
grid_space$ntrees <- c(100,200,300,400,500,1000,5000)
grid_space$max_depth <- c(1,5,10,15,20)
# grid_space$nbins <- c(6, 4, 3)
# grid_space$nbins_cats <- c(370, 449)
grid_space$mtries <- c(2, 4, 3,5,6)
# grid_space$sample_rate <- c(0.3, 0.7, 0.4,0.5,0.2)
# rm(grid_space)

rfh= h2o.grid('randomForest', grid_id = 'rfh', 
              x= predictors,y=response,
              training_frame = train,
              seed=1234,
              hyper_params = grid_space,
              search_criteria = search_criteria)

params=rfh@summary_table
print(params[order(params$logloss,decreasing=F)[1:10],])
# Hyper-Parameter Search Summary: ordered by increasing logloss
# max_depth mtries ntrees    model_ids            logloss
# 1         15      6   5000  rfh_model_7 0.8350148571792178
# 2         20      6   5000 rfh_model_11 0.8350148571792178
# 3         20      6   1000 rfh_model_13 0.8450038824390577
# 4         20      6    200  rfh_model_9 0.8457678893375086
# 5          5      5   1000 rfh_model_12  0.879558705116483
# 6          5      5    400  rfh_model_8 0.8844726628585194
# 7          5      4    200  rfh_model_4 0.9163894775721108
# 8         20      3    400  rfh_model_3 0.9757873559302014
# 9         20      3    100  rfh_model_2 0.9841496910476207
# 10         5      2    100  rfh_model_5 1.0898570721140641


# Grab the model_id for the top DL model, chosen by validation AUC
best_rf <- h2o.getModel('rfh_model_7')
# h2o.saveModel(best_rf, path = "/Users/ashle/Desktop/LAbbook/R files/", force = T)
drf= h2o.randomForest(predictors,response,training_frame = train, max_depth = 15,
                      mtries = 6, ntrees = 5000,
                      nfolds = 10)

print(best_rf)
print(drf)
# Cross-Validation Metrics Summary: 
#   mean           sd cv_1_valid cv_2_valid cv_3_valid cv_4_valid
# accuracy                0.77107143   0.07793124        0.8      0.625 0.85714287 0.71428573
# err                     0.22892857   0.07793124        0.2      0.375 0.14285715  0.2857143
# err_count                      1.8    0.6164414        2.0        3.0        1.0        2.0
# logloss                  1.2999768   0.15934344  1.2877597   1.496299  1.1067747  1.1504567
# max_per_class_error            0.9   0.21213204        1.0        1.0        1.0        1.0
# mean_per_class_accuracy 0.90588236  0.036970016 0.88235295  0.8235294  0.9411765  0.9117647
# mean_per_class_error    0.09411765  0.036970016 0.11764706  0.1764706 0.05882353  0.0882353
# mse                      0.4831165  0.061429206 0.47560742 0.58458245  0.3665536  0.4308117
# r2                       0.9735242 0.0067787957  0.9678861   0.978694 0.97886926  0.9589304
# rmse                    0.69221836  0.044442356 0.68964297 0.76457995  0.6054367  0.6563625
# cv_5_valid cv_6_valid cv_7_valid cv_8_valid cv_9_valid cv_10_valid
# accuracy                     0.875       0.75 0.71428573      0.625        1.0        0.75
# err                          0.125       0.25  0.2857143      0.375        0.0        0.25
# err_count                      1.0        2.0        2.0        3.0        0.0         2.0
# logloss                   1.021785  1.5960964  1.2949201  1.5169468 0.95234245   1.5763876
# max_per_class_error            1.0        1.0        1.0        1.0        0.0         1.0
# mean_per_class_accuracy  0.9411765 0.88235295  0.9117647  0.8235294        1.0   0.9411765
# mean_per_class_error    0.05882353 0.11764706  0.0882353  0.1764706        0.0  0.05882353
# mse                      0.3844085 0.62136775 0.50848025  0.5570955 0.37439162  0.52786624
# r2                       0.9802867 0.97728866  0.9799714 0.96861434 0.98801345   0.9566879
# rmse                    0.62000686  0.7882688   0.713078  0.7463883  0.6118755    0.726544
pred= predict(drf, test)            
labels <- as.data.frame(test[,c(51)])[,1]
predict= as.data.frame(pred)
results= predict$predict

print(labels)
print(predict$predict)
nrow(test)
# > print(labels)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  9  9  10 10 13 14 14 17
# Levels: 1 10 13 14 17 2 3 4 5 6 8 9
# > print(predict$predict)
# [1] 1  2  3  3  3  3  3  3  4  4  5  5  5  6  8  9  11 9  10 3  13 13 14 17
# Levels: 1 10 11 13 14 17 2 3 4 5 6 8 9

test_accHb4hp= 21/24*100#87.5
