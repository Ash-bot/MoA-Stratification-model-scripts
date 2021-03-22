################################################175 genes
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

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.6,0.2), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
valid <- h2o.assign(splits[[2]], "valid.hex") # 20%
test <- h2o.assign(splits[[3]], "test.hex") # 20%
response <- "Class"
response <- "Compound"
predictors <- setdiff(names(dataset), response)
print(response)

# Deeplearning hyperparamters
activation_opt <- list("Rectifier", "RectifierWithDropout", "Maxout", "MaxoutWithDropout")
l1_opt <- list(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
l2_opt <- list(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
input_dropout_ratio = list(0, 0.15, 0.3)
hidden_dropout_ratios = list(c(0,0), c(0.15,0.15),c(0.3,0.3))
hidden = list(c(32,32))
hyper_params <- list(activation = activation_opt,
                     l1 = l1_opt,
                     l2 = l2_opt,
                     input_dropout_ratio = input_dropout_ratio,
                     hidden_dropout_ratios=hidden_dropout_ratios,
                     hidden= hidden)
search_criteria <- list(strategy = "RandomDiscrete", 
                        max_runtime_secs = 120)


dl_grid <- h2o.grid("deeplearning", x = predictors, y =response,
                    grid_id = "dl_grid",
                    training_frame = train,
                    validation_frame = valid,
                    seed = 123,
                    hyper_params = hyper_params,
                    search_criteria = search_criteria)
params=dl_grid@summary_table
print(params[order(params$logloss,decreasing=T)[179:189],])

dl_gridperf <- h2o.getGrid(grid_id = "dl_grid", 
                           sort_by = "logloss", 
                           decreasing = T)


# Note that that these results are not reproducible since we are not using a single core H2O cluster
# H2O's DL requires a single core to be used in order to get reproducible results

# Grab the model_id for the top DL model, chosen by validation AUC
best_dl_model_id <- dl_gridperf@model_ids[[1]]
best_dl <- h2o.getModel(best_dl_model_id)

# Now let's evaluate the model performance on a test set
# so we get an honest estimate of top model performance
best_dl_perf <- h2o.performance(model = best_dl, 
                                newdata = test)

print(best_dl_perf)
h2o.auc(best_dl_perf)
# H2OMultinomialMetrics: deeplearning
# 
# Test Set Metrics: 
#   =====================
#   
#   MSE: (Extract with `h2o.mse`) 0.3588112
# RMSE: (Extract with `h2o.rmse`) 0.5990085
# Logloss: (Extract with `h2o.logloss`) 1.430108
# Mean Per-Class Error: 0.2058824
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>, <data>)`)
# =========================================================================
#   Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
#        1 10 11 12 13 14 15 16 17 2 3 4 5 6 7 8 9  Error     Rate
# 1      1  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 1
# 10     0  1  0  0  0  0  0  0  0 0 1 0 0 0 0 0 0 0.5000 =  1 / 2
# 11     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 12     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 13     0  0  0  0  1  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =  0 / 1
# 14     0  0  0  0  1  0  1  0  0 0 0 0 0 0 0 0 0 1.0000 =  2 / 2
# 15     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 16     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 17     0  0  0  0  0  0  0  0  1 0 0 0 0 0 0 0 0 0.0000 =  0 / 1
# 2      0  0  0  0  0  0  0  0  0 1 0 0 0 0 0 0 0 0.0000 =  0 / 1
# 3      0  0  0  0  0  0  0  0  0 0 6 0 0 0 0 0 0 0.0000 =  0 / 6
# 4      0  0  0  0  0  0  0  0  0 0 0 2 0 0 0 0 0 0.0000 =  0 / 2
# 5      0  0  0  0  0  0  0  0  0 0 0 0 3 0 0 0 0 0.0000 =  0 / 3
# 6      0  0  0  0  0  0  0  0  0 0 0 0 0 1 0 0 0 0.0000 =  0 / 1
# 7      0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =  0 / 0
# 8      0  0  0  0  0  1  0  0  0 0 0 0 0 0 0 0 0 1.0000 =  1 / 1
# 9      0  0  0  0  2  0  0  0  0 1 0 0 0 0 0 0 0 1.0000 =  3 / 3
# Totals 1  1  0  0  4  1  1  0  1 2 7 2 3 1 0 0 0 0.2917 = 7 / 24
# 
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>, <data>)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.708333
# 2   2  0.708333
# 3   3  0.750000
# 4   4  0.791667
# 5   5  0.833333
# 6   6  0.833333
# 7   7  0.916667
# 8   8  0.916667
# 9   9  0.916667
# 10 10  0.916667
average= mean(c(0.708333,0.708333,0.750000,0.791667,0.833333,0.833333,0.916667,0.916667,0.916667,0.916667))
std= sd(c(0.708333,0.708333,0.750000,0.791667,0.833333,0.833333,0.916667,0.916667,0.916667,0.916667))

#########################All genes
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

# Splitting the dataset into the Training set and Test set
# install.packages('caTools')

# install.packages("h2o")
library(h2o)
h2o.init(nthreads = -1)#
library(caTools)
splits <- h2o.splitFrame(as.h2o(dataset), c(0.6,0.2), seed=1234)
train <- h2o.assign(splits[[1]], "train.hex") # 60%
valid <- h2o.assign(splits[[2]], "valid.hex") # 20%
test <- h2o.assign(splits[[3]], "test.hex") # 20%
response <- "Class"
response <- "Compound"
predictors <- setdiff(names(dataset), response)
print(response)
validd= as.data.frame(valid)
###No HPT
dl_fit3 <- h2o.deeplearning(x = predictors, y =response,
                            training_frame = train,
                            model_id = "dl_fit3",
                            validation_frame = valid,  #in DL, early stopping is on by default
                            epochs = 20,
                            hidden = c(10,10),
                            score_interval = 1,           #used for early stopping
                            stopping_rounds = 3,          #used for early stopping
                            stopping_metric = "logloss",      #used for early stopping
                            stopping_tolerance = 0.0005,  #used for early stopping
                            seed = 123)

dl_perf3 <- h2o.performance(model = dl_fit3,
                            newdata = test)

print(dl_perf3)
print(dl_perf3@metrics$predictions)
pred= h2o.predict(dl_fit3,
                  newdata = test)

# Deeplearning hyperparamters
activation_opt <- list("Rectifier", "RectifierWithDropout", "Maxout", "MaxoutWithDropout")
l1_opt <- list(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
l2_opt <- list(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
input_dropout_ratio = list(0, 0.15, 0.3)
hidden_dropout_ratios = list(c(0,0), c(0.15,0.15),c(0.3,0.3))
hidden = list(c(32,32))
hyper_params <- list(activation = activation_opt,
                     l1 = l1_opt,
                     l2 = l2_opt,
                     input_dropout_ratio = input_dropout_ratio,
                     hidden_dropout_ratios=hidden_dropout_ratios,
                     hidden= hidden)
search_criteria <- list(strategy = "RandomDiscrete", 
                        max_runtime_secs = 120)


dl_grid <- h2o.grid("deeplearning", x = predictors, y =response,
                    grid_id = "dl_grid",
                    training_frame = train,
                    validation_frame = valid,
                    seed = 123,
                    hyper_params = hyper_params,
                    search_criteria = search_criteria)

params=dl_grid@summary_table
print(params[order(params$logloss,decreasing=T)[293:303],]) #change this to the last 10 in params

dl_gridperf <- h2o.getGrid(grid_id = "dl_grid", 
                           sort_by = "logloss", 
                           decreasing = T)


# Note that that these results are not reproducible since we are not using a single core H2O cluster
# H2O's DL requires a single core to be used in order to get reproducible results

# Grab the model_id for the top DL model, chosen by validation AUC
best_dl_model_id <- dl_gridperf@model_ids[[1]]
best_dl <- h2o.getModel(best_dl_model_id)

# Now let's evaluate the model performance on a test set
# so we get an honest estimate of top model performance
best_dl_perf <- h2o.performance(model = best_dl, 
                                newdata = test)

print(best_dl_perf)
# H2OMultinomialMetrics: deeplearning
# 
# Test Set Metrics: 
#   =====================
#   
#   MSE: (Extract with `h2o.mse`) 0.5551371
# RMSE: (Extract with `h2o.rmse`) 0.7450752
# Logloss: (Extract with `h2o.logloss`) 4.397392
# Mean Per-Class Error: 0.4705882
# Confusion Matrix: Extract with `h2o.confusionMatrix(<model>, <data>)`)
# =========================================================================
#   Confusion Matrix: Row labels: Actual class; Column labels: Predicted class
# 1 10 11 12 13 14 15 16 17 2 3 4 5 6 7 8 9  Error      Rate
# 1      1  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =   0 / 1
# 10     1  1  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0.5000 =   1 / 2
# 11     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =   0 / 0
# 12     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =   0 / 0
# 13     0  0  0  0  1  0  0  0  0 0 0 0 0 0 0 0 0 0.0000 =   0 / 1
# 14     1  0  0  0  1  0  0  0  0 0 0 0 0 0 0 0 0 1.0000 =   2 / 2
# 15     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =   0 / 0
# 16     0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =   0 / 0
# 17     1  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 1.0000 =   1 / 1
# 2      0  0  0  0  0  0  1  0  0 0 0 0 0 0 0 0 0 1.0000 =   1 / 1
# 3      0  0  0  0  0  0  0  0  0 0 5 0 0 0 0 1 0 0.1667 =   1 / 6
# 4      0  0  0  0  0  0  0  0  0 0 2 0 0 0 0 0 0 1.0000 =   2 / 2
# 5      0  0  0  0  0  0  0  0  0 0 0 0 2 0 1 0 0 0.3333 =   1 / 3
# 6      0  0  0  0  0  0  0  0  0 0 1 0 0 0 0 0 0 1.0000 =   1 / 1
# 7      0  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0     NA =   0 / 0
# 8      1  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 1.0000 =   1 / 1
# 9      0  1  0  0  2  0  0  0  0 0 0 0 0 0 0 0 0 1.0000 =   3 / 3
# Totals 5  2  0  0  4  0  1  0  0 0 8 0 2 0 1 1 0 0.5833 = 14 / 24
# 
# Hit Ratio Table: Extract with `h2o.hit_ratio_table(<model>, <data>)`
# =======================================================================
#   Top-10 Hit Ratios: 
#   k hit_ratio
# 1   1  0.416667
# 2   2  0.541667
# 3   3  0.666667
# 4   4  0.708333
# 5   5  0.750000
# 6   6  0.791667
# 7   7  0.791667
# 8   8  0.791667
# 9   9  0.833333
# 10 10  0.833333


average= mean(c(0.416667,0.541667,0.666667,0.708333,0.750000,0.791667,0.791667,0.791667,0.833333,0.833333))
std= sd(c(0.416667,0.541667,0.666667,0.708333,0.750000,0.791667,0.791667,0.791667,0.833333,0.833333))
