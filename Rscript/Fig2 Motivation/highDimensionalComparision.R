simDataset.plr = function(n = 100,
                          dims = 10,
                          shiftParts = sample(1:choose(dims,2),num_shiftedComponents),
                          varDiag = 1,
                          maxCov = 1,
                          baseNoise = 1e-11,
                          num_shiftedComponents = 5,
                          offset_c1  = 2,
                          offset_c2 = -2,
                          num_noiseComponents = 2,
                          noiseParts = NULL,
                          noiseCenter1 = -1
){


  ## function
  ## compuye number of logratios
  nfeats  = choose(dims,2)

  ## Define mean vector
  m = rep(baseNoise,nfeats)

  ### Class 1
  ## addshifted mean
  m[shiftParts] = offset_c1

  ## add noise components
  if(is.null(noiseParts)){
    pts = 1:nfeats
    pts = pts[!pts%in%shiftParts]
    noiseParts = sample(pts,num_noiseComponents)
  }
  m[noiseParts] = rnorm(length(noiseParts),noiseCenter1,log(abs(noiseCenter1))/3)

  ## Define indv logratio variance
  s = diag(varDiag,nfeats)

  ## add random covaraince
  U = matrix(stats::runif(nfeats * nfeats, 0, maxCov), nrow = nfeats)
  diag(U) = varDiag
  U_ = s + U
  U_ = Matrix::nearPD(U_)$mat

  ## simulate from mvn
  y1 = MASS::mvrnorm(n, mu = m, Sigma = U_, tol = 1e-06,
                     empirical = FALSE)
  c1 = data.frame(Status = "C1",compositions::pwlrInv(y1))


  ### Class 2
  ## Define mean vector
  m = rep(baseNoise,nfeats)

  ## addshifted mean
  m[shiftParts] = offset_c2

  ## add noise components () same between classes no need to re compute noise location; use the one specififed for class 1
  # pts = 1:nfeats
  # pts = pts[!pts%in%shiftParts]
  # noiseParts = sample(pts,num_noiseComponents)
  m[noiseParts] = noiseCenter1

  ### Define indv logratio variance
  s = diag(varDiag,nfeats)

  ## add random covaraince
  U = matrix(stats::runif(nfeats * nfeats, 0, maxCov), nrow = nfeats)
  diag(U) = varDiag
  U_ = s + U
  U_ = Matrix::nearPD(U_)$mat

  y1 = MASS::mvrnorm(n, mu = m, Sigma = U_, tol = 1e-06,
                     empirical = FALSE)
  c2 = data.frame(Status = "C2",compositions::pwlrInv(y1))


  mat= rbind(c1,c2)
  return(mat)
}




## Parms
## note -- increaseing variance decreases separation. i.e. lower variance easier classification default sould be 0.1
varDiag = .1
maxCov = .1
n = 100
baseNoise = 1e-11


sd = 6
set.seed(sd)

dims1 = 30
nrs = choose(dims1,2) # num of ratios

signalPercent = .1
numShiftedComponents = round(signalPercent*nrs)
offset_c1  = 10; offset_c2 = -10



noisePercent = .9
numNoiseComponents = round(nrs*noisePercent)
noiseCenter1 = 80
noiseCenter2 = -80


## use needs to define shifts parts to be consisitnet acors both groups
nf = choose(dims1,2)
shiftParts1 = sample(1:nf,numShiftedComponents)
pts = 1:nf
pts = pts[!pts%in%shiftParts1]
noiseParts1 = sample(pts,numNoiseComponents)


mat = simDataset.plr(n = 100,offset_c1 = offset_c1,offset_c2 = offset_c2,
                     dims = dims1 ,noiseParts = noiseParts1,
                     num_shiftedComponents = numShiftedComponents,
                     shiftParts = shiftParts1,
                     num_noiseComponents = numNoiseComponents,
                     noiseCenter1 = noiseCenter1
                     )
mat = data.frame(Dataset = "Train",mat)

mat2 = simDataset.plr(n = 100,offset_c1 = offset_c1,offset_c2 = offset_c2,
                      dims = dims1 ,noiseParts = noiseParts1,
                      num_shiftedComponents = numShiftedComponents,
                      shiftParts = shiftParts1,
                      num_noiseComponents = numNoiseComponents,
                      noiseCenter1 = noiseCenter2
                      )
mat2 = data.frame(Dataset = "Test",mat2)


mat = rbind(mat,mat2)

pc = prcomp((mat[,-2:-1]))
pc.df = data.frame(mat[,1:2],pc$x)
ggplot(pc.df,aes(PC1,PC2,col = Status,shape = Dataset))+
  geom_point()

mat1 = mat %>%
  filter(Dataset=="Train")
pc = prcomp((mat1[,-2:-1]))
pc.df = data.frame(mat1[,1:2],pc$x)
ggplot(pc.df,aes(PC1,PC2,col = Status,shape = Dataset))+
  geom_point()


mat1 = mat %>%
  filter(Dataset=="Test")
pc = prcomp((mat1[,-2:-1]))
pc.df = data.frame(mat1[,1:2],pc$x)
ggplot(pc.df,aes(PC1,PC2,col = Status,shape = Dataset))+
  geom_point()





#
# Partition Data
df = mat %>%
  filter(Dataset == "Train") %>%
  select(-Dataset)
df$Status = factor(df$Status)

## Daatset 2
df2 = mat %>%
  filter(Dataset == "Test") %>%
  select(-Dataset)
df2$Status = factor(df2$Status)
classes = as.character(unique(df$Status))


ovr.perf = data.frame()


# Relative Abundance --------------------------------------------------------------

## with feature selection
testx = (df2[,-1])
trainx = (df[,-1])


## with feature selection
cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=1,family="binomial")
features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
features = features[-1,]
features = features[abs(features)>0]
trainx = clo(subset(df,select = names(features)))
testx = clo(subset(df2,select = names(features)))

## with fes
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glmnet",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

varImp(m1$models$glmnet1)
train.all.CLR = m1$performance
mroc.all.CLR = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])


## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc.glm.RA = pROC::roc(df2$Status,p)
## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)

## train cv metrics
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glm",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
train.glm.RA = m1$performance
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::roc(df2$Status,p)


perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "LASSO",Type = "Relative_Abudance",
                  train_auc = train.glm.RA$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc)))
ovr.perf = rbind(ovr.perf,perf)


## without feature selection
testx = (df2[,-1])
trainx = (df[,-1])

m1 = trainML_Models(trainLRs = df[,-1],testLRs = df2[,-1],
                    ytrain = df[,1],y_test = df2[,1],models = "ranger",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

train.all.RA = m1$performance
mroc.all.RA = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])

perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "RF",Type = "Relative_Abudance",
                  train_auc =  m1$performance$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc.all.RA)))
ovr.perf = rbind(ovr.perf,perf)








# CLR with feature selection --------------------------------------------------------------

## with feature selection
testx = clr(df2[,-1])
trainx = clr(df[,-1])

cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=1,family="binomial")
cv.clrlasso$cvm
features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
features = features[-1,]
features = features[abs(features)>0]
trainx = clr(subset(df,select = names(features)))
testx = clr(subset(df2,select = names(features)))

p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")

coef(cv.clrlasso)

## with fes
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glmnet",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

varImp(m1$models$glmnet1)
train.all.CLR = m1$performance
mroc.all.CLR = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])

## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glmnet",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
train.glm.RA = m1$performance
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])


perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "LASSO",Type = "CLR",
                  train_auc = m1$performance$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc)))
ovr.perf = rbind(ovr.perf,perf)


## without feature selection
testx = data.frame(clr(df2[,-1]))
trainx = data.frame(clr(df[,-1]))
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "ranger",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

train.all.CLR = m1$performance
mroc.all.CLR = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])


perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "RF",Type = "CLR",
                  train_auc =  m1$performance$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc.all.CLR)))
ovr.perf = rbind(ovr.perf,perf)




# ALR with feature selection --------------------------------------------------------------
library(easyCODA)
ff = FINDALR(df[,-1])


## with feature selection
ref = ff$procrust.ref
testx = data.frame(alr(df2[,-1],ivar = ref))
trainx = data.frame(alr(df[,-1],ivar = ref))


cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=1,family="binomial")
cv.clrlasso$cvm
features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
features = features[-1,]
features = features[abs(features)>0]
trainx = subset(trainx,select = names(features))
testx = (subset(testx,select = names(features)))

## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glmnet",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
train.glm.RA = m1$performance
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])


perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "LASSO",Type = "ALR",
                  train_auc = m1$performance$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc)))
ovr.perf = rbind(ovr.perf,perf)


## without feature selection
testx = data.frame(alr(df2[,-1],ivar = ref))
trainx = data.frame(alr(df[,-1],ivar = ref))
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "ranger",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

train.all.ALR = m1$performance
mroc.all.ALR = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])

perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "RF",Type = "ALR",
                  train_auc =  m1$performance$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc.all.ALR)))
ovr.perf = rbind(ovr.perf,perf)






# Log Ratios --------------------------------------------------------------

## with feature selection
plr = calcLogRatio(df)
plr2 = calcLogRatio(df2)


testx = (plr2[,-1])
trainx = (plr[,-1])

cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=1,family="binomial")
features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
features = features[-1,]
features = features[abs(features)>0]
trainx = (subset(plr,select = names(features)))
testx = (subset(plr2,select = names(features)))


## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glmnet",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
train.glm.RA = m1$performance
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])


perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "LASSO",Type = "PLR",
                  train_auc = m1$performance$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc)))
ovr.perf = rbind(ovr.perf,perf)




## without feature selection
library(igraph)
### PLR Approach
plr = calcLogRatio(df)
plr2 = calcLogRatio(df2)
## with feature selection
testx = (plr2[,-1])
trainx = (plr[,-1])

cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = plr,
                                    includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                    rankOrder = F)
trainx = diffCompVarRcpp::mstAll(featMatrix = plr,dcvRanking = cc.dcv$lrs)
testx = subset(testx,select = colnames(trainx))



ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
ensemble = c("ranger")
m1.plr = trainML_Models(trainLRs = trainx,testLRs = testx,testIDs = data.frame(ID = 1:nrow(df2),Status =df2$Status),
                        ytrain = df[,1],y_test = df2[,1],models = ensemble,
                        mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                        numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
train.all.PLR = m1.plr$performance
pmat = m1.plr$predictionMatrix
pmat = pmat %>%
  group_by(ID,Status) %>%
  dplyr::select(-model) %>%
  summarise_all(.funs = mean)
pmat = data.frame(pmat)
classes = as.character(unique(df$Status))
mroc.all.PLR = pROC::multiclass.roc(pmat$Status,pmat[,classes])

perf = data.frame(Seed = sd,
                  Dims = dims1,signal_percent = signalPercent,noise_percent = noisePercent,
                  Model = "RF",Type = "PLR",
                  train_auc =  m1$performance$TrainAUC,
                  test_auc = as.numeric(pROC::auc(mroc.all.PLR)))
ovr.perf = rbind(ovr.perf,perf)

