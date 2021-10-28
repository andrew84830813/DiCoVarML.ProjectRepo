library(diffCompVarRcpp)
library(selEnergyPermR)
library(simplexDataAugmentation)
library(DiCoVarML)
library(caret)
library(dplyr)
library(compositions)
library(foreach)
library(parallel)
library(readr)
library(tidyr)
library(stringr)
library(glmnet) # glmnet
library(selbal) # selbal



#setwd("/nas/longleaf/home/andrew84/rProjects/DiCoVarFS_project")


# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}
source(file = 'CODA_Functions/functions_coda_penalized_regression.R')



# Setup Cluster ---------------------------------------------------------------

# # Detect Cores
# clus <- parallel::makeCluster(10)
# doParallel::registerDoParallel(clus)




# Read External Args ---------------------------------------------------------------

args = c(2,4,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
dataset = as.numeric(args[2])
permute_labels = as.logical(as.numeric(args[3])) #should be 0(False) or 1(True)





# Load Data ---------------------------------------------------------------

load("Output/microbiomeHD_cdiSchubert.Rda"); # a go
df = expResults$taxa.df;
df = data.frame(sample_id = rownames(df),df);
df$Status = str_replace_all(df$Status,pattern = "\\/","_");
df = df[!is.na(df$Status),];
keep = unique(df$Status);
keep = c("H","CDI");
df = df[df$Status%in%keep,];
df = df[,-1];
f_name = "mbiomeHD_cdiSchubert_16s";

f_name = paste0(f_name,"_seed",sd)
mxSparsePercent = .7
benchmark = data.frame()

## Partition Data
sd = 1
k_fold = 5
overll_folds = caret::createFolds(df[,1],k =k_fold,list = F)
allData = lodo_partition(df,dataset_labels = overll_folds,sd)

## Extract Test/Train Split
ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                          fold = 2,
                                          permLabels = F,
                                          maxSparisty = mxSparsePercent,
                                          extractTelAbunance = F)

## Pre-Process
trainx = data.frame(fastImputeZeroes(ttData$train_Data,impFactor = ttData$imp_factor))
testx = data.frame(fastImputeZeroes(ttData$test_data,impFactor = ttData$imp_factor))

df = data.frame(Status = ttData$y_train,trainx) %>%
  arrange(desc(Status))
trainx = df[,-1]

## compute log ratios
lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = df$Status,trainx))
lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_test,testx))



test = data.frame(Status = ttData$test_ids,testx)


## DCV Parms
scale_data = T
performRFE = F
useRidgeWeight = F
min_connected = F
ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")



# Count/RA Table ----------------------------------------------------------

## data
dd = data.frame(ID = rownames(df),clo(df[,-1]))
#dd = data.frame(ID = dd[,1],apply(dd[,-1], 2, scale))

dd =   gather(dd,key = "taxa",value = "count",2:ncol(dd))
pdf(file = "Figures/fi1_raMatrix.pdf",width = 2.25 ,height = 1.75)
ggplot(dd,aes(ID,taxa,fill = count))+
  geom_tile()+
  coord_flip()+
  #scale_fill_distiller(palette = "Purples")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()


## data
dd = data.frame(ID = rownames(testx),clo(testx))
dd =   gather(dd,key = "taxa",value = "count",2:ncol(dd))
pdf(file = "Figures/fi1_raMatrix_val.pdf",width = 2.25 ,height = 1.75)
ggplot(dd,aes(ID,taxa,fill = count))+
  geom_tile()+
  coord_flip()+
  #scale_fill_distiller(palette = "Purples")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()

dd = data.frame(ID = rownames(df),Taxa  = "fef",Status = if_else(df$Status=="H",100,-100) )
pdf(file = "Figures/fi1_raMatrixLabels.pdf",width = .4 ,height = 1.75)
ggplot(dd,aes(ID,Taxa,fill =(Status)))+
  geom_tile()+
  coord_flip()+
  scale_fill_distiller(palette = "RdBu",direction = 1)+
  #scale_fill_viridis_c()+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()

# log ratio matrix -------------------------------------------------------

#lr = data.frame(Status = lr$Status,apply(lr[,-1], 2, scale))
dd = data.frame(ID = rownames(lrs.train),lrs.train[,-1])
dd =   gather(dd,key = "taxa",value = "count",2:ncol(dd))

pdf(file = "Figures/fi1_lrMatrixLabels.pdf",width = 4 ,height = 1.75)
ggplot(dd,aes(ID,taxa,fill = count))+
  geom_tile()+
  coord_flip()+
  scale_fill_viridis_c(option = "E")+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()







# log ratio network -------------------------------------------------------

library(igraph)
cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                    includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                    rankOrder = F)
dcv = cc.dcv$lrs %>%
  separate(Ratio,into = c("num","den"),sep = "___",remove = F) %>%
  filter(rowmean>0)

g = igraph::graph_from_edgelist(as.matrix(dcv[,2:3]),directed = T)
E(g)$weight = dcv$rowmean
cc = dcv_strength(dcv)
ccc = top_n(cc,n = 10,wt = Str)

dcv$keep = (dcv$num%in%ccc$ID+dcv$den%in%ccc$ID)

cc$rank=1
cc$rank[1:10] = 2
colnames(cc)[1] = "ID"

toGephi(g,Name = "fig1",attribute_metadata = cc,edge_metadata = dcv)


# targe selection network -------------------------------------------------

## get subcomposition
train_subcomp = subset(trainx,select = cc$ID[1:10])
test_subcomp = subset(testx,select = cc$ID[1:10])

dd = data.frame(ID = rownames(trainx),train_subcomp)
dd =   gather(dd,key = "taxa",value = "count",2:ncol(dd))
pdf(file = "Figures/fig1_raMatrix_targeted.pdf",width = .75 ,height = 1.75)
ggplot(dd,aes(ID,taxa,fill = count))+
  geom_tile()+
  coord_flip()+
  #scale_fill_distiller(palette = "Purples")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()

dd = data.frame(ID = rownames(test_subcomp),test_subcomp)
dd =   gather(dd,key = "taxa",value = "count",2:ncol(dd))
pdf(file = "Figures/fig1_raMatrix_targeted_testset.pdf",width = .75 ,height =1.75 )
ggplot(dd,aes(ID,taxa,fill = count))+
  geom_tile()+
  coord_flip()+
  #scale_fill_distiller(palette = "Purples")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()




## get pair logratio from subcomp
lrs.train = selEnergyPermR::calcLogRatio(df = data.frame(Status = ttData$y_train,train_subcomp))[,-1]
lrs.test = selEnergyPermR::calcLogRatio(df = data.frame(Status = ttData$y_test,test_subcomp))[,-1]

cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = data.frame(Status = ttData$y_train,lrs.train),
                                    includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                    rankOrder = F)
dcv = cc.dcv$lrs %>%
  separate(Ratio,into = c("num","den"),sep = "___",remove = F) %>%
  filter(rowmean>0)

g = igraph::graph_from_edgelist(as.matrix(dcv[,2:3]),directed = T)
E(g)$weight = dcv$rowmean

toGephi(g,Name = "fig1_targeted")



# Targeted LR MAtix -------------------------------------------------------

lr = data.frame(Status = ttData$y_train,lrs.train)
lr = lr %>%
  arrange(desc(Status))
dd = data.frame(ID = rownames(lr),lr[,-1])
dd =   gather(dd,key = "taxa",value = "count",2:ncol(dd))

pdf(file = "Figures/fi1_lrMatrixTargeted.pdf",width = 2 ,height = 1.75)
ggplot(dd,aes(ID,taxa,fill = count))+
  geom_tile()+
  coord_flip()+
  scale_fill_viridis_c(option = "E")+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()

lr = data.frame(Status = ttData$y_test,lrs.test)
dd = data.frame(ID = rownames(lr),lr[,-1])
dd =   gather(dd,key = "taxa",value = "count",2:ncol(dd))

pdf(file = "Figures/fi1_lrMatrixTargeted_test.pdf",width = 2 ,height = 1.75)
ggplot(dd,aes(ID,taxa,fill = count))+
  geom_tile()+
  coord_flip()+
  scale_fill_viridis_c(option = "E")+
  theme_classic()+
  theme(axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.line = element_blank())
dev.off()




# Targeted Performance Estimation -----------------------------------------

perc_totalParts2Keep = .75
num_sets = 5

base_dims = ncol(ttData$train_Data)
max_parts = round(perc_totalParts2Keep*base_dims)
sets = round(seq(10,max_parts,length.out = num_sets))


# Run K-Fold Cross Validation ---------------------------------------------
inner_perf = data.frame()
classes = unique(ttData$y_train)


#ensemble = c("ranger","xgbTree","xgbLinear")
max_sparsity = .9
train_auc = data.frame()

## Tune target features
compTime1 = system.time({

  for(sd1 in 1:1){
    set.seed(sd1)
    k_fold = 2
    overll_folds = caret::createFolds(ttData$y_train,k = k_fold,list = F)
    innerfold_data = lodo_partition(data.frame(Status = ttData$y_train,ttData$train_Data),
                                    dataset_labels = overll_folds,
                                    sd1)


    ## Get within fold cross validated performance

    for(nf in 1:k_fold){

      ## Partition inner fold
      innerFold = DiCoVarML::extractTrainTestSplit(foldDataList = innerfold_data,
                                                   fold = nf,
                                                   maxSparisty = max_sparsity,
                                                   extractTelAbunance = F)

      suppressMessages(suppressWarnings({

        ## Pre-Process
        trainx = data.frame(fastImputeZeroes(innerFold$train_Data,impFactor = innerFold$imp_factor))
        testx = data.frame(fastImputeZeroes(innerFold$test_data,impFactor = innerFold$imp_factor))

        ## compute log ratios
        lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_train,trainx))
        lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_test,testx))

        cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                            includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                            rankOrder = F)

      }))



      for(tar_Features in sets){

        suppressMessages(suppressWarnings({

          tar_dcvInner = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,
                                               useRidgeWeights = useRidgeWeight,use_rfe = performRFE,scaledata = scale_data,
                                               testx = testx,
                                               dcv = cc.dcv$lrs,lrs.train = lrs.train,lrs.test = lrs.test,
                                               y_label = innerFold$y_train,
                                               seed = sd1,
                                               ensemble = ensemble,
                                               y_test = innerFold$y_test,
                                               tarFeatures = tar_Features,
                                               ts.id = innerFold$test_ids,
                                               max_sparsity = max_sparsity
          )

          perf = data.frame(Seed = sd1,Fold = nf,tar_Features ,tar_dcvInner$Performance)
          inner_perf = rbind(inner_perf,perf)

          pmat = tar_dcvInner$all_model_preds
          pmat = pmat %>%
            group_by(model) %>%
            summarise(train_auc = as.numeric(pROC::auc(Status,!!as.name(classes[1])) ))
          pmat = data.frame(pmat)
          pmat = rbind(pmat,data.frame(model =tar_dcvInner$Performance$Approach,train_auc =  as.numeric(tar_dcvInner$Performance$AUC)))
          pmat$targetFeatures = tar_Features
          train_auc = rbind(train_auc,pmat)

        }))

        message(tar_Features)

      }

    }

  }

})

train_auc1 = train_auc %>%
  group_by(model,targetFeatures) %>%
  summarise_all(.funs = mean)
ggplot(train_auc1,aes(targetFeatures,train_auc,col = model))+
  geom_point()+
  geom_line()

train_auc2 = train_auc %>%
  group_by(targetFeatures) %>%
  summarise_all(.funs = mean)
ggplot(train_auc2,aes(targetFeatures,train_auc,col = model))+
  geom_point()+
  geom_line()





# Train Models ------------------------------------------------------------

classes. as.character(uniqTsparse(ttData$y_train))

## Extract Test/Train Split
ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                          fold = 1,
                                          permLabels = F,
                                          maxSparisty = mxSparsePercent,
                                          extractTelAbunance = F)

## Pre-Process
trainx = data.frame(fastImputeZeroes(ttData$train_Data,impFactor = ttData$imp_factor))
testx = data.frame(fastImputeZeroes(ttData$test_data,impFactor = ttData$imp_factor))


## compute log ratios
lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_train,trainx))
lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_test,testx))


perc_totalParts2Keep = .75
num_sets = 5

base_dims = ncol(ttData$train_Data)
max_parts = round(perc_totalParts2Keep*base_dims)
sets = round(seq(10,max_parts,length.out = num_sets))


# Run K-Fold Cross Validation ---------------------------------------------
inner_perf = data.frame()


#ensemble = c("ranger","xgbTree","xgbLinear")
max_sparsity = .9
train_auc = data.frame()

## Tune target features
compTime1 = system.time({

  for(sd1 in 1:1){
    set.seed(sd1)
    k_fold = 2
    overll_folds = caret::createFolds(ttData$y_train,k = k_fold,list = F)
    innerfold_data = lodo_partition(data.frame(Status = ttData$y_train,ttData$train_Data),
                                    dataset_labels = overll_folds,
                                    sd1)


    ## Get within fold cross validated performance

    for(nf in 1:k_fold){

      ## Partition inner fold
      innerFold = DiCoVarML::extractTrainTestSplit(foldDataList = innerfold_data,
                                                   fold = nf,
                                                   maxSparisty = max_sparsity,
                                                   extractTelAbunance = F)

      suppressMessages(suppressWarnings({

        ## Pre-Process
        trainx = data.frame(fastImputeZeroes(innerFold$train_Data,impFactor = innerFold$imp_factor))
        testx = data.frame(fastImputeZeroes(innerFold$test_data,impFactor = innerFold$imp_factor))

        ## compute log ratios
        lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_train,trainx))
        lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_test,testx))

        cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                            includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                            rankOrder = F)

      }))



      for(tar_Features in sets){

        suppressMessages(suppressWarnings({

          tar_dcvInner = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,
                                               useRidgeWeights = useRidgeWeight,use_rfe = performRFE,scaledata = scale_data,
                                               testx = testx,
                                               dcv = cc.dcv$lrs,lrs.train = lrs.train,lrs.test = lrs.test,
                                               y_label = innerFold$y_train,
                                               seed = sd1,
                                               ensemble = ensemble,
                                               y_test = innerFold$y_test,
                                               tarFeatures = tar_Features,
                                               ts.id = innerFold$test_ids,
                                               max_sparsity = max_sparsity
          )

          perf = data.frame(Seed = sd1,Fold = nf,tar_Features ,tar_dcvInner$Performance)
          inner_perf = rbind(inner_perf,perf)

          pmat = tar_dcvInner$all_model_preds
          pmat = pmat %>%
            group_by(model) %>%
            summarise(train_auc = as.numeric(pROC::auc(Status,!!as.name(classes[1])) ))
          pmat = data.frame(pmat)
          pmat = rbind(pmat,data.frame(model =tar_dcvInner$Performance$Approach,train_auc =  as.numeric(tar_dcvInner$Performance$AUC)))
          pmat$targetFeatures = tar_Features
          train_auc = rbind(train_auc,pmat)

        }))

        message(tar_Features)

      }

    }

  }

})

train_auc1 = train_auc %>%
  group_by(model,targetFeatures) %>%
  summarise_all(.funs = mean) %>%
  filter(str_detect(model,"DCV"))
ggplot(train_auc1,aes(targetFeatures,train_auc,col = model))+
  geom_point(size = 2)+
  theme_classic()+
  geom_line()+
  theme(legend.position = "top")+
  theme(axis.text = element_blank())




## Pre-Process
trainx = data.frame(fastImputeZeroes(ttData$train_Data,impFactor = ttData$imp_factor))
testx = data.frame(fastImputeZeroes(ttData$test_data,impFactor = ttData$imp_factor))

## compute log ratios
lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_train,trainx))
lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = ttData$y_test,testx))


## Apply targted feature selection method
compTime2 = system.time({
  cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                      includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                      rankOrder = F)
  tar_Features = 24
  tar_dcv = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,
                                  useRidgeWeights = useRidgeWeight,use_rfe = performRFE,scaledata = scale_data,
                                  testx = testx,
                                  dcv = cc.dcv$lrs,lrs.train = lrs.train,lrs.test = lrs.test,
                                  y_label = ttData$y_train,
                                  seed = sd,
                                  #ensemble = ensemble,
                                  y_test = ttData$y_test,
                                  tarFeatures = tar_Features,
                                  ts.id = ttData$test_ids,
                                  max_sparsity = max_sparsity
  )
})

pmat = tar_dcv$all_model_preds
pmat = pmat %>%
  group_by(model) %>%
  summarise(auc = as.numeric(pROC::auc(Status,!!as.name(classes[1])) ))
pmat = rbind(pmat,data.frame(model =tar_dcv$Performance$Approach,auc =  as.numeric(tar_dcv$Performance$AUC)))
tt = train_auc1 %>%
  filter(targetFeatures==tar_Features)
pmat = left_join(pmat,tt)

pmat1 = pmat %>%
  filter(model != "svmRadial3") %>%
  top_n(1,train_auc)

tt = data.frame(model = "DCV",auc =  pmat1$auc[1],targetFeatures = tar_Features,train_auc = pmat1$train_auc[1])
pmat = rbind(pmat,tt)

perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                  Dataset = f_name,Seed = sd,Fold = f,Approach = pmat$model,
                  AUC = pmat$auc,train_auc = pmat$train_auc,
                  number_parts = nparts,number_ratios = nparts ,comp_time = compTime[3],
                  base_dims = ncol(train.data))
benchmark = rbind(benchmark,perf)
