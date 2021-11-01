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
library(easyCODA)



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

args = c(4,1,2,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
dataset = as.numeric(args[2])
featSetSize = as.numeric(args[3]) ## should be 1 or 2
permute_labels = as.logical(as.numeric(args[4])) #should be 0(False) or 1(True)





# Load Data ---------------------------------------------------------------

switch(dataset,

       {
         load("Output/HIV.rda");
         df = data.frame(Status = HIV$HIV_Status,HIV[,1:60]);
         f_name = "selbal_HIV_16s";
       },

       {
         load("Output/Crohn.rda")
         df = data.frame(Status = Crohn$y,Crohn[,1:48]);
         df$Status = dplyr::if_else(df$Status=="CD","pos","neg");
         #preprocessing
         #remove g__ and f__ in taxa names
         taxa = data.frame(Taxa = colnames(df[,-1]));
         taxa$taxa = stringr::str_replace(taxa$Taxa,pattern = "\\.",replace = "");
         taxa = tidyr::separate(taxa,col = 1,into = c("des","taxa"),sep = "__");
         colnames(df) = c("Status",taxa$taxa);
         rs= rowSums(df[,-1])
         bool = rs>1e3;
         df = df[bool,];
         f_name = "selbal_Crohns_16s";
       },

       {
         ## nfald
         df = read.csv(file = "Output/16S-11635_NAFLD.csv") ## good;
         rs = rowSums(df[,-1]);
         bool = rs>5e3;
         df = df[bool,];
         f_name = "qitta_NAFLD_16s";

       },

       {
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

       },

       {
         load("Output/curatedMetaGenome_RubelMA_2020.Rda"); ## good
         keep = c("control","STH");
         f_name = "cmg_RubelMA-2020_STH"
       },

       {
         load("Output/curatedMetaGenome_ZhuF_2020.Rda") ## good
         keep = c("control","schizofrenia");
         f_name = "cmg_ZhuF-2020_schizo"

       },

       {
         load("Output/curatedMetaGenome_QinN_2014.Rda")
         keep = c("control","cirrhosis");
         f_name = "cmg_QinN-2014_cirr"
       },

       {
         load("Output/curatedMetaGenome_FengQ_2015.Rda");
         keep = c("control","CRC");
         f_name = "cmg_FengQ-2015_crc"

       },



       {
         load("Output/curatedMetaGenome_WirbelJ_2018.Rda");
         keep = c("control","CRC");
         f_name = "cmg_WirbelJ-2018_crc"
       },


       {
         load("Output/curatedMetaGenome_ZellerG_2014.Rda");
         f_name = "cmg_ZellerG_2014_crc";
         keep = c("control","CRC")
       }

)






if(str_detect(f_name,"cmg")){
  ## Process Data
  df = expResults$taxa.df
  df = data.frame(sample_id = rownames(df),df)
  df$Status = str_replace_all(df$Status,pattern = "\\/","_")
  rs = rowSums(df[,-2:-1])
  bool = rs>10e6
  df = df[bool,]
  df = df[df$Status %in% keep,]
  md = expResults$metadata
  ## by taxa level
  cname.table = expResults$taxa_info [,-1]
  taxa_levels =colnames(cname.table)
  cname.table = data.frame(ID = cname.table$species,cname.table)
  rdphf = tidyr::gather(df,"ID","Count",3:ncol(df))
  cname.table$ID = stringr::str_replace_all(cname.table$ID,pattern = " ",replacement = "\\.")
  rdphf$ID = stringr::str_replace_all(rdphf$ID,pattern = " ",replacement = "\\.")
  rdphf = left_join(rdphf,cname.table)
  level.list = list()
  for(taxa_level in c("species")){
    rdf = subset(rdphf,select = c("sample_id","Status",taxa_level,"Count"))
    rdf = rdf %>%
      dplyr::group_by(sample_id,get(taxa_level),Status) %>%
      dplyr::summarise(count = sum(Count) ) %>%
      rename(tlevel = "get(taxa_level)" ) %>%
      spread("tlevel","count",fill = 0)
    rdf = data.frame(rdf[,-1],row.names = rdf$sample_id )
    # rdf = rdf[rdf$Status %in% c("control","CRC"),]
    colnames(rdf)[-1] = paste0(str_sub(taxa_level,1,1),"__",colnames(rdf)[-1])
    level.list[[taxa_level]] = rdf
  }
  df = level.list$species
  df = df[!is.na(df$Status),]
}


f_name = paste0(f_name,"_seed",sd,"_featSize",featSetSize)
message("\n",f_name,"\n")
mxSparsePercent = .9
benchmark = data.frame()
seed_ = sd

#for(sd in 1:5){
set.seed(sd)

## Partition Data
k_fold = 2
overll_folds = caret::createFolds(df[,1],k =k_fold,list = F)
allData = lodo_partition(df,dataset_labels = overll_folds,sd)

## DCV Parms
scale_data = T
performRFE = F
useRidgeWeight = F
min_connected = F
ensemble = c("ranger","xgbTree","xgbLinear")
ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
lasso_alpha = 0 # normally lasso_alpha = 0


## Define Sets
perc_totalParts2Keep = .75
num_sets = 8
got = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                          fold = 1,
                                          permLabels = permute_labels,
                                          maxSparisty = .9,
                                          extractTelAbunance = F)
##get train test partitions
got1 = got$train_Data
got = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                       fold = 2,
                                       permLabels = permute_labels,
                                       maxSparisty = .9,
                                       extractTelAbunance = F)
##get train test partitions
got2 = got$train_Data
base_dims = min(c(ncol(got1),ncol(got2)))
max_parts = round(perc_totalParts2Keep*base_dims)
sets = unique(round(seq(5,max_parts,length.out = num_sets)))

if(featSetSize == 1){
  feat_s = 1:4
}else{
  feat_s = 5:8
}


for(tarfeatureSet in feat_s){
  target_features = sets[tarfeatureSet]

  message("\n\n``````` Start Seed:  ",sd,"````````````\n\n")

  #  Perform 2-fold cross validation -------------------------------

  for(f in 1:k_fold){


    ## Extract Test/Train Split
    ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                              fold = f,
                                              permLabels = permute_labels,
                                              maxSparisty = .9,
                                              extractTelAbunance = F)
    ##get train test partitions
    train.data = ttData$train_Data
    test.data = ttData$test_data
    ## Compute Total Parts
    number_parts = ncol(train.data);number_parts
    nsamps = nrow(train.data)
    table(ttData$y_test)
    classes = as.character(unique(ttData$y_train))




    # ALR/PLR ---------------------------------------------------------------------

    library(easyCODA)



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

            ## find best ref
            ff = FINDALR(trainx)

            ## with feature selection
            ref = ff$procrust.ref
            rt =  paste0(colnames(trainx)[-ref],"___",colnames(trainx)[ref])

            lrs.train = data.frame(Status = innerFold$y_train,getLogratioFromList(Ratio = rt,trainx,"test"))
            lrs.test = data.frame(Status = innerFold$y_test, getLogratioFromList(Ratio = rt,testx,"test"))


            cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                                includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                                rankOrder = F)

          }))



          for(tar_Features in target_features){

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




    ## Pre-Process
    trainx = data.frame(fastImputeZeroes(ttData$train_Data,impFactor = ttData$imp_factor))
    testx = data.frame(fastImputeZeroes(ttData$test_data,impFactor = ttData$imp_factor))

    ## find best ref
    library(easyCODA)
    ff = FINDALR(trainx)

    ## with feature selection
    ref = ff$procrust.ref
    rt =  paste0(colnames(trainx)[-ref],"___",colnames(trainx)[ref])

    lrs.train = data.frame(Status = ttData$y_train,getLogratioFromList(Ratio = rt,trainx,"test"))
    lrs.test = data.frame(Status = ttData$y_test, getLogratioFromList(Ratio = rt,testx,"test"))

    ## Apply targted feature selection method
    compTime2 = system.time({
      cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                          includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                          rankOrder = F)
      tar_Features = train_auc2$targetFeatures[which.max(train_auc2$train_auc)]
      tar_dcv = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,alpha_ = lasso_alpha,
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


    ## test either ensemble or regression model
    pmat1 = pmat %>%
      filter(str_detect(model,"DCV")) %>%
      top_n(1,train_auc)

    tt = data.frame(model = "ALR/PLR_DCV",auc =  pmat1$auc[1],targetFeatures = tar_Features,train_auc = pmat1$train_auc[1])
    pmat = rbind(pmat,tt)

    cn  =colnames(tar_dcv$weighted_features$train)
    nparts = (str_split(cn,pattern = "___",simplify = T))
    nparts = n_distinct(c(nparts[,1],nparts[,2]))

    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = paste0("ALR/PLR_",pmat$model),
                      AUC = pmat$auc,train_auc = pmat$train_auc,
                      number_parts = nparts,number_ratios = length(cn) ,comp_time = compTime2[3],
                      base_dims = ncol(train.data))
    perf$feat_set =tarfeatureSet

    benchmark = rbind(benchmark,perf)



    # ALR ---------------------------------------------------------------------



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

            ## find best ref
            ff = FINDALR(trainx)

            ## with feature selection
            ref = ff$procrust.ref
            rt =  paste0(colnames(trainx)[-ref],"___",colnames(trainx)[ref])

            lrs.train = data.frame(Status = innerFold$y_train,getLogratioFromList(Ratio = rt,trainx,"test"))
            lrs.test = data.frame(Status = innerFold$y_test, getLogratioFromList(Ratio = rt,testx,"test"))


            cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                                includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                                rankOrder = F)

          }))



          for(tar_Features in target_features){

            suppressMessages(suppressWarnings({

              tar_dcvInner = targeted_dcvSelection.alr(trainx = trainx,minConnected = min_connected,
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




    ## Pre-Process
    trainx = data.frame(fastImputeZeroes(ttData$train_Data,impFactor = ttData$imp_factor))
    testx = data.frame(fastImputeZeroes(ttData$test_data,impFactor = ttData$imp_factor))

    ## find best ref
    library(easyCODA)
    ff = FINDALR(trainx)

    ## with feature selection
    ref = ff$procrust.ref
    rt =  paste0(colnames(trainx)[-ref],"___",colnames(trainx)[ref])

    lrs.train = data.frame(Status = ttData$y_train,getLogratioFromList(Ratio = rt,trainx,"test"))
    lrs.test = data.frame(Status = ttData$y_test, getLogratioFromList(Ratio = rt,testx,"test"))

    ## Apply targted feature selection method
    compTime2 = system.time({
      cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                          includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                          rankOrder = F)
      tar_Features = train_auc2$targetFeatures[which.max(train_auc2$train_auc)]
      tar_dcv = targeted_dcvSelection.alr(trainx = trainx,minConnected = min_connected,alpha_ = lasso_alpha,
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

    ## test either ensemble or regression model
    pmat1 = pmat %>%
      filter(str_detect(model,"DCV")) %>%
      top_n(1,train_auc)

    tt = data.frame(model = "ALR_DCV",auc =  pmat1$auc[1],targetFeatures = tar_Features,train_auc = pmat1$train_auc[1])
    pmat = rbind(pmat,tt)

    cn  =colnames(tar_dcv$weighted_features$train)
    nparts = (str_split(cn,pattern = "___",simplify = T))
    nparts = n_distinct(c(nparts[,1],nparts[,2]))

    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = paste0("ALR_",pmat$model),
                      AUC = pmat$auc,train_auc = pmat$train_auc,
                      number_parts = nparts,number_ratios = length(cn) ,comp_time = compTime2[3],
                      base_dims = ncol(train.data))
    perf$feat_set =tarfeatureSet

    benchmark = rbind(benchmark,perf)




    # DCV --------------------------------------------------------------



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



          for(tar_Features in target_features){

            suppressMessages(suppressWarnings({

              tar_dcvInner = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,alpha_ = lasso_alpha,
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
      tar_Features = train_auc2$targetFeatures[which.max(train_auc2$train_auc)]
      tar_dcv = targeted_dcvSelection(trainx = trainx,minConnected = min_connected,alpha_ = lasso_alpha,
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

    ## test either ensemble or regression model
    pmat1 = pmat %>%
      filter(str_detect(model,"DCV")) %>%
      top_n(1,train_auc)

    tt = data.frame(model = "PLR_DCV",auc =  pmat1$auc[1],targetFeatures = tar_Features,train_auc = pmat1$train_auc[1])
    pmat = rbind(pmat,tt)

    cn  =colnames(tar_dcv$weighted_features$train)
    nparts = (str_split(cn,pattern = "___",simplify = T))
    nparts = n_distinct(c(nparts[,1],nparts[,2]))

    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = paste0("PLR_",pmat$model),
                      AUC = pmat$auc,train_auc = pmat$train_auc,
                      number_parts = nparts,number_ratios = length(cn) ,comp_time = compTime2[3],
                      base_dims = ncol(train.data))
    perf$feat_set =tarfeatureSet
    benchmark = rbind(benchmark,perf)



  }
}





#}



write_csv(x = benchmark,file = paste0("Results/tarSelectionComp/",f_name,".csv"))

