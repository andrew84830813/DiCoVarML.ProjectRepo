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



# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}

# Setup Cluster ---------------------------------------------------------------

## Detect Cores
# clus <- parallel::makeCluster(10)
# doParallel::registerDoParallel(clus)




# Read External Args ---------------------------------------------------------------

args = c(2,3,0)
args = commandArgs(trailingOnly = TRUE)
sd = as.numeric(args[1]) # random seed selection
f = as.numeric(args[2]) # random seed selection
permute_labels = as.logical(as.numeric(args[3])) #should be 0(False) or 1(True)


## set random seed
set.seed(sd)

## Name Iteration
f_name = paste0("NEC_seed",sd,"_permute",permute_labels,"_fold",f)

## read data
load("Output/microbiomeDB_dataAtLeast7DaysBeforeNEC.Rda")
df = exp$mbiome
md1 = exp$metadata

## Get fold matrix
folds_matrix = data.frame(read_csv(file = "Output/commonFolds_NEC_crossValidation.csv"))


## Partition data
allData = lodo_partition(df,dataset_labels = folds_matrix[,sd],sd)

##Define Out Data frame
benchmark = data.frame()

## DCV Parms
scale_data = T
performRFE = F
useRidgeWeight = F
min_connected = F
ensemble = c("ranger","xgbTree","xgbLinear")
ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")


message("\n\n``````` Start Seed:  ",sd,"````````````\n\n","fold",f)

## Compute Performance for the f - th fold

  ## Extract Test/Train Split
  ttData = DiCoVarML::extractTrainTestSplit(foldDataList = allData,
                                            fold = f,
                                            maxSparisty = .9,
                                            permLabels = permute_labels,
                                            extractTelAbunance = F)
  ##get train test partitions
  train.data = ttData$train_Data
  test.data = ttData$test_data
  ## Compute Total Parts
  number_parts = ncol(train.data);number_parts
  nsamps = nrow(train.data)
  table(ttData$y_test)
  classes = as.character(unique(ttData$y_train))

  ##get metadata
  train.md = data.frame(X1 = rownames(train.data))
  train.md = left_join(train.md,md1)
  test.md = data.frame(X1 = rownames(test.data))
  test.md = left_join(test.md,md1)
  y_label = ttData$y_train


 ## Pre-Process Metadata
  mm.char = train.md %>%
    select(Sex,
           Chorioamnionitis,
           Maternal.antepartum.antibiotics.administered,
           Baby.s.delivery.mode,
           Born.stunted,
           Host.diet)
  mm.char.test = subset(test.md,select = colnames(mm.char))

  ## convert variables to factors
  newdata.train <- mm.char %>% mutate_if(is.character,as.factor)
  newdata.test <- mm.char.test %>% mutate_if(is.character,as.factor)

  ##-----------------------------
  ## Add covairate
  # ## one hot features
  # dummy <- dummyVars(" ~ .", data=rbind(mm.char,mm.char.test))
  # newdata <- data.frame(predict(dummy, newdata = rbind(mm.char,mm.char.test)))
  # newdata.train = newdata[1:nrow(mm.char),]
  # newdata.test =  newdata[-1:-nrow(mm.char),]


  ## cont features
  train_metadata = data.frame(Age = train.md$Age.at.sample.collection..days.,
                             gesAge = train.md$Gestational.age.at.birth..days.,
                             host_weight = train.md$Host.weight..g.,
                             newdata.train
  )


  test.metadata = data.frame(Age = test.md$Age.at.sample.collection..days.,
                             gesAge = test.md$Gestational.age.at.birth..days.,
                             host_weight = test.md$Host.weight..g.,
                             newdata.test
  )



  # RFE - Mbiome --------------------------------------------------------------------
  message("Compute Performance - Mbiome Target DCV")
  suppressMessages(suppressWarnings({

    perc_totalParts2Keep = .75
    num_sets = 5

    base_dims = ncol(ttData$train_Data)
    max_parts = round(perc_totalParts2Keep*base_dims)
    sets = round(seq(5,max_parts,length.out = num_sets))


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

    pmat = data.frame(model =tar_dcv$Performance$Approach,auc =  as.numeric(tar_dcv$Performance$AUC))
    # pmat  =pmat %>%
    #   file(model=="DCV-ridgeEnsemble")

    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = pmat$model,
                      AUC = pmat$auc,#train_auc = pmat$train_auc,
                      number_parts = tar_Features,number_ratios = ncol(tar_dcv$weighted_features$train) ,comp_time = NA,
                      base_dims = ncol(train.data))


    benchmark = rbind(benchmark,perf)

  }))


  # GLM -  Meta --------------------------------------------------------------------
  message("Compute Performance - Metadata Alone GLM")
  suppressMessages(suppressWarnings({


    ## retrieve test and train data
    train.data = cbind(train_metadata)
    test.data = cbind(test.metadata)
    y_label = ttData$y_train
    y_test = ttData$y_test



    ## Apply Penalized Regression
    ## Tune Alpha
    type_family = if_else(length(classes)>2,"multinomial","binomial")
    mdl = glm(formula = y~x,data = data.frame(x =  train_metadata$Age,y = y_label),family = type_family)

    ## make predictions
    p = predict(mdl, newdata= data.frame(x =  test.metadata$Age),type = "response")
    if(type_family=="binomial"){
      mroc = pROC::roc(ttData$y_test,p)
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }else{
      ## multiclass
      mroc = pROC::multiclass.roc(ttData$y_test,p[,,1])
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }


    ## Save Performance
    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = "meta_dataGLM",
                      AUC = as.numeric(pROC::auc(mroc)),
                      number_parts = 1,number_ratios =NA ,comp_time = NA,
                      base_dims = 1)
    benchmark = rbind(benchmark,perf)


  }))



  # GLM - Mbiome + Meta --------------------------------------------------------------------
  message("Compute Performance - Mbiome + Metadata GLM")
  suppressMessages(suppressWarnings({


    ## retrieve test and train data
    train.data = cbind(tar_dcv$weighted_features$train ,train_metadata)
    test.data = cbind(tar_dcv$weighted_features$test,test.metadata)
    y_label = ttData$y_train
    y_test = ttData$y_test


    ## Apply Penalized Regression
    ## Tune Alpha
    type_family = if_else(length(classes)>2,"multinomial","binomial")


    ## Train GLM
    compTime2 = system.time({
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(train.data),y_label, standardize=T, alpha = 0,family=type_family)
    })
    if(type_family=="binomial"){
      features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
      features = features[-1,]
      features = features[abs(features)>0]
      length(features)
      c = as.matrix(coef(cv.clrlasso, s = "lambda.min"))[-1,]
    }else{
      features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
      feat.df = data.frame()
      for(o in 1:length(features)){
        ph = as.matrix(features[[o]])
        feat = ph[-1,]
        keep = feat[abs(feat)>0]
        feat.df = rbind(feat.df,data.frame(Ratio = names(keep),coef = as.numeric(keep)))
      }
      feat.df =feat.df %>%
        group_by(Ratio) %>%
        summarise(coef = sum(coef)) %>%
        filter(coef!=0)
    }

    ## make predictions
    p = predict(cv.clrlasso, newx = as.matrix(test.data), s = "lambda.min",type = "response")
    if(type_family=="binomial"){
      mroc = pROC::roc(ttData$y_test,p)
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }else{
      ## multiclass
      mroc = pROC::multiclass.roc(ttData$y_test,p[,,1])
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    }

    ## Compute Number of Part
    cn =names(c[abs(c)>0])
    n_ratios = length(cn)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)


    # ## Apply Penalized Regression
    # ## Tune Alpha
    # type_family = if_else(length(classes)>2,"multinomial","binomial")
    # mdl = glm(formula = y~.,data = data.frame(train.data,y = y_label),family = type_family)
    #
    # ## make predictions
    # p = predict(mdl, newdata= data.frame(test.data),type = "response")
    # if(type_family=="binomial"){
    #   mroc = pROC::roc(ttData$y_test,p)
    #   mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    # }else{
    #   ## multiclass
    #   mroc = pROC::multiclass.roc(ttData$y_test,p[,,1])
    #   mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
    # }



    ## Save Performance
    perf = data.frame(Scenario = if_else(permute_labels,"Permuted","Empirical"),
                      Dataset = f_name,Seed = sd,Fold = f,Approach = "mbiome + meta_dataGLM",
                      AUC = as.numeric(pROC::auc(mroc)),
                      number_parts = n_parts,number_ratios = ncol(train.data) ,comp_time = NA,
                      base_dims = ncol(train.data))
    benchmark = rbind(benchmark,perf)


  }))






## Write Performance Estimates
write_csv(x = benchmark,file = paste0("Results/NEC_caseStudyFinal/",f_name,".csv"))



