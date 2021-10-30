
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
library(ggplot2)




# Load Helper Functions  ---------------------------------------------------------------
fnames = dir("Helper_Functions/")
for(f in fnames){
  source(paste0("Helper_Functions/",f))
}


clus <- parallel::makeCluster(10)
doParallel::registerDoParallel(clus)


d = 50
s = .1
noise = .9
signal_center = 2
noise_center = 15
nsamps = 100

performance = data.frame()
for(d in c(50,100,150)){
  for(s in c(.05,.1,.2)){
    for(noise in seq(0,1,length.out = 10)){
      alr_dimensions = d
      signalPercent = s

      shiftparts = ceiling(alr_dimensions*signalPercent)
      ## sim parms
      if(signalPercent==0){
        cl_center = 0
        c2_center = 0
      }else{
        cl_center = signal_center
        c2_center = -signal_center
      }


      percentNoise = noise
      noisePercent = percentNoise *((1-signalPercent))
      noi_com = round(alr_dimensions*noisePercent)
      #percentNoise = noisePercent/(1-signalPercent)
      ## sim parms
      if(percentNoise==0){
        noiseC1 = 0
        noiseC2 = 0
      }else{
        noiseC1 = noise_center
        noiseC2 = -noise_center
      }


      ## SIm Dataset
      s1 = simDataset(n = nsamps,dims = alr_dimensions,num_shiftedComponents = shiftparts,
                      num_noiseComponents = noi_com,
                      offset_c1 = cl_center,
                      offset_c2 = c2_center,
                      noiseCenter1 = noiseC1,convertCounts = F )
      # s22 = s1[,-1]
      # s22[s22<1e-11] = 0
      # s22 = fastImputeZeroes(s22)
      # s1[,-1]=s22
      s1 = data.frame(Dataset = "Train",s1)
      s2 = simDataset(n = nsamps,dims = alr_dimensions,num_shiftedComponents = shiftparts,
                      num_noiseComponents = noi_com,
                      offset_c1 = cl_center,
                      offset_c2 = c2_center,
                      noiseCenter1 = noiseC2 ,convertCounts = F)
      # s22 = s2[,-1]
      # s22[s22<1e-11] = 0
      # s22 = fastImputeZeroes(s22)
      # s2[,-1]=s22
      s2 = data.frame(Dataset = "Test",s2)
      df = rbind(s1,s2)
      tbl = data.frame(df)
      tbl$Status = str_replace_all(string = df$Status,pattern = "S",replacement = "C" )



      #
      # Partition Data
      df = tbl %>%
        filter(Dataset == "Train") %>%
        select(-Dataset)
      df$Status = factor(df$Status)

      ## Daatset 2
      df2 = tbl %>%
        filter(Dataset == "Test") %>%
        select(-Dataset)
      df2$Status = factor(df2$Status)
      classes = as.character(unique(df$Status))



      pc = prcomp((tbl[,-2:-1]))
      pc.df = data.frame(tbl[,1:2],pc$x)
      ggplot(pc.df,aes(PC1,PC2,col =Status,shape = Dataset))+
        geom_point()


      classes = unique(df$Status)
      parts = ncol(df[,-1])
      ovr.perf = data.frame()

      # Relative Abundance --------------------------------------------------------------

      ## with feature selection
      testx = (df2[,-1])
      trainx = (df[,-1])


      ## Train Model
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha = 1,family=type_family)
      features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
      features = features[-1,]
      features = features[abs(features)>0]
      trainx = data.frame((subset(trainx,select = names(features))))
      testx = data.frame((subset(testx,select = names(features))))

      ## if no features are selecting the AUC is .5 since nothing could be evaluated
      ## this implied LASSO failed on the data
      if(length(features)<=1){
        accc = .5
        taccc = .5
      }else{
        type_family = if_else(length(classes)>2,"multinomial","binomial")
        infld = 2
        flds = caret::createFolds(y = df[,1],k = infld,list = F)
        compTime2 = system.time({
          aseq = seq(0,1,length.out = 10)
          min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{

            aucc = c()
            for(f in 1:infld){
              bool  = flds==f
              compTime2 = system.time({
                cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx[bool,]),df$Status[bool],
                                                 standardize=T, alpha=a,family=type_family)
              })

              ## make predictions
              p = predict(cv.clrlasso, newx = as.matrix(trainx[!bool,]), s = "lambda.min",type = "response")
              if(type_family=="binomial"){
                mroc = pROC::roc(df$Status[!bool],p)
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }else{
                ## multiclass
                mroc = pROC::multiclass.roc(df$Status[!bool],p[,,1])
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }
              aucc = c(aucc,as.numeric(mroc.dcvlasso))
            }
            data.frame(a,auc = mean(aucc))
          }
        })
        min_dev = min_dev %>%
          arrange(desc(auc))
        cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=min_dev$a[1],family="binomial")
        p = predict(cv.clrlasso, newx = as.matrix(testx), s = "lambda.min",type = "response")
        mroc = pROC::roc(df2$Status,p)
        accc =as.numeric(pROC::auc(mroc))
        taccc = min_dev$auc[1]
      }


      perf = data.frame(Seed = sd,
                        Dims = parts,
                        Model = "LASSO",Type = "Relative_Abudance",
                        train_auc = taccc,
                        test_auc = accc)
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
                        Dims = parts,
                        Model = "RF",Type = "Relative_Abudance",
                        train_auc =  m1$performance$TrainAUC,
                        test_auc = as.numeric(pROC::auc(mroc.all.RA)))
      ovr.perf = rbind(ovr.perf,perf)








      # CLR with feature selection --------------------------------------------------------------

      ## with feature selection
      testx = data.frame(clr(df2[,-1]))
      trainx = data.frame(clr(df[,-1]))



      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=1,family="binomial")
      features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
      features = features[-1,]
      features = features[abs(features)>0]
      trainx = data.frame(clr(subset(df,select = names(features))))
      testx = data.frame(clr(subset(df2,select = names(features))))


      ## if no features are selecting the AUC is .5 since nothing could be evaluated
      ## this implied LASSO failed on the data
      if(length(features)<=1){
        accc = .5
        taccc = .5
      }else{
        type_family = if_else(length(classes)>2,"multinomial","binomial")
        infld = 2
        flds = caret::createFolds(y = df[,1],k = infld,list = F)
        compTime2 = system.time({
          aseq = seq(0,1,length.out = 10)
          min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{

            aucc = c()
            for(f in 1:infld){
              bool  = flds==f
              compTime2 = system.time({
                cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx[bool,]),df$Status[bool],
                                                 standardize=T, alpha=a,family=type_family)
              })

              ## make predictions
              p = predict(cv.clrlasso, newx = as.matrix(trainx[!bool,]), s = "lambda.min",type = "response")
              if(type_family=="binomial"){
                mroc = pROC::roc(df$Status[!bool],p)
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }else{
                ## multiclass
                mroc = pROC::multiclass.roc(df$Status[!bool],p[,,1])
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }
              aucc = c(aucc,as.numeric(mroc.dcvlasso))
            }
            data.frame(a,auc = mean(aucc))
          }
        })
        min_dev = min_dev %>%
          arrange(desc(auc))
        ## Train Model
        cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha = min_dev$a[1],family=type_family)
        p = predict(cv.clrlasso, newx = as.matrix(testx), s = "lambda.min",type = "response")
        mroc = pROC::roc(df2$Status,p)
        accc =as.numeric(pROC::auc(mroc))
        taccc = min_dev$auc[1]
      }


      perf = data.frame(Seed = sd,
                        Dims = parts,
                        Model = "LASSO",Type = "CLR",
                        train_auc = taccc,
                        test_auc = accc)
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
                        Dims = parts,
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


      ## Train Model
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha = 1,family=type_family)
      features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
      features = features[-1,]
      features = features[abs(features)>0]
      trainx = data.frame((subset(trainx,select = names(features))))
      testx = data.frame((subset(testx,select = names(features))))





      if(length(features)<=1){
        accc = .5
        taccc = .5
      }else{
        type_family = if_else(length(classes)>2,"multinomial","binomial")
        infld = 2
        flds = caret::createFolds(y = df[,1],k = infld,list = F)
        compTime2 = system.time({
          aseq = seq(0,1,length.out = 10)
          min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{

            aucc = c()
            for(f in 1:infld){
              bool  = flds==f
              compTime2 = system.time({
                cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx[bool,]),df$Status[bool],
                                                 standardize=T, alpha=a,family=type_family)
              })

              ## make predictions
              p = predict(cv.clrlasso, newx = as.matrix(trainx[!bool,]), s = "lambda.min",type = "response")
              if(type_family=="binomial"){
                mroc = pROC::roc(df$Status[!bool],p)
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }else{
                ## multiclass
                mroc = pROC::multiclass.roc(df$Status[!bool],p[,,1])
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }
              aucc = c(aucc,as.numeric(mroc.dcvlasso))
            }
            data.frame(a,auc = mean(aucc))
          }
        })
        min_dev = min_dev %>%
          arrange(desc(auc))
        cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=min_dev$a[1],family="binomial")
        p = predict(cv.clrlasso, newx = as.matrix(testx), s = "lambda.min",type = "response")
        mroc = pROC::roc(df2$Status,p)
        accc =as.numeric(pROC::auc(mroc))
        taccc = min_dev$auc[1]
      }




      perf = data.frame(Seed = sd,
                        Dims = parts,
                        Model = "LASSO",Type = "ALR",
                        train_auc = taccc,
                        test_auc =accc)
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
                        Dims = parts,
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


      ## Train Model
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha = 1,family=type_family)
      features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
      features = features[-1,]
      features = features[abs(features)>0]
      trainx = data.frame((subset(trainx,select = names(features))))
      testx = data.frame((subset(testx,select = names(features))))


      if(length(features)<=1){
        accc = .5
        taccc = .5
      }else{
        type_family = if_else(length(classes)>2,"multinomial","binomial")
        infld = 2
        flds = caret::createFolds(y = df[,1],k = infld,list = F)
        compTime2 = system.time({
          aseq = seq(0,1,length.out = 10)
          min_dev =  foreach(a = aseq,.combine = rbind)%dopar%{

            aucc = c()
            for(f in 1:infld){
              bool  = flds==f
              compTime2 = system.time({
                cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx[bool,]),df$Status[bool],
                                                 standardize=T, alpha=a,family=type_family)
              })

              ## make predictions
              p = predict(cv.clrlasso, newx = as.matrix(trainx[!bool,]), s = "lambda.min",type = "response")
              if(type_family=="binomial"){
                mroc = pROC::roc(df$Status[!bool],p)
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }else{
                ## multiclass
                mroc = pROC::multiclass.roc(df$Status[!bool],p[,,1])
                mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
              }
              aucc = c(aucc,as.numeric(mroc.dcvlasso))
            }
            data.frame(a,auc = mean(aucc))
          }
        })
        min_dev = min_dev %>%
          arrange(desc(auc))
        cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=min_dev$a[1],family="binomial")
        p = predict(cv.clrlasso, newx = as.matrix(testx), s = "lambda.min",type = "response")
        mroc = pROC::roc(df2$Status,p)
        accc =as.numeric(pROC::auc(mroc))
        taccc = min_dev$auc[1]
      }








      perf = data.frame(Seed = sd,
                        Dims = parts,
                        Model = "LASSO",Type = "PLR",
                        train_auc = taccc,
                        test_auc = accc)
      ovr.perf = rbind(ovr.perf,perf)




      ## without feature selection
      library(igraph)
      ### PLR Approach
      plr = calcLogRatio(df)
      plr2 = calcLogRatio(df2)
      ## with feature selection
      testx = (plr2[,-1])
      trainx = (plr[,-1])

      # cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = plr,
      #                                     includeInfoGain = T, nfolds = 1, numRepeats = 1,
      #                                     rankOrder = F)
      # trainx = diffCompVarRcpp::mstAll(featMatrix = plr,dcvRanking = cc.dcv$lrs)
      # testx = subset(testx,select = colnames(trainx))


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
                        Dims = parts,
                        Model = "RF",Type = "PLR",
                        train_auc =  m1$performance$TrainAUC,
                        test_auc = as.numeric(pROC::auc(mroc.all.PLR)))
      ovr.perf = rbind(ovr.perf,perf)




      ovr.perf$percentNoise = percentNoise
      ovr.perf$signalPercent = signalPercent


      performance = rbind(performance,ovr.perf)
      View(performance)

    }
  }
}

ggplot(performance,aes(percentNoise,test_auc,col =Type ))+
  geom_point()+
  geom_line()+
  facet_grid(.~Model)+
  theme_bw()



