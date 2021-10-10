
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





# Sim Dataset -------------------------------------------------------------

## Shift Data Within Simplex
## D1
set.seed(08272008)
comp =  compositions::rDirichlet.acomp(2,alpha = 200*c(4,2,1));plot(comp);comp = data.frame(comp)[1,]


s1 = data.frame(Dataset = "D1",Status = "C1",comp)

cnames = colnames(comp)
rwnames = rownames(comp)
madison = data.frame()
nsamps = nrow(comp)
alphaRow = which(a==1.2)
lrs.df = data.frame()
a = round(seq(-10,10,by = .1),1)  #<===== Input =======


## Define Decision Boundary
for(j in 1:nsamps){
  x1 = comp[j,]
  c.df = data.frame()
  label = c()
  i = 1
  for(A in a){
    ph = data.frame((clo(x1^A)))
    c.df = rbind(c.df,ph)
    label[i] = "black"
    i = i+1
  }
  # lrs = calcLogRatio(data.frame(Status = paste0("p",j),c.df))
  # lrs.df = rbind(lrs.df,lrs)#

  #madison = rbind(madison,data.frame(c.df[alphaRow,]))
  madison = rbind(madison,data.frame(Dataset = j,Status = a,c.df))
  #ADD Indv Points
  colnames(c.df) = cnames
  rownames(c.df) = paste("A",as.character(round(a,digits = 3)),sep = "_") #with alpha values
  labels = str_split(rownames(c.df),pattern = "_",simplify = T,n = 2)
  labels = as.factor(labels[,1])
  #plot
  ac = acomp(c.df)
  if(j==1){
    plot(ac,col = labels,"l")
    plot(acomp(comp),add = T,pch = 18,col = "red")
    #plot(compMean,col = "steelblue",add = T,pch = 8,cex = 2)
    #plot(acomp(baryCenter),col = "green",add = T,cex= 2,pch = 16)
  }else{
    plot(ac,col = labels,"l",add = T)
  }
}


md.df = calcLogRatio(df = madison[,-1])

x = .5*( (2*madison$V2+madison$V3) / (rowSums(madison[,-2:-1])) )
y = (sqrt(3)/2)*( (madison$V3) / (rowSums(madison[,-2:-1])) )
dec_bound = data.frame(x,y,prob = 0)
dec_bound1 = rbind(dec_bound[-1,],dec_bound[nrow(dec_bound),])
colnames(dec_bound1) = c("x1","y1","prob1")
dec_bound = cbind.data.frame(dec_bound,dec_bound1)


alpha = runif(n = 10,min = 0.05,max = .125)
c = simplexDataAugmentation::linearShift(data.frame(madison[,-2:-1]),
                                         a = alpha,alpha_n = alpha,pertubationDirection = clo(c(100,1,1)),
                                         directionType = "centroid",evComp = 1)
c = c$linAdjustedData
s2 = data.frame(Dataset = "d3",Status = "C2",sample_n(c,size = 200,replace = F))



alpha = .1
c = simplexDataAugmentation::linearShift(data.frame(s2[,-2:-1]),
                                         a = alpha,alpha_n = alpha,pertubationDirection = clo(c(1,100,1)),
                                         directionType = "centroid",evComp = 1)
c = c$linAdjustedData
s3 = data.frame(Dataset = "d3",Status = "C1",c)



tbl = rbind(s2,s3)
tbl$Dataset = if_else(tbl$V1>.33,"Test","Train")

## test points
x = .5*( (2*tbl$V2+tbl$V3) / (rowSums(tbl[,-2:-1])) )
y = (sqrt(3)/2)*( (tbl$V3) / (rowSums(tbl[,-2:-1])) )
points = data.frame(x,y,prob = 0,Class = tbl$Status,Partition = tbl$Dataset)





## define grid

# Define Test Grid
ph = data.frame()
incr  = 0.005
test_grid = foreach(a =seq(0,1,by = incr),.combine = rbind)%dopar%{
  leftover = 1-a
  b = seq(0,leftover,by =incr)
  c = 1-(b+a)
  data.frame(a,b,c)
}
##convert to cartesian
x = .5*( (2*test_grid$b+test_grid$c) / (rowSums(test_grid)) )
y = (sqrt(3)/2)*( (test_grid$c) / (rowSums(test_grid)) )




## Visulaize Syn Data
border = data.frame(start_x = c(0,.5,1),start_y = c(0,sqrt(3)/2,0),
                    end_x = c(.5,1,0),
                    end_yy = c(sqrt(3)/2,0,0),prob = 0)
pdf(file = "Figures/fig2_simplexlinear_data.pdf",width = 2.25 ,height = 2.25)
  ggplot(dec_bound,aes(x,y))+
    geom_point(data = points,aes(x,y,col = Class,shape = Partition),size = 2,alpha = .5)+
    scale_shape_manual(values = c(3,16))+
    scale_color_manual(values = c("red","blue"))+
    geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
    geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
    theme_classic()+
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          legend.margin=margin(-5,-10,-10,-10),
          #axis.text = element_text(size = 8),
          #panel.grid = element_blank(),
          legend.key.size = unit(.15,units = "in"),
          legend.text = element_text(size = 8))

dev.off()


## Partition Data
df = tbl %>%
  filter(Dataset == "Train") %>%
  select(-Dataset)
df$Status = factor(df$Status)

## Daatset 2
df2 = tbl %>%
  filter(Dataset == "Test") %>%
  select(-Dataset)
df2$Status = factor(df2$Status)







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


## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
## train cv metrics
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glm",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::roc(df2$Status,p)


##compute probs on test grid
test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
test_ran = clo(subset(test_ran,select = names(features)))
p = predict.glm(cv.clrlasso, newdata = data.frame(test_ran), type = "response")
cc = data.frame(x,y,prob = p)


pdf(file = "Figures/fig2_simplexlinear_glmnet.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())+
  theme(legend.position = "right",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 8,face = "bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        #axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
  )
dev.off()



## without feature selection
testx = (df2[,-1])
trainx = (df[,-1])

m1 = trainML_Models(trainLRs = df[,-1],testLRs = df2[,-1],
                    ytrain = df[,1],y_test = df2[,1],models = "ranger",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])
r = pROC::roc(df2$Status,m1$predictionMatrix$C1)
tt = data.frame(pROC::coords(r)) %>%
  filter(!is.infinite(threshold))

test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )
cc = data.frame(x,y,prob = preds$C1)

pdf(file = "Figures/fig2_simplexLinear_ranger.pdf",width = 3 ,height = 3)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())
dev.off()





# CLR with feature selection --------------------------------------------------------------

## with feature selection
testx = clr(df2[,-1])
trainx = clr(df[,-1])

cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=F, alpha=1,family="binomial",type.measure = "auc")
features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
features = features[-1,]
features = features[abs(features)>0]
trainx = clr(subset(df,select = names(features)))
testx = clr(subset(df2,select = names(features)))

## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
## Train CV metic
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glm",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::roc(df2$Status,p)



##compute probs on test grid
test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
test_ran = clr(subset(test_ran,select = names(features)))
p = predict.glm(cv.clrlasso, newdata = data.frame(test_ran), type = "response")
cc = data.frame(x,y,prob = p)


pdf(file = "Figures/fig2_simplexLinear_clrGLMNET.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())+
  theme(legend.position = "right",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 8,face = "bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        #axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
  )
dev.off()




## without feature selection
testx = data.frame(clr(df2[,-1]))
trainx = data.frame(clr(df[,-1]))
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "ranger",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)

m1$performance
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])


test_ran = data.frame(clr(test_grid))
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )
cc = data.frame(x,y,prob = preds$C1)


pdf(file = "Figures/fig2_simplexLinear_ranger.pdf",width = 3 ,height = 3)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())
dev.off()




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

## Train CV metric
m1 = trainML_Models(trainLRs = trainx,testLRs = testx,
                    ytrain = df[,1],y_test = df2[,1],models = "glm",
                    mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                    numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
m1$performance


## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::roc(df2$Status,p)

##compute probs on test grid
test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
pp  =calcLogRatio(data.frame(Status = "Test",fastImputeZeroes(test_ran)))
test_ran =subset(pp,select = names(features))
p = predict.glm(cv.clrlasso, newdata = data.frame(test_ran), type = "response")
cc = data.frame(x,y,prob = p)


pdf(file = "Figures/fig2_simplexLinear_plrGLMNET.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())+
  theme(legend.position = "right",panel.grid = element_blank(),
        plot.title = element_text(size = 8,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        #axis.title = element_text(size = 8,face = "bold"),
        #axis.title.y = element_blank(),
        #axis.title.x = element_text(size = 8,face = "bold"),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        #axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")+
  )
dev.off()



## without feature selection

### PLR Approach
plr = calcLogRatio(df)
plr2 = calcLogRatio(df2)
## with feature selection
testx = (plr2[,-1])
trainx = (plr[,-1])

ensemble = c("ranger","pls","svmRadial","glmnet","rangerE")
m1.plr = trainML_Models(trainLRs = trainx,testLRs = testx,testIDs = data.frame(ID = 1:nrow(df2),Status =df2$Status),
                        ytrain = df[,1],y_test = df2[,1],models = ensemble,
                        mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                        numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
pmat = m1.plr$predictionMatrix
pmat = pmat %>%
  group_by(ID,Status) %>%
  dplyr::select(-model) %>%
  summarise_all(.funs = mean)
pmat = data.frame(pmat)
classes = as.character(unique(df$Status))
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])
r = pROC::roc(pmat$Status,pmat$C1)
tt = data.frame(pROC::coords(r)) %>%
  filter(!is.infinite(threshold))


## compute ratios
test_ran = test_grid ## <-------------- Input
colnames(test_ran) = colnames(df[,-1])
tt.lr = calcLogRatio(data.frame(Status = "test",test_ran))

pmat = data.frame()
for(i in 1:length( m1.plr$models)){
  testh = m1.plr$models[[i]]
  ph = data.frame(ID = 1:nrow(tt.lr),
                  predict.train(testh,tt.lr[,-1],type = "prob" ),model = names(m1.plr$models)[i])
  pmat = rbind(pmat,ph)
}

pmat = pmat %>%
  group_by(ID) %>%
  dplyr::select(-model) %>%
  summarise_all(.funs = mean)




cc = data.frame(x,y,prob = preds$C1)
ggplot(cc,aes(x,y,col = prob))+
  geom_point()+
  #scale_colour_gradient2(low = "red",high = "blue",midpoint = .5)+
  scale_color_distiller(palette = "RdBu")+
  theme_classic()+
  geom_segment(aes(x = x, y = y, xend = x1, yend = y1), data = dec_bound,size = 1,col  ="black")+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())

