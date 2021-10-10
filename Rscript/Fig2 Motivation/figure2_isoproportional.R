
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


## sim parms
cl_center = -.5
c2_center = .5
noiseC1 = 2
noiseC2 = -2
nsamps = 100
alr_dimensions = 2

## SIm Dataset
s1 = simDataset(n = nsamps,dims = alr_dimensions,
                offset_c1 = cl_center,
                offset_c2 = c2_center,
                noiseCenter1 = noiseC1,convertCounts = F )
s1 = data.frame(Dataset = "Train",s1)
s2 = simDataset(n = nsamps,dims = alr_dimensions,
                offset_c1 = cl_center,
                offset_c2 = c2_center,
                noiseCenter1 = noiseC2 ,convertCounts = F)
s2 = data.frame(Dataset = "Test",s2)
df = rbind(s1,s2)
tbl = data.frame(df)
tbl$Status = str_replace_all(string = df$Status,pattern = "S",replacement = "C" )
x = .5*( (2*tbl$V2+tbl$V3) / (rowSums(tbl[,-2:-1])) )
y = (sqrt(3)/2)*( (tbl$V3) / (rowSums(tbl[,-2:-1])) )
points = data.frame(x,y,prob = 0,Class = tbl$Status,Partition = tbl$Dataset)


## Define Decision Boundary
LR = mean(log(tbl$V1/tbl$V3))
xx = data.frame(Dataset = "",Status = "Expected Decison Bound",a = seq(.01,1,length.out = 1000) )
xx$c = xx$a/exp(LR)
xx = xx %>%
  dplyr::mutate(S = a+c) %>%
  filter(S<=1) %>%
  dplyr::select(-S) %>%
  dplyr::mutate(b = 1 - (a + c)) %>%
  dplyr::select(a,b,c)
x = .5*( (2*xx$b+xx$c) / (rowSums(xx)) )
y = (sqrt(3)/2)*( (xx$c) / (rowSums(xx)) )
dec_bound = data.frame(x,y,prob = 0)



# Define Test Grid
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

pdf(file = "Figures/fig2_isoProportional_data.pdf",width = 2.25 ,height = 2.25)
ggplot(dec_bound,aes(x,y))+
  geom_point(data = points,aes(x,y,col = Class,shape = Partition),size = 2,alpha = .5)+
  geom_line(aes(x,y),col = "black")+
  scale_shape_manual(values = c(3,16))+
  scale_color_manual(values = c("red","blue"))+
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
classes = as.character(unique(df$Status))





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
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::roc(df2$Status,p)


##compute probs on test grid
test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
test_ran = clo(subset(test_ran,select = names(features)))
p = predict.glm(cv.clrlasso, newdata = data.frame(test_ran), type = "response")
cc = data.frame(x,y,prob = p)



pdf(file = "Figures/fig2_isoproportional_glmnet.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_line(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  #scale_color_gradient2(mid = "white",low = gplots::redblue(3)[1],high = gplots::redblue(3)[3],midpoint = .5)+
  theme_classic()+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
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

m1$performance
mroc = pROC::multiclass.roc(df2$Status,m1$predictionMatrix[,classes])
r = pROC::roc(df2$Status,m1$predictionMatrix$C1)
tt = data.frame(pROC::coords(r)) %>%
  filter(!is.infinite(threshold))

test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
testh =m1$models[[1]]
preds = predict.train(testh,test_ran,type = "prob" )
cc = data.frame(x,y,prob = preds$C1)



pdf(file = "Figures/fig2_isoproportional_ranger.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_line(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  #scale_color_gradient2(mid = "white",low = gplots::redblue(3)[1],high = gplots::redblue(3)[3],midpoint = .5)+
  theme_classic()+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
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





# CLR with feature selection --------------------------------------------------------------

## with feature selection
testx = clr(df2[,-1])
trainx = clr(df[,-1])

cv.clrlasso <- glmnet::cv.glmnet(as.matrix(trainx),df$Status, standardize=T, alpha=1,family="binomial",type.measure = "auc")
cv.clrlasso$cvm
features = as.matrix(coef(cv.clrlasso, s = "lambda.min"))
features = features[-1,]
features = features[abs(features)>0]
trainx = clr(subset(df,select = names(features)))
testx = clr(subset(df2,select = names(features)))

## compute AUC
glm.test = data.frame(Status = df$Status,trainx)
cv.clrlasso <-glm(formula = Status~.,data = glm.test,family = binomial)
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::roc(df2$Status,p)

##compute probs on test grid
test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
test_ran = clr(subset(test_ran,select = names(features)))
p = predict.glm(cv.clrlasso, newdata = data.frame(test_ran), type = "response")
cc = data.frame(x,y,prob = p)



pdf(file = "Figures/fig2_isoproportional_clrGLMNET.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_line(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  #scale_color_gradient2(mid = "white",low = gplots::redblue(3)[1],high = gplots::redblue(3)[3],midpoint = .5)+
  theme_classic()+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
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



pdf(file = "Figures/fig2_isoproportional_clrRanger.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_line(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  #scale_color_gradient2(mid = "white",low = gplots::redblue(3)[1],high = gplots::redblue(3)[3],midpoint = .5)+
  theme_classic()+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
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
p = predict.glm(cv.clrlasso, newdata = data.frame(testx), type = "response")
mroc = pROC::roc(df2$Status,p)

##compute probs on test grid
test_ran = test_grid
colnames(test_ran) = colnames(df[,-1])
pp  =calcLogRatio(data.frame(Status = "Test",fastImputeZeroes(test_ran)))
test_ran =subset(pp,select = names(features))
p = predict.glm(cv.clrlasso, newdata = data.frame(test_ran), type = "response")
cc = data.frame(x,y,prob = p)



pdf(file = "Figures/fig2_isoproportional_plrGLMNET.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_line(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  #scale_color_gradient2(mid = "white",low = gplots::redblue(3)[1],high = gplots::redblue(3)[3],midpoint = .5)+
  theme_classic()+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
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

ensemble = c("ranger")
m1.plr = trainML_Models(trainLRs = trainx,testLRs = testx,testIDs = data.frame(ID = 1:nrow(df2),Status =df2$Status),
                        ytrain = df[,1],y_test = df2[,1],models = ensemble,
                        mtry_ = round(sqrt(ncol(df[,-1]))),ntrees = 500,
                        numFolds = 10,cvMethod = "repeatedcv",numRepeats = 5
)
m1.plr$performance
pmat = m1.plr$predictionMatrix
pmat = pmat %>%
  group_by(ID,Status) %>%
  dplyr::select(-model) %>%
  summarise_all(.funs = mean)
pmat = data.frame(pmat)
classes = as.character(unique(df$Status))
mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes])
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
                  predict.train(testh,tt.lr[,-1],type = "prob" ),
                  model = names(m1.plr$models)[i])
  pmat = rbind(pmat,ph)
}

pmat = pmat %>%
  group_by(ID) %>%
  dplyr::select(-model) %>%
  summarise_all(.funs = mean)
cc = data.frame(x,y,prob = pmat$C1)


pdf(file = "Figures/fig2_isoproportional_plrEnsemble.pdf",width = 2.25 ,height = 2.25)
ggplot(cc,aes(x,y,col = prob))+
  geom_point(size = 1)+
  geom_line(data = dec_bound,aes(x,y),col = "black")+
  scale_color_distiller(palette = "RdBu")+
  #scale_color_gradient2(mid = "white",low = gplots::redblue(3)[1],high = gplots::redblue(3)[3],midpoint = .5)+
  theme_classic()+
  #geom_segment(aes(x = start_x, y = start_y, xend = end_x, yend = end_yy), data = border,size = 1,col  ="black")+
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



