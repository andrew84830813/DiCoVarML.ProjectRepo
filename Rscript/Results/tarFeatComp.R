library(dplyr)
library(ggplot2)
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


fnames = dir("Results/tarSelectionComp/")


results = data.frame()
for(i in fnames){
  ph = read_csv(file = paste0("Results/tarSelectionComp/",i))
  results = rbind(results,ph)
}

results1 = results%>%
  filter(Approach %in%  c("PLR_PLR_DCV","ALR_ALR_DCV","ALR/PLR_ALR/PLR_DCV"))%>%
  group_by(Approach,Seed,Dataset,feat_set,Scenario)%>%
  summarise_all(.funs = mean)
results1 = tidyr::separate(results1,col = 3,into = c("Type","Dataset","Disease"),sep = "_")
results1$Type = if_else(results1$Type=="cmg","WGS","16S")

results1$feat_set
ggplot(results1,aes(feat_set,AUC,col = Approach))+
  theme_bw()+
  #geom_line(aes(group =Seed),col = "gray",alpha = .5)+
  stat_summary(fun.y = mean, geom = "line",size = .75,aes(group = Approach))+
  stat_summary(fun.y =mean,geom = "point",width = .1)+
  #geom_point(size = 1,col = "gray")+
  ggsci::scale_fill_d3()+
  ggsci::scale_color_d3()+
  facet_wrap(.~Type,nrow = 1,
             scales = "free_y"
  )
  #scale_shape_manual(values = c(21:26,8))+
  stat_summary(fun.y = mean, geom = "point",size = 3,aes(fill = Approach,col = Approach))+
  geom_point(data = res.df,aes(Approach,AUC),fill = res.df$col,col = res.df$col,size = if_else(is.na(res.df$col),2,5))+
  facet_wrap(.~Dataset,nrow = 2,
             scales = "free_y"
  )+
  theme(legend.position = "top",
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        strip.background = element_blank(),strip.text = element_text(face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        #panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )

