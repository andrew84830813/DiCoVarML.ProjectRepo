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


fnames = dir("Results/SimData2/")




## Final Data Sets to include
nm = "16S_meanShift"
nm[2] = "WGS_meanShift"


results_all = data.frame()
for(f_name in nm){
  bool = str_detect(fnames,f_name)
  f = fnames[bool]
  results = data.frame()
  for(i in f){
    ph = read_csv(file = paste0("Results/SimData2/",i))
    ph$corrected_fold = rep(c(rep(1,17),rep(2,17)),n_distinct(ph$Seed))
    results = rbind(results,ph)
  }
  results_all= rbind(results_all,results)
}

results_all1 = results_all %>%
  group_by(Scenario,Dataset,Seed,corrected_fold,Approach,permuteLabel,shift_parm,sparsity) %>%
  summarise_all(.funs = mean)

results_all1 = tidyr::separate(results_all1,col = 2,into = c("Dataset","shift","Permute","S","signalSparsity"),sep = "_")


res = results_all1 %>%
  filter(Approach %in% c("SELBAL","CLR-LASSO","DCV-ridgeEnsemble","DCV-ridgeRegression","Coda-LASSO")) %>%
  filter(shift_parm %in% c(4:6))

res$Approach = factor(res$Approach,levels = c("DCV-ridgeEnsemble","DCV-ridgeRegression","SELBAL","CLR-LASSO","Coda-LASSO"))



# res = results_all1 %>%
#   filter(!is.na(train_auc))
# seeds = unique(res$Seed)
# opt = data.frame()
# ens.df = data.frame(Approach = c("DCV-ridgeEnsemble","DCV-ridgeRegression",paste0(ensemble,1:length(ensemble))))
# for(s in seeds){
#   for(f in 1:2){
#     ph = res %>%
#       filter(Seed==s & corrected_fold == f) %>%
#       group_by(Scenario,Dataset,Seed,corrected_fold,permuteLabel,shift_parm,sparsity) %>%
#       summarise(train_auc = max(train_auc))
#     ph1 = res %>%
#       filter(Seed==s & corrected_fold == f) %>%
#       dplyr::select(Approach,AUC,train_auc,Scenario,Dataset,Seed,corrected_fold,permuteLabel,shift_parm,sparsity)
#     ph = left_join(ph,ph1)
#
#     ph.tie = ph %>%
#       group_by(Scenario,Dataset,Seed,corrected_fold,permuteLabel,shift_parm,sparsity) %>%
#       summarise(n = n()) %>%
#       filter(n>1) %>%
#       dplyr::select(-n)
#     ph2 = ph %>%
#       group_by(Scenario,Dataset,Seed,corrected_fold,permuteLabel,shift_parm,sparsity) %>%
#       summarise(n = n()) %>%
#       filter(n==1) %>%
#       dplyr::select(-n)
#     ph2 = left_join(ph2,ph)
#     ## handle ties by order/priority of ens.df
#     ph1 = data.frame()
#     for( i in 1:nrow(ph.tie)){
#       pp = ph.tie[i,]
#       pp = left_join(pp,ph)
#       pp = right_join(pp,ens.df)
#       ph1 = rbind(ph1,pp[1,])
#     }
#     ph = rbind(ph2,ph1)
#     opt = rbind(opt,ph)
#   }
# }

# opt$Approach = "DCV"
# res = results_all1 %>%
#   filter(is.na(train_auc))
# res = rbind(res,opt)


# res = res %>%
#   group_by(Scenario,Dataset,shift_parm) %>%
#   summarise_all(mean)
# ds = unique(res$Dataset)
# res.df = data.frame()
# for(d in ds){
#   ph = res %>%
#     filter(Dataset==d)
#   ph$col = "black"
#   i = which.max(ph$AUC)
#   ph$col[i] = "red"
#   res.df = rbind(res.df,ph)
# }
# res = data.frame(res)

#tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)

results_all1 = res %>%
  filter(Scenario=="Empirical")
# results_all2 = results_all %>%
#   filter(Scenario!="Empirical")
sets = data.frame(shift_parm = 1:6,shift_per = 100*(seq(1,1.3,length.out = 6)-1))
results_all1 = left_join(results_all1,sets)
dw = 0
results_all1$Signal_Density = (1-results_all1$sparsity)*100

pdf(file = "Figures/benchamrk_SimData.pdf",width = 7 ,height = 4)
ggplot(results_all1,aes(Signal_Density,AUC,col = Approach))+
  theme_bw()+
  ggsci::scale_color_d3()+
  stat_summary(fun.y = mean, geom = "line",size = .75,aes(group =Approach,col = Approach),
               #position = position_dodge2(dw)
               )+
  stat_summary(fun.y = mean, geom = "point",size = 3,aes(group =Approach,fill = Approach),
               #position = position_dodge2(dw)
               )+
  #stat_summary(fun.data = mean_cl_normal,geom = "pointrange",aes(group =Approach,col = Approach),width = .25,position = position_dodge2(dw))+
  facet_grid(Dataset~shift_per)+
  scale_shape_manual(values = 21:25)+
  xlab("% Mean Shifted Features")+
  scale_x_continuous(breaks = results_all1$Signal_Density)+
  theme(legend.position = "top",panel.grid = element_blank(),
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 8),
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8),
        #axis.title.y = element_blank(),
        #axis.text.y = element_text(size = 7),
        #axis.text.y = element_blank(),
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

dev.off()




#
#
#
# f <- function(x) {
#   r <- quantile(x, probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
#   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
#   r
# }
#
# sets = data.frame(shift_parm = 1:6,shift_per = (seq(1,1.3,length.out = 6)-1)*100)
# results_all = left_join(results_all,sets)
#
#
# pdf(file = "Figures/benchamrk_SimData.pdf",width = 7 ,height = 4)
# ggplot(results_all,aes(shift_per,AUC,col = Scenario,shape =Scenario))+
#   theme_bw()+
#   #geom_point()+
#   ggsci::scale_color_d3()+
#   xlab("% Mean Shift")+
#   #scale_x_continuous(breaks = unique(results_all$shift_parm))+
#   stat_summary(fun.y = mean, geom = "line",size = .75)+
#   stat_summary(fun.y = mean, geom = "point",size = 2)+
#   stat_summary(fun.data = ggplot2::mean_cl_normal,geom = "errorbar",width = .25)+
#   facet_grid(Dataset~Approach)+
#   theme(legend.position = "top",panel.grid = element_blank(),
#         strip.background = element_blank(),strip.text = element_text(face = "bold",size = 8),
#         plot.title = element_text(size = 7,hjust = .5,face = "bold"),
#         #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
#         axis.title = element_text(size = 12),
#         #axis.title.y = element_blank(),
#         #axis.text.y = element_text(size = 7),
#         #axis.text.y = element_blank(),
#         #legend.margin=margin(-1,-1,-1,-1),
#         strip.switch.pad.wrap = margin(0,0,0,0),
#         legend.margin=margin(-5,-10,-10,-10),
#         axis.text = element_text(size = 12),
#         #panel.grid = element_blank(),
#         legend.key.size = unit(.15,units = "in"),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         #legend.background = element_rect(colour = "black")
#   )
# dev.off()
#
#
#
#
#
#
#
# ## Pairwsie comparsion
# cx = list( c("DCV-rfRFE","CLR-LASSO"), c("DCV-rfRFE","Coda-LASSO"),
#            c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"),
#            c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO")  )
#
#
# dd = results_all1 %>%
#   dplyr::group_by(Scenario,Dataset,shift_parm) %>%
#   rstatix::wilcox_test(data =., AUC ~ Approach,paired = T,comparisons =  cx) %>%
#   rstatix::adjust_pvalue(method = "BH") %>%
#   rstatix::add_significance("p.adj") %>%
#   dplyr::arrange(dplyr::desc(-p))
#
#

