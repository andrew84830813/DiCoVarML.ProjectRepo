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


fnames = dir("Results/Exp_Benchamark3/")


# nm = "cmg_RubelMA-2020_STH"
# nm[2] = "cmg_ZhuF-2020_schizo"
# nm[3] = "cmg_QinN-2014_cirr"
# #nm[4] = "cmg_NielsenHB-2014_ibd"
# nm[4] = "cmg_FengQ-2015_crc"
# #nm[6] = "cmg_ThomasAM_2019_crc" ## sim
#  nm[7] = "cmg_WirbelJ-2018_crc"
# # nm[8] = "cmg_YachidaS-2019_crc" ## sim
#  nm[9] = "cmg_ZellerG_2014_crc"
# # #f_name = "cmg_ZVogtmannE_2016_crc" ## sim
#  nm[10] = "cmg_YuJ_2015_crc"
#
#

## 16S
nm = "selbal_Crohns_16s"
nm[2] = "selbal_HIV_16s"
nm[3] = "qitta_NAFLD_16s"
nm[4] = "mbiomeHD_cdiSchubert_16s"
## WGS
nm[5] = "cmg_FengQ-2015_crc"
nm[6] = "cmg_WirbelJ-2018_crc"
nm[7] = "cmg_YuJ_2015_crc"
nm[8] = "cmg_ZellerG_2014_crc"
## WGS
nm[5] = "cmg_QinN-2014_cirr"
nm[6] = "cmg_RubelMA-2020_STH"
nm[7] = "cmg_ZhuF-2020_schizo"
nm[8] = "cmg_ZellerG_2014_crc"



results_all = data.frame()
for(f_name in nm){
  bool = str_detect(fnames,paste0(f_name,"_seed"))
  f = fnames[bool]
  results = data.frame()
  for(i in f){
    ph = read_csv(file = paste0("Results/Exp_Benchamark3/",i))
    results = rbind(results,ph)
  }
  results = separate(results,col = 2,into = c("Dataset","s"),sep = "_seed") %>%
    dplyr::select(-s)
  ## correct fold mislabeling
  #results$corrected_fold = rep(c(rep(1,5),rep(2,5)),n_distinct(results$Seed))
  results$data_type = if_else(str_detect(f_name,"cmg"),"WGS","16S")
  results_all = rbind(results_all,results)
}
results_all$seed_fold = paste0(results_all$Seed,"_",results_all$Fold)


# seeds = sample(1:15,size = 5,replace = F)
# results_all = results_all %>%
#   filter(Seed %in% seeds)


#
# ## Assign Rank
# seed_fold = unique(results_all$seed_fold)
# ds = unique(results_all$Dataset)
# results_all1 = data.frame()
# for(sf in seed_fold){
#   ph = results_all %>%
#     filter(seed_fold==sf)
#   for(d in ds){
#     phh = ph %>%
#       filter(Dataset == d)
#     phh$rank = rank(-phh$AUC)
#     results_all1 = rbind(results_all1,phh)
#   }
# }


unique(results_all$Approach)


results_all = results_all %>%
  filter(Approach %in% c("SELBAL","CLR-LASSO","DCV-ridgeEnsemble","DCV-ridgeRegression","Coda-LASSO","ALR-GLMNET"))
results_all$Approach = factor(results_all$Approach,levels = c("DCV-ridgeEnsemble","DCV-ridgeRegression","ALR-GLMNET","CLR-LASSO","Coda-LASSO","SELBAL"))
results_all = results_all %>%
  group_by(Approach,Dataset,data_type,Seed) %>%
  summarise_all(mean)


res = results_all %>%
  group_by(Approach,Dataset,data_type) %>%
  summarise_all(mean)
ds = unique(res$Dataset)
res.df = data.frame()
for(d in ds){
  ph = res %>%
    filter(Dataset==d)
  ph$col = NA
  i = which.max(ph$AUC)
  ph$col[i] = "red"
  res.df = rbind(res.df,ph)
}
res = data.frame(res)

## Pairwsie comparsion
cx = list( c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"),c("DCV-ridgeEnsemble","SELBAL"),c("DCV-ridgeEnsemble","ALR-GLMNET"),
           c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO"),c("DCV-ridgeRegression","SELBAL"),c("DCV-ridgeRegression","ALR-GLMNET"))










#
#
#
#
# results_all = results_all %>%
#   filter(Approach %in% c("SELBAL","CLR-LASSO","DCV-ridgeEnsemble","DCV-ridgeRegression","Coda-LASSO","ALR-LASSO"))
# results_all$Approach = factor(results_all$Approach,levels = c("DCV-ridgeEnsemble","DCV-ridgeRegression","ALR-LASSO","CLR-LASSO","Coda-LASSO","SELBAL"))
# results_all = results_all %>%
#   group_by(Approach,Dataset,data_type,Seed) %>%
#   summarise_all(mean)
#
#
# res = results_all %>%
#   group_by(Approach,Dataset,data_type) %>%
#   summarise_all(mean)
# ds = unique(res$Dataset)
# res.df = data.frame()
# for(d in ds){
#   ph = res %>%
#     filter(Dataset==d)
#   ph$col = NA
#   i = which.max(ph$AUC)
#   ph$col[i] = "red"
#   res.df = rbind(res.df,ph)
# }
# res = data.frame(res)
#
# ## Pairwsie comparsion
# cx = list( c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"),c("DCV-ridgeEnsemble","SELBAL"),c("DCV-ridgeEnsemble","ALR-LASSO"),
#            c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO"),c("DCV-ridgeRegression","SELBAL"),c("DCV-ridgeRegression","ALR-LASSO"))

# cx = list( c("DCV","CLR-LASSO"), c("DCV","Coda-LASSO"),c("DCV","SELBAL"))

# res1 = results_all %>%
#   group_by(Dataset,data_type) %>%
#   summarise_all(mean) %>%
#   select(Dataset,AUC)
#
# res2 = results_all %>%
#   select(Dataset,data_type,Approach,seed_fold,AUC) %>%
#   spread("Approach","AUC")
#
#
#
#
#
#
#
#
#
# res1 = results_all %>%
#   group_by(Dataset,data_type,Approach) %>%
#   summarise_all(mean)

res2 = results_all %>%
  dplyr::select(Dataset,data_type,Approach,Seed,AUC) %>%
  spread("Approach","AUC")

# res2 = results_all %>%
#   dplyr::select(Dataset,data_type,Approach,seed_fold,AUC) %>%
#   spread("Approach","AUC")

wt.df = data.frame()

for(d in ds){
  ph = res2 %>%
    filter(Dataset==d)
  phh = data.frame()
  for(i in 1:length(cx)){
    wt = wilcox.test(pull(ph,cx[[i]][1]),pull(ph,cx[[i]][2]),paired = T)
    # st = coin::wilcoxsign_test( pull(ph,cx[[i]][1])~pull(ph,cx[[i]][2]) ,alternative = "greater")
    # st@statistic@teststatistic/sqrt(length(pull(ph,cx[[i]][1])))
    dd = data.frame(dif = pull(ph,cx[[i]][1]) - pull(ph,cx[[i]][2]))
    mean(dd$dif)
    n = length(pull(ph,cx[[i]][1]))
    dtt = data.frame(auc = c(pull(ph,cx[[i]][1]),pull(ph,cx[[i]][2])),
                     label = c(rep(cx[[i]][1],n), rep(cx[[i]][2],n)))
    ef = rstatix::wilcox_effsize(data = dtt,formula = auc~label,alternative = "greater",ref.group = cx[[i]][1],paired = T)

    phh = rbind(phh,data.frame(Data_type = unique(ph$data_type),Dataset = d,g1 = cx[[i]][1],g2 = cx[[i]][2],p = wt$p.value,stat = wt$statistic,
                               mean_diff = mean(dd$dif),sd_diff =sd(dd$dif) ,ef ))
  }
  phh$p.adj = p.adjust(phh$p,method = "BH")
  wt.df = rbind(wt.df,phh)
}


wt.df$signf = if_else(wt.df$p.adj<0.05,T,F)
wt.df$fill_p = if_else(wt.df$p.adj<0.05,wt.df$p.adj,NULL)
wt.df$star = gtools::stars.pval(wt.df$p.adj)
wt.df$star = if_else(wt.df$signf,wt.df$star,"NS")
wt.df$p.adj1 = if_else(wt.df$p.adj>0.0001,as.character(sprintf("%.4f",round(wt.df$p.adj,digits = 4))),"<0.0001")

#sprintf("%.3f", p.adj)

w1 = wt.df
w1$p.adj1 = if_else(str_detect(w1$p.adj1,"<"),paste0("p",w1$p.adj1),paste0("p=",w1$p.adj1))
w1$effsize = as.numeric(w1$effsize)
w1 = w1 %>%
  mutate(effsize = if_else(signf,effsize,NULL))
w1$p.adj1 = if_else(w1$signf,w1$p.adj1,"N.S.")

pdf(file = "Figures/benchamrk_ExpData_effectSize.pdf",width = 7 ,height = 3)
ggplot(w1,aes(mean_diff,g2,label =p.adj1,fill  = g1 ))+
  geom_col(position = position_dodge2(width = .1),col = "black") +
  geom_text(size = 2,position = position_dodge2(width = 1,padding = 1),hjust = -.1)+
  ggsci::scale_fill_d3()+
  facet_wrap(.~Dataset,nrow = 2)+
  geom_vline(xintercept = 0)+
  #scale_x_continuous(breaks = c(0,0.01,0.05))+
  theme_bw()+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",size = 8),
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8),
        #axis.text.y = element_text(size = 7),
        #axis.text.x =element_text(hjust = 1),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        #panel.grid.minor.x   = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        #legend.background = element_rect(colour = "black")
  )
dev.off()

pdf(file = "Figures/benchamrk_ExpData_signf.pdf",width = 7 ,height = 3)
ggplot(wt.df,aes(g1,g2,fill = fill_p,label =p.adj1 ))+
  geom_tile(col="white")+
  geom_text(size = 2)+
  scale_fill_distiller(direction = -1,palette = "Reds")+
  facet_wrap(.~Dataset,nrow = 2)+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),strip.text = element_text(face = "bold",size = 8),
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 8),
        #axis.text.y = element_text(size = 7),
        #axis.text.x =element_text(hjust = 1),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )

dev.off()

library(rstatix)
colnames(res1)[2] = "y.position"
r = results_all %>%
  dplyr::select(Dataset,Approach,seed_fold,AUC) %>%
  filter(Dataset == results_all$Dataset[1])
dd = r %>%
  dplyr::group_by(Dataset) %>%
  rstatix::wilcox_test(data =., AUC ~ Approach,paired = T,alternative = "greater",comparisons = cx) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::arrange(dplyr::desc(-p)) %>%
  filter(p.adj<0.05)




#tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
pdf(file = "Figures/benchamrk_ExpData.pdf",width = 7 ,height = 3)
ggplot(results_all,aes(Approach,AUC,shape = Approach))+
  theme_bw()+
  geom_line(aes(group =Seed),col = "gray",alpha = .5)+
  stat_summary(fun.y = mean, geom = "line",size = .75,col = "black",aes(group =1))+
  #stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
  geom_point(size = 1,col = "gray")+
  ggsci::scale_fill_d3()+
  ggsci::scale_color_d3()+
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
dev.off()








#tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
pdf(file = "Figures/benchamrk_ExpData_numberParts.pdf",width = 7 ,height = 4)
ggplot(results_all,aes(Approach,number_parts,fill = Approach))+
  theme_bw()+
  stat_summary(fun.y = mean,geom = "col")+
  stat_summary(fun.y = mean, geom = "point")+
  stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
  scale_shape_manual(values = 21:26)+
  facet_wrap(.~Dataset,nrow = 2,
             scales = "free_y"
  )+
  theme(legend.position = "top",
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        strip.background = element_blank(),strip.text = element_text(face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )
dev.off()


#
#
#
# ggplot(results_all,aes(Approach,comp_time,fill = Approach))+
#   theme_bw()+
#   stat_summary(fun.y = mean,geom = "col")+
#   stat_summary(fun.y = mean, geom = "point")+
#   stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
#   scale_shape_manual(values = 21:26)+
#   facet_wrap(.~Dataset,nrow = 2,
#              scales = "free_y"
#   )+
#   theme(legend.position = "top",
#         plot.title = element_text(size = 7,hjust = .5,face = "bold"),
#         strip.background = element_blank(),strip.text = element_text(face = "bold"),
#         #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
#         axis.title = element_text(size = 8),
#         axis.title.y = element_blank(),
#         #axis.text.y = element_text(size = 7),
#         axis.text.x = element_blank(),
#         #legend.margin=margin(-1,-1,-1,-1),
#         strip.switch.pad.wrap = margin(0,0,0,0),
#         legend.margin=margin(-5,-10,-10,-10),
#         axis.text = element_text(size = 8),
#         panel.grid = element_blank(),
#         legend.key.size = unit(.15,units = "in"),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         #legend.background = element_rect(colour = "black")
#   )
#
#
#
# #
# ggplot(results_all1,aes(Approach,rank,shape = Approach))+
#   theme_bw()+
#   geom_line(aes(group =seed_fold),col = "gray",alpha = .2)+
#   scale_y_reverse()+
#   #geom_boxplot()+
#   #coord_flip()+
stat_summary(fun.y library(diffCompVarRcpp)
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


             fnames = dir("Results/Exp_Benchamark/")


             # nm = "cmg_RubelMA-2020_STH"
             # nm[2] = "cmg_ZhuF-2020_schizo"
             # nm[3] = "cmg_QinN-2014_cirr"
             # #nm[4] = "cmg_NielsenHB-2014_ibd"
             # nm[4] = "cmg_FengQ-2015_crc"
             # #nm[6] = "cmg_ThomasAM_2019_crc" ## sim
             #  nm[7] = "cmg_WirbelJ-2018_crc"
             # # nm[8] = "cmg_YachidaS-2019_crc" ## sim
             #  nm[9] = "cmg_ZellerG_2014_crc"
             # # #f_name = "cmg_ZVogtmannE_2016_crc" ## sim
             #  nm[10] = "cmg_YuJ_2015_crc"
             #
             #

             ## 16S
             nm = "selbal_Crohns_16s"
             nm[2] = "selbal_HIV_16s"
             nm[3] = "qitta_NAFLD_16s"
             nm[4] = "mbiomeHD_cdiSchubert_16s"
             ## WGS
             nm[5] = "cmg_QinN-2014_cirr"
             nm[6] = "cmg_RubelMA-2020_STH"
             nm[7] = "cmg_ZhuF-2020_schizo"
             nm[8] = "cmg_ZellerG_2014_crc"



             results_all = data.frame()
             for(f_name in nm){
               bool = str_detect(fnames,paste0(f_name,"_seed"))
               f = fnames[bool]
               results = data.frame()
               for(i in f){
                 ph = read_csv(file = paste0("Results/Exp_Benchamark/",i))
                 results = rbind(results,ph)
               }
               results = separate(results,col = 2,into = c("Dataset","s"),sep = "_seed") %>%
                 dplyr::select(-s)
               ## correct fold mislabeling
               results$corrected_fold = rep(c(rep(1,5),rep(2,5)),n_distinct(results$Seed))
               results$data_type = if_else(str_detect(f_name,"cmg"),"WGS","16S")
               results_all = rbind(results_all,results)
             }
             results_all$seed_fold = paste0(results_all$Seed,"_",results_all$corrected_fold)



             ## Assign Rank
             seed_fold = unique(results_all$seed_fold)
             ds = unique(results_all$Dataset)
             results_all1 = data.frame()
             for(sf in seed_fold){
               ph = results_all %>%
                 filter(seed_fold==sf)
               for(d in ds){
                 phh = ph %>%
                   filter(Dataset == d)
                 phh$rank = rank(-phh$AUC)
                 results_all1 = rbind(results_all1,phh)
               }
             }






             res = results_all %>%
               group_by(Approach,Dataset,data_type) %>%
               summarise_all(mean)
             ds = unique(res$Dataset)
             res.df = data.frame()
             for(d in ds){
               ph = res %>%
                 filter(Dataset==d)
               ph$col = "lightblue"
               i = which.max(ph$AUC)
               ph$col[i] = "red"
               res.df = rbind(res.df,ph)
             }
             res = data.frame(res)

             ## Pairwsie comparsion
             cx = list( c("DCV-ridgeEnsemble","CLR-LASSO"), c("DCV-ridgeEnsemble","Coda-LASSO"),c("DCV-ridgeEnsemble","SELBAL"),
                        c("DCV-ridgeRegression","CLR-LASSO"),c("DCV-ridgeRegression","Coda-LASSO"),c("DCV-ridgeRegression","SELBAL")  )



             res1 = results_all %>%
               group_by(Dataset,data_type) %>%
               summarise_all(mean) %>%
               select(Dataset,AUC)

             res2 = results_all %>%
               select(Dataset,data_type,Approach,seed_fold,AUC) %>%
               spread("Approach","AUC")

             wt.df = data.frame()
             dt = unique(results_all$data_type)
             for(s in dt){
               for(d in ds){
                 ph = res2 %>%
                   filter(Dataset==d)
                 phh = data.frame()
                 for(i in 1:length(cx)){
                   wt = wilcox.test(pull(ph,cx[[i]][1]),pull(ph,cx[[i]][2]),paired = T)
                   phh = rbind(phh,data.frame(Data_type = s,Dataset = d,g1 = cx[[i]][1],g2 = cx[[i]][2],p = wt$p.value,stat = wt$statistic ))
                 }
                 phh$p.adj = p.adjust(phh$p,method = "BH")
                 wt.df = rbind(wt.df,phh)
               }
             }

             wt.df$signf = if_else(wt.df$p.adj<0.05,T,F)
             wt.df$fill_p = if_else(wt.df$p.adj<0.05,wt.df$p.adj,NULL)
             pdf(file = "Figures/benchamrk_ExpData_signf.pdf",width = 7 ,height = 3)
             ggplot(wt.df,aes(g1,g2,fill = fill_p,label = round(p.adj,6)))+
               geom_tile(col="white")+
               geom_text(size = 2)+
               scale_fill_distiller(direction = -1,palette = "Reds")+
               facet_wrap(.~Dataset,nrow = 2)+
               theme_bw()+
               theme(legend.position = "none",
                     strip.background = element_blank(),strip.text = element_text(face = "bold",size = 8),
                     plot.title = element_text(size = 7,hjust = .5,face = "bold"),
                     #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
                     axis.title = element_text(size = 8),
                     axis.title.y = element_blank(),
                     axis.text = element_text(size = 8),
                     #axis.text.y = element_text(size = 7),
                     axis.text.x =element_text(hjust = 1),
                     #legend.margin=margin(-1,-1,-1,-1),
                     strip.switch.pad.wrap = margin(0,0,0,0),
                     legend.margin=margin(-5,-10,-10,-10),
                     panel.grid = element_blank(),
                     legend.key.size = unit(.15,units = "in"),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 8),
                     #legend.background = element_rect(colour = "black")
               )

             dev.off()

             colnames(res1)[2] = "y.position"
             dd = results_all %>%
               dplyr::group_by(Dataset) %>%
               rstatix::wilcox_test(data =., AUC ~ Approach,paired = T,alternative = "greater",comparisons = cx) %>%
               rstatix::adjust_pvalue(method = "BH") %>%
               rstatix::add_significance("p.adj") %>%
               dplyr::arrange(dplyr::desc(-p)) %>%
               filter(p.adj<0.05) %>%
               left_join(res1)


             #tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
             pdf(file = "Figures/benchamrk_ExpData.pdf",width = 7 ,height = 3)
             ggplot(results_all,aes(Approach,AUC,shape = Approach))+
               theme_bw()+
               geom_line(aes(group =seed_fold),col = "gray",alpha = .2)+
               stat_summary(fun.y = mean, geom = "line",size = .75,col = "black",aes(group =1))+
               stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
               geom_point(col  ="gray",size = 1,alpha = .5)+
               scale_shape_manual(values = 21:26)+
               geom_point(data = res.df,aes(Approach,AUC),fill = res.df$col,col = "black",size = 2)+
               facet_wrap(.~Dataset,nrow = 2,
                          scales = "free_y"
               )+
               theme(legend.position = "top",
                     plot.title = element_text(size = 7,hjust = .5,face = "bold"),
                     strip.background = element_blank(),strip.text = element_text(face = "bold"),
                     #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
                     axis.title = element_text(size = 8),
                     axis.title.y = element_blank(),
                     #axis.text.y = element_text(size = 7),
                     axis.text.x = element_blank(),
                     #legend.margin=margin(-1,-1,-1,-1),
                     strip.switch.pad.wrap = margin(0,0,0,0),
                     legend.margin=margin(-5,-10,-10,-10),
                     axis.text = element_text(size = 8),
                     panel.grid = element_blank(),
                     legend.key.size = unit(.15,units = "in"),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 8),
                     #legend.background = element_rect(colour = "black")
               )
             dev.off()








             #tiff(filename =paste0(f_name,".tiff"),width = 4.5,height = 5.5,units = "in",res = 300)
             pdf(file = "Figures/benchamrk_ExpData_numberParts.pdf",width = 7 ,height = 4)
             ggplot(results_all,aes(Approach,number_parts,fill = Approach))+
               theme_bw()+
               stat_summary(fun.y = mean,geom = "col")+
               stat_summary(fun.y = mean, geom = "point")+
               stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
               scale_shape_manual(values = 21:26)+
               facet_wrap(.~Dataset,nrow = 2,
                          scales = "free_y"
               )+
               theme(legend.position = "top",
                     plot.title = element_text(size = 7,hjust = .5,face = "bold"),
                     strip.background = element_blank(),strip.text = element_text(face = "bold"),
                     #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
                     axis.title = element_text(size = 8),
                     axis.title.y = element_blank(),
                     #axis.text.y = element_text(size = 7),
                     axis.text.x = element_blank(),
                     #legend.margin=margin(-1,-1,-1,-1),
                     strip.switch.pad.wrap = margin(0,0,0,0),
                     legend.margin=margin(-5,-10,-10,-10),
                     axis.text = element_text(size = 8),
                     panel.grid = element_blank(),
                     legend.key.size = unit(.15,units = "in"),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 8),
                     #legend.background = element_rect(colour = "black")
               )
             dev.off()





             ggplot(results_all,aes(Approach,comp_time,fill = Approach))+
               theme_bw()+
               stat_summary(fun.y = mean,geom = "col")+
               stat_summary(fun.y = mean, geom = "point")+
               stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
               scale_shape_manual(values = 21:26)+
               facet_wrap(.~Dataset,nrow = 2,
                          scales = "free_y"
               )+
               theme(legend.position = "top",
                     plot.title = element_text(size = 7,hjust = .5,face = "bold"),
                     strip.background = element_blank(),strip.text = element_text(face = "bold"),
                     #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
                     axis.title = element_text(size = 8),
                     axis.title.y = element_blank(),
                     #axis.text.y = element_text(size = 7),
                     axis.text.x = element_blank(),
                     #legend.margin=margin(-1,-1,-1,-1),
                     strip.switch.pad.wrap = margin(0,0,0,0),
                     legend.margin=margin(-5,-10,-10,-10),
                     axis.text = element_text(size = 8),
                     panel.grid = element_blank(),
                     legend.key.size = unit(.15,units = "in"),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 8),
                     #legend.background = element_rect(colour = "black")
               )



             #
             ggplot(results_all1,aes(Approach,rank,shape = Approach))+
               theme_bw()+
               geom_line(aes(group =seed_fold),col = "gray",alpha = .2)+
               scale_y_reverse()+
               #geom_boxplot()+
               #coord_flip()+
               stat_summary(fun.y = mean, geom = "line",size = 1,col = "black",aes(group =1))+
               stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
               stat_summary(fun.y = mean,geom = "point",width = .1,aes(shape = Approach),size = 2,pch = 17)+

               # geom_signif(comparisons = cx,tip_length = .05,
               #             test = "wilcox.test",hjust = -.5,vjust = -1,textsize = 3,
               #             test.args = list(paired = T),y_position = seq(1,1.25,length.out = 6),
               #             map_signif_level = function(p) sprintf("p = %.2g", p)
               #            )+
               geom_point(col  ="gray",size = 3,alpha = .2,pch = 16)+
               scale_shape_manual(values = 21:26)+
               #geom_point(data = res.df,aes(Approach,AUC),fill = res.df$col,col = "black",size = 1)+
               facet_wrap(.~Dataset,nrow = 2,
                          #scales = "free_y"
               )+
               theme(legend.position = "top",
                     plot.title = element_text(size = 7,hjust = .5,face = "bold"),
                     #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
                     axis.title = element_text(size = 8),
                     axis.title.y = element_blank(),
                     #axis.text.y = element_text(size = 7),
                     axis.text.x = element_blank(),
                     #legend.margin=margin(-1,-1,-1,-1),
                     strip.switch.pad.wrap = margin(0,0,0,0),
                     legend.margin=margin(-5,-10,-10,-10),
                     axis.text = element_text(size = 8),
                     panel.grid = element_blank(),
                     legend.key.size = unit(.15,units = "in"),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 8),
                     #legend.background = element_rect(colour = "black")
               )

             = mean, geom = "line",size = 1,col = "black",aes(group =1))+
  stat_summary(fun.data = mean_cl_normal,geom = "errorbar",width = .1)+
  stat_summary(fun.y = mean,geom = "point",width = .1,aes(shape = Approach),size = 2,pch = 17)+

  # geom_signif(comparisons = cx,tip_length = .05,
  #             test = "wilcox.test",hjust = -.5,vjust = -1,textsize = 3,
  #             test.args = list(paired = T),y_position = seq(1,1.25,length.out = 6),
  #             map_signif_level = function(p) sprintf("p = %.2g", p)
  #            )+
  geom_point(col  ="gray",size = 3,alpha = .2,pch = 16)+
  scale_shape_manual(values = 21:26)+
  #geom_point(data = res.df,aes(Approach,AUC),fill = res.df$col,col = "black",size = 1)+
  facet_wrap(.~Dataset,nrow = 2,
             #scales = "free_y"
  )+
  theme(legend.position = "top",
        plot.title = element_text(size = 7,hjust = .5,face = "bold"),
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        #legend.margin=margin(-1,-1,-1,-1),
        strip.switch.pad.wrap = margin(0,0,0,0),
        legend.margin=margin(-5,-10,-10,-10),
        axis.text = element_text(size = 8),
        panel.grid = element_blank(),
        legend.key.size = unit(.15,units = "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.background = element_rect(colour = "black")
  )

