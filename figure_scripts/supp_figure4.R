#Loading the required packages
library(tidyverse, lib.loc = "/home/myildiz/RLibs/")
library(ggforce, lib.loc = "/home/myildiz/RLibs/")
library(ggpubr, lib.loc = "/home/myildiz/RLibs/")
library(ggridges, lib.loc = "/home/myildiz/RLibs/")
library(ggrepel, lib.loc = "/home/myildiz/RLibs/")
library(RColorBrewer)
library(scales, lib.loc = "/home/myildiz/RLibs/")
library(igraph, lib.loc = "/home/myildiz/RLibs/")
library(patchwork)
library(unikn, lib.loc = "/home/myildiz/RLibs/")     
library(plyr)
library(gplots, lib.loc = "/home/myildiz/RLibs/")  # Visual plotting of tables
library(ggplot2, lib.loc = "/home/myildiz/RLibs/")
library(dplyr)
library(reshape2)
library(rlist)

#Loading up the Tabula Muris Results 
load("mus_musculus/tabula_muris_senis/R/results/Results_tabula_muris_senis.RData")#Tabula Muris Senis

exp_cons_cor_18m <- lapply(names(Results_tabula_muris_senis), function(tissue){
  
  exp_cons_cor_18m <- lapply(names(Results_tabula_muris_senis[[tissue]]), function(cell_type){
    
    if(sum(Results_tabula_muris_senis[[tissue]][[cell_type]][["age"]] == "18") > 2){
      
      complete_results = t(apply(Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]][["All"]][["all_genes_exp"]][,Results_tabula_muris_senis[[tissue]][[cell_type]][["age"]] == "18"], 2 , function(x){
        cor(c(-log(Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]][["All"]][["all_genes_dnds"]]$dNdS)), x,  method= "spearman",
            use = "everything")
      }))
      
    }
    
  })
  
  names(exp_cons_cor_18m) <- names(Results_tabula_muris_senis[[tissue]])
  
  exp_cons_cor_18m
  
})

names(exp_cons_cor_18m) <- names(Results_tabula_muris_senis)
#removing cell types that failed to meet indv. cutoff
exp_cons_cor_18m <- list.clean(exp_cons_cor_18m, fun = is.null, recursive = TRUE)
exp_cons_cor_18m <- exp_cons_cor_18m[-5] #removing skin, as no cell type passes our indv. limit


for(tissue in names(exp_cons_cor_18m)){
  
  exp_cons_cor_df <- unlist(exp_cons_cor_18m[[tissue]])
  exp_cons_cor_df <- data.frame(cor = exp_cons_cor_df, cell_type = gsub('.{1}$', '', names(exp_cons_cor_df)), Tissue = tissue)
  exp_cons_cor_18m[[tissue]] <- exp_cons_cor_df 
  
}

exp_cons_cor_18m_long <- rbind(exp_cons_cor_18m[["lung"]], exp_cons_cor_18m[["liver"]], 
                           exp_cons_cor_18m[["muscle"]], exp_cons_cor_18m[["brain"]], 
                           exp_cons_cor_18m[["kidney"]])
exp_cons_cor_18m_long$cell_type_per_tissue <- paste0(exp_cons_cor_18m_long$cell_type,".", exp_cons_cor_18m_long$Tissue)



##Cell types present in more than one tissue
aaa=  table(exp_cons_cor_18m_long[,c('cell_type','Tissue')])>0
aaa2= apply(aaa, 1, sum)
aaa3= aaa2[aaa2>1]
names(aaa3)

exp_cons_cor_18m_longx= exp_cons_cor_18m_long[!(exp_cons_cor_18m_long$cell_type%in%names(aaa3)),]
tissuesx <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#9e1b42", "#027cd9", "#afb370"), c('brain','liver', 'skin','kidney', 'lung', 'muscle'))
max(exp_cons_cor_18m_longx$cor)


p1=  ggplot(exp_cons_cor_18m_longx, aes(x = reorder(cell_type, cor, median), y = cor)) +
  geom_boxplot(aes(col=Tissue, fill=Tissue), lwd=0.4, width=0.3, outlier.shape = NA, alpha=0.6) + 
  scale_color_manual(name = "", values = tissuesx)+
  scale_fill_manual(name = "",  values = tissuesx)+
  scale_y_continuous(limits = c(0,0.32))+
  theme(legend.text = element_text(size = 7))+
  ylab("Expression-conservation correlation")+
  xlab("Cell type")+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_text(vjust=0.5, hjust=1, angle=90))+
  ggtitle("18 months old")

exp_cons_cor_18m_longy= exp_cons_cor_18m_long[(exp_cons_cor_18m_long$cell_type%in%names(aaa3)),]

p2= ggplot(exp_cons_cor_18m_longy, aes(x = cell_type, y = cor)) +
  geom_boxplot(aes(col=Tissue, fill=Tissue), lwd=0.4, width=0.3, outlier.shape= NA, alpha=0.6) + 
  scale_color_manual(name = "", values= tissuesx)+
  scale_fill_manual(name = "",  values= tissuesx)+
  scale_y_continuous(limits = c(0,0.32))+
  theme(legend.text= element_text(size = 7))+
  ylab("")+
  xlab("Cell type")+
  theme_pubr(base_size = 12)+
  guides(fill=F)+
  guides(color=F)+
  theme(axis.text.x = element_text(vjust=0.5, hjust=1, angle=90))



exp_cons_cor_24m <- lapply(names(Results_tabula_muris_senis), function(tissue){
  
  exp_cons_cor_24m <- lapply(names(Results_tabula_muris_senis[[tissue]]), function(cell_type){
    
    if(sum(Results_tabula_muris_senis[[tissue]][[cell_type]][["age"]] == "24") > 2){
      
      complete_results = t(apply(Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]][["All"]][["all_genes_exp"]][,Results_tabula_muris_senis[[tissue]][[cell_type]][["age"]] == "24"], 2 , function(x){
        cor(c(-log(Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]][["All"]][["all_genes_dnds"]]$dNdS)), x,  method= "spearman",
            use = "everything")
      }))
      
    }
    
  })
  
  names(exp_cons_cor_24m) <- names(Results_tabula_muris_senis[[tissue]])
  
  exp_cons_cor_24m
  
})

names(exp_cons_cor_24m) <- names(Results_tabula_muris_senis)
#removing cell types that failed to meet indv. cutoff
exp_cons_cor_24m <- list.clean(exp_cons_cor_24m, fun = is.null, recursive = TRUE)
exp_cons_cor_24m <- exp_cons_cor_24m[-c(2,4,5)] #removing liver, brain, andskin, as no cell type passes our indv. limit


for(tissue in names(exp_cons_cor_24m)){
  
  exp_cons_cor_df <- unlist(exp_cons_cor_24m[[tissue]])
  exp_cons_cor_df <- data.frame(cor = exp_cons_cor_df, cell_type = gsub('.{1}$', '', names(exp_cons_cor_df)), Tissue = tissue)
  exp_cons_cor_24m[[tissue]] <- exp_cons_cor_df 
  
}

exp_cons_cor_24m_long <- rbind(exp_cons_cor_24m[["lung"]], 
                               exp_cons_cor_24m[["muscle"]], 
                               exp_cons_cor_24m[["kidney"]])
exp_cons_cor_24m_long$cell_type_per_tissue <- paste0(exp_cons_cor_24m_long$cell_type,".", exp_cons_cor_24m_long$Tissue)



##Cell types present in more than one tissue
aaa=  table(exp_cons_cor_24m_long[,c('cell_type','Tissue')])>0
aaa2= apply(aaa, 1, sum)
aaa3= aaa2[aaa2>1]
names(aaa3)

exp_cons_cor_24m_longx= exp_cons_cor_24m_long[!(exp_cons_cor_24m_long$cell_type%in%names(aaa3)),]
tissuesx <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#9e1b42", "#027cd9", "#afb370"), c('brain','liver', 'skin','kidney', 'lung', 'muscle'))
max(exp_cons_cor_24m_longx$cor)


p3=  ggplot(exp_cons_cor_24m_longx, aes(x = reorder(cell_type, cor, median), y = cor)) +
  geom_boxplot(aes(col=Tissue, fill=Tissue), lwd=0.4, width=0.3, outlier.shape = NA, alpha=0.6) + 
  scale_color_manual(name = "", values = tissuesx)+
  scale_fill_manual(name = "",  values = tissuesx)+
  scale_y_continuous(limits = c(0,0.32))+
  theme(legend.text = element_text(size = 7))+
  ylab("Expression-conservation correlation")+
  xlab("Cell type")+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_text(vjust=0.5, hjust=1, angle=90))+
  ggtitle("24 months old")

exp_cons_cor_24m_longy= exp_cons_cor_24m_long[(exp_cons_cor_24m_long$cell_type%in%names(aaa3)),]

p4= ggplot(exp_cons_cor_24m_longy, aes(x = cell_type, y = cor)) +
  geom_boxplot(aes(col=Tissue, fill=Tissue), lwd=0.4, width=0.3, outlier.shape= NA, alpha=0.6) + 
  scale_color_manual(name = "", values= tissuesx)+
  scale_fill_manual(name = "",  values= tissuesx)+
  scale_y_continuous(limits = c(0,0.32))+
  theme(legend.text= element_text(size = 7))+
  ylab("")+
  xlab("Cell type")+
  theme_pubr(base_size = 12)+
  guides(fill=F)+
  guides(color=F)+
  theme(axis.text.x = element_text(vjust=0.5, hjust=1, angle=90))

supp_fig_4 <- p1 + p2 +  p3 + p4 + plot_layout(ncol = 2, width = c(3,1), guides = "collect")& 
  #  theme(legend.position = 'top')&
  plot_annotation(tag_levels = c('A', 'B'))& 
  theme(plot.tag = element_text(size = 14, face = "bold"))+
  theme(plot.margin = unit(c(0.6, 0.6, 0.1, 0.1), "cm"))

ggsave("suppFig4-celltype-diff-conservation.png", plot = supp_fig_4, device = png,
       path = "results_graphs/", type = "cairo",
       scale = 1, units = "cm", width = 42, height = 38, dpi = 300)
ggsave("suppFig4-celltype-diff-conservation.tiff", plot = supp_fig_4, device = "tiff",
       path = "results_graphs/", type = "cairo",
       scale = 1, units = "cm", width = 42, height = 38, dpi = 300)
ggsave("suppFig4-celltype-diff-conservation.pdf", plot = supp_fig_4, device = "pdf",
       path = "results_graphs/", 
       scale = 1, units = "cm", width = 42, height = 38, dpi = 300)

### Two-way ANOVA of tissue and cell type difference
summary(aov(cor ~ Tissue + cell_type, data = exp_cons_cor_18m_long))
summary(aov(cor ~ Tissue + cell_type, data = exp_cons_cor_24m_long))

### Nested ANOVA for immune cell type difference
immune_status <- read.csv("supplements/immune_status.csv")
immune_cells <- as.character(immune_status$cellTypes[immune_status$immuneStatus == "Immune"])
exp_cons_cor_18m_long$immune <- ifelse(exp_cons_cor_18m_long$cell_type %in% immune_cells, "Immune", "Non-immune")
exp_cons_cor_24m_long$immune <- ifelse(exp_cons_cor_24m_long$cell_type %in% immune_cells, "Immune", "Non-immune")

library(nlme)
anova(lme(cor~immune, random= ~1|Tissue, data= exp_cons_cor_18m_long))
anova(lme(cor~immune, random= ~1|Tissue, data= exp_cons_cor_24m_long))

