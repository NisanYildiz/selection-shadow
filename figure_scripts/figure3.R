# Code used to produce Figure 3

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

# Setting up the data

## Transcriptome conservation of 3-month old samples
exp_cons_cor_3m <- lapply(names(Results_tabula_muris_senis), function(tissue){
  
  exp_cons_cor_3m <- lapply(names(Results_tabula_muris_senis[[tissue]]), function(cell_type){
    
    if(sum(Results_tabula_muris_senis[[tissue]][[cell_type]][["age"]] == "3") > 2){
      
      complete_results = t(apply(Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]][["All"]][["all_genes_exp"]][,Results_tabula_muris_senis[[tissue]][[cell_type]][["age"]] == "3"], 2 , function(x){
        cor(c(-log(Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]][["All"]][["all_genes_dnds"]]$dNdS)), x,  method= "spearman",
            use = "everything")
      }))
      
    }
    
  })
  
  names(exp_cons_cor_3m) <- names(Results_tabula_muris_senis[[tissue]])
  
  exp_cons_cor_3m
  
})

names(exp_cons_cor_3m) <- names(Results_tabula_muris_senis)
#removing "NULL" cell types that failed to meet the 3 indv. cut-off
exp_cons_cor_3m <- list.clean(exp_cons_cor_3m, fun = is.null, recursive = TRUE) 
exp_cons_cor_3m <- exp_cons_cor_3m[-2] #removing liver, as no cell type passes our indv. limit



for(tissue in names(exp_cons_cor_3m)){
  
  exp_cons_cor_df <- unlist(exp_cons_cor_3m[[tissue]])
  exp_cons_cor_df <- data.frame(cor = exp_cons_cor_df, cell_type = gsub('.{1}$', '', names(exp_cons_cor_df)), Tissue = tissue)
  exp_cons_cor_3m[[tissue]] <- exp_cons_cor_df 
  
}

exp_cons_cor_long <- rbind(exp_cons_cor_3m[["lung"]], exp_cons_cor_3m[["muscle"]], exp_cons_cor_3m[["brain"]],
                           exp_cons_cor_3m[["skin"]], exp_cons_cor_3m[["kidney"]])
exp_cons_cor_long$cell_type_per_tissue <- paste0(exp_cons_cor_long$cell_type,".", exp_cons_cor_long$Tissue)



##Cell types present in more than one tissue
aaa=  table(exp_cons_cor_long[,c('cell_type','Tissue')])>0
aaa2= apply(aaa, 1, sum)
aaa3= aaa2[aaa2>1]
names(aaa3)

exp_cons_cor_longx= exp_cons_cor_long[!(exp_cons_cor_long$cell_type%in%names(aaa3)),]
tissuesx <- setNames(c("#1B9E77" ,"#7570B3" ,"#9E1B42", "#027CD9" ,"#AFB370" ), c('brain','skin', 'kidney','lung', 'muscle'))

max(exp_cons_cor_longx$cor)


p1=  ggplot(exp_cons_cor_longx, aes(x = reorder(cell_type, cor, median), y = cor)) +
  geom_boxplot(aes(col=Tissue, fill=Tissue), lwd=0.4, width=0.42, outlier.shape = NA, alpha=0.6) + 
  scale_color_manual(name = "", values = tissuesx)+
  scale_fill_manual(name = "",  values = tissuesx)+
  scale_y_continuous(limits = c(0,0.32))+
  theme(legend.text = element_text(size = 7))+
  ylab("Expression-conservation correlation")+
  xlab("Cell type")+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_text(vjust=0.5, hjust=1, angle=90))

exp_cons_cor_longy= exp_cons_cor_long[(exp_cons_cor_long$cell_type%in%names(aaa3)),]
p2= ggplot(exp_cons_cor_longy, aes(x = cell_type, y = cor)) +
  geom_boxplot(aes(col=Tissue, fill=Tissue), lwd=0.4, width=0.42, outlier.shape= NA, alpha=0.6) + 
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

figure3= p1 + p2 + plot_layout(ncol = 2, width = c(3,1))& 
  #  theme(legend.position = 'top')&
  plot_annotation(tag_levels = c('A', 'B'))& 
  theme(plot.tag = element_text(size = 14, face = "bold"))+
  theme(plot.margin = unit(c(0.6, 0.6, 0.1, 0.1), "cm"))

ggsave("figure3-celltype-diff-conservation.png", plot = figure3, device = "png",
       path = "results_graphs/", type = "cairo",
       scale = 1, units = "cm", width = 42, height = 22, dpi = 300)
ggsave("figure3-celltype-diff-conservation.tiff", plot = figure3, device = "tiff",
       path = "results_graphs/",
       scale = 1, units = "cm", width = 42, height = 22, dpi = 300)
ggsave("figure3-celltype-diff-conservation.pdf", plot = figure3, device = "pdf",
       path = "results_graphs/", 
       scale = 1, units = "cm", width = 42, height = 22, dpi = 300)


### Two-way ANOVA of tissue and cell type difference
summary(aov(cor ~ Tissue + cell_type, data = exp_cons_cor_long))

### Nested ANOVA for immune cell type difference
immune_status <- read.csv("supplements/immune_status.csv")
immune_cells <- as.character(immune_status$cellTypes[immune_status$immuneStatus == "Immune"])
exp_cons_cor_long$immune <- ifelse(exp_cons_cor_long$cell_type %in% immune_cells, "Immune", "Non-immune")

library(nlme)
anova(lme(cor~immune, random= ~1|Tissue, data= exp_cons_cor_long))

ggplot(data = exp_cons_cor_long)+
  geom_boxplot(aes(y = cor, x = immune))


