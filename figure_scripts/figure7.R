#Script used to produce figure 7

#Loading required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

#Loading results 
load("mus_musculus/tabula_muris_senis/R/results/Results_tabula_muris_senis.RData")#Tabula Muris Senis

#loading misc. functions
source("functions.R")

## Tests and calculations

#conservation of increasing genes relative to no change
inc_cons_tabula <- lapply(Results_tabula_muris_senis, function(x){
  
  inc_cons <- sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]])) - 
      mean(c(-log(tissue[["deg_list"]][["Non_sig"]][["nonsig_genes_dnds"]][["dNdS"]])))
    
  })
})

#conservation of decreasing genes relative to no change
dec_cons_tabula <- lapply(Results_tabula_muris_senis, function(x){
  
  dec_cons <-sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]])) - 
      mean(c(-log(tissue[["deg_list"]][["Non_sig"]][["nonsig_genes_dnds"]][["dNdS"]])))
    
  })
})

#Creating mean&confInt dataframes for relative increasing and decreasing conservation scores
inc_cons_means_tabula <- c(brain = sapply(inc_cons_tabula[["brain"]], mean), 
                           kidney = sapply(inc_cons_tabula[["kidney"]], mean), 
                           liver = sapply(inc_cons_tabula[["liver"]], mean), 
                           lung = sapply(inc_cons_tabula[["lung"]], mean), 
                           muscle = sapply(inc_cons_tabula[["muscle"]], mean), 
                           skin = sapply(inc_cons_tabula[["skin"]], mean))
names_inc_cons_means_tabula <- names(inc_cons_means_tabula)
inc_cons_means_tabula <- as.data.frame(inc_cons_means_tabula)
inc_cons_means_tabula$names <-  names_inc_cons_means_tabula
inc_cons_means_tabula$upper <- c(brain = sapply(inc_cons_tabula[["brain"]], function(x){conf_int(x)[1]}), 
                                 kidney = sapply(inc_cons_tabula[["kidney"]], function(x){conf_int(x)[1]}), 
                                 liver = sapply(inc_cons_tabula[["liver"]], function(x){conf_int(x)[1]}), 
                                 lung = sapply(inc_cons_tabula[["lung"]], function(x){conf_int(x)[1]}), 
                                 muscle = sapply(inc_cons_tabula[["muscle"]], function(x){conf_int(x)[1]}), 
                                 skin = sapply(inc_cons_tabula[["skin"]], function(x){conf_int(x)[1]}))

inc_cons_means_tabula$lower <- c(brain = sapply(inc_cons_tabula[["brain"]], function(x){conf_int(x)[2]}), 
                                 kidney = sapply(inc_cons_tabula[["kidney"]], function(x){conf_int(x)[2]}), 
                                 liver = sapply(inc_cons_tabula[["liver"]], function(x){conf_int(x)[2]}), 
                                 lung = sapply(inc_cons_tabula[["lung"]], function(x){conf_int(x)[2]}), 
                                 muscle = sapply(inc_cons_tabula[["muscle"]], function(x){conf_int(x)[2]}), 
                                 skin = sapply(inc_cons_tabula[["skin"]], function(x){conf_int(x)[2]}))

inc_cons_means_tabula$Tissue <- sapply(strsplit(inc_cons_means_tabula$names, "\\."), function(x){x[[1]]})

inc_cons_means_tabula$names <- sapply(strsplit(inc_cons_means_tabula$names, "\\."), function(x){x[[2]]})
inc_cons_means_tabula$title <- c("Genes Increasing with Age")
colnames(inc_cons_means_tabula)[1] <- "cons_means"

dec_cons_means_tabula <- c(brain = sapply(dec_cons_tabula[["brain"]], mean), 
                           kidney = sapply(dec_cons_tabula[["kidney"]], mean), 
                           liver = sapply(dec_cons_tabula[["liver"]], mean), 
                           lung = sapply(dec_cons_tabula[["lung"]], mean), 
                           muscle = sapply(dec_cons_tabula[["muscle"]], mean), 
                           skin = sapply(dec_cons_tabula[["skin"]], mean))
names_dec_cons_means_tabula <- names(dec_cons_means_tabula)
dec_cons_means_tabula <- as.data.frame(dec_cons_means_tabula)
dec_cons_means_tabula$names <-  names_dec_cons_means_tabula
dec_cons_means_tabula$upper <- c(brain = sapply(dec_cons_tabula[["brain"]], function(x){conf_int(x)[1]}), 
                                 kidney = sapply(dec_cons_tabula[["kidney"]], function(x){conf_int(x)[1]}), 
                                 liver = sapply(dec_cons_tabula[["liver"]], function(x){conf_int(x)[1]}), 
                                 lung = sapply(dec_cons_tabula[["lung"]], function(x){conf_int(x)[1]}), 
                                 muscle = sapply(dec_cons_tabula[["muscle"]], function(x){conf_int(x)[1]}), 
                                 skin = sapply(dec_cons_tabula[["skin"]], function(x){conf_int(x)[1]}))

dec_cons_means_tabula$lower <- c(brain = sapply(dec_cons_tabula[["brain"]], function(x){conf_int(x)[2]}), 
                                 kidney = sapply(dec_cons_tabula[["kidney"]], function(x){conf_int(x)[2]}), 
                                 liver = sapply(dec_cons_tabula[["liver"]], function(x){conf_int(x)[2]}), 
                                 lung = sapply(dec_cons_tabula[["lung"]], function(x){conf_int(x)[2]}), 
                                 muscle = sapply(dec_cons_tabula[["muscle"]], function(x){conf_int(x)[2]}), 
                                 skin = sapply(dec_cons_tabula[["skin"]], function(x){conf_int(x)[2]}))

dec_cons_means_tabula$Tissue <- sapply(strsplit(dec_cons_means_tabula$names, "\\."), function(x){x[[1]]})

dec_cons_means_tabula$names <- sapply(strsplit(dec_cons_means_tabula$names, "\\."), function(x){x[[2]]})
dec_cons_means_tabula$title <- c("Genes Decreasing with Age")
colnames(dec_cons_means_tabula)[1] <- "cons_means"

cons_means_tabula <- rbind(inc_cons_means_tabula, dec_cons_means_tabula)
cons_means_tabula$title <- factor(cons_means_tabula$title, levels= c("Genes Increasing with Age", "Genes Decreasing with Age"))


#Spearman's correlation

corr_tabula <- lapply(Results_tabula_muris_senis, function(results){
  
  corr_sig <- sapply(results, function(cell){
    #corr values for differentially expressed genes.
    corr_sig = cor(cell[["all_results"]][1,], cell[["age"]], 
                   method = "spearman")
    corr_sig = signif(corr_sig,2)
    
    corr_sig
  })
  
  
  p_all <- sapply(results, function(cell){
    p_all= cor.test(cell[["all_results"]][1,], cell[["age"]],
                    method = "spearman", exact = F)$p.val
    p_all
  })
  
  p_adjust <- p.adjust(p_all, "BH")
  
  cell_type <- as.factor(as.character(names(corr_sig)))
  corr <- data.frame(corr_sig, cell_type, p_all, p_adjust)
  
})

for(names in names(corr_tabula)){
  corr_tabula[[names]]$Tissue <- names
  
}


corr_tabula_df <- bind_rows(corr_tabula)
corr_tabula_df$p_adjust <- p.adjust(corr_tabula_df$p_all, "BH")


#Loading immune status of celll types
immune_status <- read.csv("supplements/immune_status.csv")
immune_cells <- as.character(immune_status$cellTypes[immune_status$immuneStatus == "Immune"])

corr_tabula_df$immune <- ifelse(corr_tabula_df$cell_type %in% immune_cells, "Immune", "Non-immune")

cons_means_tabula$immune <- ifelse(cons_means_tabula$names %in% immune_cells, "Immune", "Non-immune")

cons_diff_tabula <- sapply(1, function(x){
  
  inc_cons <- cons_means_tabula[cons_means_tabula$title == "Genes Increasing with Age",]
  dec_cons <- cons_means_tabula[cons_means_tabula$title == "Genes Decreasing with Age",]
  
  if (identical(inc_cons$names, dec_cons$names)){
    
      dec_cons$cons_means - inc_cons$cons_means
    
  } else{print("ERROR: reorder cell types to match")}
      
})

#using "Genes inc. with age" to extract names, they repeat for e dec. ones.
cons_diff_tabula <- data.frame(cell_type = cons_means_tabula[cons_means_tabula$title == "Genes Increasing with Age",]$names,
                               tissue = cons_means_tabula[cons_means_tabula$title == "Genes Increasing with Age",]$Tissue,
                               cons_diff = cons_diff_tabula)
cons_diff_tabula$immune <- ifelse(cons_diff_tabula$cell_type %in% immune_cells, "Immune", "Non-immune")

#Effect size
cohens_d(cons_diff_tabula[cons_diff_tabula$immune == "Non-immune",]$cons_diff, 
         cons_diff_tabula[cons_diff_tabula$immune == "Immune",]$cons_diff)

cons_diff_anova <- aov(cons_diff_tabula$cons_diff ~ 
                         cons_diff_tabula$immune / cons_diff_tabula$cell_type)

summary(cons_diff_anova)

#taking the mean of repeating cell types
mean_cons_diff_tabula <- as.data.frame(t(sapply(unique(cons_diff_tabula$cell_type), function(cell_type){
  
  mean_cons_diff = mean(cons_diff_tabula[cons_diff_tabula$cell_type == cell_type,]$cons_diff)
  cbind(as.character(cell_type), mean_cons_diff)
})))
colnames(mean_cons_diff_tabula) <- c("cell_type", "mean_cons_diff")
mean_cons_diff_tabula$mean_cons_diff <- 
  as.numeric(as.character(mean_cons_diff_tabula$mean_cons_diff))
mean_cons_diff_tabula$immune <- ifelse(mean_cons_diff_tabula$cell_type %in% immune_cells, "Immune", "Non-immune")


mean_corr_tabula_df <- as.data.frame(t(sapply(unique(corr_tabula_df$cell_type), function(cell_type){
  
  mean_corr = mean(corr_tabula_df[corr_tabula_df$cell_type == cell_type,]$corr_sig)
  cbind(as.character(cell_type), mean_corr)

  })))
colnames(mean_corr_tabula_df) <- c("cell_type", "mean_corr")
mean_corr_tabula_df$mean_corr <- 
  as.numeric(as.character(mean_corr_tabula_df$mean_corr))
mean_corr_tabula_df$immune <- ifelse(mean_corr_tabula_df$cell_type %in% immune_cells, "Immune", "Non-immune")


#Effect sizes.
cohens_d(mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Non-immune",]$mean_cons_diff, 
         mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Immune",]$mean_cons_diff)

cohens_d(mean_corr_tabula_df[mean_corr_tabula_df$immune == "Non-immune",]$mean_corr, 
         mean_corr_tabula_df[mean_corr_tabula_df$immune == "Immune",]$mean_corr)


immune_non_immune_diff_p <- data.frame(p = 
                                         paste0("p = ", 
                                                signif(wilcox.test(
                                                  mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Immune",]$mean_cons_diff, 
                                                  mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Non-immune",]$mean_cons_diff)$p.val,3)))


immune_non_immune_rho_p <- data.frame(p = 
                                        paste0("p = ", 
                                               signif(wilcox.test(
                                                 mean_corr_tabula_df$mean_corr[mean_corr_tabula_df$immune == "Immune"], 
                                                 mean_corr_tabula_df$mean_corr[mean_corr_tabula_df$immune == "Non-immune"])$p.val,3)))



first= 
  ggplot(mean_cons_diff_tabula, aes(x = immune, y = mean_cons_diff), fill = immune) +
  theme_pubr(base_size = 12)+
  geom_boxplot(width = 0.3, lwd=0.4, aes(fill = immune), alpha = 0.4) + 
  # , aes(fill = type), alpha = 0.2
  # geom_boxplot(lwd=0.4, color=c("#BEBADA","#8DD3C7"), outlier.shape = NA) + 
  ylab("Young-biased - Old-biased MRCS")+
  xlab("")+
  geom_jitter(width = 0.1, color ='gray36' , cex=1.7, alpha = 0.6)+
  # geom_jitter(width = 0.3, aes(color = type), alpha = 0.99,  cex=1.7)+
  scale_color_manual(values = c("#BEBADA","#8DD3C7"), )+
  scale_fill_manual(values = c("#BEBADA","#8DD3C7"))+
  geom_text(data    =immune_non_immune_diff_p,
            mapping = aes(x = Inf, y = -Inf, label = p),vjust = -1.1,hjust = 1.1,inherit.aes = FALSE, size = 3.4)+
  theme(legend.position = "none")


second= ggplot(mean_corr_tabula_df, aes(x = immune, y = mean_corr)) +
  theme_pubr(base_size = 12)+
  geom_boxplot(width = 0.3, lwd=0.4, aes(fill = immune), alpha = 0.4) + 
  ylab("rho")+
  xlab("")+
  geom_jitter(width = 0.1, color ='gray36' , cex=1.7, alpha = 0.6)+
  scale_color_manual(values = c("#BEBADA","#8DD3C7"), )+
  scale_fill_manual(values = c("#BEBADA","#8DD3C7"))+
  geom_text(data    =immune_non_immune_rho_p,
            mapping = aes(x = Inf, y = -Inf, label = p),vjust = -1.1,hjust = 1.1,inherit.aes = FALSE,size = 3.4)+
  theme(legend.position = "none") 

third = ggarrange(first,  second , labels = c('A','B'), 
                 nrow = 1, ncol =2)+
  theme(plot.margin = unit(c(0.55, 0.2, 0.001, 0.1), "cm"))

ggsave(third, 
       filename = "results_graphs/figure7-immune-diff.pdf",
       unit="cm", width = 24, height = 10, useDingbats = F,limitsize = FALSE)

ggsave(third, 
       filename = "results_graphs/figure7-immune-diff.png",
       device = png, type="cairo", unit="cm", width = 24, height = 10)

ggsave(third, 
       filename = "results_graphs/figure7-immune-diff.tiff",
       type='cairo',
       unit="cm", width = 24, height = 10)
