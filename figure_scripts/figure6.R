# Code used to produce Figure 6

#Loading the required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Loading up the Results 
load("mus_musculus/tabula_muris_senis/R/results/Results_tabula_muris_senis.RData")#Tabula Muris Senis

source("functions.R")
#Setting up the data
#conservation of old-biased genes relative to no change
inc_cons_tabula <- lapply(Results_tabula_muris_senis, function(x){
  
  inc_cons <- sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]])) - 
      mean(c(-log(tissue[["deg_list"]][["Non_sig"]][["nonsig_genes_dnds"]][["dNdS"]])))
    
  })
})

#conservation of young-biased genes relative to no change
dec_cons_tabula <- lapply(Results_tabula_muris_senis, function(x){
  
  dec_cons <-sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]])) - 
      mean(c(-log(tissue[["deg_list"]][["Non_sig"]][["nonsig_genes_dnds"]][["dNdS"]])))
    
  })
})

#Creating mean&confInt dataframes for relative old and young-biased conservation scores

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
inc_cons_means_tabula$title <- c("Old-biased genes")
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
dec_cons_means_tabula$title <- c("Young-biased genes")
colnames(dec_cons_means_tabula)[1] <- "cons_means"

cons_means_tabula <- rbind(inc_cons_means_tabula, dec_cons_means_tabula)
cons_means_tabula$title <- factor(cons_means_tabula$title, levels= c("Old-biased genes", "Young-biased genes"))

unique(cons_means_tabula$names)

inc_tissue_dist_bulk <- cons_means_tabula[cons_means_tabula$title == "Old-biased genes",]$cons_means

dec_tissue_dist_bulk <- cons_means_tabula[cons_means_tabula$title == "Young-biased genes",]$cons_means

inc_tissue_dist <- sapply(unique(cons_means_tabula$Tissue), function(tissue){
  
  data <- cons_means_tabula[(cons_means_tabula$Tissue == tissue),]
  inc_cons <- data[data$title == "Old-biased genes",]$cons_means
  
})

dec_tissue_dist <- sapply(unique(cons_means_tabula$Tissue), function(tissue){
  
  data <- cons_means_tabula[(cons_means_tabula$Tissue == tissue),]
  inc_cons <- data[data$title == "Young-biased genes",]$cons_means
  
})

inc_dec_diff_test_bulk_p <- wilcox.test(inc_tissue_dist_bulk, dec_tissue_dist_bulk, paired = T)

inc_dec_diff_test_per_tissue_p <- sapply(names(inc_tissue_dist), function(tissue){
  
  wilcox.test(inc_tissue_dist[[tissue]], dec_tissue_dist[[tissue]], paired = T)$p.val
  
})

tabula_violin_text_p <- data.frame(Tissue = names(inc_dec_diff_test_per_tissue_p), 
                                   p = paste0("p = ", signif(inc_dec_diff_test_per_tissue_p, 3)))

tabula_violin_text_n <- data.frame(Tissue = names(inc_dec_diff_test_per_tissue_p), 
                                   n = paste0(summary(as.factor(cons_means_tabula$Tissue))/2, " cell types"))

#Plotting 
figure6 <- ggplot(cons_means_tabula, aes(x = title, y = cons_means)) +
  #geom_violin(aes(fill = Tissue), alpha = 0.6, show.legend = F) +
  geom_line(aes(group = names), color = 'gray36',
            alpha = 0.4, size = 0.3, show.legend = F) +
  facet_wrap(~Tissue, scale = "free_x") +
  theme_bw()+
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = 10,),
        axis.title.y = element_text(size=14), 
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size=14),
        strip.background = element_rect(fill='white', size = 1))+
  geom_hline(yintercept = 0)+
  scale_color_manual(values = c("#1B9E77", "#9e1b42", "#D95F02", "#027cd9", "#afb370", "#7570B3"))+
  scale_fill_manual(values = c("#1B9E77", "#9e1b42", "#D95F02", "#027cd9", "#afb370", "#7570B3"))+
  ylab("Mean Relative Conservation Score per Cell Type")+
  xlab("Gene Class.")+
#  geom_point(aes(color = Tissue), alpha = 1, show.legend = F) +
  geom_boxplot(width = 0.3, lwd=0.4, aes(fill= Tissue), alpha = 0.4, show.legend = F)+
  geom_jitter(width = 0.01, color ='gray36' , cex=1.4, alpha = 0.6, show.legend = F) +
  scale_x_discrete(labels=c("Old-biased genes" = "Old-biased", "Young-biased genes" = "Young-biased"))+
  geom_text(data    = tabula_violin_text_p,
            mapping = aes(x = Inf, y = -Inf, label = p),
            vjust = -1.1,
            hjust = 1.1,
            inherit.aes = FALSE,
            size = 4.5)+
  geom_text(data    = tabula_violin_text_n,
            mapping = aes(x = -Inf, y = -Inf, label = n),
            vjust = -1.1,
            hjust = -0.1,
            inherit.aes = FALSE,
            size = 4.5)

ggsave("figure6-tabula-old-vs-young-cons-diff.png", plot = figure6, device = "png",
       path = "results_graphs/",  type = "cairo",
       scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("figure6-tabula-old-vs-young-cons-diff.tiff", plot = figure6, device = "tiff",
       path = "results_graphs/",  type = "cairo",
       scale = 1, width = 12, height = 8, units = "in", dpi = 300)
ggsave("figure6-tabula-old-vs-young-cons-diff.pdf", plot = figure6, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 12, height = 8, units = "in", dpi = 300)
