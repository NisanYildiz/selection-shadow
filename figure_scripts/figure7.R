#Script used to produce figure 7

#Loading required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

library(org.Mm.eg.db)
library(GO.db)
library(topGO)
library(ALL)
library(ggh4x)


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

corr_tabula_df$immune <- ifelse(corr_tabula_df$cell_type %in% immune_cells, "Immune cells", "Non-immune cells")

cons_means_tabula$immune <- ifelse(cons_means_tabula$names %in% immune_cells, "Immune cells", "Non-immune cells")

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
cons_diff_tabula$immune <- ifelse(cons_diff_tabula$cell_type %in% immune_cells, "Immune cells", "Non-immune cells")

#Effect size
cohens_d(cons_diff_tabula[cons_diff_tabula$immune == "Non-immune cells",]$cons_diff, 
         cons_diff_tabula[cons_diff_tabula$immune == "Immune cells",]$cons_diff)

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
mean_cons_diff_tabula$immune <- ifelse(mean_cons_diff_tabula$cell_type %in% immune_cells, "Immune cells", "Non-immune cells")


mean_corr_tabula_df <- as.data.frame(t(sapply(unique(corr_tabula_df$cell_type), function(cell_type){
  
  mean_corr = mean(corr_tabula_df[corr_tabula_df$cell_type == cell_type,]$corr_sig)
  cbind(as.character(cell_type), mean_corr)

  })))
colnames(mean_corr_tabula_df) <- c("cell_type", "mean_corr")
mean_corr_tabula_df$mean_corr <- 
  as.numeric(as.character(mean_corr_tabula_df$mean_corr))
mean_corr_tabula_df$immune <- ifelse(mean_corr_tabula_df$cell_type %in% immune_cells, "Immune cells", "Non-immune cells")


#Effect sizes.
cohens_d(mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Non-immune cells",]$mean_cons_diff, 
         mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Immune cells",]$mean_cons_diff)

cohens_d(mean_corr_tabula_df[mean_corr_tabula_df$immune == "Non-immune cells",]$mean_corr, 
         mean_corr_tabula_df[mean_corr_tabula_df$immune == "Immune cells",]$mean_corr)


immune_non_immune_diff_p <- data.frame(p = 
                                         paste0("p = ", 
                                                signif(wilcox.test(
                                                  mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Immune cells",]$mean_cons_diff, 
                                                  mean_cons_diff_tabula[mean_cons_diff_tabula$immune == "Non-immune cells",]$mean_cons_diff)$p.val,3)))


immune_non_immune_rho_p <- data.frame(p = 
                                        paste0("p = ", 
                                               signif(wilcox.test(
                                                 mean_corr_tabula_df$mean_corr[mean_corr_tabula_df$immune == "Immune cells"], 
                                                 mean_corr_tabula_df$mean_corr[mean_corr_tabula_df$immune == "Non-immune cells"])$p.val,3)))


###----- Gene-wise

immune_gene_results <- AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0002376"), columns = c('ENSEMBL'), keytype = "GOALL")
immune_gene_ids <- unique(immune_gene_results$ENSEMBL)

###----- Average dnds of increasing and decreasing immune genes

immune_cons <- data.frame(cons = c(), tissue = c(), cell_type = c(), bias = c())

for (tissue in names(Results_tabula_muris_senis)){
  
  for (cell_type in names(Results_tabula_muris_senis[[tissue]])){
    
    
    deg_list <- Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]]
    
    inc_immune_genes_dnds <- deg_list$inc_sig$inc_genes_dnds
    inc_immune_genes_dnds <- inc_immune_genes_dnds[inc_immune_genes_dnds$ensembl_gene_id %in% immune_gene_ids,]$dNdS
    
    dec_immune_genes_dnds <- deg_list$dec_sig$dec_genes_dnds
    dec_immune_genes_dnds <- dec_immune_genes_dnds[dec_immune_genes_dnds$ensembl_gene_id %in% immune_gene_ids,]$dNdS
    
    
    old_biased_cons <- mean(-log(inc_immune_genes_dnds))
    young_biased_cons <- mean(-log(dec_immune_genes_dnds))
    
    table <- rbind(data.frame(cons = old_biased_cons, tissue, cell_type, bias = "old-biased"), 
          data.frame(cons = young_biased_cons, tissue, cell_type, bias = "young-biased"))
    
    immune_cons <- rbind(immune_cons, table)
    
    
  }
  
}



immune_cons$immune_status <- "Immune genes"

#----

nonimmune_cons <-data.frame(old_biased_cons = c(), young_biased_cons = c(), tissue = c(), cell_type = c(), bias =c())

for (tissue in names(Results_tabula_muris_senis)){
  
  for (cell_type in names(Results_tabula_muris_senis[[tissue]])){
    
    
    deg_list <- Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]]
    
    inc_nonimmune_genes_dnds <- deg_list$inc_sig$inc_genes_dnds
    inc_nonimmune_genes_dnds <- inc_nonimmune_genes_dnds[!(inc_nonimmune_genes_dnds$ensembl_gene_id %in% immune_gene_ids),]$dNdS
    
    dec_nonimmune_genes_dnds <- deg_list$dec_sig$dec_genes_dnds
    dec_nonimmune_genes_dnds <- dec_nonimmune_genes_dnds[!(dec_nonimmune_genes_dnds$ensembl_gene_id %in% immune_gene_ids),]$dNdS
    
    old_biased_cons <- mean(-log(inc_nonimmune_genes_dnds))
    young_biased_cons <- mean(-log(dec_nonimmune_genes_dnds))
    
    table <- rbind(data.frame(cons = old_biased_cons, tissue, cell_type, bias = "old-biased"), 
                   data.frame(cons = young_biased_cons, tissue, cell_type, bias = "young-biased"))
    
    nonimmune_cons <- rbind(nonimmune_cons, table)
    
  }
  
}

nonimmune_cons$immune_status <- "Non-immune genes"


cons <- rbind(nonimmune_cons, immune_cons)

cons_p <- data.frame(
  
  p = c( paste0("p = ", signif(wilcox.test(
  cons$cons[cons$bias == "young-biased" & cons$immune_status == "Non-immune genes"],
  cons$cons[cons$bias == "young-biased" & cons$immune_status == "Immune genes"])$p.val,3)),
  
  paste0("p = ", signif(wilcox.test(
  cons$cons[cons$bias == "old-biased" & cons$immune_status == "Non-immune genes"],
  cons$cons[cons$bias == "old-biased" & cons$immune_status == "Immune genes"])$p.val,3))),
  
  bias = c("young-biased", "old-biased")
)




#### odds ratio


dnds <- Results_tabula_muris_senis[["lung"]][["bronchial smooth muscle cell"]][["dnds_exp_list"]][["dnds"]]
low_conserv_genes = dnds$ensembl_gene_id[dnds$dNdS > quantile(dnds$dNdS, 0.75)] 
high_conserv_genes = dnds$ensembl_gene_id[dnds$dNdS < quantile(dnds$dNdS, 0.25)]


### Fisher's test for overrepresentation of immune genes in old-biased low-conserv genes



fishers_test_results_inc <- data.frame()

for (tissue in names(Results_tabula_muris_senis)){
  
  for (cell_type in names(Results_tabula_muris_senis[[tissue]])){
    
    deg_list <- Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]]
    inc_genes <- deg_list[["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]]
    
    inc_low_conserv <- intersect(inc_genes, low_conserv_genes)
    inc_remain <- setdiff(inc_genes, low_conserv_genes)
    
    
    inc_low_conserv_immune <- intersect(inc_low_conserv, immune_gene_ids)
    inc_low_conserv_non_immune <- setdiff(inc_low_conserv, immune_gene_ids)
    
    inc_remain_immune <- intersect(inc_remain, immune_gene_ids)
    inc_remain_non_immune <- setdiff(inc_remain, immune_gene_ids)
    
    cont_table <- matrix(c(length(inc_low_conserv_immune), length(inc_low_conserv_non_immune),
                           length(inc_remain_immune),length(inc_remain_non_immune)),
                         nrow = 2,
                         dimnames =
                           list(c("immune", "non_immune"),
                                c("inc_low_conserv", "inc_remain")))
    
    
    fisher_result <- fisher.test(cont_table)
    values <- data.frame(tissue, cell_type, fisher_result$estimate, fisher_result$p.value)
    colnames(values) <- c("tissue", "cell_type", "odds_ratio", "p_val")
    
    fishers_test_results_inc <- rbind(fishers_test_results_inc, values)
    
    
  }
  
}

fishers_test_results_inc$p_adj <- p.adjust(fishers_test_results_inc$p_val, "BH")
fishers_test_results_inc$bias <- "old-biased"

fishers_test_results_dec <- data.frame()

for (tissue in names(Results_tabula_muris_senis)){
  
  for (cell_type in names(Results_tabula_muris_senis[[tissue]])){
    
    deg_list <- Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]]
    dec_genes <- deg_list[["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]]
    
    dec_high_conserv <- intersect(dec_genes, high_conserv_genes)
    dec_remain <- setdiff(dec_genes, high_conserv_genes)
    
    
    dec_high_conserv_immune <- intersect(dec_high_conserv, immune_gene_ids)
    dec_high_conserv_non_immune <- setdiff(dec_high_conserv, immune_gene_ids)
    
    dec_remain_immune <- intersect(dec_remain, immune_gene_ids)
    dec_remain_non_immune <- setdiff(dec_remain, immune_gene_ids)
    
    cont_table <- matrix(c(length(dec_high_conserv_immune), length(dec_high_conserv_non_immune),
                           length(dec_remain_immune),length(dec_remain_non_immune)),
                           nrow = 2,
                         dimnames =
                           list(c("immune", "non_immune"),
                                c("dec_high_conserv", "dec_remain")))
    
    
    fisher_result <- fisher.test(cont_table)
    values <- data.frame(tissue, cell_type, fisher_result$estimate, fisher_result$p.value)
    colnames(values) <- c("tissue", "cell_type", "odds_ratio", "p_val")
    
    fishers_test_results_dec <- rbind(fishers_test_results_dec, values)
    
    
  }
  
}
fishers_test_results_dec$p_adj <- p.adjust(fishers_test_results_dec$p_val, "BH")
fishers_test_results_dec$bias <- "young-biased"



fishers_test_results <- rbind(fishers_test_results_inc, fishers_test_results_dec)

mean(fishers_test_results_inc$odds_ratio)
mean(fishers_test_results_inc$odds_ratio[fishers_test_results_inc$p_adj < 0.1])

mean(fishers_test_results_dec$odds_ratio[fishers_test_results_dec$odds_ratio != Inf])





### mosaic plot example

tissue <- "brain"
cell_type <- "neuron"

deg_list <- Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]]
inc_genes <- deg_list[["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]]

inc_low_conserv <- intersect(inc_genes, low_conserv_genes)
inc_remain <- setdiff(inc_genes, low_conserv_genes)


inc_low_conserv_immune <- intersect(inc_low_conserv, immune_gene_ids)
inc_low_conserv_non_immune <- setdiff(inc_low_conserv, immune_gene_ids)

inc_remain_immune <- intersect(inc_remain, immune_gene_ids)
inc_remain_non_immune <- setdiff(inc_remain, immune_gene_ids)

cont_table <- matrix(c(length(inc_low_conserv_immune), length(inc_low_conserv_non_immune),
                       length(inc_remain_immune),length(inc_remain_non_immune)),
                     nrow = 2,
                     dimnames =
                       list(c("immune", "non_immune"),
                            c("inc_low_conserv", "inc_remain")))

cont_df <- rbind(data.frame(gene_id = inc_low_conserv_immune, immune_status = "Immune", category= "low conservation"),
data.frame(gene_id =inc_remain_immune, immune_status = "Immune", category= "non-low conservation"),
data.frame(gene_id =inc_low_conserv_non_immune, immune_status = "Non-immune", category= "low conservation"),
data.frame(gene_id =inc_remain_non_immune, immune_status = "Non-immune", category= "non-low conservation"))

cont_df_p <- data.frame(
  p = paste0("p = ", 
             signif(fishers_test_results_inc$p_adj[fishers_test_results_inc$tissue == "brain" & 
                                          fishers_test_results_inc$cell_type == "neuron"], 3)),
  odds_ratio = paste0("OR = ",
                      signif(fishers_test_results_inc$odds_ratio[fishers_test_results_inc$tissue == "brain" & 
                                                              fishers_test_results_inc$cell_type == "neuron"], 3)))



write.csv(fishers_test_results, file = "supplements/immune_odds_ratios.csv", quote = T, row.names = F)

###---- PLOTS


first= 
  ggplot(mean_cons_diff_tabula, aes(x = immune, y = mean_cons_diff), fill = immune) +
  theme_pubr(base_size = 12)+
  geom_boxplot(width = 0.3, lwd=0.4, aes(fill = immune), alpha = 0.4) + 
  # , aes(fill = type), alpha = 0.2
  # geom_boxplot(lwd=0.4, color=c("#BEBADA","#8DD3C7"), outlier.shape = NA) + 
  ylab("Young-biased - Old-biased MRCS")+
  xlab("")+
  ggtitle("Dist. of MRCS differences for different cell types")+
  geom_hline(yintercept = 0)+
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
  ggtitle("Dist. of expression-conservation correlation vs. age rho values")+
  xlab("")+
  geom_hline(yintercept = 0)+
  geom_jitter(width = 0.1, color ='gray36' , cex=1.7, alpha = 0.6)+
  scale_color_manual(values = c("#BEBADA","#8DD3C7"), )+
  scale_fill_manual(values = c("#BEBADA","#8DD3C7"))+
  geom_text(data    =immune_non_immune_rho_p,
            mapping = aes(x = Inf, y = -Inf, label = p),vjust = -1.1,hjust = 1.1,inherit.aes = FALSE,size = 3.4)+
  theme(legend.position = "none") 

third <- ggplot(data= cons, 
                       aes(x = immune_status, cons)) +
  geom_boxplot(width = 0.3, lwd=0.4, aes(fill = immune_status), alpha = 0.4)+
  scale_y_continuous(name="Mean Conservation")+
  theme_bw()+
  xlab("")+
  ggtitle("Mean conservation distribution of gene sets across all cell types")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  #  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_fill_manual(values = c("#BEBADA","#8DD3C7"))+
  theme_pubr(base_size = 12)+
  #facet_grid2(~ tissue, scales = "free_x", independent = "x")+
  facet_grid(~ bias)+
  guides(fill = guide_legend(override.aes = list(size = 7, shape = 22)))+
  geom_text(data    =cons_p,
           mapping = aes(x = Inf, y = -Inf, label = p),vjust = -1.1,hjust = 1.1,inherit.aes = FALSE, size = 3.4)+
  theme(legend.position = "none")

fourth <- ggplot(data = cont_df) +
  geom_mosaic(aes(x = product(immune_status, category), fill=immune_status))+
  theme_bw()+
  xlab("")+
  ggtitle("Neuron old-biased genes")+
  ylab("Proportion of genes")+
  theme_pubr(base_size = 12)+
  scale_fill_manual(values = c("#BEBADA","#8DD3C7"))+
  geom_text(data    =cont_df_p,
            mapping = aes(x = Inf, y = Inf, label = odds_ratio),vjust = 3,hjust = 1.5,inherit.aes = FALSE,size = 3.4)+
  guides(fill = guide_legend(override.aes = list(size = 7, shape = 22))) +
  geom_text(data    =cont_df_p,
            mapping = aes(x = Inf, y = Inf, label = p),vjust = 4.5,hjust = 1.45,inherit.aes = FALSE,size = 3.4)+
  guides(fill = guide_legend(override.aes = list(size = 7, shape = 22))) +
  theme(legend.position = "none")


fourth
fifth <- ggplot(data= fishers_test_results, 
                aes(x = bias, y = odds_ratio)) +
  geom_boxplot(width = 0.3, lwd=0.4, aes(fill = bias), alpha = 0.4)+
  scale_y_continuous(name="Odds ratio for immune genes")+
  theme_bw()+
  xlab("")+
  ggtitle("Odds ratio distribution across all cell types")+
  geom_hline(yintercept = 1)+
  scale_x_discrete(guide = guide_axis(n.dodge=2),
                   labels=c("old-biased" = "Old-biased (low conserv. vs. others)",
                            "young-biased" = "Young-biased (high conserv. vs. others)"))+
  #  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_fill_manual(values=c("#888888", "#6699CC"))+
  theme_pubr(base_size = 12)+
  #facet_grid2(~ tissue, scales = "free_x", independent = "x")+
  guides(fill = guide_legend(override.aes = list(size = 7, shape = 22)))+
  #geom_text(data    =cons_p,
  #          mapping = aes(x = Inf, y = -Inf, label = p),vjust = -1.1,hjust = 1.1,inherit.aes = FALSE, size = 3.4)+
  theme(legend.position = "none")

#layout_matrix <- matrix(c(1, 1, 2, 2, 4, 3, 3, 4), nrow = 2, byrow = TRUE)

#grid.arrange(first,  second, third, layout_matrix = layout_matrix)

top_row = ggarrange(first,  second, labels = c('A','B'), 
                 nrow = 1, ncol =2)+
  theme(plot.margin = unit(c(0.55, 0.2, 0.001, 0.1), "cm"))

empty <- ggplot() + theme_void()

mid_row <- ggarrange(third,fourth, labels = c('C','D'),
                        nrow = 1, ncol =2)+
  theme(plot.margin = unit(c(0.55, 0.2, 0.001, 0.1), "cm"))

bottom_row <- ggarrange(empty, fifth, empty, labels = c("","E",""), nrow = 1, ncol =3, widths = c(0.25, 0.5, 0.25))+
  theme(plot.margin = unit(c(0.55, 0.2, 0.001, 0.1), "cm"))
  

final <- ggarrange(top_row, mid_row, bottom_row,ncol=1,nrow=3) + bgcolor("white")

final
ggsave(final, 
       filename = "results_graphs/figure7-immune-diff.pdf",
       unit="cm", width = 33, height = 33, useDingbats = F,limitsize = FALSE)

ggsave(final, 
       filename = "results_graphs/figure7-immune-diff.png",
       device = png, type="cairo", unit="cm", width = 33, height = 33)

ggsave(final, 
       filename = "results_graphs/figure7-immune-diff.tiff",
       type='cairo',
       unit="cm", width = 33, height = 33)
