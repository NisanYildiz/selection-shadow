#supp. fig x

#Loading the required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(latex2exp)

#Loading up the Results 
load("gallus_gallus/R/results/Results_GSE114129_rho_03.RData") # chicken
load("drosophila_melanogaster/R/results/Results_pacifico_rho_03.RData") # drosophila
load("n_furzeri/GSE66712/R/results/Results_GSE66712_rho_03.RData") # killifish

#loading general functions
source("functions.R")

Results_rho_03 <- list(Results_GSE114129 = Results_GSE114129, 
                       Results_GSE66712 = Results_GSE66712, 
                       Results_pacifico = Results_pacifico)

load("gallus_gallus/R/results/Results_GSE114129_rho_07.RData") # chicken
load("drosophila_melanogaster/R/results/Results_pacifico_rho_07.RData") # drosophila
load("n_furzeri/GSE66712/R/results/Results_GSE66712_rho_07.RData") # killifish


Results_rho_07 <- list(Results_GSE114129 = Results_GSE114129, 
                       Results_GSE66712 = Results_GSE66712, 
                       Results_pacifico = Results_pacifico)

## A 


rho_chick = c(signif( cor(as.numeric(Results_rho_03[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results_rho_03[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ), 
              signif( cor(as.numeric(Results_rho_03[["Results_GSE114129"]][["brain"]][["all_results"]]), Results_rho_03[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ) )
p_chick = c(signif( cor.test(as.numeric(Results_rho_03[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results_rho_03[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
            signif( cor.test(as.numeric(Results_rho_03[["Results_GSE114129"]][["brain"]][["all_results"]]), Results_rho_03[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_chick = c(length(Results_rho_03[["Results_GSE114129"]][["brain"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
            length(Results_rho_03[["Results_GSE114129"]][["brain"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fish_liv = c(signif( cor(as.numeric(Results_rho_03[["Results_GSE66712"]][["liver"]][["sig_results"]]), Results_rho_03[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman"), 2 ), 
                 signif( cor(as.numeric(Results_rho_03[["Results_GSE66712"]][["liver"]][["all_results"]]), Results_rho_03[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman"), 2 ) )
p_fish_liv = c(signif( cor.test(as.numeric(Results_rho_03[["Results_GSE66712"]][["liver"]][["sig_results"]]), Results_rho_03[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
               signif( cor.test(as.numeric(Results_rho_03[["Results_GSE66712"]][["liver"]][["all_results"]]), Results_rho_03[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fish_liv = c(length(Results_rho_03[["Results_GSE66712"]][["liver"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
               length(Results_rho_03[["Results_GSE66712"]][["liver"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fish_skin = c(signif( cor(as.numeric(Results_rho_03[["Results_GSE66712"]][["skin"]][["sig_results"]]), Results_rho_03[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman"), 2 ) , 
                  signif( cor(as.numeric(Results_rho_03[["Results_GSE66712"]][["skin"]][["all_results"]]), Results_rho_03[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman"), 2 ) )
p_fish_skin = c(signif( cor.test(as.numeric(Results_rho_03[["Results_GSE66712"]][["skin"]][["sig_results"]]), Results_rho_03[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
                signif( cor.test(as.numeric(Results_rho_03[["Results_GSE66712"]][["skin"]][["all_results"]]), Results_rho_03[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fish_skin = c(length(Results_rho_03[["Results_GSE66712"]][["skin"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
                length(Results_rho_03[["Results_GSE66712"]][["skin"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fly_brain = c(signif( cor(as.numeric(Results_rho_03[["Results_pacifico"]][["head"]][["sig_results"]]), Results_rho_03[["Results_pacifico"]][["head"]][["age"]], method = "spearman"), 2 ), 
                  signif( cor(as.numeric(Results_rho_03[["Results_pacifico"]][["head"]][["all_results"]]), Results_rho_03[["Results_pacifico"]][["head"]][["age"]], method = "spearman"), 2 ) )
p_fly_brain = c(signif( cor.test(as.numeric(Results_rho_03[["Results_pacifico"]][["head"]][["sig_results"]]), Results_rho_03[["Results_pacifico"]][["head"]][["age"]], method = "spearman", exact = F)$p.val, 2 ),
                signif( cor.test(as.numeric(Results_rho_03[["Results_pacifico"]][["head"]][["all_results"]]), Results_rho_03[["Results_pacifico"]][["head"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fly_brain = c(length(Results_rho_03[["Results_pacifico"]][["head"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
                length(Results_rho_03[["Results_pacifico"]][["head"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_vals_03 <- data.frame(organism = c("G. gallus", "G. gallus", "N. furzeri", "N. furzeri", 
                                    " N. furzeri ", " N. furzeri ", "D. melanogaster", "D. melanogaster"),
                       rho = c(rho_chick, rho_fish_liv, rho_fish_skin, rho_fly_brain),
                       Tissue = rep(c("Brain", "Liver", "Skin", "Brain"), each = 2),
                       type = factor(rep(c("Differentially Expressed Genes", "All Genes")), levels = c("All Genes", "Differentially Expressed Genes")),
                       p = c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain),
                       p_adj = p.adjust(c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain), "BH"),
                       sig = sig_symbols(p.adjust(c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain), "BH"), ""))
rho_vals_03$organism = factor(rho_vals_03$organism, levels = c("G. gallus", "D. melanogaster", "N. furzeri", " N. furzeri "))


## B 


rho_chick = c(signif( cor(as.numeric(Results_rho_07[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results_rho_07[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ), 
              signif( cor(as.numeric(Results_rho_07[["Results_GSE114129"]][["brain"]][["all_results"]]), Results_rho_07[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ) )
p_chick = c(signif( cor.test(as.numeric(Results_rho_07[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results_rho_07[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
            signif( cor.test(as.numeric(Results_rho_07[["Results_GSE114129"]][["brain"]][["all_results"]]), Results_rho_07[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_chick = c(length(Results_rho_07[["Results_GSE114129"]][["brain"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
            length(Results_rho_07[["Results_GSE114129"]][["brain"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fish_liv = c(signif( cor(as.numeric(Results_rho_07[["Results_GSE66712"]][["liver"]][["sig_results"]]), Results_rho_07[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman"), 2 ), 
                 signif( cor(as.numeric(Results_rho_07[["Results_GSE66712"]][["liver"]][["all_results"]]), Results_rho_07[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman"), 2 ) )
p_fish_liv = c(signif( cor.test(as.numeric(Results_rho_07[["Results_GSE66712"]][["liver"]][["sig_results"]]), Results_rho_07[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
               signif( cor.test(as.numeric(Results_rho_07[["Results_GSE66712"]][["liver"]][["all_results"]]), Results_rho_07[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fish_liv = c(length(Results_rho_07[["Results_GSE66712"]][["liver"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
               length(Results_rho_07[["Results_GSE66712"]][["liver"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fish_skin = c(signif( cor(as.numeric(Results_rho_07[["Results_GSE66712"]][["skin"]][["sig_results"]]), Results_rho_07[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman"), 2 ) , 
                  signif( cor(as.numeric(Results_rho_07[["Results_GSE66712"]][["skin"]][["all_results"]]), Results_rho_07[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman"), 2 ) )
p_fish_skin = c(signif( cor.test(as.numeric(Results_rho_07[["Results_GSE66712"]][["skin"]][["sig_results"]]), Results_rho_07[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
                signif( cor.test(as.numeric(Results_rho_07[["Results_GSE66712"]][["skin"]][["all_results"]]), Results_rho_07[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fish_skin = c(length(Results_rho_07[["Results_GSE66712"]][["skin"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
                length(Results_rho_07[["Results_GSE66712"]][["skin"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fly_brain = c(signif( cor(as.numeric(Results_rho_07[["Results_pacifico"]][["head"]][["sig_results"]]), Results_rho_07[["Results_pacifico"]][["head"]][["age"]], method = "spearman"), 2 ), 
                  signif( cor(as.numeric(Results_rho_07[["Results_pacifico"]][["head"]][["all_results"]]), Results_rho_07[["Results_pacifico"]][["head"]][["age"]], method = "spearman"), 2 ) )
p_fly_brain = c(signif( cor.test(as.numeric(Results_rho_07[["Results_pacifico"]][["head"]][["sig_results"]]), Results_rho_07[["Results_pacifico"]][["head"]][["age"]], method = "spearman", exact = F)$p.val, 2 ),
                signif( cor.test(as.numeric(Results_rho_07[["Results_pacifico"]][["head"]][["all_results"]]), Results_rho_07[["Results_pacifico"]][["head"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fly_brain = c(length(Results_rho_07[["Results_pacifico"]][["head"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
                length(Results_rho_07[["Results_pacifico"]][["head"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_vals_07 <- data.frame(organism = c("G. gallus", "G. gallus", "N. furzeri", "N. furzeri", 
                                       " N. furzeri ", " N. furzeri ", "D. melanogaster", "D. melanogaster"),
                          rho = c(rho_chick, rho_fish_liv, rho_fish_skin, rho_fly_brain),
                          Tissue = rep(c("Brain", "Liver", "Skin", "Brain"), each = 2),
                          type = factor(rep(c("Differentially Expressed Genes", "All Genes")), levels = c("All Genes", "Differentially Expressed Genes")),
                          p = c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain),
                          p_adj = p.adjust(c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain), "BH"),
                          sig = sig_symbols(p.adjust(c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain), "BH"), ""))
rho_vals_07$organism = factor(rho_vals_03$organism, levels = c("G. gallus", "D. melanogaster", "N. furzeri", " N. furzeri "))


df_03 <- rho_vals_03[rho_vals_03$type == "Differentially Expressed Genes",]
df_03$type = "DEG rho cutoff = 0.3"
df_07 <- rho_vals_07[rho_vals_07$type == "Differentially Expressed Genes",]
df_07$type = "DEG rho cutoff = 0.7"

df_combined = rbind(df_03, df_07)

figure_s3 <- ggplot(data = df_combined, aes(x = organism, y = rho)) +
  geom_bar(stat = "identity",alpha = 0.8, aes(fill = Tissue), colour = "black", show.legend = T)+
  theme_bw()+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_y_continuous(limits = c(-1, 0.25))+
  xlab("")+
  ggtitle("")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  ylab("Expression-Conservation Correlation vs. Age rho")+
  geom_hline(yintercept = 0)+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 13,  face = "italic"),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size = 13),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~type)+
  geom_text(aes(label=sig), position=position_dodge(width=0.9), size=6, vjust = 1.6)


figure_s3

ggsave("suppFig3-deg-rho-cutoffs.png", plot = figure_s3, device = png,
       path = "results_graphs/", scale = 1, 
       width = 10, height = 4.5, units = "in", dpi = 300)
ggsave("suppFig3-deg-rho-cutoffs.tiff", plot = figure_s3, device = "tiff",
       path = "results_graphs/",
       scale = 1, width = 1, height = 4.5, units = "in", dpi = 300)
ggsave("suppFig3-deg-rho-cutoffs.pdf", plot = figure_s3, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 4.5, units = "in", dpi = 300)


