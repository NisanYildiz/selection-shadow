# Code used to produce Figure 1

#Loading the required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)

#Loading up the Results 
load("gallus_gallus/R/results/Results_GSE114129.RData") # chicken
load("drosophila_melanogaster/R/results/Results_pacifico.RData") # drosophila
load("n_furzeri/GSE66712/R/results/Results_GSE66712.RData") # killifish
load("mus_musculus/GSE99791/R/results/Results_GSE99791.RData") #mouse astrocyte
load("naked_mole_rat/GSE30337/R/results/Results_GSE30337.RData") #NMR 

Results <- list(Results_GSE114129 = Results_GSE114129, Results_GSE30337 = Results_GSE30337,
                Results_GSE66712 = Results_GSE66712, Results_GSE99791 = Results_GSE99791,
                Results_pacifico = Results_pacifico)

#loading general functions
source("functions.R")

#Setting up the data

#For Figure 1A

gallus_gene_names <- rownames(Results_GSE114129[["brain"]][["dnds_exp_list"]][["exp_matrix"]])
gallus_DE_gene_names <- rownames(Results_GSE114129[["brain"]][["deg_list"]][["All_sig"]][["sig_genes_exp"]])

GSE114129_exp_dnds_cor <- data.frame(exp = c(log(Results_GSE114129[["brain"]][["dnds_exp_list"]][["exp_matrix"]][(gallus_gene_names %in% gallus_DE_gene_names),1] + 10),
                                             log(Results_GSE114129[["brain"]][["dnds_exp_list"]][["exp_matrix"]][(gallus_gene_names %in% gallus_DE_gene_names),13] + 10)),
                                     cons= c(-log(Results_GSE114129[["brain"]][["dnds_exp_list"]][["dnds"]][gallus_gene_names %in% gallus_DE_gene_names,][["dNdS"]]),
                                             -log(Results_GSE114129[["brain"]][["dnds_exp_list"]][["dnds"]][gallus_gene_names %in% gallus_DE_gene_names,][["dNdS"]])),
                                     type= rep(c("100 days old", 
                                                 "1825 days old"), each = length(gallus_DE_gene_names)))


GSE114129_exp_dnds_cor_text <- data.frame(rho = c(paste0("rho = ", signif( cor(GSE114129_exp_dnds_cor$exp[GSE114129_exp_dnds_cor$type == "100 days old"], GSE114129_exp_dnds_cor$cons[GSE114129_exp_dnds_cor$type == "100 days old"], method = "spearman"), 2 ) ),
                                                  paste0("rho = ", signif( cor(GSE114129_exp_dnds_cor$exp[GSE114129_exp_dnds_cor$type == "1825 days old"], GSE114129_exp_dnds_cor$cons[GSE114129_exp_dnds_cor$type == "1825 days old"], method = "spearman"), 2 ) )),
                                          p = c(signif(cor.test(GSE114129_exp_dnds_cor$exp[GSE114129_exp_dnds_cor$type == "100 days old"], GSE114129_exp_dnds_cor$cons[GSE114129_exp_dnds_cor$type == "100 days old"], method = "spearman", exact = F)$p.val, 2),
                                                signif(cor.test(GSE114129_exp_dnds_cor$exp[GSE114129_exp_dnds_cor$type == "1825 days old"], GSE114129_exp_dnds_cor$cons[GSE114129_exp_dnds_cor$type == "1825 days old"], method = "spearman", exact = F)$p.val, 2)),
                                          gene_no = c(paste0(length(GSE114129_exp_dnds_cor$exp[GSE114129_exp_dnds_cor$type == "100 days old"]), " genes"),
                                                      paste0(length(GSE114129_exp_dnds_cor$exp[GSE114129_exp_dnds_cor$type == "1825 days old"]), " genes")),
                                          type = c("100 days old", "1825 days old"))

GSE114129_exp_dnds_cor_text$p_sci <- paste0("p == ", sub("e", "%*%10^", as.character(GSE114129_exp_dnds_cor_text$p)))


##For Figure 1B

GSE114129_df <- data.frame(results = c(as.numeric(Results[["Results_GSE114129"]][["brain"]][["sig_results"]]),
                                           as.numeric(Results[["Results_GSE114129"]][["brain"]][["all_results"]])), 
                               age = rep(as.numeric(Results[["Results_GSE114129"]][["brain"]][["age"]]), 2),
                               type = rep(c("Differentially Expressed Genes", 
                                            "All Genes"), each = 13))
GSE114129_df$type <- factor(GSE114129_df$type, 
                                levels = c("Differentially Expressed Genes",
                                           "All Genes"))

### rho and p-values of transcriptome conservation vs. age
GSE114129_df_text <- data.frame(
  rho = c(paste0("rho = ", signif( cor(as.numeric(Results[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ) ), 
          paste0("rho = ", signif( cor(as.numeric(Results[["Results_GSE114129"]][["brain"]][["all_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ) )),
  p = c(signif(cor.test(as.numeric(Results[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2),
        signif( cor.test(as.numeric(Results[["Results_GSE114129"]][["brain"]][["all_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2)),
  gene_no = c(paste0(length(Results[["Results_GSE114129"]][["brain"]][["deg_list"]][["All_sig"]][["sig_genes_exp"]][,1]), " genes"),
              paste0(length(Results[["Results_GSE114129"]][["brain"]][["deg_list"]][["All"]][["all_genes_exp"]][,1]), " genes")),
  type = c("Differentially Expressed Genes",
           "All Genes")
)

GSE114129_df_text$p_sci <- paste0("p == ", sub("e", "%*%10^", as.character(GSE114129_df_text$p)))

##For Figure 1B

rho_chick = c(signif( cor(as.numeric(Results[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ), 
              signif( cor(as.numeric(Results[["Results_GSE114129"]][["brain"]][["all_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman"), 2 ) )
p_chick = c(signif( cor.test(as.numeric(Results[["Results_GSE114129"]][["brain"]][["sig_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
            signif( cor.test(as.numeric(Results[["Results_GSE114129"]][["brain"]][["all_results"]]), Results[["Results_GSE114129"]][["brain"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_chick = c(length(Results[["Results_GSE114129"]][["brain"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
            length(Results[["Results_GSE114129"]][["brain"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fish_liv = c(signif( cor(as.numeric(Results[["Results_GSE66712"]][["liver"]][["sig_results"]]), Results[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman"), 2 ), 
                 signif( cor(as.numeric(Results[["Results_GSE66712"]][["liver"]][["all_results"]]), Results[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman"), 2 ) )
p_fish_liv = c(signif( cor.test(as.numeric(Results[["Results_GSE66712"]][["liver"]][["sig_results"]]), Results[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
               signif( cor.test(as.numeric(Results[["Results_GSE66712"]][["liver"]][["all_results"]]), Results[["Results_GSE66712"]][["liver"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fish_liv = c(length(Results[["Results_GSE66712"]][["liver"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
               length(Results[["Results_GSE66712"]][["liver"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fish_skin = c(signif( cor(as.numeric(Results[["Results_GSE66712"]][["skin"]][["sig_results"]]), Results[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman"), 2 ) , 
                  signif( cor(as.numeric(Results[["Results_GSE66712"]][["skin"]][["all_results"]]), Results[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman"), 2 ) )
p_fish_skin = c(signif( cor.test(as.numeric(Results[["Results_GSE66712"]][["skin"]][["sig_results"]]), Results[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ,
                signif( cor.test(as.numeric(Results[["Results_GSE66712"]][["skin"]][["all_results"]]), Results[["Results_GSE66712"]][["skin"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fish_skin = c(length(Results[["Results_GSE66712"]][["skin"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
                length(Results[["Results_GSE66712"]][["skin"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_fly_brain = c(signif( cor(as.numeric(Results[["Results_pacifico"]][["head"]][["sig_results"]]), Results[["Results_pacifico"]][["head"]][["age"]], method = "spearman"), 2 ), 
                  signif( cor(as.numeric(Results[["Results_pacifico"]][["head"]][["all_results"]]), Results[["Results_pacifico"]][["head"]][["age"]], method = "spearman"), 2 ) )
p_fly_brain = c(signif( cor.test(as.numeric(Results[["Results_pacifico"]][["head"]][["sig_results"]]), Results[["Results_pacifico"]][["head"]][["age"]], method = "spearman", exact = F)$p.val, 2 ),
                signif( cor.test(as.numeric(Results[["Results_pacifico"]][["head"]][["all_results"]]), Results[["Results_pacifico"]][["head"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )
n_fly_brain = c(length(Results[["Results_pacifico"]][["head"]][["deg_list"]][["All_sig"]][["sig_genes_cor"]][,1]),
                length(Results[["Results_pacifico"]][["head"]][["deg_list"]][["All"]][["all_genes_cor"]][,1]))

rho_vals <- data.frame(organism = c("G. gallus", "G. gallus", "N. furzeri", "N. furzeri", 
                                    " N. furzeri ", " N. furzeri ", "D. melanogaster", "D. melanogaster"),
                       rho = c(rho_chick, rho_fish_liv, rho_fish_skin, rho_fly_brain),
                       Tissue = rep(c("Brain", "Liver", "Skin", "Brain"), each = 2),
                       type = factor(rep(c("Differentially Expressed Genes", "All Genes")), levels = c("All Genes", "Differentially Expressed Genes")),
                       p = c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain),
                       p_adj = p.adjust(c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain), "BH"),
                       sig = sig_symbols(p.adjust(c(p_chick, p_fish_liv, p_fish_skin, p_fly_brain), "BH"), ""))
rho_vals$organism = factor(rho_vals$organism, levels = c("G. gallus", "D. melanogaster", "N. furzeri", " N. furzeri "))

#Plotting 

## Figure 1A
exp_dnds <- 
  ggplot(data = GSE114129_exp_dnds_cor, aes(cons, exp)) +
  geom_point(size = 2, alpha = 0.8, aes(colour = "#1B9E77"), show.legend = F)+
  theme_bw()+
  stat_smooth(method = "lm", col = "#000000")+
  scale_color_manual(values = c("#1B9E77"))+
  ggtitle(expression(paste(italic("G. gallus")," brain, DE genes")))+
  xlab("Conservation Score")+
  ylab("log(Expression + 10)")+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size = 13),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~type)+
  geom_text(data    = GSE114129_exp_dnds_cor_text,
            mapping = aes(x = Inf, y = -Inf, label = rho),
            hjust   = 1.1,
            vjust   = -1,
            size = 4.4)
#  geom_text(data    = GSE114129_exp_dnds_cor_text,
#            mapping = aes(x = -Inf, y = -Inf, label = p_sci),
#            hjust   = -0.11,
#            vjust   = -1.6,
#            size = 4.4,
#            parse = T)+
#  geom_text(data    = GSE114129_exp_dnds_cor_text,
#            mapping = aes(x = Inf, y = -Inf, label = gene_no),
#            hjust = 1.1,
#            vjust = -1,
#            size = 4.4)

## Figure 1B
genes_cor <- ggplot(data = GSE114129_df, aes(age, results)) +
  geom_point(size = 4, alpha = 0.8, aes(colour = "#1B9E77"), show.legend = F)+
  theme_bw()+
  stat_smooth(method = "lm", col = "#1B9E77")+
  scale_color_manual(values = c("#1B9E77"))+
  ggtitle(expression(paste(italic("G. gallus")," brain")))+
  xlab("Age (days)")+
  ylab("Expression-Conservation Correlation Score")+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size = 13),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~type)+
  geom_text(data    = GSE114129_df_text,
            mapping = aes(x = -Inf, y = -Inf, label = rho),
            hjust   = -0.1,
            vjust   = -1,
            size = 4.4)+
  geom_text(data    = GSE114129_df_text,
            mapping = aes(x = -Inf, y = -Inf, label = p_sci),
            hjust   = -0.14,
            vjust   = -1.6,
            size = 4.4,
            parse = T)+
  geom_text(data    = GSE114129_df_text,
            mapping = aes(x = Inf, y = -Inf, label = gene_no),
            hjust = 1.1,
            vjust = -1,
            size = 4.4)


## Figure 1C
rho_barplot <- ggplot(data = rho_vals, aes(x = organism, y = rho)) +
  geom_bar(stat = "identity",alpha = 0.8, aes(fill = Tissue), colour = "black", show.legend = T)+
  theme_bw()+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_y_continuous(limits = c(-1, 0.25))+
  xlab("")+
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

## Combined Figure 1

figure1 <- exp_dnds / genes_cor / rho_barplot + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16))

ggsave("figure1-transcriptome-cons-vs-age.png", plot = figure1, device = png,
       path = "results_graphs/", scale = 1, 
       width = 10, height = 13, units = "in", dpi = 300)
ggsave("figure1-transcriptome-cons-vs-age.tiff", plot = figure1, device = "tiff",
       path = "results_graphs/",
       scale = 1, width = 1, height = 13, units = "in", dpi = 300)
ggsave("figure1-transcriptome-cons-vs-age.pdf", plot = figure1, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 13, units = "in", dpi = 300)

