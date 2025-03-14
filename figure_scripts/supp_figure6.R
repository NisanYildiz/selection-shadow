# Code to produce Supp Figure 6 

#Loading the required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)


#Loading up the Results 
load("mus_musculus/GSE99791/R/results/Results_GSE99791.RData") #mouse astrocyte

#Setting up the data
GSE99791_df <- data.frame(results = c(as.numeric(Results_GSE99791[["hypothalamus"]][["sig_results"]]),
                                      as.numeric(Results_GSE99791[["hypothalamus"]][["all_results"]])), 
                          age = rep(as.numeric(Results_GSE99791[["hypothalamus"]][["age"]]), 2),
                          type = rep(c("Differentially Expressed Genes", 
                                       "All Genes"), each = 6))
GSE99791_df$type <- factor(GSE99791_df$type, 
                           levels = c("Differentially Expressed Genes",
                                      "All Genes"))
GSE99791_df_text <- data.frame(
  rho = c(paste0("rho = ", signif( cor(as.numeric(Results_GSE99791[["hypothalamus"]][["sig_results"]]), Results_GSE99791[["hypothalamus"]][["age"]], method = "spearman"), 2 ) ), 
          paste0("rho = ", signif( cor(as.numeric(Results_GSE99791[["hypothalamus"]][["all_results"]]), Results_GSE99791[["hypothalamus"]][["age"]], method = "spearman"), 2 ) )),
  p = c(paste0("p = ", signif( cor.test(as.numeric(Results_GSE99791[["hypothalamus"]][["sig_results"]]), Results_GSE99791[["hypothalamus"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) ),
        paste0("p = ", signif( cor.test(as.numeric(Results_GSE99791[["hypothalamus"]][["all_results"]]), Results_GSE99791[["hypothalamus"]][["age"]], method = "spearman", exact = F)$p.val, 2 ) )),
  gene_no = c(paste0(length(Results_GSE99791[["hypothalamus"]][["deg_list"]][["All_sig"]][["sig_genes_exp"]][,1]), " genes"),
              paste0(length(Results_GSE99791[["hypothalamus"]][["deg_list"]][["All"]][["all_genes_exp"]][,1]), " genes")),
  type = c("Differentially Expressed Genes",
           "All Genes")
)
GSE99791_df_text_wilcox <- data.frame(
  rho = c(paste0("p = ", signif( wilcox.test(as.numeric(Results_GSE99791[["hypothalamus"]][["sig_results"]][Results_GSE99791[["hypothalamus"]][["age"]] == 4]), as.numeric(Results_GSE99791[["hypothalamus"]][["sig_results"]][Results_GSE99791[["hypothalamus"]][["age"]] == 24]), alternative = "greater")$p.val, 2 ) ), 
          paste0("p = ", signif( wilcox.test(as.numeric(Results_GSE99791[["hypothalamus"]][["all_results"]][Results_GSE99791[["hypothalamus"]][["age"]] == 4]), as.numeric(Results_GSE99791[["hypothalamus"]][["all_results"]][Results_GSE99791[["hypothalamus"]][["age"]] == 24]), alternative = "greater")$p.val, 2 ) )),
  gene_no = c(paste0(length(Results_GSE99791[["hypothalamus"]][["deg_list"]][["All_sig"]][["sig_genes_exp"]][,1]), " genes"),
              paste0(length(Results_GSE99791[["hypothalamus"]][["deg_list"]][["All"]][["all_genes_exp"]][,1]), " genes")),
  type = c("Differentially Expressed Genes",
           "All Genes")
)

sup_figure6 <- ggplot(data = GSE99791_df, aes(age, results)) +
  geom_point(size = 4, alpha = 0.8, aes(colour = "#1B9E77"), show.legend = F)+
  theme_bw()+
  stat_smooth(method = "lm", col = "#1B9E77")+
  scale_color_manual(values = c("#1B9E77"))+
  ggtitle(expression(paste(italic("M. musculus")," astrocytes from hypothalamus")))+
  xlab("Age (months)")+
  ylab("Expression-Conservation Correlation")+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~type)+
  geom_text(data    = GSE99791_df_text,
            mapping = aes(x = -Inf, y = -Inf, label = p),
            hjust   = -0.12,
            vjust   = -3,
            size = 4.4)+
  geom_text(data    = GSE99791_df_text,
            mapping = aes(x = -Inf, y = -Inf, label = rho),
            hjust   = -0.1,
            vjust   = -1,
            size = 4.4)+
  geom_text(data    = GSE99791_df_text,
            mapping = aes(x = Inf, y = -Inf, label = gene_no),
            vjust = -1,
            hjust = 1.1,
            size = 4.4)

ggsave("suppFig6-astrocyte-hypothalamus-transcriptome-cons-vs-age.png", plot = sup_figure6, device = png,
       path = "results_graphs/",
       scale = 1, width = 10, height = 5, units = "in", dpi = 300)
ggsave("suppFig6-astrocyte-hypothalamus-transcriptome-cons-vs-age.tiff", plot = sup_figure6, device = "tiff",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 10, height = 5, units = "in", dpi = 300)
ggsave("suppFig6-astrocyte-hypothalamus-transcriptome-cons-vs-age.pdf", plot = sup_figure6, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 5, units = "in", dpi = 300)
