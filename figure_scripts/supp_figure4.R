
#Loading required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

#Loading up Results 

load("naked_mole_rat/GSE30337/R/results/Results_GSE30337.RData") #NMR 

#dataframes
GSE30337_brain_cons <- data.frame(cons_val = 
                                    c(-log(Results_GSE30337[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]),
                                      -log(Results_GSE30337[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]])),
                                  type = 
                                    c(rep("dec", length(Results_GSE30337[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]])), 
                                      rep("inc", length(Results_GSE30337[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]))),
                                  dataset = "H. glaber brain")

GSE30337_kidney_cons <- data.frame(cons_val = 
                                     c(-log(Results_GSE30337[["kidney"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]),
                                       -log(Results_GSE30337[["kidney"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]])),
                                   type = 
                                     c(rep("dec", length(Results_GSE30337[["kidney"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]])), 
                                       rep("inc", length(Results_GSE30337[["kidney"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]))),
                                   dataset = "H. glaber kidney")

GSE30337_liver_cons <- data.frame(cons_val = 
                                    c(-log(Results_GSE30337[["liver"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]),
                                      -log(Results_GSE30337[["liver"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]])),
                                  type = 
                                    c(rep("dec", length(Results_GSE30337[["liver"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]])), 
                                      rep("inc", length(Results_GSE30337[["liver"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]))),
                                  dataset = "H. glaber liver")
NMR_cons <- rbind(GSE30337_brain_cons, GSE30337_kidney_cons, GSE30337_liver_cons)

NMR_t_p <- data.frame(dataset = c("H. glaber brain", "H. glaber kidney", "H. glaber liver"), 
                      p = c(signif(t.test(-log(Results_GSE30337[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]),
                                                             -log(Results_GSE30337[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]))$p.val, 3),
                                signif(t.test(-log(Results_GSE30337[["kidney"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]),
                                                             -log(Results_GSE30337[["kidney"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]))$p.val, 3),
                                signif(t.test(-log(Results_GSE30337[["liver"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]),
                                                             -log(Results_GSE30337[["liver"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]))$p.val, 3)))
NMR_t_p$label <- paste0("p == ", sub("e", "%*%10^", as.character(NMR_t_p$p)))

NMR_cons_plot <- ggplot(NMR_cons, aes(cons_val, fill = type)) + 
  geom_density(alpha = 0.7)+
  scale_x_continuous(limits = c(-2, 7))+
  ylab("Density")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white', size = 1))+
  scale_fill_manual(values=c("#888888", "#6699CC"), labels = c("Young-biased","Old-biased"))+
  xlab("Conservation Score")+
  facet_wrap(~dataset, ncol = 1)+
  geom_text(data    = NMR_t_p,
            mapping = aes(x = Inf, y = Inf, label = label),
            vjust = 2,
            hjust = 1.1,
            inherit.aes = FALSE,
            size = 4.4,
            parse =T)

ggsave("suppFig4-NMR-cons-density.png", plot = NMR_cons_plot, device = png,
       path = "results_graphs/",
       scale = 1, width = 10, height = 12, units = "in", dpi = 300)
ggsave("suppFig4-NMR-cons-density.tiff", plot = NMR_cons_plot, device = "tiff",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 10, height = 12, units = "in", dpi = 300)
ggsave("suppFig4-NMR-cons-density.pdf", plot = NMR_cons_plot, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 12, units = "in", dpi = 300)
