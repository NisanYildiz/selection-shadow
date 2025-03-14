# Code used to produce Supp. Figure 3

#Loading the required packages
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

#Loading up the Results 
load("gallus_gallus/R/results/Results_GSE114129.RData") # chicken
load("drosophila_melanogaster/R/results/Results_pacifico.RData") # drosophila
load("n_furzeri/GSE66712/R/results/Results_GSE66712.RData") # killifish
load("mus_musculus/GSE99791/R/results/Results_GSE99791.RData") #mouse astrocyte
load("naked_mole_rat/GSE30337/R/results/Results_GSE30337.RData") #NMR 

Results <- list(Results_GSE114129 = Results_GSE114129, Results_GSE30337 = Results_GSE30337,
                Results_GSE66712 = Results_GSE66712, Results_GSE99791 = Results_GSE99791,
                Results_pacifico = Results_pacifico)

inc_dnds <- lapply(Results, function(x){
  
  inc_dnds <- sapply(x, function(tissue){
    tissue[["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]
  })
})


dec_dnds <- lapply(Results, function(x){
  dec_dnds <-sapply(x, function(tissue){
    tissue[["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]
  })
})

inc_dec_df <- rbind(
  
  data.frame(dnds = c(inc_dnds$Results_GSE114129, dec_dnds$Results_GSE114129),
             class = c(rep("old-biased", length(inc_dnds$Results_GSE114129)),
                       rep("young-biased", length(dec_dnds$Results_GSE114129))),
             organism = "G. gallus brain",
             tissue = "brain"),
  
  data.frame(dnds = c(inc_dnds$Results_GSE66712$liver, dec_dnds$Results_GSE66712$liver),
             class = c(rep("old-biased", length(inc_dnds$Results_GSE66712$liver)),
                       rep("young-biased", length(dec_dnds$Results_GSE66712$liver))),
             organism = "N. furzeri liver",
             tissue = "liver"),
  
  data.frame(dnds = c(inc_dnds$Results_GSE66712$skin, dec_dnds$Results_GSE66712$skin),
             class = c(rep("old-biased", length(inc_dnds$Results_GSE66712$skin)),
                       rep("young-biased", length(dec_dnds$Results_GSE66712$skin))),
             organism = "N. furzeri skin",
             tissue = "skin"),
  
  data.frame(dnds = c(inc_dnds$Results_pacifico, dec_dnds$Results_pacifico),
             class = c(rep("old-biased", length(inc_dnds$Results_pacifico)),
                       rep("young-biased", length(dec_dnds$Results_pacifico))),
             organism = "D. melanogaster brain",
             tissue = "brain"), 
  
  data.frame(dnds = c(inc_dnds$Results_GSE99791$cerebellum, dec_dnds$Results_GSE99791$cerebellum),
             class = c(rep("old-biased", length(inc_dnds$Results_GSE99791$cerebellum)),
                       rep("young-biased", length(dec_dnds$Results_GSE99791$cerebellum))),
             organism = "M. musculus astrocyte (cerebellum)",
             tissue = "astrocyte"), 
  
  data.frame(dnds = c(inc_dnds$Results_GSE99791$hypothalamus, dec_dnds$Results_GSE99791$hypothalamus),
             class = c(rep("old-biased", length(inc_dnds$Results_GSE99791$hypothalamus)),
                       rep("young-biased", length(dec_dnds$Results_GSE99791$hypothalamus))),
             organism = "M. musculus astrocyte (hypothalamus)",
             tissue = "astrocyte")
  
)


supp_figure3 <- ggplot(data= inc_dec_df, 
                       aes(x = class, y = dnds)) +
  geom_boxplot(aes(fill = class), shape = 21, show.legend = T)+
  scale_y_continuous(name="dN/dS")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  theme_bw()+
  xlab("")+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_fill_manual(values=c("#6699CC", "#888888"))+
  theme(legend.title = element_blank(), 
        legend.text=element_text(size=11),
        axis.text.x = element_text(size = 0), 
        axis.title.y = element_text(size=15), 
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size=11, face = "italic"),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~ organism) +
  guides(fill = guide_legend(override.aes = list(size = 7, shape = 22)))+ 
  stat_compare_means(method = "wilcox.test",
                      label.x = 1.3, 
                      label.y = -0.05)

ggsave("suppFig3-dnds-difference-bulk.png", plot = supp_figure3, device = png,
       path = "results_graphs/",
       scale = 1, width = 10, height = 5, units = "in", dpi = 300)
ggsave("suppFig3-dnds-difference-bulk.tiff", plot = supp_figure3, device = "tiff",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 10, height = 5, units = "in", dpi = 300)
ggsave("suppFig3-dnds-difference-bulk.pdf", plot = supp_figure3, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 5, units = "in", dpi = 300)