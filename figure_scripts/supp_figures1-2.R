## PCA supplement
library(ggplot2)
library(patchwork)
library(viridis)
library(wesanderson)
library(ggsci)
library(ggthemes)
library(scico)
library(gridExtra)

load("drosophila_melanogaster/R/results/pca_dros.RData")
load("gallus_gallus/R/results/pca_gallus.RData")
load("n_furzeri/GSE66712/R/results/pca_killifish.RData")
load("mus_musculus/GSE99791/R/results/pca_mus.RData")

gallus <- pca_gallus_1_2 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (days)") +
  ggtitle(expression(paste("A  ", italic("Gallus gallus "), "brain" ))) +
  pca_gallus_3_4 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (days)")  +
  ggtitle("") + plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_gallus_brain

dros <- pca_dros_1_2 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (days)") + 
  ggtitle(expression(paste("B  ", italic("Drosophila melanogaster "), "brain" )))+
  pca_dros_3_4 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (days)")+
  ggtitle("")+  plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_dros_brain

killifish_liver <- pca_killifish_liver_1_2 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (weeks)") + 
  ggtitle(expression(paste("C  ", italic("Nothobranchius furzeri "), "liver" )))+
  pca_killifish_liver_3_4 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (weeks)") +
  ggtitle("")+  plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_furzeri_liver

killifish_skin <- pca_killifish_skin_1_2 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (weeks)") + 
  ggtitle(expression(paste("D  ",italic("Nothobranchius furzeri "), "skin" )))+
  pca_killifish_skin_3_4 + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (weeks)") +
  ggtitle("")+  plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_furzeri_skin

mouse_cereb <- pca_mus_1_2$cerebellum + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") + 
  ggtitle(expression(paste("A  ",italic("Mus musculus "), "astrocyte from cerebellum" )))+
  pca_mus_3_4$cerebellum + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") +
  ggtitle("") +  plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_musculus_cerebellum

mouse_hypo <- pca_mus_1_2$hypothalamus + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") + 
     ggtitle(expression(paste("B  ",italic("Mus musculus "), "astrocyte from hypothalamus" )))+
     pca_mus_3_4$hypothalamus + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") +
     ggtitle("")  +  plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_musculus_hypothalamus

mouse_motor <- pca_mus_1_2$`motor cortex` + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") + 
      ggtitle(expression(paste("C  ", italic("Mus musculus "), "astrocyte from motor cortex" )))+
      pca_mus_3_4$`motor cortex` + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") +
      ggtitle("") +  plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_musculus_motorcortex

mouse_visual <- pca_mus_1_2$`visual cortex` + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") + 
        ggtitle(expression(paste("D  ",italic("Mus musculus "), "astrocyte from visual cortex" )))+
        pca_mus_3_4$`visual cortex` + scale_color_scico_d(palette = "acton", direction = -1, name = "Age (months)") +
        ggtitle("")  +  plot_layout(guides = "collect")
#prop 5x10, 1000x500, name = supp_pca_musculus_visualcortex

supp_fig_1_pca <- gallus / dros / killifish_liver / killifish_skin + 
  theme(plot.tag = element_text(size = 16))

ggsave("suppFig1-bulkPCA.png", plot = supp_fig_1_pca, device = "png",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 10, height = 12, units = "in", dpi = 300)
ggsave("suppFig1-bulkPCA.tiff", plot = supp_fig_1_pca, device = "tiff",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 10, height = 12, units = "in", dpi = 300)
ggsave("suppFig1-bulkPCA.pdf", plot = supp_fig_1_pca, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 12, units = "in", dpi = 300)

supp_fig_2_pca <- mouse_cereb / mouse_hypo / mouse_motor / mouse_visual + 
  theme(plot.tag = element_text(size = 16))

ggsave("suppFig2-mousePCA.png", plot = supp_fig_2_pca, device = "png",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 11.3, height = 12, units = "in", dpi = 300)
ggsave("suppFig2-mousePCA.tiff", plot = supp_fig_2_pca, device = "tiff",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 11.3, height = 12, units = "in", dpi = 300)
ggsave("suppFig2-mousePCA.pdf", plot = supp_fig_2_pca, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 11.3, height = 12, units = "in", dpi = 300)
