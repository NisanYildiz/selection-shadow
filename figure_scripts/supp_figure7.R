#Loading required packages
library(tidyr)
library(patchwork)

all_table <- data.frame(transcript_id = "", mean_tajimasD = 0, gene_id = "", 
                        gene_class = "", tissue = "", pop = "")

for (file in list.files("/mnt/NEOGENE4/projects/melih_2020/mus_musculus/", pattern = "mean_tajimasD.*csv")){
  table <- read.csv(paste0("/mnt/NEOGENE4/projects/melih_2020/mus_musculus/",file), row.names = 1)
  tissue_and_pop <- gsub("mean_tajimasD_scores_per_transcript_GSE99791_", "", file)
  tissue_and_pop <- gsub(".csv", "", tissue_and_pop)
  tissue <- strsplit(tissue_and_pop, "_")[[1]][1]
  pop <- gsub(paste0(tissue,"_"), "", tissue_and_pop)
  table$tissue <- tissue
  table$pop <- pop
  all_table <- rbind(all_table, table)
}

all_table <- all_table[-1,]
hypothalamus_table <- all_table[all_table$tissue == "hypothalamus",]
cerebellum_table <- all_table[all_table$tissue == "cerebellum",]

supp_fig_7_a <- ggplot(data = cerebellum_table[cerebellum_table$gene_class != "non-biased",], aes(y = mean_tajimasD, x = gene_class)) +
  geom_boxplot(aes(fill = gene_class)) + 
  theme(legend.title = element_blank(), 
        legend.text=element_text(size=11),
        axis.text.x = element_text(size = 0), 
        axis.title.y = element_text(size=15), 
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size=15, face = "italic"),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~pop, ncol = 4)+
  theme_bw()+
  ggtitle("Cerebellum")+
  ylab("Gene-wise Mean Tajima's D")+
  xlab("Gene Class")+
  scale_fill_manual(values=c("#6699CC", "#888888"), guide = "none")+
  # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",
                     label.x = 1.4, label.y = -2)  


supp_fig_7_b <- ggplot(data = hypothalamus_table[hypothalamus_table$gene_class != "non-biased",], aes(y = mean_tajimasD, x = gene_class)) +
  geom_boxplot(aes(fill = gene_class)) + 
  theme(legend.title = element_blank(), 
        legend.text=element_text(size=11),
        axis.text.x = element_text(size = 0), 
        axis.title.y = element_text(size=15), 
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size=15, face = "italic"),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~pop, ncol = 4)+
  theme_bw()+
  ggtitle("Hypothalamus")+
  ylab("Gene-wise Mean Tajima's D")+
  xlab("Gene Class")+
  scale_fill_manual(values=c("#6699CC", "#888888"), guide = "none")+
  # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",
                     label.x = 1.4, label.y = -2)  

supp_fig_7 <- supp_fig_7_a / supp_fig_7_b + plot_annotation(tag_levels = 'A') 
supp_fig_7 
ggsave("suppFig7-tajimasD.png", plot = supp_fig_7, device = png,
       path = "results_graphs/",
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)
ggsave("suppFig7-tajimasD.tiff", plot = supp_fig_7, device = "tiff",
       path = "results_graphs/", type = "cairo",
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)
ggsave("suppFig7-tajimasD.pdf", plot = supp_fig_7, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)
