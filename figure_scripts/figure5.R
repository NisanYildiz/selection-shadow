library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(tidytext)
library(ggpubr)

#source('./turan/scripts/00_Source.R')

#Loading up the Results 
load("mus_musculus/tabula_muris_senis/R/results/Results_tabula_muris_senis.RData")#Tabula Muris Senis

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
corr_tabula <- split(corr_tabula_df, corr_tabula_df$Tissue)
head(corr_tabula)

corr_tabulax=rbind()
for (i in 1:6){
  aa= corr_tabula[[i]]
  aa2= data.frame(rbindlist(list(aa)))
  corr_tabulax= rbind(corr_tabulax,aa2)}
class(corr_tabulax)
head(corr_tabulax)
corr_tabulax[corr_tabulax$Tissue%in%'muscle',]

corr_tabulax$p_adjust= ifelse(corr_tabulax$p_adjust < 0.1, "*", "")
corr_tabulaxa= corr_tabulax[,c('corr_sig','cell_type', 'Tissue','p_adjust')]
head(corr_tabulaxa)

corr_tabulax2= reshape2::melt(corr_tabulaxa)
head(corr_tabulax2)
corr_tabulax2$Tissue= factor(corr_tabulax2$Tissue, levels= c('lung','brain','liver','kidney', 'muscle', 'skin'))

# the item we want to reorder
# what we want to reorder by
# the groups or categories we want to reorder within
tissuesx <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#9e1b42", "#027cd9", "#afb370"), c('brain','liver', 'skin','kidney', 'lung', 'muscle'))
figure5= ggplot(corr_tabulax2, aes(x= reorder_within(cell_type, value, Tissue), y=value, fill=Tissue))+
  geom_bar(position = 'dodge',stat = 'identity')+
  geom_text(aes(label=p_adjust), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(~Tissue, drop=T, scales = "free", space = "free")+
  scale_fill_manual(name='', values = tissuesx)+
  scale_color_manual(name='', values = tissuesx)+
  theme_pubr(base_size = 12) +
  scale_x_reordered() +
  ylab(expression(rho))+
  xlab("Cell type")+
  guides(fill=F)+
  guides(color=F)+
  theme(axis.title=element_text(size=14))+
  theme(axis.text.x = element_text(vjust= 0.5, hjust=1, angle=90))+
  theme(strip.text= element_text(size=14),
        strip.background = element_rect(fill='white', color = "black", size = 1.3)) 

ggsave("figure5-tabula-ADICT.pdf", plot = figure5, device = cairo_pdf,
       path = "results_graphs/",
       units="cm", width = 40, height = 20)
ggsave("figure5-tabula-ADICT.png", plot = figure5, device = png,
       path = "results_graphs/", type = "cairo",
       units="cm", width = 40, height = 20)
ggsave("figure5-tabula-ADICT.tiff", plot = figure5, device = "tiff",
       path = "results_graphs/", type = "cairo",
       units="cm", width = 40, height = 20)

write.table(data.frame(cell_type= corr_tabula_df$cell_type, 
                     tissue = corr_tabula_df$Tissue,
                     p = corr_tabula_df$p_all, 
                     p_adj = corr_tabula_df$p_adjust,
                     rho = corr_tabula_df$corr_sig),
          file = "supplements/tabula_muris_rho_vals.tsv",
          quote = F, row.names = F, sep = "\t"
          )
