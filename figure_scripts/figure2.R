# Code used to produce Figure 2 

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

#loading the confidence interval function
source("functions.R")

#Setting up the data

##conservation of increasing genes
inc_cons <- lapply(Results, function(x){
  
  inc_cons <- sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]]))
    
  })
})

##conservation of decreasing genes
dec_cons <- lapply(Results, function(x){
  
  dec_cons <-sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]]))
    
  })
})

##conservation of increasing genes relative to no change
inc_cons_relative <- lapply(Results, function(x){
  
  inc_cons <- sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]])) - 
      mean(c(-log(tissue[["deg_list"]][["Non_sig"]][["nonsig_genes_dnds"]][["dNdS"]])))
    
  })
})

##conservation of decreasing genes relative to no change
dec_cons_relative <- lapply(Results, function(x){
  
  dec_cons <-sapply(x, function(tissue){
    c(-log(tissue[["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]])) - 
      mean(c(-log(tissue[["deg_list"]][["Non_sig"]][["nonsig_genes_dnds"]][["dNdS"]])))
    
  })
})
##creating the data frames for Figure 2B

inc_cons_relative_means <- c(chicken_brain = mean(inc_cons_relative[["Results_GSE114129"]]), killifish_liver = mean(inc_cons_relative[["Results_GSE66712"]][["liver"]]), killifish_skin = mean(inc_cons_relative[["Results_GSE66712"]][["skin"]]), drosophila_brain = mean(inc_cons_relative[["Results_pacifico"]]))
names_inc_cons_relative_means <- names(inc_cons_relative_means)
inc_cons_relative_means <- as.data.frame(inc_cons_relative_means)
inc_cons_relative_means$names <-  names_inc_cons_relative_means
inc_cons_relative_means$upper <- c(chicken_brain = conf_int(inc_cons_relative[["Results_GSE114129"]])[1], killifish_liver = conf_int(inc_cons_relative[["Results_GSE66712"]][["liver"]])[1], killifish_skin = conf_int(inc_cons_relative[["Results_GSE66712"]][["skin"]])[1], drosophila_brain = conf_int(inc_cons_relative[["Results_pacifico"]])[1])
inc_cons_relative_means$lower <- c(chicken_brain = conf_int(inc_cons_relative[["Results_GSE114129"]])[2], killifish_liver = conf_int(inc_cons_relative[["Results_GSE66712"]][["liver"]])[2], killifish_skin = conf_int(inc_cons_relative[["Results_GSE66712"]][["skin"]])[2], drosophila_brain = conf_int(inc_cons_relative[["Results_pacifico"]])[2])
inc_cons_relative_means$species <- c("Chicken","Turquoise Killifish", "Turquoise Killifish", "Fruit Fly")
inc_cons_relative_means$Tissue <- c("Brain", "Liver", "Skin", "Brain")
inc_cons_relative_means$tag <- c("G. gallus", " N. furzeri ", 
                                 "N. furzeri", "D. melanogaster")
inc_cons_relative_means$title <- c("Old-biased genes")
colnames(inc_cons_relative_means)[1] <- "cons_means"

dec_cons_relative_means <- c(chicken_brain = mean(dec_cons_relative[["Results_GSE114129"]]), killifish_liver = mean(dec_cons_relative[["Results_GSE66712"]][["liver"]]), killifish_skin = mean(dec_cons_relative[["Results_GSE66712"]][["skin"]]), drosophila_brain = mean(dec_cons_relative[["Results_pacifico"]]))
names_dec_cons_relative_means <- names(dec_cons_relative_means)
dec_cons_relative_means <- as.data.frame(dec_cons_relative_means)
dec_cons_relative_means$names <-  names_dec_cons_relative_means
dec_cons_relative_means$upper <- c(chicken_brain = conf_int(dec_cons_relative[["Results_GSE114129"]])[1], killifish_liver = conf_int(dec_cons_relative[["Results_GSE66712"]][["liver"]])[1], killifish_skin = conf_int(dec_cons_relative[["Results_GSE66712"]][["skin"]])[1], drosophila_brain = conf_int(dec_cons_relative[["Results_pacifico"]])[1])
dec_cons_relative_means$lower <- c(chicken_brain = conf_int(dec_cons_relative[["Results_GSE114129"]])[2], killifish_liver = conf_int(dec_cons_relative[["Results_GSE66712"]][["liver"]])[2], killifish_skin = conf_int(dec_cons_relative[["Results_GSE66712"]][["skin"]])[2], drosophila_brain = conf_int(dec_cons_relative[["Results_pacifico"]])[2])
dec_cons_relative_means$species <- c("Chicken","Turquoise Killifish", "Turquoise Killifish", "Fruit Fly")
dec_cons_relative_means$Tissue <- c("Brain", "Liver", "Skin", "Brain")
dec_cons_relative_means$tag <- c("G. gallus", " N. furzeri ", 
                                 "N. furzeri", "D. melanogaster")
dec_cons_relative_means$title <- c("Young-biased genes")
colnames(dec_cons_relative_means)[1] <- "cons_means"

cons_means <- rbind(inc_cons_relative_means, dec_cons_relative_means)
cons_means$title <- factor(cons_means$title, levels= c("Old-biased genes", "Young-biased genes"))

##creating the dataframes for figure 2A

GSE66712_liver_cons <- data.frame(cons_val = c(dec_cons[["Results_GSE66712"]][["liver"]],
                                               inc_cons[["Results_GSE66712"]][["liver"]]), 
                                  type = c(rep("dec", length(dec_cons[["Results_GSE66712"]][["liver"]])), 
                                           rep("inc", length(inc_cons[["Results_GSE66712"]][["liver"]]))),
                                  dataset = "N. furzeri liver",
                                  t_p = sig_symbols(t.test(dec_cons[["Results_GSE66712"]][["liver"]],
                                                                     inc_cons[["Results_GSE66712"]][["liver"]])$p.val))

GSE66712_skin_cons <- data.frame(cons_val = c(dec_cons[["Results_GSE66712"]][["skin"]],
                                              inc_cons[["Results_GSE66712"]][["skin"]]), 
                                 type = c(rep("dec", length(dec_cons[["Results_GSE66712"]][["skin"]])), 
                                          rep("inc", length(inc_cons[["Results_GSE66712"]][["skin"]]))),
                                 dataset = "N. furzeri skin",
                                 t_p = sig_symbols(t.test(dec_cons[["Results_GSE66712"]][["skin"]],
                                                                    inc_cons[["Results_GSE66712"]][["skin"]])$p.val))


GSE114129_brain_cons <- data.frame(cons_val = c(dec_cons[["Results_GSE114129"]],
                                                inc_cons[["Results_GSE114129"]]), 
                                   type = c(rep("dec", length(dec_cons[["Results_GSE114129"]])), 
                                            rep("inc", length(inc_cons[["Results_GSE114129"]]))),
                                   dataset = "G. gallus brain",
                                   t_p = sig_symbols(t.test(dec_cons[["Results_GSE114129"]],
                                                                      inc_cons[["Results_GSE114129"]])$p.val))

Pacifico_brain_cons <- data.frame(cons_val = c(dec_cons[["Results_pacifico"]],
                                               inc_cons[["Results_pacifico"]]), 
                                  type = c(rep("dec", length(dec_cons[["Results_pacifico"]])), 
                                           rep("inc", length(inc_cons[["Results_pacifico"]]))),
                                  dataset = "D. melanogaster brain",
                                  t_p = sig_symbols(t.test(dec_cons[["Results_pacifico"]],
                                                                     inc_cons[["Results_pacifico"]])$p.val))

cons <- rbind(GSE114129_brain_cons, Pacifico_brain_cons, GSE66712_liver_cons, GSE66712_skin_cons)

cons_t_p <- data.frame(dataset = c("G. gallus brain", "D. melanogaster brain", "N. furzeri liver", "N. furzeri skin"),
                       p = c(signif(t.test(dec_cons[["Results_GSE114129"]],
                                           inc_cons[["Results_GSE114129"]])$p.val, 3),
                             signif(t.test(dec_cons[["Results_pacifico"]],
                                           inc_cons[["Results_pacifico"]])$p.val,3),
                             signif(t.test(dec_cons[["Results_GSE66712"]][["liver"]],
                                           inc_cons[["Results_GSE66712"]][["liver"]])$p.val, 3),
                             signif(t.test(dec_cons[["Results_GSE66712"]][["skin"]],
                                           inc_cons[["Results_GSE66712"]][["skin"]])$p.val, 3)))
cons_t_p$t_p <- paste0("p == ", sub("e", "%*%10^", as.character(cons_t_p$p)))
cons_t_p$p_adj <- p.adjust(cons_t_p$p, "BH")
#Plotting

## Figure 2A
cons_plot <- ggplot(cons, aes(cons_val, fill = type)) + 
  geom_density(alpha = 0.7, aes(fill = type))+
  scale_x_continuous(limits = c(-1, 7))+
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
  facet_wrap(~dataset)+
  geom_text(data    = cons_t_p,
            mapping = aes(x = Inf, y = Inf, label = t_p),
            vjust = 1.3,
            hjust = 1.1,
            inherit.aes = FALSE,
            size = 4.4,
            parse = T)

## Figure 2B
inc_dec_plot <- ggplot(data= cons_means, 
                       aes(x = factor(tag, levels = c("G. gallus",
                                                      "D. melanogaster", 
                                                      " N. furzeri ", 
                                                      "N. furzeri")), cons_means)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.14, show.legend = F) +
  geom_point(size = 8, aes(fill = Tissue), shape = 21, show.legend = T)+
  scale_y_continuous(name="Conservation relative to no change", limits=c(-0.44, 0.44))+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  xlab("")+
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))+
  theme(legend.title = element_blank(), 
        legend.text=element_text(size=11),
        axis.text.x = element_text(size = 13, face = "italic"), 
        axis.title.y = element_text(size=15), 
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white', size = 1))+
  facet_wrap(~ title)+
  guides(fill = guide_legend(override.aes = list(size = 7, shape = 22)))

## Combined Figure 2
figure2 <- cons_plot / inc_dec_plot + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(3, 2)) & 
  theme(plot.tag = element_text(size = 16))

ggsave("figure2-cons_point-density-4samples.png", plot = figure2, device = png,
       path = "results_graphs/",
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)
ggsave("figure2-cons_point-density-4samples.tiff", plot = figure2, device = "tiff",
       path = "results_graphs/",
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)
ggsave("figure2-cons_point-density-4samples.pdf", plot = figure2, device = "pdf",
       path = "results_graphs/",
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)
