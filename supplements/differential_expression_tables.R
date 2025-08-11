### Differential Expression gene class tables.

library(tidyverse)

load("drosophila_melanogaster/R/results/Results_pacifico.RData")

pacifico_inc_genes <- Results_pacifico[["head"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
pacifico_inc_genes$class <- "Old-Biased"

pacifico_dec_genes <- Results_pacifico[["head"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
pacifico_dec_genes$class <- "Young-Biased"

pacifico_nonsig <- Results_pacifico[["head"]][["deg_list"]][["Non_sig"]][["nonsig_genes_cor"]] %>%
  as.data.frame()
pacifico_nonsig$class <- "Not-significant"

pacifico_gene_class <- bind_rows(pacifico_inc_genes, pacifico_dec_genes, pacifico_nonsig)

write.csv(pacifico_gene_class, "supplements/differential_expression_tables/pacifico_gene_class.csv",
          quote = F)

load("gallus_gallus/R/results/Results_GSE114129.RData") # chicken

chicken_inc_genes <-Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
chicken_inc_genes$class <- "Old-Biased"

chicken_dec_genes <-Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
chicken_dec_genes$class <- "Young-Biased"

chicken_nonsig <-Results_GSE114129[["brain"]][["deg_list"]][["Non_sig"]][["nonsig_genes_cor"]] %>%
  as.data.frame()
chicken_nonsig$class <- "Not-significant"

chicken_gene_class <- bind_rows(chicken_inc_genes, chicken_dec_genes, chicken_nonsig)

write.csv(chicken_gene_class, "supplements/differential_expression_tables/GSE114129_gene_class.csv",
          quote = F)


load("n_furzeri/GSE66712/R/results/Results_GSE66712.RData") # killifish

#liver
killifish_liver_inc_genes <-Results_GSE66712[["liver"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
killifish_liver_inc_genes$class <- "Old-Biased"

killifish_liver_dec_genes <-Results_GSE66712[["liver"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
killifish_liver_dec_genes$class <- "Young-Biased"

killifish_liver_nonsig <-Results_GSE66712[["liver"]][["deg_list"]][["Non_sig"]][["nonsig_genes_cor"]] %>%
  as.data.frame()
killifish_liver_nonsig$class <- "Not-significant"

killifish_liver_gene_class <- bind_rows(killifish_liver_inc_genes, killifish_liver_dec_genes, killifish_liver_nonsig)

write.csv(killifish_liver_gene_class, "supplements/differential_expression_tables/GSE66712_liver_gene_class.csv",
          quote = F)

#skin
killifish_skin_inc_genes <-Results_GSE66712[["skin"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
killifish_skin_inc_genes$class <- "Old-Biased"

killifish_skin_dec_genes <-Results_GSE66712[["skin"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
killifish_skin_dec_genes$class <- "Young-Biased"

killifish_skin_nonsig <-Results_GSE66712[["skin"]][["deg_list"]][["Non_sig"]][["nonsig_genes_cor"]] %>%
  as.data.frame()
killifish_skin_nonsig$class <- "Not-significant"

killifish_skin_gene_class <- bind_rows(killifish_skin_inc_genes, killifish_skin_dec_genes, killifish_skin_nonsig)

write.csv(killifish_skin_gene_class, "supplements/differential_expression_tables/GSE66712_skin_gene_class.csv",
          quote = F)



load("mus_musculus/GSE99791/R/results/Results_GSE99791.RData") #mouse astrocyte

#cerebellum
mouse_cerebellum_inc_genes <-Results_GSE99791[["cerebellum"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
mouse_cerebellum_inc_genes$class <- "Old-Biased"

mouse_cerebellum_dec_genes <-Results_GSE99791[["cerebellum"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
mouse_cerebellum_dec_genes$class <- "Young-Biased"

mouse_cerebellum_nonsig <-Results_GSE99791[["cerebellum"]][["deg_list"]][["Non_sig"]][["nonsig_genes_cor"]] %>%
  as.data.frame()
mouse_cerebellum_nonsig$class <- "Not-significant"

mouse_cerebellum_gene_class <- bind_rows(mouse_cerebellum_inc_genes, mouse_cerebellum_dec_genes, mouse_cerebellum_nonsig)

write.csv(mouse_cerebellum_gene_class, "supplements/differential_expression_tables/GSE99791_cerebellum_gene_class.csv",
          quote = F)

#hypothalamus
mouse_hypothalamus_inc_genes <-Results_GSE99791[["hypothalamus"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
mouse_hypothalamus_inc_genes$class <- "Old-Biased"

mouse_hypothalamus_dec_genes <-Results_GSE99791[["hypothalamus"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
mouse_hypothalamus_dec_genes$class <- "Young-Biased"

mouse_hypothalamus_nonsig <-Results_GSE99791[["hypothalamus"]][["deg_list"]][["Non_sig"]][["nonsig_genes_cor"]] %>%
  as.data.frame()
mouse_hypothalamus_nonsig$class <- "Not-significant"

mouse_hypothalamus_gene_class <- bind_rows(mouse_hypothalamus_inc_genes, mouse_hypothalamus_dec_genes, mouse_hypothalamus_nonsig)

write.csv(mouse_hypothalamus_gene_class, "supplements/differential_expression_tables/GSE99791_hypothalamus_gene_class.csv",
          quote = F)


load("naked_mole_rat/GSE30337/R/results/Results_GSE30337.RData") #NMR 

#brain

nmr_brain_inc_genes <-Results_GSE30337[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
nmr_brain_inc_genes$class <- "Old-Biased"

nmr_brain_dec_genes <-Results_GSE30337[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
nmr_brain_dec_genes$class <- "Young-Biased"

nmr_brain_gene_class <- bind_rows(nmr_brain_inc_genes, nmr_brain_dec_genes)

write.csv(nmr_brain_gene_class, "supplements/differential_expression_tables/GSE30337_brain_gene_class.csv",
          quote = F)

#kidney

nmr_kidney_inc_genes <-Results_GSE30337[["kidney"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
nmr_kidney_inc_genes$class <- "Old-Biased"

nmr_kidney_dec_genes <-Results_GSE30337[["kidney"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
nmr_kidney_dec_genes$class <- "Young-Biased"

nmr_kidney_gene_class <- bind_rows(nmr_kidney_inc_genes, nmr_kidney_dec_genes)

write.csv(nmr_kidney_gene_class, "supplements/differential_expression_tables/GSE30337_kidney_gene_class.csv",
          quote = F)


#liver

nmr_liver_inc_genes <-Results_GSE30337[["liver"]][["deg_list"]][["inc_sig"]][["inc_genes_cor"]] %>%
  as.data.frame()
nmr_liver_inc_genes$class <- "Old-Biased"

nmr_liver_dec_genes <-Results_GSE30337[["liver"]][["deg_list"]][["dec_sig"]][["dec_genes_cor"]] %>% 
  as.data.frame()
nmr_liver_dec_genes$class <- "Young-Biased"

nmr_liver_gene_class <- bind_rows(nmr_liver_inc_genes, nmr_liver_dec_genes)

write.csv(nmr_liver_gene_class, "supplements/differential_expression_tables/GSE30337_liver_gene_class.csv",
          quote = F)

