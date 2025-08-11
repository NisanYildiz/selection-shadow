#gene age, dN/dS analysis

musculus_gene_age <- read.csv("/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/Mus_musculus.csv")
musculus_gene_age$gene_age <- as.numeric(gsub(">", "", musculus_gene_age$gene_age))


load("mus_musculus/GSE99791/R/results/Results_GSE99791.RData") #mouse astrocyte

mus_cerebellum_common_inc_dnds <- Results_GSE99791[["cerebellum"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]][Results_GSE99791[["cerebellum"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_cerebellum_common_inc_gene_names <- Results_GSE99791[["cerebellum"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]][Results_GSE99791[["cerebellum"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_cerebellum_common_inc_gene_age <- musculus_gene_age[musculus_gene_age$ensembl_gene_id %in% mus_cerebellum_common_inc_gene_names,]
mus_cerebellum_common_inc_gene_age <- mus_cerebellum_common_inc_gene_age[order(mus_cerebellum_common_inc_gene_age$ensembl_gene_id),]
mus_cerebellum_common_inc_max_exp <- apply(X= Results_GSE99791[["cerebellum"]][["deg_list"]][["inc_sig"]][["inc_genes_exp"]][Results_GSE99791[["cerebellum"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id,],
                            MARGIN = 1, function(x){max(x)})

identical(as.character(mus_cerebellum_common_inc_gene_age$ensembl_gene_id), mus_cerebellum_common_inc_gene_names)

mus_cerebellum_common_dec_dnds <- Results_GSE99791[["cerebellum"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]][Results_GSE99791[["cerebellum"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_cerebellum_common_dec_gene_names <- Results_GSE99791[["cerebellum"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]][Results_GSE99791[["cerebellum"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_cerebellum_common_dec_gene_age <- musculus_gene_age[musculus_gene_age$ensembl_gene_id %in% mus_cerebellum_common_dec_gene_names,]
mus_cerebellum_common_dec_gene_age <- mus_cerebellum_common_dec_gene_age[order(mus_cerebellum_common_dec_gene_age$ensembl_gene_id),]
mus_cerebellum_common_dec_max_exp <- apply(X= Results_GSE99791[["cerebellum"]][["deg_list"]][["dec_sig"]][["dec_genes_exp"]][Results_GSE99791[["cerebellum"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id,],
                            MARGIN = 1, function(x){max(x)})

identical(as.character(mus_cerebellum_common_dec_gene_age$ensembl_gene_id), mus_cerebellum_common_dec_gene_names)


mus_astrocyte_cerebellum_dnds_gene_age <- rbind(
  
  data.frame(gene_age = mus_cerebellum_common_inc_gene_age$gene_age ,
             gene_class = "old-biased",
             dnds = mus_cerebellum_common_inc_dnds,
             max_exp = mus_cerebellum_common_inc_max_exp),
  
  data.frame(gene_age = mus_cerebellum_common_dec_gene_age$gene_age ,
             gene_class = "young-biased",
             dnds = mus_cerebellum_common_dec_dnds,
             max_exp = mus_cerebellum_common_dec_max_exp)
)

summary(lm(dnds ~ max_exp + gene_age + gene_class, data= mus_astrocyte_cerebellum_dnds_gene_age))


mus_hypothalamus_common_inc_dnds <- Results_GSE99791[["hypothalamus"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]][Results_GSE99791[["hypothalamus"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_hypothalamus_common_inc_gene_names <- Results_GSE99791[["hypothalamus"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]][Results_GSE99791[["hypothalamus"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_hypothalamus_common_inc_gene_age <- musculus_gene_age[musculus_gene_age$ensembl_gene_id %in% mus_hypothalamus_common_inc_gene_names,]
mus_hypothalamus_common_inc_gene_age <- mus_hypothalamus_common_inc_gene_age[order(mus_hypothalamus_common_inc_gene_age$ensembl_gene_id),]
mus_hypothalamus_common_inc_max_exp <- apply(X= Results_GSE99791[["hypothalamus"]][["deg_list"]][["inc_sig"]][["inc_genes_exp"]][Results_GSE99791[["hypothalamus"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id,],
                                           MARGIN = 1, function(x){max(x)})

identical(as.character(mus_hypothalamus_common_inc_gene_age$ensembl_gene_id), mus_hypothalamus_common_inc_gene_names)

mus_hypothalamus_common_dec_dnds <- Results_GSE99791[["hypothalamus"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]][Results_GSE99791[["hypothalamus"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_hypothalamus_common_dec_gene_names <- Results_GSE99791[["hypothalamus"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]][Results_GSE99791[["hypothalamus"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id]
mus_hypothalamus_common_dec_gene_age <- musculus_gene_age[musculus_gene_age$ensembl_gene_id %in% mus_hypothalamus_common_dec_gene_names,]
mus_hypothalamus_common_dec_gene_age <- mus_hypothalamus_common_dec_gene_age[order(mus_hypothalamus_common_dec_gene_age$ensembl_gene_id),]
mus_hypothalamus_common_dec_max_exp <- apply(X= Results_GSE99791[["hypothalamus"]][["deg_list"]][["dec_sig"]][["dec_genes_exp"]][Results_GSE99791[["hypothalamus"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id,],
                                           MARGIN = 1, function(x){max(x)})

identical(as.character(mus_hypothalamus_common_dec_gene_age$ensembl_gene_id), mus_hypothalamus_common_dec_gene_names)

mus_astrocyte_hypothalamus_dnds_gene_age <- rbind(
  
  data.frame(gene_age = mus_hypothalamus_common_inc_gene_age$gene_age ,
             gene_class = "old-biased",
             dnds = mus_hypothalamus_common_inc_dnds,
             max_exp = mus_hypothalamus_common_inc_max_exp),
  
  data.frame(gene_age = mus_hypothalamus_common_dec_gene_age$gene_age ,
             gene_class = "young-biased",
             dnds = mus_hypothalamus_common_dec_dnds,
             max_exp = mus_hypothalamus_common_dec_max_exp)
)

summary(lm(dnds ~ max_exp + gene_age + gene_class, data= mus_astrocyte_hypothalamus_dnds_gene_age))


gallus_gene_age <- read.csv("/mnt/NEOGENE4/projects/melih_2020/data/gallus_gallus/Gallus_gallus.csv")
gallus_gene_age$gene_age <- as.numeric(gsub(">", "", gallus_gene_age$gene_age))


load("gallus_gallus/R/results/Results_GSE114129.RData") # chicken

gallus_brain_common_inc_dnds <- Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]][Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id]
gallus_brain_common_inc_gene_names <- Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]][Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id]
gallus_brain_common_inc_gene_age <- gallus_gene_age[gallus_gene_age$ensembl_gene_id %in% gallus_brain_common_inc_gene_names,]
gallus_brain_common_inc_gene_age <- gallus_brain_common_inc_gene_age[order(gallus_brain_common_inc_gene_age$ensembl_gene_id),]
gallus_brain_common_inc_max_exp <- apply(X= Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_exp"]][Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id,],
                                           MARGIN = 1, function(x){max(x)})

identical(as.character(gallus_brain_common_inc_gene_age$ensembl_gene_id), gallus_brain_common_inc_gene_names)

gallus_brain_common_dec_dnds <- Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]][Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id]
gallus_brain_common_dec_gene_names <- Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]][Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id]
gallus_brain_common_dec_gene_age <- gallus_gene_age[gallus_gene_age$ensembl_gene_id %in% gallus_brain_common_dec_gene_names,]
gallus_brain_common_dec_gene_age <- gallus_brain_common_dec_gene_age[order(gallus_brain_common_dec_gene_age$ensembl_gene_id),]
gallus_brain_common_dec_max_exp <- apply(X= Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_exp"]][Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id,],
                                           MARGIN = 1, function(x){max(x)})

identical(as.character(gallus_brain_common_dec_gene_age$ensembl_gene_id), gallus_brain_common_dec_gene_names)

gallus_brain_dnds_gene_age <- rbind(
  
  data.frame(gene_age = gallus_brain_common_inc_gene_age$gene_age ,
             gene_class = "old-biased",
             dnds = gallus_brain_common_inc_dnds,
             max_exp = gallus_brain_common_inc_max_exp),
  
  data.frame(gene_age = gallus_brain_common_dec_gene_age$gene_age ,
             gene_class = "young-biased",
             dnds = gallus_brain_common_dec_dnds,
             max_exp = gallus_brain_common_dec_max_exp)
)

summary(lm(dnds ~ max_exp + gene_age + gene_class, data= gallus_brain_dnds_gene_age))


melanogaster_gene_age <- read.csv("/mnt/NEOGENE4/projects/melih_2020/data/drosophila_melanogaster/Drosophila_melanogaster.csv")
melanogaster_gene_age$gene_age <- as.numeric(gsub(">", "", melanogaster_gene_age$gene_age))
melanogaster_gene_age$ensembl_gene_id <- gsub("FBGN", "FBgn", melanogaster_gene_age$ensembl_gene_id)

load("drosophila_melanogaster/R/results/Results_pacifico.RData") # drosophila

melanogaster_brain_common_inc_dnds <- Results_pacifico[["head"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]][Results_pacifico[["head"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% melanogaster_gene_age$ensembl_gene_id]
melanogaster_brain_common_inc_gene_names <- Results_pacifico[["head"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]][Results_pacifico[["head"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% melanogaster_gene_age$ensembl_gene_id]
melanogaster_brain_common_inc_gene_age <- melanogaster_gene_age[melanogaster_gene_age$ensembl_gene_id %in% melanogaster_brain_common_inc_gene_names,]
melanogaster_brain_common_inc_gene_age <- melanogaster_brain_common_inc_gene_age[order(melanogaster_brain_common_inc_gene_age$ensembl_gene_id),]
melanogaster_brain_common_inc_max_exp <- apply(X= Results_pacifico[["head"]][["deg_list"]][["inc_sig"]][["inc_genes_exp"]][Results_pacifico[["head"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% melanogaster_gene_age$ensembl_gene_id,],
                                         MARGIN = 1, function(x){max(x)})

identical(as.character(melanogaster_brain_common_inc_gene_age$ensembl_gene_id), melanogaster_brain_common_inc_gene_names)

melanogaster_brain_common_dec_dnds <- Results_pacifico[["head"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]][Results_pacifico[["head"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% melanogaster_gene_age$ensembl_gene_id]
melanogaster_brain_common_dec_gene_names <- Results_pacifico[["head"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]][Results_pacifico[["head"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% melanogaster_gene_age$ensembl_gene_id]
melanogaster_brain_common_dec_gene_age <- melanogaster_gene_age[melanogaster_gene_age$ensembl_gene_id %in% melanogaster_brain_common_dec_gene_names,]
melanogaster_brain_common_dec_gene_age <- melanogaster_brain_common_dec_gene_age[order(melanogaster_brain_common_dec_gene_age$ensembl_gene_id),]
melanogaster_brain_common_dec_max_exp <- apply(X= Results_pacifico[["head"]][["deg_list"]][["dec_sig"]][["dec_genes_exp"]][Results_pacifico[["head"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% melanogaster_gene_age$ensembl_gene_id,],
                                         MARGIN = 1, function(x){max(x)})

identical(as.character(melanogaster_brain_common_dec_gene_age$ensembl_gene_id), melanogaster_brain_common_dec_gene_names)

melanogaster_brain_dnds_gene_age <- rbind(
  
  data.frame(gene_age = melanogaster_brain_common_inc_gene_age$gene_age ,
             gene_class = "old-biased",
             dnds = melanogaster_brain_common_inc_dnds,
             max_exp = melanogaster_brain_common_inc_max_exp),
  
  data.frame(gene_age = melanogaster_brain_common_dec_gene_age$gene_age ,
             gene_class = "young-biased",
             dnds = melanogaster_brain_common_dec_dnds,
             max_exp = melanogaster_brain_common_dec_max_exp)
)

summary(lm(dnds ~ max_exp + gene_age + gene_class, data= melanogaster_brain_dnds_gene_age))
