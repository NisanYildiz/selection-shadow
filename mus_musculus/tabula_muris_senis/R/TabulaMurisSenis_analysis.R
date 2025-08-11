### Script for the analysis of the Tabula Muris Senis data 

#Loading the libraries.
library("biomaRt")
library(ggplot2)
library(tidyr)
library(grid)
library(dplyr)

#Loading the dataset
celltype_perind_3m <- readRDS("mus_musculus/tabula_muris_senis/R/results/3m.celltype.per.ind.rds")
celltype_perind_18m <- readRDS("mus_musculus/tabula_muris_senis/R/results/18m.celltype.per.ind.rds")
celltype_perind_24m <- readRDS("mus_musculus/tabula_muris_senis/R/results/24m.celltype.per.ind.rds")

## For each tissue present in the lists, find the cell types that are common in all 3 age groups.

#Getting the list of tissues present in all samples
tissues <- Reduce(intersect, list(names(celltype_perind_3m),names(celltype_perind_18m),
                                  names(celltype_perind_24m)))

unique_cell_types_per_tissue <- 
  sapply(tissues, function(x){
    cell_types_3m <- names(celltype_perind_3m[[x]])
    cell_types_18m <- names(celltype_perind_18m[[x]])
    cell_types_24m <- names(celltype_perind_24m[[x]])
    cell_types_all <- unique(c(cell_types_3m, cell_types_18m, cell_types_24m))
    cell_types_all
  })


#for each tissue, finding cell types that are present in all age groups
#We'll continue the analyses with these cell types
common_cell_types_per_tissue <- 
  sapply(tissues, function(x){
    cell_types_3m <- names(celltype_perind_3m[[x]])
    cell_types_18m <- names(celltype_perind_18m[[x]])
    cell_types_24m <- names(celltype_perind_24m[[x]])
    cell_types_common <- Reduce(intersect, list(cell_types_3m, cell_types_18m, cell_types_24m))
    cell_types_common
  })

##Getting all the genes present in the dataset
#since all the genes within the tissues are the common genes,
#using the gene set for a single cell type per tissue is sufficent 

genes_3m <- lapply(tissues, function(x){
  rownames(celltype_perind_3m[[x]][[1]])
})
genes_3m <- unique(Reduce(c,genes_3m))

genes_18m <- lapply(tissues, function(x){
  rownames(celltype_perind_18m[[x]][[1]])
})
genes_18m <- unique(Reduce(c,genes_18m))

genes_24m <- lapply(tissues, function(x){
  rownames(celltype_perind_24m[[x]][[1]])
})
genes_24m <- unique(Reduce(c,genes_24m))

genes <- unique(Reduce(c,list(genes_3m, genes_18m, genes_24m)))

##Gene name to ENSEMBL ID conversion

ensembl99 <- useEnsembl("ensembl", version = 99, 
                        dataset = "mmusculus_gene_ensembl")

gene_ids <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), 
                  filters = 'external_gene_name', values = genes, 
                  mart = ensembl99, useCache = F)
gene_ids <- gene_ids[!(gene_ids$external_gene_name %in% gene_ids$external_gene_name[which(duplicated(gene_ids$external_gene_name))]),]#getting rid of duplicated genes

##Creating a list of expression matrices

expression_list <- lapply(tissues, function(x){
  lapply(common_cell_types_per_tissue[[x]], function(y){
    exp_3m <- as.data.frame(celltype_perind_3m[[x]][[y]])
    exp_3m$rownames <- rownames(exp_3m)
    exp_18m <- as.data.frame(celltype_perind_18m[[x]][[y]])
    exp_18m$rownames <- rownames(exp_18m)
    exp_24m <-  as.data.frame(celltype_perind_24m[[x]][[y]])
    exp_24m$rownames <- rownames(exp_24m)
    
    exp_all <- full_join(exp_3m, exp_18m, by = "rownames")
    rownames(exp_all) <- exp_all$rownames 
    exp_all <- full_join(exp_all, exp_24m, by = "rownames")
    rownames(exp_all) <- exp_all$rownames 
    exp_all <- subset(exp_all, select = -c(rownames))
    
    exp_all[setdiff(genes, rownames(exp_all)), ] <- NA
    
    #setting NA values to zero
    exp_all[is.na(exp_all)] <- 0
    
    exp_all
  })
  
  
})

expression_list_sep <- lapply(tissues, function(x){
  lapply(common_cell_types_per_tissue[[x]], function(y){
    exp_3m <- as.data.frame(celltype_perind_3m[[x]][[y]])
    exp_3m$rownames <- rownames(exp_3m)
    exp_18m <- as.data.frame(celltype_perind_18m[[x]][[y]])
    exp_18m$rownames <- rownames(exp_18m)
    exp_24m <-  as.data.frame(celltype_perind_24m[[x]][[y]])
    exp_24m$rownames <- rownames(exp_24m)
    
    exp_all <- full_join(exp_3m, exp_18m, by = "rownames")
    rownames(exp_all) <- exp_all$rownames 
    exp_all <- full_join(exp_all, exp_24m, by = "rownames")
    rownames(exp_all) <- exp_all$rownames 
    exp_all <- subset(exp_all, select = -c(rownames))
    
    #setting NA values to zero
    exp_all[is.na(exp_all)] <- 0
    
    exp_all
  })
  
  
})

names(expression_list) <- tissues
names(expression_list_sep) <- tissues

for(tissue in tissues){
  names(expression_list[[tissue]]) <- common_cell_types_per_tissue[[tissue]]
  names(expression_list_sep[[tissue]]) <- common_cell_types_per_tissue[[tissue]]
}



## Loading dnds values between rat and mice

#loading the values

dnds <- getBM(attributes=c("ensembl_gene_id", "rnorvegicus_homolog_ensembl_gene", 
                           "rnorvegicus_homolog_dn", "rnorvegicus_homolog_ds"),
              filters = "ensembl_gene_id", values = gene_ids$ensembl_gene_id,
              mart = ensembl99,
              useCache = F)

#cleaning the dnds data
dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_ds)),] #removing NA dS values
dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_dn)),] #removing NA dN values
dnds <- dnds[!(dnds$rnorvegicus_homolog_dn == 0),] #removing 0 dN values
dnds <- dnds[!(dnds$rnorvegicus_homolog_ds == 0),] #removing 0 dS values
dNdS <- (dnds$rnorvegicus_homolog_dn / dnds$rnorvegicus_homolog_ds) #calculating dn/ds
dnds <- cbind(dnds, dNdS)#adding dn/ds vector to the dataframe 
###removing duplicate gene ids, leaves only 1to1 orth. genes
dnds <- dnds[!((duplicated(dnds$ensembl_gene_id)) | (duplicated(dnds$ensembl_gene_id, fromLast = T))) ,]
dnds <- dnds[!((duplicated(dnds$rnorvegicus_homolog_ensembl_gene)) | (duplicated(dnds$rnorvegicus_homolog_ensembl_gene, fromLast = T))) ,]
###removing genes with dnds > 0.8
dnds <- dnds[dnds$dNdS < 0.8,]

#saveRDS(dnds, file = "/mnt/NAS/projects/2020_melih/mus_musculus/tabula_muris_senis/R/dnds.rds")

## Gene ID conversion and getting dnds values for expression matrix

#Converts gene names in rownames to ENSEMBL gene ids and outputs the matrix
IDconv <- function(expMatrix, gene_ids){
  exp_matrix <- expMatrix[rownames(expMatrix) %in% gene_ids$external_gene_name,]
  gene_names <- gene_ids[gene_ids$external_gene_name %in% rownames(exp_matrix),]
  exp_matrix <- exp_matrix[order(rownames(exp_matrix)),]
  gene_names <- gene_names[order(gene_names$external_gene_name),]
  
  rownames(exp_matrix) <- gene_names$ensembl_gene_id
  
  exp_matrix
}

#Subsets the genes that are present both in dnds data.frame and matrix
#outputs the subset expression matrix and dnds data.frame as a list
#Both have the same order of genes
dndsCommon <- function(exp_matrix, dnds){
  exp_matrix_subset <- exp_matrix[rownames(exp_matrix) %in% dnds$ensembl_gene_id,]
  dnds_subset <- dnds[dnds$ensembl_gene_id %in% rownames(exp_matrix_subset),]
  #ordering
  dnds_subset <- dnds_subset[(order(dnds_subset$ensembl_gene_id)),]
  exp_matrix_subset <- exp_matrix_subset[order(rownames(exp_matrix_subset)),]
  
  a_list <- list(dnds = dnds_subset, exp_matrix = exp_matrix_subset)
  a_list
}


##Differential expression analysis

DEG_genes <- function(expMatrix, ageVector, rhoCutoff, dnds){
  Aging_exp_cor <- t(apply(expMatrix,1,function(x){
    test <- cor.test(x, ageVector, method="spearman", exact = F, use= "everything")
    c(test$est,test$p.val)
  }))
  
  adj_p.val <- p.adjust(Aging_exp_cor[,2], "BH") #BH p adjustment
  Aging_exp_cor <- cbind(Aging_exp_cor, adj_p.val)
  colnames(Aging_exp_cor) = c("rho", "p", "adjusted_p")
  
  all_genes_cor <- Aging_exp_cor
  all_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(all_genes_cor),]
  all_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(all_genes_cor),]
  
  #removing genes not expressed at all
  Aging_exp_cor <- Aging_exp_cor[complete.cases(Aging_exp_cor), ]
  
  selected_genes_cor <- Aging_exp_cor[Aging_exp_cor[,1] > rhoCutoff | Aging_exp_cor [,1] < -(rhoCutoff),]
  selected_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(selected_genes_cor),]
  selected_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(selected_genes_cor),]
  
  inc_genes_cor <- Aging_exp_cor[Aging_exp_cor[,1] > rhoCutoff ,]
  inc_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(inc_genes_cor),]
  inc_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(inc_genes_cor),]
  
  dec_genes_cor <- Aging_exp_cor[Aging_exp_cor [,1] < -(rhoCutoff) ,]
  dec_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(dec_genes_cor),]
  dec_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(dec_genes_cor),]
  
  nonsig_genes_cor <- Aging_exp_cor[Aging_exp_cor[,1] < rhoCutoff | Aging_exp_cor [,1] > -(rhoCutoff),]
  nonsig_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(nonsig_genes_cor),]
  nonsig_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(nonsig_genes_cor),]
  
  listy <- list(Aging_exp_cor = all_genes_cor, 
                All = list(all_genes_cor = all_genes_cor, all_genes_exp = all_genes_exp, all_genes_dnds = all_genes_dnds),
                Non_sig = list(nonsig_genes_cor = nonsig_genes_cor, nonsig_genes_exp = nonsig_genes_exp, nonsig_genes_dnds = nonsig_genes_dnds), 
                All_sig = list(sig_genes_cor = selected_genes_cor, sig_genes_exp = selected_genes_exp, sig_genes_dnds = selected_genes_dnds), 
                inc_sig = list(inc_genes_cor = inc_genes_cor, inc_genes_exp = inc_genes_exp, inc_genes_dnds = inc_genes_dnds), 
                dec_sig = list(dec_genes_cor = dec_genes_cor, dec_genes_exp = dec_genes_exp, dec_genes_dnds = dec_genes_dnds))
  listy
}



Results_tabula <- lapply(names(expression_list), function(tissues){
  
  
  Results_rho05 <- lapply(expression_list[[tissues]], function(expression){
    exp_matrix <- as.matrix(expression)
    exp_matrix <- IDconv(expMatrix = exp_matrix, gene_ids = gene_ids) #converting ids
    
    ##Finding a large enough subset of shared individuals
    
    individuals_3m <- lapply(celltype_perind_3m[[tissues]], FUN = colnames)
    
    presence <- sapply(unique(unlist(individuals_3m)), function(indv){
      
      sapply(individuals_3m, function(cell_type){
        
        indv %in% cell_type
        
      })
      
    })
    
    common_ind_3m <- unique(unlist(individuals_3m))[apply(presence, 2, function(indv){
      
      sum(indv) / length(indv)
      
    }) > 0.7]
    
    # common_ind_3m <- Reduce(intersect, individuals_3m[sapply(individuals_3m, function(x){
    #   length(x) > 4})])
    
    individuals_18m <- lapply(celltype_perind_18m[[tissues]], FUN = colnames)
    
    presence <- sapply(unique(unlist(individuals_18m)), function(indv){
      
      sapply(individuals_18m, function(cell_type){
        
        indv %in% cell_type
        
      })
      
    })
    
    common_ind_18m <- unique(unlist(individuals_18m))[apply(presence, 2, function(indv){
      
      sum(indv) / length(indv)
      
    }) > 0.7]
    # common_ind_18m <- Reduce(intersect, individuals_18m[sapply(individuals_18m, function(x){
    #   length(x) > 3})])
    
    individuals_24m <- lapply(celltype_perind_24m[[tissues]], FUN = colnames)
    
    presence <- sapply(unique(unlist(individuals_24m)), function(indv){
      
      sapply(individuals_24m, function(cell_type){
        
        indv %in% cell_type
        
      })
      
    })
    
    common_ind_24m <- unique(unlist(individuals_24m))[apply(presence, 2, function(indv){
      
      sum(indv) / length(indv)
      
    }) > 0.7]
    
    
    
    # common_ind_24m <- Reduce(intersect, individuals_24m[sapply(individuals_24m, function(x){
    #   length(x) > 1})])
    
    #END
    
    exp_matrix <- exp_matrix[,(colnames(exp_matrix) %in% common_ind_3m | colnames(exp_matrix) %in% common_ind_18m | colnames(exp_matrix) %in% common_ind_24m)]
    
    dnds_exp_list <- dndsCommon(exp_matrix = exp_matrix, dnds = dnds) 
    
    ## Getting age vector 
    
    age <- colnames(exp_matrix)
    age <- sub("(.*)_.*", "\\1", age)
    age <- sub("(.*)_.*", "\\1", age)
    age <- as.numeric(age)
    
    deg_list <- DEG_genes(expMatrix = dnds_exp_list[["exp_matrix"]], ageVector = age, 
                          rhoCutoff = 0.5, dnds = dnds_exp_list[["dnds"]])
    
    complete_results = t(apply(deg_list[["All"]][["all_genes_exp"]], 2 , function(x){
      cor(c(-log(deg_list[["All"]][["all_genes_dnds"]]$dNdS)), x,  method= "spearman",
          use = "everything")
    }))
    
    all_results = t(apply(deg_list[["All_sig"]][["sig_genes_exp"]], 2 , function(x){
      cor(c(-log(deg_list[["All_sig"]][["sig_genes_dnds"]]$dNdS)), x,  method= "spearman",
          use = "everything")
    }))
    
    inc_results = t(apply(deg_list[["inc_sig"]][["inc_genes_exp"]], 2 , function(x){
      cor(c(-log(deg_list[["inc_sig"]][["inc_genes_dnds"]]$dNdS)), x,  method= "spearman",
          use = "everything")
    }))
    
    dec_results <- t(apply(deg_list[["dec_sig"]][["dec_genes_exp"]], 2 , function(x){
      cor(c(-log(deg_list[["dec_sig"]][["dec_genes_dnds"]]$dNdS)), x,  method= "spearman",
          use = "everything")
    }))
    
    nonsig_results <- t(apply(deg_list[["Non_sig"]][["nonsig_genes_exp"]], 2 , function(x){
      cor(c(-log(deg_list[["Non_sig"]][["nonsig_genes_dnds"]]$dNdS)), x,  method= "spearman",
          use = "everything")
    }))
    
    list(age=age , complete_results = complete_results, all_results = all_results, 
         inc_results = inc_results, 
         dec_results = dec_results, nonsig_results = nonsig_results, deg_list = deg_list,
         dnds_exp_list = dnds_exp_list)
    
  })
  
  
})

names(Results_tabula) <- names(expression_list)

Results_tabula_muris_senis <- Results_tabula

save(Results_tabula_muris_senis, file = "mus_musculus/tabula_muris_senis/R/results/Results_tabula_muris_senis.RData")

##Conservation metric statistics

dnds <- getBM(attributes=c("ensembl_gene_id", "rnorvegicus_homolog_ensembl_gene", 
                           "rnorvegicus_homolog_dn", "rnorvegicus_homolog_ds"),
              filters = "ensembl_gene_id", values = gene_ids$ensembl_gene_id,
              mart = ensembl99,
              useCache = F)

dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_ds)),] #removing NA dS values
dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_dn)),] #removing NA dN values

table <- sapply(tissues, function(tissue){
  
  dnds = dndsCommon(IDconv(expression_list_sep[[tissue]][[1]], gene_ids), dnds)$dnds #tissue level intersection
  ###removing duplicate gene ids, leaves only 1to1 orth. genes
  dnds <- dnds[!((duplicated(dnds$ensembl_gene_id)) | (duplicated(dnds$ensembl_gene_id, fromLast = T))) ,]
  dnds <- dnds <- dnds[!((duplicated(dnds$rnorvegicus_homolog_ensembl_gene)) | (duplicated(dnds$rnorvegicus_homolog_ensembl_gene, fromLast = T))) ,]
  n_1to1 <- nrow(dnds)
  
  dnds <- dnds[!(dnds$rnorvegicus_homolog_dn == 0),] #removing 0 dN values
  dnds <- dnds[!(dnds$rnorvegicus_homolog_ds == 0),] #removing 0 dS values
  dNdS <- (dnds$rnorvegicus_homolog_dn / dnds$rnorvegicus_homolog_ds) #calculating dn/ds
  dnds <- cbind(dnds, dNdS)#adding dn/ds vector to the dataframe 
  dnds_higher_08 <- sum(dnds$dNdS >= 0.8)
  dnds <- dnds[dnds$dNdS < 0.8,]
  mean <- signif(mean(dnds$dNdS), 3)
  median <- signif(median(dnds$dNdS), 3)
  final <- nrow(dnds)
  
  cbind(n_1to1 = n_1to1, dnds_higher_08 = dnds_higher_08, final_n = final, mean = mean, median = median)
})

cons_metric_table <- t(table)
colnames(cons_metric_table) = c("n_1to1", "dnds >= 0.8", "n_final", "mean", "median")



### Calculating sample size, number of cells etc. 

#calculating number of cells and cell type

metadata <- readRDS("mus_musculus/tabula_muris_senis/R/results/metadata.rds")


n_cells <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    no_of_cells <- sum(metadata[[tissue]][["mouse.id"]] %in% colnames(Results_tabula[[tissue]][[cell_type]][["complete_results"]]) & 
                         metadata[[tissue]][["cell_ontology_class"]] == cell_type )
    
    
  })
})
names(n_cells) <- tissues

n_cells_3m <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    no_of_cells <- sum(metadata[[tissue]][["mouse.id"]] %in% colnames(Results_tabula[[tissue]][[cell_type]][["complete_results"]]) & 
                         (metadata[[tissue]][["cell_ontology_class"]] == cell_type & metadata[[tissue]][["age"]] == "3m") )
    
    
  })
})
names(n_cells_3m) <- tissues

n_cells_18m <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    no_of_cells <- sum(metadata[[tissue]][["mouse.id"]] %in% colnames(Results_tabula[[tissue]][[cell_type]][["complete_results"]]) & 
                         (metadata[[tissue]][["cell_ontology_class"]] == cell_type & metadata[[tissue]][["age"]] == "18m") )
    
    
  })
})
names(n_cells_18m) <- tissues

n_cells_24m <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    no_of_cells <- sum(metadata[[tissue]][["mouse.id"]] %in% colnames(Results_tabula[[tissue]][[cell_type]][["complete_results"]]) & 
                         (metadata[[tissue]][["cell_ontology_class"]] == cell_type & metadata[[tissue]][["age"]] == "24m") )
    
    
  })
})
names(n_cells_24m) <- tissues

n_indv <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    
    no_of_indv <- length(Results_tabula[[tissue]][[cell_type]][["age"]])
    
  })
})
names(n_indv) <- tissues

n_indv_3m <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    
    no_of_indv <- sum(Results_tabula[[tissue]][[cell_type]][["age"]] == 3)
    
  })
})
names(n_indv_3m) <- tissues

n_indv_18m <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    
    no_of_indv <- sum(Results_tabula[[tissue]][[cell_type]][["age"]] == 18)
    
  })
})
names(n_indv_18m) <- tissues

n_indv_24m <- lapply(tissues, function(tissue){
  
  
  sapply(common_cell_types_per_tissue[[tissue]], function(cell_type){
    
    
    no_of_indv <- sum(Results_tabula[[tissue]][[cell_type]][["age"]] == 24)
    
  })
})
names(n_indv_24m) <- tissues


n_table <- data.frame(n_cells = unlist(n_cells),
                      n_cells_3m = unlist(n_cells_3m),
                      n_cells_18m = unlist(n_cells_18m),
                      n_cells_24m = unlist(n_cells_24m),
                      n_indv = unlist(n_indv), 
                      n_indv_3m = unlist(n_indv_3m),
                      n_indv_18m = unlist(n_indv_18m),
                      n_indv_24m = unlist(n_indv_24m))

write.csv(n_table, file = "supplements/tabula_afterfilter_samplesize.csv",quote = T)

####add n_cells_3m, n_cells_18m, n_cells_24m and n_indv_3m, n_indv_18m, n_indv_24m


### Fisher's test for overrepresentation of immune genes in old-biased low-conserv 
# and young-biased high-conserv genes

#using the dnds of a random cell, bcs all the same set
dnds <- Results_tabula_muris_senis[["lung"]][["bronchial smooth muscle cell"]][["dnds_exp_list"]][["dnds"]]
low_conserv_genes = dnds$ensembl_gene_id[dnds$dNdS > quantile(dnds$dNdS, 0.75)] 
high_conserv_genes = dnds$ensembl_gene_id[dnds$dNdS < quantile(dnds$dNdS, 0.25)]


immune_gene_results <- AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0002376"), columns = c('ENSEMBL'), keytype = "GOALL")
immune_gene_ids <- unique(immune_gene_results$ENSEMBL)




fishers_test_results_inc <- data.frame()

for (tissue in names(Results_tabula_muris_senis)){
  
  for (cell_type in names(Results_tabula_muris_senis[[tissue]])){
    
    deg_list <- Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]]
    inc_genes <- deg_list[["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]]
    
    inc_low_conserv <- intersect(inc_genes, low_conserv_genes)
    inc_remain <- setdiff(inc_genes, low_conserv_genes)
    
    
    inc_low_conserv_immune <- intersect(inc_low_conserv, immune_gene_ids)
    inc_low_conserv_non_immune <- setdiff(inc_low_conserv, immune_gene_ids)
    
    inc_remain_immune <- intersect(inc_remain, immune_gene_ids)
    inc_remain_non_immune <- setdiff(inc_remain, immune_gene_ids)
    
    cont_table <- matrix(c(length(inc_remain_immune), length(inc_remain_non_immune)
                           , length(inc_low_conserv_immune), length(inc_low_conserv_non_immune))
                         , nrow = 2,
                         dimnames =
                           list(c("immune", "non_immune"),
                                c("inc_rest", "inc_low_conserv")))
    
    
    fisher_result <- fisher.test(cont_table)
    values <- data.frame(tissue, cell_type, fisher_result$estimate, fisher_result$p.value)
    colnames(values) <- c("tissue", "cell_type", "odds_ratio", "p_val")
    
    fishers_test_results_inc <- rbind(fishers_test_results_inc, values)
    
    
  }
  
}

fishers_test_results_inc$p_adj <- p.adjust(fishers_test_results_inc$p_val, "BH")


fishers_test_results_dec <- data.frame()

for (tissue in names(Results_tabula_muris_senis)){
  
  for (cell_type in names(Results_tabula_muris_senis[[tissue]])){
    
    deg_list <- Results_tabula_muris_senis[[tissue]][[cell_type]][["deg_list"]]
    dec_genes <- deg_list[["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]]
    
    dec_high_conserv <- intersect(dec_genes, high_conserv_genes)
    dec_remain <- setdiff(dec_genes, high_conserv_genes)
    
    
    dec_high_conserv_immune <- intersect(dec_high_conserv, immune_gene_ids)
    dec_high_conserv_non_immune <- setdiff(dec_high_conserv, immune_gene_ids)
    
    dec_remain_immune <- intersect(dec_remain, immune_gene_ids)
    dec_remain_non_immune <- setdiff(dec_remain, immune_gene_ids)
    
    cont_table <- matrix(c(length(dec_remain_immune), length(dec_remain_non_immune)
                           , length(dec_high_conserv_immune), length(dec_high_conserv_non_immune))
                         , nrow = 2,
                         dimnames =
                           list(c("immune", "non_immune"),
                                c("dec_rest", "dec_high_conserv")))
    
    
    fisher_result <- fisher.test(cont_table)
    values <- data.frame(tissue, cell_type, fisher_result$estimate, fisher_result$p.value)
    colnames(values) <- c("tissue", "cell_type", "odds_ratio", "p_val")
    
    fishers_test_results_dec <- rbind(fishers_test_results_dec, values)
    
    
  }
  
}
fishers_test_results_dec$p_adj <- p.adjust(fishers_test_results_dec$p_val, "BH")




mean(fishers_test_results_inc$odds_ratio)
mean(fishers_test_results_inc$odds_ratio[fishers_test_results_inc$p_adj < 0.1])

mean(fishers_test_results_dec$odds_ratio[fishers_test_results_dec$odds_ratio != Inf])


