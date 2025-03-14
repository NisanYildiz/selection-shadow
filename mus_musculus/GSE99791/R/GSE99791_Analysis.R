# Analysis of GSE99791

#Loading the libraries
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(plyr)
library("biomaRt")
library(grid)


#Loading up the expression matrix
expmat <- read.table(header = T, row.names = 1,
                     file = "data/mus_musculus/GSE99791/all_counts.txt")


expmat <- as.matrix(expmat)
#Loading the metadata 
metadata <- read.csv(file = "mus_musculus/GSE99791/SraRunTable_GSE99791.txt",
                     row.names = 1)

input_samples <- c("GSM2652319", "GSM2652323", "GSM2652327", "GSM2652331", "GSM2652335",
                   "GSM2652339", "GSM2652343", "GSM2652347", "GSM2698018")



##reordering the metadata 
metadata <- metadata[order(row.names(metadata)),]

## removing input samples 
metadata <- metadata[!(metadata$GEO_Accession..exp. %in% input_samples),]
expmat <- expmat[,colnames(expmat) %in% rownames(metadata)]

##Checking if the order is correct
identical(row.names(metadata), colnames(expmat)) ###TRUE

#Converting age to numeric values
metadata$AGE <- plyr::revalue(metadata$AGE, c("adult (4 month old)"=4, "aged (2 year old)"=24))

#removing somatosensory cortex data, as it only has 4-month-old individuals 
expmat <- exp <- expmat[,!(metadata$Tissue == "somatosensory cortex")]
metadata <- metadata[!(metadata$Tissue == "somatosensory cortex"),]
metadata <- droplevels(metadata)


#separating exp matrix and metadata into different tissues
tissue_exp <- lapply(levels(metadata$Tissue), function(tissue){
  exp <- expmat[,metadata$Tissue == tissue]
  meta <- metadata[metadata$Tissue == tissue,]
  a_list <- list(expmat = exp, metadata = meta)
  a_list
})
names(tissue_exp) <- levels(metadata$Tissue)

##Normalization

### Creating the DESeq object
tissue_dds <- lapply(levels(metadata$Tissue), function(tissue){
  dds <- DESeqDataSetFromMatrix(countData = tissue_exp[[tissue]][["expmat"]], 
                                colData = tissue_exp[[tissue]][["metadata"]], 
                                design = ~AGE)
  
})
names(tissue_dds) <- levels(metadata$Tissue)

### Removing 0 count genes

tissue_dds <- lapply(tissue_dds, function(dds){
  
  nonzero_genes <- rowSums(counts(dds)) > 0
  dds <- dds[ nonzero_genes, ]
  
})

### Size factor estimatiıon for DESeq normalization

tissue_dds <- lapply(tissue_dds, function(dds){

  dds <- estimateSizeFactors(dds)
  
})

### Normalizing

normalized_matrices <- lapply(tissue_dds, function(dds){
  
  counts(dds, normalized=TRUE)
  
})

## Basic exploratory visualizations 

for( tissue in names(normalized_matrices) ){
  
  boxplot(log2(normalized_matrices[[tissue]] + 1), notch=TRUE ,main = paste0("Size-factor-normalized read counts
        Mus musculus astrocyte from ", tissue, "(GSE114129)"), col = "lightblue3", las = 2, cex.axis = .8)
  hist(log2(normalized_matrices[[tissue]] + 1), col = "lightblue3", main = paste0("Distribution of size-factor-normalized counts
     Mus musculus astrocyte from ", tissue, "(GSE114129)"),
       xlab = "log2(NormalizedCounts + 1)")  
  
  
}

pca_mus_1_2 <- lapply(names(normalized_matrices), function(tissue){
  
  pc_mus = prcomp(t(normalized_matrices[[tissue]]), scale = T)
  percent = paste0(colnames(pc_mus$x[,1:4])," (", round(summary(pc_mus)$importance[2,1:4] * 100, 2), "%) ")
  
  data.frame(pc_mus$x[,1:4], Age= factor(tissue_exp[[tissue]][["metadata"]][["AGE"]], levels = c("4", "24"))) %>% 
    ggplot(aes(PC1, PC2, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent[1]) + ylab(percent[2])
  
})
names(pca_mus_1_2) <- names(normalized_matrices)

pca_mus_3_4 <- lapply(names(normalized_matrices), function(tissue){
  
  pc_mus = prcomp(t(normalized_matrices[[tissue]]), scale = T)
  percent = paste0(colnames(pc_mus$x[,1:4])," (", round(summary(pc_mus)$importance[2,1:4] * 100, 2), "%) ")
  
  
  data.frame(pc_mus$x[,1:4], Age= factor(tissue_exp[[tissue]][["metadata"]][["AGE"]], levels = c("4", "24"))) %>% 
    ggplot(aes(PC3, PC4, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent[3]) + ylab(percent[4])
  
})
names(pca_mus_3_4) <- names(normalized_matrices)

save(pca_mus_1_2, pca_mus_3_4, file = "mus_musculus/GSE99791/R/results/pca_mus.RData")


## Downloading the ENSEMBL99 dn/ds values between mouse and rat

ensembl99 <- useEnsembl("ensembl", version = 99, 
                        dataset = "mmusculus_gene_ensembl")

tissue_dnds <- lapply(normalized_matrices, function(tissue){
  
  dnds <- getBM(attributes=c("ensembl_gene_id", "rnorvegicus_homolog_ensembl_gene", 
                             "rnorvegicus_homolog_dn", "rnorvegicus_homolog_ds"),
                filters = "ensembl_gene_id", values = rownames(tissue),
                mart = ensembl99,
                useCache = FALSE)
  
})

##cleaning the dnds data

tissue_dnds <- lapply(tissue_dnds, function(dnds){
  
  dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_ds)),] #removing NA dS values
  dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_dn)),] #removing NA dN values
  dnds <- dnds[!(dnds$rnorvegicus_homolog_dn == 0),] #removing 0 dN values
  dnds <- dnds[!(dnds$rnorvegicus_homolog_ds == 0),] #removing 0 dS values
  dNdS <- (dnds$rnorvegicus_homolog_dn / dnds$rnorvegicus_homolog_ds) #calculating dn/ds
  dnds <- cbind(dnds, dNdS)#adding dn/ds vector to the dataframe 
  ###removing duplicate gene ids, leaves only 1to1 orth. genes
  dnds <- dnds[!((duplicated(dnds$ensembl_gene_id)) | (duplicated(dnds$ensembl_gene_id, fromLast = T))) ,]
  dnds <- dnds[!((duplicated(dnds$rnorvegicus_homolog_ensembl_gene)) | (duplicated(dnds$rnorvegicus_homolog_ensembl_gene, fromLast = T))) ,]
  ###removing genes with dnds >= 0.8
  dnds <- dnds[dnds$dNdS < 0.8,]
  
})

###--- Setting up functions ----


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

DEG_genes <- function(expMatrix, ageVector, pAdjCutoff, rhoCutoff, dnds){
  Aging_exp_cor <- t(apply(expMatrix,1,function(x){
    test <- cor.test(x, ageVector, method="spearman", exact = F, use= "everything")
    c(test$est,test$p.val)
  }))
  
  adj_p.val <- p.adjust(Aging_exp_cor[,2], "BH") #BH p adjustment
  Aging_exp_cor <- cbind(Aging_exp_cor, adj_p.val)
  colnames(Aging_exp_cor) = c("rho", "p", "adjusted_p")
  
  Aging_exp_cor <- Aging_exp_cor[complete.cases(Aging_exp_cor), ]
  
  all_genes_cor <- Aging_exp_cor
  all_genes_exp <- expMatrix
  all_genes_dnds <- dnds
  
  DEG_genes_cor <- Aging_exp_cor[Aging_exp_cor[,"adjusted_p"] < pAdjCutoff & abs(Aging_exp_cor[,1]) > rhoCutoff ,]
  DEG_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(DEG_genes_cor),]
  DEG_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(DEG_genes_cor),]
  
  inc_genes_cor <- DEG_genes_cor[DEG_genes_cor[,1] > rhoCutoff ,]
  inc_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(inc_genes_cor),]
  inc_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(inc_genes_cor),]
  
  dec_genes_cor <- DEG_genes_cor[DEG_genes_cor [,1] < -(rhoCutoff) ,]
  dec_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(dec_genes_cor),]
  dec_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(dec_genes_cor),]
  
  nonsig_genes_cor <- Aging_exp_cor[Aging_exp_cor[,"adjusted_p"] > pAdjCutoff,]
  nonsig_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(nonsig_genes_cor),]
  nonsig_genes_dnds <- dnds[dnds$ensembl_gene_id %in% rownames(nonsig_genes_cor),]
  
  listy <- list(Aging_exp_cor = Aging_exp_cor, 
                All = list(all_genes_cor = all_genes_cor, all_genes_exp = all_genes_exp, all_genes_dnds = all_genes_dnds),
                Non_sig = list(nonsig_genes_cor = nonsig_genes_cor, nonsig_genes_exp = nonsig_genes_exp, nonsig_genes_dnds = nonsig_genes_dnds), 
                All_sig = list(sig_genes_cor = DEG_genes_cor, sig_genes_exp = DEG_genes_exp, sig_genes_dnds = DEG_genes_dnds), 
                inc_sig = list(inc_genes_cor = inc_genes_cor, inc_genes_exp = inc_genes_exp, inc_genes_dnds = inc_genes_dnds), 
                dec_sig = list(dec_genes_cor = dec_genes_cor, dec_genes_exp = dec_genes_exp, dec_genes_dnds = dec_genes_dnds))
  listy
}

###---- DONE ---- 


###DEG analysis

Results_GSE99791 <- lapply(levels(metadata$Tissue), function(tissue){
  
  
  dnds_exp_list <- dndsCommon(exp_matrix = normalized_matrices[[tissue]], 
                              dnds = tissue_dnds[[tissue]]) 
  
  ## Getting age vector 
  age <- as.numeric(as.character(tissue_exp[[tissue]][["metadata"]][["AGE"]]))
  
  deg_list <- DEG_genes(expMatrix = dnds_exp_list[["exp_matrix"]], ageVector = age,
                        pAdjCutoff = 0.1, rhoCutoff = 0.5, dnds = dnds_exp_list[["dnds"]])
  
  all_results = t(apply(deg_list[["All"]][["all_genes_exp"]], 2 , function(x){
    cor(c(-log(deg_list[["All"]][["all_genes_dnds"]]$dNdS)), x,  method= "spearman",
        use = "everything")
  }))
  
  sig_results = t(apply(deg_list[["All_sig"]][["sig_genes_exp"]], 2 , function(x){
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
  
  Results <- list(age=age , all_results = all_results, sig_results = sig_results, inc_results = inc_results, 
                  dec_results = dec_results, nonsig_results = nonsig_results, deg_list = deg_list,
                  dnds_exp_list = dnds_exp_list)
  
  
})
names(Results_GSE99791) <- levels(metadata$Tissue)

save(Results_GSE99791, file = "mus_musculus/GSE99791/R/results/Results_GSE99791.RData")

### Exp-conservation correlation vs. Age


for (tissue in names(Results_GSE99791)){
  
  p = cor.test(as.numeric(Results_GSE99791[[tissue]][["sig_results"]]), Results_GSE99791[[tissue]][["age"]], method = "spearman", exact = F)$p.val
  p = signif(p,2)
  corr = cor(as.numeric(Results_GSE99791[[tissue]][["sig_results"]]), Results_GSE99791[[tissue]][["age"]], method = "spearman")
  corr = signif(corr,2)
  plot(Results_GSE99791[[tissue]][["age"]], Results_GSE99791[[tissue]][["sig_results"]], main = paste0("Differentially Expressed Genes in Aging
     Mus musculus astrocyte from ",tissue, "(GSE99791)"), xlab = "Age(days)",
       ylab = "Expression-Conservation rho", col = "lightblue3",
       pch = 19, cex = 1.6)
  abline(lm(as.numeric(Results_GSE99791[[tissue]][["sig_results"]]) ~ Results_GSE99791[[tissue]][["age"]]),
         col = "#2B4C6F")
  grid.text(x = unit(0.18, "npc"), y = unit(0.25, "npc"), just = c(1,1), 
            label = paste("p =",p))
  grid.text(x = unit(0.18, "npc"), y = unit(0.22, "npc"), just = c(1,1),
            label = paste("corr =",corr))
  grid.text(x = unit(0.20, "npc"), y = unit(0.19, "npc"), just = c(1,1),
            label = paste("# of genes =",length(Results_GSE99791[[tissue]][["deg_list"]][["All_sig"]][["sig_genes_exp"]][,1])))
  
  
}


