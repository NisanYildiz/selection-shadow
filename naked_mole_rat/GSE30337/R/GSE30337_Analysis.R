#Analysis of GSE30337

#Loading the libraries
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(plyr)
library("biomaRt")
library(grid)

#misc. functions
source("functions.R")

#Loading up the expression matrix
expmat <- read.table(header = T, row.names = 1,
                     file = "data/naked_mole_rat/GSE30337/all_counts.txt")

#Loading the metadata 
metadata <- read.csv(file = "naked_mole_rat/GSE30337/SraRunTable-GSE30337.txt",
                     row.names = 1)

##Checking if the order is correct
identical(row.names(metadata), colnames(expmat)) ###TRUE

##removing hypoxia treated samples, mixed tissue samples and newborn samples
metadata <- metadata[-grep("hypoxia", metadata$source_name),]
metadata <- metadata[-grep("mixed", metadata$Tissue),]
metadata <- metadata[-grep("new born", metadata$AGE),]
expmat <- expmat[,colnames(expmat) %in% rownames(metadata)]

metadata <- droplevels(metadata)

#Converting age to numeric values
metadata$AGE <- plyr::revalue(metadata$AGE, c("4 years"=4, "20 years"=20))

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
  
  boxplot(log2(normalized_matrices[[tissue]] + 1), notch=TRUE ,main ="Size-factor-normalized read counts
        Heterocephalus glaber (GSE30377)", col = "lightblue3", las = 2, cex.axis = .8)
  hist(log2(normalized_matrices[[tissue]] + 1), col = "lightblue3", main = "Distribution of size-factor-normalized counts
     Heterocephalus glaber (GSE30377)",
       xlab = "log2(NormalizedCounts + 1)")  
  
}

## Downloading the ENSEMBL99 dn/ds values between NMR and Guinea Pig

ensembl99 <- useEnsembl("ensembl", version = 99,
                        dataset = "hgfemale_gene_ensembl")

tissue_dnds <- lapply(normalized_matrices, function(tissue){
  
  dnds <- getBM(attributes=c("ensembl_gene_id", "cporcellus_homolog_ensembl_gene", 
                             "cporcellus_homolog_dn", "cporcellus_homolog_ds"),
                filters = "ensembl_gene_id", values = rownames(tissue),
                mart = ensembl99,
                useCache = FALSE)
  
})

##cleaning the dnds data

tissue_dnds <- lapply(tissue_dnds, function(dnds){
  
  dnds <- dnds[!(is.na(dnds$cporcellus_homolog_ds)),] #removing NA dS values
  dnds <- dnds[!(is.na(dnds$cporcellus_homolog_dn)),] #removing NA dN values
  dnds <- dnds[!(dnds$cporcellus_homolog_dn == 0),] #removing 0 dN values
  dnds <- dnds[!(dnds$cporcellus_homolog_ds == 0),] #removing 0 dS values
  dNdS <- (dnds$cporcellus_homolog_dn / dnds$cporcellus_homolog_ds) #calculating dn/ds
  dnds <- cbind(dnds, dNdS)#adding dn/ds vector to the dataframe 
  ###removing duplicate gene ids, leaves only 1to1 orth. genes
  dnds <- dnds[!((duplicated(dnds$ensembl_gene_id)) | (duplicated(dnds$ensembl_gene_id, fromLast = T))) ,]
  ###removing genes with dnds > 0.8
  dnds <- dnds[dnds$dNdS < 0.8,]
  
})

###--- Setting up functions ----

##Differential expression analysis

###EDITED TO REMOVE P CUTOFF AS THERE ARE ONLY TWO POINTS,
###
DEG_genes <- function(expMatrix, ageVector, rhoCutoff, dnds){
  Aging_exp_cor <- t(apply(expMatrix,1,function(x){
    test <- cor.test(x, ageVector, method="spearman", exact = F, use= "everything")
    c(test$est)
  }))
  Aging_exp_cor <- t(Aging_exp_cor)
  colnames(Aging_exp_cor) <- c("rho")
  
  Aging_exp_cor <- as.matrix(Aging_exp_cor[complete.cases(Aging_exp_cor),])
  colnames(Aging_exp_cor) <- c("rho")
  
  all_genes_cor <- Aging_exp_cor
  all_genes_exp <- expMatrix
  all_genes_dnds <- dnds
  
  DEG_genes_cor <- Aging_exp_cor[abs(Aging_exp_cor) > rhoCutoff,]
  DEG_genes_exp <- expMatrix[rownames(expMatrix) %in% names(DEG_genes_cor),]
  DEG_genes_dnds <- dnds[dnds$ensembl_gene_id %in% names(DEG_genes_cor),]
  
  inc_genes_cor <- DEG_genes_cor[DEG_genes_cor > rhoCutoff]
  inc_genes_exp <- expMatrix[rownames(expMatrix) %in% names(inc_genes_cor),]
  inc_genes_dnds <- dnds[dnds$ensembl_gene_id %in% names(inc_genes_cor),]
  
  dec_genes_cor <- DEG_genes_cor[DEG_genes_cor < -(rhoCutoff)]
  dec_genes_exp <- expMatrix[rownames(expMatrix) %in% names(dec_genes_cor),]
  dec_genes_dnds <- dnds[dnds$ensembl_gene_id %in% names(dec_genes_cor),]
  
  nonsig_genes_cor <- Aging_exp_cor[abs(Aging_exp_cor) < rhoCutoff]
  nonsig_genes_exp <- expMatrix[rownames(expMatrix) %in% names(nonsig_genes_cor),]
  nonsig_genes_dnds <- dnds[dnds$ensembl_gene_id %in% names(nonsig_genes_cor),]
  
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

Results_GSE30337 <- lapply(levels(metadata$Tissue), function(tissue){
  
  
  dnds_exp_list <- subsetCommon(matrix1 = tissue_dnds[[tissue]],
                                colname1 = "ensembl_gene_id", 
                                matrix2 = normalized_matrices[[tissue]]) 
  names(dnds_exp_list) <- c("dnds", "exp_matrix")
  ## Getting age vector 
  age <- as.numeric(as.character(tissue_exp[[tissue]][["metadata"]][["AGE"]]))
  
  deg_list <- DEG_genes(expMatrix = dnds_exp_list[["exp_matrix"]], ageVector = age,
                        rhoCutoff = 0.5, dnds = dnds_exp_list[["dnds"]])
  
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
names(Results_GSE30337) <- levels(metadata$Tissue)

save(Results_GSE30337, file = "naked_mole_rat/GSE30337/R/results/Results_GSE30337.RData")
