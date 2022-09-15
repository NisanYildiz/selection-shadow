# Analysis of GSE114129

#Loading the libraries
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(plyr)
library("biomaRt")
library(grid)
library(dplyr)


#Loading up the expression matrix

expmat <- read.table(header = T, row.names = 1,
  file = "data/gallus_gallus/GSE114129/all_counts.txt")


#Loading the metadata 
metadata <- read.csv(file = "gallus_gallus/SraRunTable-GSE114129.txt",
                     row.names = 1)

##reordering the metadata 
metadata <- metadata[order(row.names(metadata)),]

##Checking if the order is correct
identical(row.names(metadata), colnames(expmat)) ###TRUE

##Removing embryonic samples
expmat <- expmat[, - grep("Embryonic", metadata$AGE) ]
metadata <- metadata[ - grep("Embryonic", metadata$AGE) ,]

#identical?
identical(row.names(metadata), colnames(expmat)) ###TRUE

##Converting age to numeric vector
metadata$AGE <- plyr::revalue(metadata$AGE, c("1 year old"=365, "100 days old"=100,
                        "3 years old"=1095, "300 days old" = 300,
                        "5 years old"= 1825))


##Normalization

### Creating the DESeq object
dds <- DESeqDataSetFromMatrix(countData = expmat, colData = metadata, design = ~AGE)

### Removing 0 count genes
nonzero_genes <-rowSums(counts(dds)) > 0
dds <- dds[ nonzero_genes, ]

### Size factor estimatiıon for DESeq normalization
dds <- estimateSizeFactors(dds)

### Normalizing
normalized_matrix <- counts(dds, normalized=TRUE)

## Basic exploratory visualizations 

boxplot(log2(counts(dds, normalize= TRUE)+1), notch=TRUE ,main = "Size-factor-normalized read counts
        Gallus gallus(GSE114129)", col = "lightblue3", las = 2, cex.axis = .8)
hist(log2(normalized_matrix + 1), col = "lightblue3", main = "Distribution of size-factor-normalized counts
     Gallus gallus(GSE113129)",
     xlab = "log2(NormalizedCounts + 1)")

pc_gallus = prcomp(t(normalized_matrix), scale = T)
percent = paste0(colnames(pc_gallus$x[,1:4])," (", round(summary(pc_gallus)$importance[2,1:4] * 100, 2), "%) ")

pca_gallus_1_2 <- data.frame(pc_gallus$x[,1:4], Age= factor(metadata$AGE, levels = c("100", "300", "365", "1095", "1825", "Embryonic day 12", "Embryonic day 16", "Embryonic day 20"))) %>% 
  ggplot(aes(PC1, PC2, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent[1]) + ylab(percent[2])
pca_gallus_3_4 <- data.frame(pc_gallus$x[,1:4], Age= factor(metadata$AGE, levels = c("100", "300", "365", "1095", "1825", "Embryonic day 12", "Embryonic day 16", "Embryonic day 20"))) %>% 
  ggplot(aes(PC3, PC4, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent[3]) + ylab(percent[4])

save(pca_gallus_1_2, pca_gallus_3_4, file = "gallus_gallus/R/results/pca_gallus.RData")

## Downloading the ENSEMBL99 dn/ds values between chicken and turkey

ensembl99 <- useEnsembl("ensembl", version = 99, 
                        dataset = "ggallus_gene_ensembl")
dnds <- getBM(attributes=c("ensembl_gene_id", "mgallopavo_homolog_ensembl_gene", 
                               "mgallopavo_homolog_dn", "mgallopavo_homolog_ds"),
                  filters = "ensembl_gene_id", values = rownames(normalized_matrix),
                  mart = ensembl99,
              useCache = F)

##cleaning the dnds data
dnds <- dnds[!(is.na(dnds$mgallopavo_homolog_ds)),] #removing NA dS values
dnds <- dnds[!(is.na(dnds$mgallopavo_homolog_dn)),] #removing NA dN values
dnds <- dnds[!(dnds$mgallopavo_homolog_dn == 0),] #removing 0 dN values
dnds <- dnds[!(dnds$mgallopavo_homolog_ds == 0),] #removing 0 dS values
dNdS <- (dnds$mgallopavo_homolog_dn / dnds$mgallopavo_homolog_ds) #calculating dn/ds
dnds <- cbind(dnds, dNdS)#adding dn/ds vector to the dataframe 
###removing duplicate gene ids, leaves only 1to1 orth. genes
dnds <- dnds[!((duplicated(dnds$ensembl_gene_id)) | (duplicated(dnds$ensembl_gene_id, fromLast = T))) ,]
###removing genes with dnds > 0.8
dnds <- dnds[dnds$dNdS < 0.8,]


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

dnds_exp_list <- dndsCommon(exp_matrix = normalized_matrix, dnds = dnds) 
  
## Getting age vector 
age <- as.numeric(as.character(metadata$AGE))
  
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
  
Results_GSE114129 <- list(age=age , all_results = all_results, sig_results = sig_results, inc_results = inc_results, 
           dec_results = dec_results, nonsig_results = nonsig_results, deg_list = deg_list,
           dnds_exp_list = dnds_exp_list)

Results_GSE114129 <- list(brain = Results_GSE114129)
save(Results_GSE114129, file = "gallus_gallus/R/results/Results_GSE114129.RData")


### Exp-conservation correlation vs. Age
p = cor.test(as.numeric(sig_results), age, method = "spearman", exact = F)$p.val
p = signif(p,2)
corr = cor(as.numeric(sig_results), age, method = "spearman")
corr = signif(corr,2)
plot(age, sig_results, main = "Differentially Expressed Genes in Aging
     Gallus gallus Brain(GSE114129)", xlab = "Age(days)",
     ylab = "Expression-Conservation rho", col = "lightblue3",
     pch = 19, cex = 1.6)
abline(lm(as.numeric(sig_results) ~ age),
       col = "#2B4C6F")
grid.text(x = unit(0.18, "npc"), y = unit(0.25, "npc"), just = c(1,1), 
          label = paste("p =",p))
grid.text(x = unit(0.18, "npc"), y = unit(0.22, "npc"), just = c(1,1),
          label = paste("corr =",corr))
grid.text(x = unit(0.20, "npc"), y = unit(0.19, "npc"), just = c(1,1),
          label = paste("# of genes =",length(deg_list[["All_sig"]][["sig_genes_exp"]][,1])))
