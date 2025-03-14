# Analysis of pacifico et al. data


#Loading the libraries
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(plyr)
library("biomaRt")
library(grid)
library(openxlsx)

#Loading up the expression matrix
expmat <- read.delim(file = "data/drosophila_melanogaster/pone.0209405.s009.txt")


##removing duplicate genes
expmat <- expmat[!((duplicated(expmat$Gene)) | (duplicated(expmat$Gene, fromLast = T))) ,]

##gene name conversion 
ensembl_metazoa <- useMart("metazoa_mart", host = "sep2019-metazoa.ensembl.org", 
                           dataset = "dmelanogaster_eg_gene")
gene_ids <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), 
                  filters = 'external_gene_name', values = expmat$Gene, mart = ensembl_metazoa,
                  useCache = F)
gene_ids <- gene_ids[!(gene_ids$external_gene_name %in% gene_ids$external_gene_name[which(duplicated(gene_ids$external_gene_name))]),]#getting rid of duplicated genes

row.names(expmat) <- expmat$Gene

##Converts gene names in rownames to ENSEMBL gene ids and outputs the matrix
IDconv <- function(expMatrix, gene_ids){
  exp_matrix <- expMatrix[rownames(expMatrix) %in% gene_ids$external_gene_name,]
  gene_names <- gene_ids[gene_ids$external_gene_name %in% rownames(exp_matrix),]
  exp_matrix <- exp_matrix[order(rownames(exp_matrix)),]
  gene_names <- gene_names[order(gene_names$external_gene_name),]
  
  rownames(exp_matrix) <- gene_names$ensembl_gene_id
  
  exp_matrix
}
expmat <- IDconv(expmat, gene_ids)

##removing the gene name column
expmat <- expmat[,-1]
expmat <- as.matrix(expmat)

# Defining the metadata 

metadata <- t(sapply(colnames(expmat), function(x){
  
  AGE <- strsplit(x, "_")[[1]][[2]]
  Sex <- strsplit(x, "_")[[1]][[1]]
  c(AGE,Sex)
  
}))

metadata <- as.data.frame(metadata)
colnames(metadata) <- c("AGE", "Sex")

##Converting age to numeric vector
metadata$AGE <- plyr::revalue(metadata$AGE, c("20d" = 20, "30d" = 30,
                                              "40d" = 40, "5d" = 5))

#identical order?
identical(row.names(metadata), colnames(expmat)) ###TRUE

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
        Drosophila melanogaster(Pacifico et al.)", col = "lightblue3", las = 2, cex.axis = .8)
hist(log2(normalized_matrix + 1), col = "lightblue3", main = "Distribution of size-factor-normalized counts
     Drosophila melanogaster(Pacifico et al.)",
     xlab = "log2(NormalizedCounts + 1)")

pc_dros = prcomp(t(normalized_matrix), scale = T)
percent = paste0(colnames(pc_dros$x[,1:4])," (", round(summary(pc_dros)$importance[2,1:4] * 100, 2), "%) ")

pca_dros_1_2 <- data.frame(pc_dros$x[,1:4], Age= factor(metadata$AGE, levels = c("5", "20", "30", "40")), Sex = metadata$Sex) %>% 
  ggplot(aes(PC1, PC2, colour=Age, shape = Sex)) + geom_point(size = 4) + theme_bw() + xlab(percent[1]) + ylab(percent[2])
pca_dros_3_4 <- data.frame(pc_dros$x[,1:4], Age= factor(metadata$AGE, levels = c("5", "20", "30", "40")), Sex = metadata$Sex) %>% 
  ggplot(aes(PC3, PC4, colour=Age, shape = Sex)) + geom_point(size = 4) + theme_bw() + xlab(percent[3]) + ylab(percent[4])

save(pca_dros_1_2, pca_dros_3_4, file = "drosophila_melanogaster/R/results/pca_dros.RData")

# Downloading the ENSEMBL Genomes 45 dn/ds values between d melanogaster and 
# d simulans
dnds <- getBM(attributes=c("ensembl_gene_id", "dsimulans_eg_homolog_ensembl_gene", 
                           "dsimulans_eg_homolog_dn", "dsimulans_eg_homolog_ds"),
              filters = "ensembl_gene_id", values = rownames(normalized_matrix),
              mart = ensembl_metazoa,
              useCache = F)

##cleaning the dnds data
dnds <- dnds[!(is.na(dnds$dsimulans_eg_homolog_ds)),] #removing NA dS values
dnds <- dnds[!(is.na(dnds$dsimulans_eg_homolog_dn)),] #removing NA dN values
dnds <- dnds[!(dnds$dsimulans_eg_homolog_dn == 0),] #removing 0 dN values
dnds <- dnds[!(dnds$dsimulans_eg_homolog_ds == 0),] #removing 0 dS values
dNdS <- (dnds$dsimulans_eg_homolog_dn / dnds$dsimulans_eg_homolog_ds) #calculating dn/ds
dnds <- cbind(dnds, dNdS)#adding dn/ds vector to the dataframe 
###removing duplicate gene ids, leaves only 1to1 orth. genes
dnds <- dnds[!((duplicated(dnds$ensembl_gene_id)) | (duplicated(dnds$ensembl_gene_id, fromLast = T))) ,]
dnds <- dnds[!((duplicated(dnds$dsimulans_eg_homolog_ensembl_gene)) | (duplicated(dnds$dsimulans_eg_homolog_ensembl_gene, fromLast = T))) ,]
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

Results_pacifico <- list(age=age , all_results = all_results, sig_results = sig_results, inc_results = inc_results, 
                          dec_results = dec_results, nonsig_results = nonsig_results, deg_list = deg_list,
                          dnds_exp_list = dnds_exp_list)

Results_pacifico <- list(head = Results_pacifico)

save(Results_pacifico, file = "drosophila_melanogaster/R/results/Results_pacifico.RData")

### Exp-conservation correlation vs. Age
p = cor.test(as.numeric(sig_results), age, method = "spearman", exact = F)$p.val
p = signif(p,2)
corr = cor(as.numeric(sig_results), age, method = "spearman")
corr = signif(corr,2)
plot(age, sig_results, main = "Differentially Expressed Genes in Aging
     Drosophila melanogaster head(Pacifico et al)", xlab = "Age(days)",
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

