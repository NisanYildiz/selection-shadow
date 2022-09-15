# Analysis of GSE66712

#Loading the libraries
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(plyr)
library("biomaRt")
library(grid)
library(openxlsx)

#Loading up the expression matrix
expmat <- read.delim(row.names = 1,
  file = "data/nothobranchius_furzeri/GSE66712/GSE66712_nfu_counts_185_samples.txt")


gene_info <- expmat[,1:6]
expmat <- expmat[,-(1:6)]

#Downloading the metadata
Sys.setenv(VROOM_CONNECTION_SIZE=500072)
geodata <- getGEO('GSE66712',GSEMatrix=TRUE)

## ordering geodata into metadata as expected by the DESeq package
## excluding the tailfin tissue
metadata <- data.frame(geodata[["GSE66712-GPL17944_series_matrix.txt.gz"]]@phenoData@data[["description"]],
                  geodata[["GSE66712-GPL17944_series_matrix.txt.gz"]]@phenoData@data[["age at death:ch1"]],
                  geodata[["GSE66712-GPL17944_series_matrix.txt.gz"]]@phenoData@data[["tissue:ch1"]],
                  geodata[["GSE66712-GPL17944_series_matrix.txt.gz"]]@phenoData@data[["treatment:ch1"]],
                  row.names = 1)
colnames(metadata) <- c("AGE", "Tissue", "Treatment")

## keeping only non-chemical-treated individuals 
metadata <- metadata[metadata$Treatment == "none",]

## limiting expmat to non-treated individuals and liver/skin tissues 
expmat <- expmat[,colnames(expmat) %in% rownames(metadata)]

## re-ordering the metadata and the expmat
metadata <- metadata[order(row.names(metadata)),]
expmat <- expmat[,order(colnames(expmat))]

##Checking if the order is correct
identical(row.names(metadata), colnames(expmat)) ###TRUE

## converting age to numeric values
metadata$AGE <- plyr::revalue(metadata$AGE, c("12 weeks" = 12, "20 weeks" = 20,
                                              "27 weeks" = 27, "39 weeks" = 39,
                                              "5 weeks" = 5, "9 weeks" = 9))

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
        Nothobranchius furzeri ", tissue, "(GSE66712)"), col = "lightblue3", las = 2, cex.axis = .8)
  hist(log2(normalized_matrices[[tissue]] + 1), col = "lightblue3", main = paste0("Distribution of size-factor-normalized counts
     Nothobranchius furzeri ", tissue, "(GSE66712)"),
       xlab = "log2(NormalizedCounts + 1)")  
  
}

pc_killifish_liver = prcomp(t(normalized_matrices[["liver"]]), scale = T)
percent_liver = paste0(colnames(pc_killifish_liver$x[,1:4])," (", round(summary(pc_killifish_liver)$importance[2,1:4] * 100, 2), "%) ")

pc_killifish_skin = prcomp(t(normalized_matrices[["skin"]]), scale = T)
percent_skin = paste0(colnames(pc_killifish_skin$x[,1:4])," (", round(summary(pc_killifish_skin)$importance[2,1:4] * 100, 2), "%) ")

pca_killifish_liver_1_2 <- data.frame(pc_killifish_liver$x[,1:4], Age= factor(tissue_exp[["liver"]][["metadata"]][["AGE"]], levels = c("5", "12", "20", "27", "39"))) %>% 
  ggplot(aes(PC1, PC2, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent_liver[1]) + ylab(percent_liver[2])
pca_killifish_liver_3_4 <- data.frame(pc_killifish_liver$x[,1:4], Age= factor(tissue_exp[["liver"]][["metadata"]][["AGE"]], levels = c("5", "12", "20", "27", "39"))) %>% 
  ggplot(aes(PC3, PC4, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent_liver[3]) + ylab(percent_liver[4])

pca_killifish_skin_1_2 <- data.frame(pc_killifish_skin$x[,1:4], Age= factor(tissue_exp[["skin"]][["metadata"]][["AGE"]], levels = c("5", "12", "20", "27", "39"))) %>% 
  ggplot(aes(PC1, PC2, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent_skin[1]) + ylab(percent_skin[2])
pca_killifish_skin_3_4 <- data.frame(pc_killifish_skin$x[,1:4], Age= factor(tissue_exp[["skin"]][["metadata"]][["AGE"]], levels = c("5", "12", "20", "27", "39"))) %>% 
  ggplot(aes(PC3, PC4, colour=Age)) + geom_point(size = 4) + theme_bw() + xlab(percent_skin[3]) + ylab(percent_skin[4])

save(pca_killifish_liver_1_2, pca_killifish_liver_3_4, pca_killifish_skin_1_2, pca_killifish_skin_3_4,
     file = "n_furzeri/GSE66712/R/results/pca_killifish.RData")

## loading the dn/ds values 

dnds <- read.xlsx("n_furzeri/acel12577-sup-0001-tabless1-13.xlsx",
                  sheet = 3, startRow = 4)

###removing duplicate gene ids, leaves only 1to1 orth. genes
dnds <- dnds[!((duplicated(dnds$Gene.ID)) | (duplicated(dnds$Gene.ID, fromLast = T))) ,]
###removing genes with dnds > 0.8
dnds <- dnds[dnds$`Ka/Ks.tested.branch` < 0.8,]
###removing genes with dnds = 0 
dnds <- dnds[dnds$`Ka/Ks.tested.branch` > 0,]
colnames(dnds)[9] <- "dNdS"


###--- Setting up functions ----


#Subsets the genes that are present both in dnds data.frame and matrix
#outputs the subset expression matrix and dnds data.frame as a list
#Both have the same order of genes
dndsCommon <- function(exp_matrix, dnds){
  exp_matrix_subset <- exp_matrix[rownames(exp_matrix) %in% dnds$Gene.ID,]
  dnds_subset <- dnds[dnds$Gene.ID %in% rownames(exp_matrix_subset),]
  #ordering
  dnds_subset <- dnds_subset[(order(dnds_subset$Gene.ID)),]
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
  DEG_genes_dnds <- dnds[dnds$Gene.ID %in% rownames(DEG_genes_cor),]
  
  inc_genes_cor <- DEG_genes_cor[DEG_genes_cor[,1] > rhoCutoff ,]
  inc_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(inc_genes_cor),]
  inc_genes_dnds <- dnds[dnds$Gene.ID %in% rownames(inc_genes_cor),]
  
  dec_genes_cor <- DEG_genes_cor[DEG_genes_cor [,1] < -(rhoCutoff) ,]
  dec_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(dec_genes_cor),]
  dec_genes_dnds <- dnds[dnds$Gene.ID %in% rownames(dec_genes_cor),]
  
  nonsig_genes_cor <- Aging_exp_cor[Aging_exp_cor[,"adjusted_p"] > pAdjCutoff,]
  nonsig_genes_exp <- expMatrix[rownames(expMatrix) %in% rownames(nonsig_genes_cor),]
  nonsig_genes_dnds <- dnds[dnds$Gene.ID %in% rownames(nonsig_genes_cor),]
  
  listy <- list(Aging_exp_cor = Aging_exp_cor, 
                All = list(all_genes_cor = all_genes_cor, all_genes_exp = all_genes_exp, all_genes_dnds = all_genes_dnds),
                Non_sig = list(nonsig_genes_cor = nonsig_genes_cor, nonsig_genes_exp = nonsig_genes_exp, nonsig_genes_dnds = nonsig_genes_dnds), 
                All_sig = list(sig_genes_cor = DEG_genes_cor, sig_genes_exp = DEG_genes_exp, sig_genes_dnds = DEG_genes_dnds), 
                inc_sig = list(inc_genes_cor = inc_genes_cor, inc_genes_exp = inc_genes_exp, inc_genes_dnds = inc_genes_dnds), 
                dec_sig = list(dec_genes_cor = dec_genes_cor, dec_genes_exp = dec_genes_exp, dec_genes_dnds = dec_genes_dnds))
  listy
}

###---- DONE ---- 

Results_GSE66712 <- lapply(levels(metadata$Tissue), function(tissue){
  
  
  dnds_exp_list <- dndsCommon(exp_matrix = normalized_matrices[[tissue]], 
                              dnds = dnds) 
  
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
names(Results_GSE66712) <- levels(metadata$Tissue)

save(Results_GSE66712, file = "n_furzeri/GSE66712/R/results/Results_GSE66712.RData")

### Exp-conservation correlation vs. Age


for (tissue in names(Results_GSE66712)){
  
  p = cor.test(as.numeric(Results_GSE66712[[tissue]][["sig_results"]]), Results_GSE66712[[tissue]][["age"]], method = "spearman", exact = F)$p.val
  p = signif(p,2)
  corr = cor(as.numeric(Results_GSE66712[[tissue]][["sig_results"]]), Results_GSE66712[[tissue]][["age"]], method = "spearman")
  corr = signif(corr,2)
  plot(Results_GSE66712[[tissue]][["age"]], Results_GSE66712[[tissue]][["sig_results"]], main = paste0("Differentially Expressed Genes in Aging
     Nothobranchius furzeri ",tissue, "(GSE66712)"), xlab = "Age(days)",
       ylab = "Expression-Conservation rho", col = "lightblue3",
       pch = 19, cex = 1.6)
  abline(lm(as.numeric(Results_GSE66712[[tissue]][["sig_results"]]) ~ Results_GSE66712[[tissue]][["age"]]),
         col = "#2B4C6F")
  grid.text(x = unit(0.18, "npc"), y = unit(0.25, "npc"), just = c(1,1), 
            label = paste("p =",p))
  grid.text(x = unit(0.18, "npc"), y = unit(0.22, "npc"), just = c(1,1),
            label = paste("corr =",corr))
  grid.text(x = unit(0.20, "npc"), y = unit(0.19, "npc"), just = c(1,1),
            label = paste("# of genes =",length(Results_GSE66712[[tissue]][["deg_list"]][["All_sig"]][["sig_genes_exp"]][,1])))
  
  
}


