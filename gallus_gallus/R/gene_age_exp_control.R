#Controlling for gene age and gene expression effect

#loading gene age data
gallus_gene_age <- read.csv("/mnt/NEOGENE4/projects/melih_2020/data/gallus_gallus/Gallus_gallus.csv")
gallus_gene_age$gene_age <- as.numeric(gsub(">", "", gallus_gene_age$gene_age))

#loading differential expression and dnds results
load("gallus_gallus/R/results/Results_GSE114129.RData")

#loading tpm results
tpm <- read.csv("kallisto/gallus_gallus/GSE114129_tpm.csv", row.names = 1)
rownames(tpm) <- tpm$transcript_id
tpm <- tpm[,-1]

#Loading the metadata 
metadata <- read.csv(file = "gallus_gallus/SraRunTable-GSE114129.txt",
                     row.names = 1)
##reordering the metadata 
metadata <- metadata[order(row.names(metadata)),]

##Checking if the order is correct
identical(row.names(metadata), colnames(tpm)) ###TRUE

##Removing embryonic samples
tpm <- tpm[, - grep("Embryonic", metadata$AGE) ]
metadata <- metadata[ - grep("Embryonic", metadata$AGE) ,]

#identical?
identical(row.names(metadata), colnames(tpm)) ###TRUE

##Converting age to numeric vector
metadata$AGE <- plyr::revalue(metadata$AGE, c("1 year old"=365, "100 days old"=100,
                                              "3 years old"=1095, "300 days old" = 300,
                                              "5 years old"= 1825))

ensembl99 <- useEnsembl("ensembl", version = 99, 
                        dataset = "ggallus_gene_ensembl")
transcript_gene_ids <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"),
                             filters = "ensembl_transcript_id", values = rownames(tpm),
                             mart = ensembl99,
                             useCache = FALSE)
transcript_gene_ids <- 
  transcript_gene_ids[order(transcript_gene_ids$ensembl_transcript_id),]
tpm <- tpm[order(rownames(tpm)),]

tpm$gene_ids <- transcript_gene_ids$ensembl_gene_id

#converting tpm values to genewise values by summing up different transcripts for genes
genewise_tpm <- sapply(unique(tpm$gene_ids), function(gene_id){
  
  subset_tpm <- tpm[(tpm$gene_ids == gene_id),]
  colSums(subset_tpm[,(1:(length(tpm)-1))])
  
})

genewise_tpm <- t(genewise_tpm)
genewise_tpm <- as.matrix(genewise_tpm)
identical(row.names(metadata), colnames(genewise_tpm)) ###TRUE

#removing all zero genes
nonzero_genes <- rowSums(genewise_tpm) > 0
genewise_tpm <- genewise_tpm[nonzero_genes,]

#basic plots

boxplot(log2(genewise_tpm + 1), notch=TRUE ,main = paste0("TPM values
        for Gallus gallus brain (GSE114129)"), col = "lightblue3", las = 2, cex.axis = .8)
hist(log2(genewise_tpm + 1), col = "lightblue3", main = paste0("Distribution of TPM values
     for Gallus gallus brain (GSE114129)"),
     xlab = "log2(NormalizedCounts + 1)")  

x_genewise_tpm <- genewise_tpm[
    rownames(genewise_tpm) %in%
    Results_GSE114129[["brain"]][["dnds_exp_list"]][["dnds"]][["ensembl_gene_id"]],]
  
x_genewise_tpm <- x_genewise_tpm[rownames(x_genewise_tpm) %in% 
                                         gallus_gene_age$ensembl_gene_id,]
  
  
gallus_common_inc_dnds <- Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]][
  (Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id) & 
    (Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_genewise_tpm))]
gallus_common_inc_gene_names <- Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]][
  (Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id) &
    (Results_GSE114129[["brain"]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_genewise_tpm))]
gallus_common_inc_gene_age <- gallus_gene_age[gallus_gene_age$ensembl_gene_id %in% gallus_common_inc_gene_names,]
gallus_common_inc_gene_age <- gallus_common_inc_gene_age[order(gallus_common_inc_gene_age$ensembl_gene_id),]
gallus_common_inc_exp <- x_genewise_tpm[rownames(x_genewise_tpm) %in% gallus_common_inc_gene_names,]
gallus_common_inc_exp <- gallus_common_inc_exp[order(rownames(gallus_common_inc_exp)),]
gallus_common_inc_max_exp <- apply(X= gallus_common_inc_exp, MARGIN = 1, function(x){max(x)})
gallus_common_inc_mean_exp <- apply(X= gallus_common_inc_exp, MARGIN = 1, function(x){mean(x)})

gallus_common_dec_dnds <- Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]][
  (Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id) & 
    (Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_genewise_tpm))]
gallus_common_dec_gene_names <- Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]][
  (Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% gallus_gene_age$ensembl_gene_id) &
    (Results_GSE114129[["brain"]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_genewise_tpm))]
gallus_common_dec_gene_age <- gallus_gene_age[gallus_gene_age$ensembl_gene_id %in% gallus_common_dec_gene_names,]
gallus_common_dec_gene_age <- gallus_common_dec_gene_age[order(gallus_common_dec_gene_age$ensembl_gene_id),]
gallus_common_dec_exp <- x_genewise_tpm[rownames(x_genewise_tpm) %in% gallus_common_dec_gene_names,]
gallus_common_dec_exp <- gallus_common_dec_exp[order(rownames(gallus_common_dec_exp)),]
gallus_common_dec_max_exp <- apply(X= gallus_common_dec_exp, MARGIN = 1, function(x){max(x)})
gallus_common_dec_mean_exp <- apply(X= gallus_common_dec_exp, MARGIN = 1, function(x){mean(x)})

if (identical(rownames(gallus_common_inc_exp), gallus_common_inc_gene_names) &
    identical(as.character(gallus_common_dec_gene_age$ensembl_gene_id), gallus_common_dec_gene_names)){
    gallus_astrocyte_dnds_gene_age <- rbind(
    
    data.frame(gene_age = gallus_common_inc_gene_age$gene_age ,
               gene_class = "old-biased",
               dnds = gallus_common_inc_dnds,
               max_exp = gallus_common_inc_max_exp,
               mean_exp = gallus_common_inc_mean_exp),
    
    data.frame(gene_age = gallus_common_dec_gene_age$gene_age ,
               gene_class = "young-biased",
               dnds = gallus_common_dec_dnds,
               max_exp = gallus_common_dec_max_exp,
               mean_exp = gallus_common_dec_mean_exp)
  )
    
  print(summary(lm(dnds ~ max_exp + gene_age + gene_class, data = gallus_astrocyte_dnds_gene_age)))
  print(summary(lm(dnds ~ mean_exp + gene_age + gene_class, data = gallus_astrocyte_dnds_gene_age)))
    
    
} else{
  print("ERROR: Gene names are not identical")
}
  

