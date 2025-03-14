#Controlling for gene age and gene expression effect

#loading gene age data
musculus_gene_age <- read.csv("/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/Mus_musculus.csv")
musculus_gene_age$gene_age <- as.numeric(gsub(">", "", musculus_gene_age$gene_age))

#loading differential expression and dnds results
load("mus_musculus/GSE99791/R/results/Results_GSE99791.RData")

#loading tpm results
tpm <- read.csv("kallisto/mus_musculus/GSE99791_tpm.csv", row.names = 1)
rownames(tpm) <- tpm$transcript_id
tpm <- tpm[,-1]
#tpm <- as.matrix(tpm)
metadata <- read.csv(file = "mus_musculus/GSE99791/SraRunTable_GSE99791.txt",
                     row.names = 1)

input_samples <- c("GSM2652319", "GSM2652323", "GSM2652327", "GSM2652331", "GSM2652335",
                   "GSM2652339", "GSM2652343", "GSM2652347", "GSM2698018")
##reordering the metadata 
metadata <- metadata[order(row.names(metadata)),]

## removing input samples 
metadata <- metadata[!(metadata$GEO_Accession..exp. %in% input_samples),]
tpm <- tpm[,colnames(tpm) %in% rownames(metadata)]

##Checking if the order is correct
identical(row.names(metadata), colnames(tpm)) ###TRUE

#Converting age to numeric values
metadata$AGE <- plyr::revalue(metadata$AGE, c("adult (4 month old)"=4, "aged (2 year old)"=24))

#removing somatosensory cortex data, as it only has 4-month-old individuals 
tpm <- exp <- tpm[,!(metadata$Tissue == "somatosensory cortex")]
metadata <- metadata[!(metadata$Tissue == "somatosensory cortex"),]
metadata <- droplevels(metadata)
metadata$Tissue <- as.factor(metadata$Tissue)

ensembl99 <- useEnsembl("ensembl", version = 99, 
                        dataset = "mmusculus_gene_ensembl")
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
  colSums(subset_tpm[,(1:24)])
  
})

genewise_tpm <- t(genewise_tpm)
genewise_tpm <- as.matrix(genewise_tpm)
identical(row.names(metadata), colnames(genewise_tpm)) ###TRUE

#separating exp matrix and metadata into different tissues
tissue_exp_tpm <- lapply(levels(metadata$Tissue), function(tissue){
  exp <- genewise_tpm[,metadata$Tissue == tissue]
  meta <- metadata[metadata$Tissue == tissue,]
  a_list <- list(genewise_tpm = exp, metadata = meta)
  a_list
})
names(tissue_exp_tpm) <- levels(metadata$Tissue)


for (tissue in names(tissue_exp_tpm)){
  
  nonzero_transcripts <- rowSums(tissue_exp_tpm[[tissue]][["genewise_tpm"]]) > 0
  
  tissue_exp_tpm[[tissue]][["genewise_tpm"]] <- 
    tissue_exp_tpm[[tissue]][["genewise_tpm"]][nonzero_transcripts, ]
  
}


for( tissue in names(tissue_exp_tpm) ){
  
  boxplot(log2(tissue_exp_tpm[[tissue]][["tpm"]] + 1), notch=TRUE ,main = paste0("TPM values
        Mus musculus astrocyte from ", tissue, "(GSE99791)"), col = "lightblue3", las = 2, cex.axis = .8)
  hist(log2(tissue_exp_tpm[[tissue]][["tpm"]] + 1), col = "lightblue3", main = paste0("Distribution of TPM values
     Mus musculus astrocyte from ", tissue, "(GSE99791)"),
       xlab = "log2(NormalizedCounts + 1)")  
  
  
}


for (tissue in c("cerebellum", "hypothalamus")){

  x_tissue_exp_tpm <- tissue_exp_tpm[[tissue]][["genewise_tpm"]][
    rownames(tissue_exp_tpm[[tissue]][["genewise_tpm"]]) %in%
    Results_GSE99791[[tissue]][["dnds_exp_list"]][["dnds"]][["ensembl_gene_id"]],]
  
  x_tissue_exp_tpm <- x_tissue_exp_tpm[rownames(x_tissue_exp_tpm) %in% 
                                         musculus_gene_age$ensembl_gene_id,]
  
  
  mus_common_inc_dnds <- Results_GSE99791[[tissue]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["dNdS"]][
    (Results_GSE99791[[tissue]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id) & 
      (Results_GSE99791[[tissue]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_tissue_exp_tpm))]
  mus_common_inc_gene_names <- Results_GSE99791[[tissue]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]][
    (Results_GSE99791[[tissue]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id) &
      (Results_GSE99791[[tissue]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_tissue_exp_tpm))]
  mus_common_inc_gene_age <- musculus_gene_age[musculus_gene_age$ensembl_gene_id %in% mus_common_inc_gene_names,]
  mus_common_inc_gene_age <- mus_common_inc_gene_age[order(mus_common_inc_gene_age$ensembl_gene_id),]
  mus_common_inc_exp <- x_tissue_exp_tpm[rownames(x_tissue_exp_tpm) %in% mus_common_inc_gene_names,]
  mus_common_inc_exp <- mus_common_inc_exp[order(rownames(mus_common_inc_exp)),]
  mus_common_inc_max_exp <- apply(X= mus_common_inc_exp, MARGIN = 1, function(x){max(x)})
  mus_common_inc_mean_exp <- apply(X= mus_common_inc_exp, MARGIN = 1, function(x){mean(x)})

  
  mus_common_dec_dnds <- Results_GSE99791[[tissue]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["dNdS"]][
    (Results_GSE99791[[tissue]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id) & 
      (Results_GSE99791[[tissue]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_tissue_exp_tpm))]
  mus_common_dec_gene_names <- Results_GSE99791[[tissue]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]][
    (Results_GSE99791[[tissue]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% musculus_gene_age$ensembl_gene_id) &
      (Results_GSE99791[[tissue]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]] %in% rownames(x_tissue_exp_tpm))]
  mus_common_dec_gene_age <- musculus_gene_age[musculus_gene_age$ensembl_gene_id %in% mus_common_dec_gene_names,]
  mus_common_dec_gene_age <- mus_common_dec_gene_age[order(mus_common_dec_gene_age$ensembl_gene_id),]
  mus_common_dec_exp <- x_tissue_exp_tpm[rownames(x_tissue_exp_tpm) %in% mus_common_dec_gene_names,]
  mus_common_dec_exp <- mus_common_dec_exp[order(rownames(mus_common_dec_exp)),]
  mus_common_dec_max_exp <- apply(X= mus_common_dec_exp, MARGIN = 1, function(x){max(x)})
  mus_common_dec_mean_exp <- apply(X= mus_common_dec_exp, MARGIN = 1, function(x){mean(x)})
  
  if (identical(rownames(mus_common_inc_exp), mus_common_inc_gene_names) &
      identical(as.character(mus_common_dec_gene_age$ensembl_gene_id), mus_common_dec_gene_names)){

    mus_astrocyte_dnds_gene_age <- rbind(
      
      data.frame(gene_age = mus_common_inc_gene_age$gene_age ,
                 gene_class = "old-biased",
                 dnds = mus_common_inc_dnds,
                 max_exp = mus_common_inc_max_exp,
                 mean_exp = mus_common_inc_mean_exp),
      
      data.frame(gene_age = mus_common_dec_gene_age$gene_age ,
                 gene_class = "young-biased",
                 dnds = mus_common_dec_dnds,
                 max_exp = mus_common_dec_max_exp,
                 mean_exp = mus_common_dec_mean_exp)
    )
    
    print(paste("Tissue:", tissue))
    print(summary(lm(dnds ~ max_exp + gene_age + gene_class, data = mus_astrocyte_dnds_gene_age)))
    print(summary(lm(dnds ~ mean_exp + gene_age + gene_class, data = mus_astrocyte_dnds_gene_age)))
    
    
  } else{
    print("ERROR: Gene names are not identical")
  }
  
  
}
