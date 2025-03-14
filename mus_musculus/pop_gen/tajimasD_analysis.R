library(spatstat.utils)
library(tidyr)

#gene id, transcript id lookup table
transcript_geneid_table <- read.csv("/mnt/NEOGENE4/projects/melih_2020/mus_musculus/gene_transcript_name_table.csv", row.names = 1)
transcript_geneid_table$gene_name <- sapply(transcript_geneid_table$gene_name, function(x){
  strsplit(x, ".", fixed = TRUE)[[1]][1]
})
transcript_geneid_table$longest_transcript <- sapply(transcript_geneid_table$longest_transcript, function(x){
  strsplit(x, ".", fixed = TRUE)[[1]][1]
})


tajimas <- lapply(list.files("/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/", "*mus_musculus_tajimasD_cds_map.bed"),
       function(tajima_file){
         read.delim(file = paste0("/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/", tajima_file), sep = "\t", header = F)
         
       })
pop <- gsub("_mus_musculus_tajimasD_cds_map.bed", "",list.files("/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/", "*mus_musculus_tajimasD_cds_map.bed"))
names(tajimas) <- pop

tajimas_nongenic <- tajimas
for (x in pop){
  tajimas_nongenic[[x]] <- tajimas_nongenic[[x]][-(union((which(is.na(tajimas_nongenic[[x]]$V5)) - 1), which(is.na(tajimas_nongenic[[x]]$V5)))),]
}

tajimas_filled <- tajimas
for (x in pop){
  tajimas_filled[[x]] <- tajimas_filled[[x]] %>% tidyr::fill(V6)
  tajimas_filled[[x]] <- tajimas_filled[[x]][is.na(tajimas_filled[[x]]$V5),]
  tajimas_filled[[x]] <- tajimas_filled[[x]][, -5]
  tajimas_filled[[x]]$exon_length <- tajimas_filled[[x]]$V3 - tajimas_filled[[x]]$V2
}

tajimas_loc <- tajimas

for (x in pop){
  tajimas_loc[[x]]$V2[is.na(tajimas_loc[[x]]$V6)] <- NA
  tajimas_loc[[x]]$V3[is.na(tajimas_loc[[x]]$V6)] <- NA
  tajimas_loc[[x]] <- tajimas_loc[[x]] %>% tidyr::fill(V2)
  tajimas_loc[[x]] <- tajimas_loc[[x]] %>% tidyr::fill(V3)
  tajimas_loc[[x]] <- tajimas_loc[[x]][is.na(tajimas_loc[[x]]$V5),]
  
  
  if(identical(tajimas_loc[[x]]$V4, tajimas_filled[[x]]$V4)){
    
    tajimas_filled[[x]]$bin_start <- tajimas_loc[[x]]$V2
    tajimas_filled[[x]]$bin_end <- tajimas_loc[[x]]$V3
    
  }else{print("not identical")}
}


for (x in pop){
  
  tajimas_filled[[x]]$intersect_length <- sapply(1:nrow(tajimas_filled[[x]]), function(i){
    row = tajimas_filled[[x]][i,]
    intersect_range <- intersect.ranges(c(row$bin_start, row$bin_end), c(row$V2, row$V3), fatal = TRUE)
    intersect_range[2]- intersect_range[1]
  })
    
  
}

mean_tajimasD <- lapply(pop, function(x){
  mean_tajimasD <- sapply(unique(tajimas_filled[[x]]$V4), function(transcript){
    subset <- tajimas_filled[[x]][tajimas_filled[[x]]$V4 == transcript,]
    mean_tajimasD_per_transcript <- 
      weighted.mean(x = subset$V6, w = subset$intersect_length)
    mean_tajimasD_per_transcript
  })
  mean_tajimasD
})
names(mean_tajimasD) <- pop

for (x in pop){
  mean_tajimasD[[x]] <- data.frame(transcript_id = names(mean_tajimasD[[x]]), mean_tajimasD = mean_tajimasD[[x]])
  mean_tajimasD[[x]] <- mean_tajimasD[[x]][mean_tajimasD[[x]]$transcript_id %in% transcript_geneid_table$longest_transcript,]
  mean_tajimasD[[x]] <- mean_tajimasD[[x]][order(mean_tajimasD[[x]]$transcript_id),]
}

transcript_geneid_table_list <- lapply(pop, function(x){
  
  transcript_geneid_table <- transcript_geneid_table[
    transcript_geneid_table$longest_transcript %in% mean_tajimasD[[x]]$transcript_id,
  ]
  transcript_geneid_table <- transcript_geneid_table[
    order(transcript_geneid_table$longest_transcript), 
  ]
  
})

names(transcript_geneid_table_list) <- pop
all(sapply(pop, function(x){identical(as.character(transcript_geneid_table_list[[x]]$longest_transcript), as.character(mean_tajimasD[[x]]$transcript_id))}))

for (x in pop){
  
  mean_tajimasD[[x]]$gene_id <- transcript_geneid_table_list[[x]]$gene_name
  
}


#gettingg dnds values to check for correlation between tajimas D and dnds
ensembl99 <- useEnsembl("ensembl", version = 99, 
                        dataset = "mmusculus_gene_ensembl")
dnds <- getBM(attributes=c("ensembl_gene_id", "rnorvegicus_homolog_ensembl_gene", 
                           "rnorvegicus_homolog_dn", "rnorvegicus_homolog_ds"),
              filters = "ensembl_gene_id", values = transcript_geneid_table$gene_name,
              mart = ensembl99,
              useCache = FALSE)

##cleaning the dnds data
dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_ds)),] #removing NA dS values
dnds <- dnds[!(is.na(dnds$rnorvegicus_homolog_dn)),] #removing NA dN values
dnds <- dnds[!(dnds$rnorvegicus_homolog_dn == 0),] #removing 0 dN values
dnds <- dnds[!(dnds$rnorvegicus_homolog_ds == 0),] #removing 0 dS values
dNdS <- (dnds$rnorvegicus_homolog_dn / dnds$rnorvegicus_homolog_ds) #calculating dn/ds
dnds <- cbind(dnds, dNdS)#adding dn/ds vector to the dataframe 
###removing duplicate gene ids, leaves only 1to1 orth. genes
dnds <- dnds[!((duplicated(dnds$ensembl_gene_id)) | (duplicated(dnds$ensembl_gene_id, fromLast = T))) ,]
dnds <- dnds[!((duplicated(dnds$rnorvegicus_homolog_ensembl_gene)) | (duplicated(dnds$rnorvegicus_homolog_ensembl_gene, fromLast = T))) ,]

dnds <- dnds[order(dnds$ensembl_gene_id),]

dnds_mean_tajimasD_table_list <- lapply(pop, function(x){
  
  dnds_mean_tajimasD_table <- mean_tajimasD[[x]]
  dnds_mean_tajimasD_table <- dnds_mean_tajimasD_table[
    dnds_mean_tajimasD_table$gene_id %in% dnds$ensembl_gene_id,
  ]
  dnds_mean_tajimasD_table <- dnds_mean_tajimasD_table[
    order(dnds_mean_tajimasD_table$gene_id),
  ]
  
  dnds <- dnds[dnds$ensembl_gene_id %in% dnds_mean_tajimasD_table$gene_id,]
  
  if (identical(dnds$ensembl_gene_id, dnds_mean_tajimasD_table$gene_id)){
    
    dnds_mean_tajimasD_table$dNdS <- dnds$dNdS
  }
  dnds_mean_tajimasD_table
})


names(dnds_mean_tajimasD_table_list) <- pop


#cor.test(dnds_mean_tajimasD_table$mean_tajimasD, 
##    dnds_mean_tajimasD_table$dNdS, method = "spearman",
##    exact = F)


#DEG resultts for mus musculus astrocyte data
load("/mnt/NEOGENE4/projects/melih_2020/mus_musculus/GSE99791/R/results/Results_GSE99791.RData")

for (tissue in c("cerebellum", "hypothalamus")){
  
  #DEG genes for hypothalamus
  inc_genes <- Results_GSE99791[[tissue]][["deg_list"]][["inc_sig"]][["inc_genes_dnds"]][["ensembl_gene_id"]]
  dec_genes <- Results_GSE99791[[tissue]][["deg_list"]][["dec_sig"]][["dec_genes_dnds"]][["ensembl_gene_id"]]
  non_biased <- Results_GSE99791[[tissue]][["deg_list"]][["Non_sig"]][["nonsig_genes_dnds"]][["ensembl_gene_id"]]
  
  inc_genes_transcriptid <- sapply(inc_genes, function(x){
    transcript_id <- transcript_geneid_table$longest_transcript[
      transcript_geneid_table$gene_name == x]
  })
  
  dec_genes_transcriptid <- sapply(dec_genes, function(x){
    transcript_id <- transcript_geneid_table$longest_transcript[
      transcript_geneid_table$gene_name == x]
    
  })
  
  non_biased_transcriptid <- sapply(non_biased, function(x){
    transcript_id <- transcript_geneid_table$longest_transcript[
      transcript_geneid_table$gene_name == x]
    
  })
  
  for (p in pop){
    
    #subsetting tajimas score table to only include genes in DEG dataset
    mean_tajimasD_subset <- mean_tajimasD[[p]][mean_tajimasD[[p]]$transcript_id %in% c(inc_genes_transcriptid, dec_genes_transcriptid, non_biased_transcriptid),]
    
    
    #adding gene class info
    mean_tajimasD_subset$gene_class <- ""
    mean_tajimasD_subset$gene_class[mean_tajimasD_subset$transcript_id %in% inc_genes_transcriptid] <- "old-biased"
    mean_tajimasD_subset$gene_class[mean_tajimasD_subset$transcript_id %in% dec_genes_transcriptid] <- "young-biased"
    mean_tajimasD_subset$gene_class[mean_tajimasD_subset$transcript_id %in% non_biased_transcriptid] <- "non-biased"
    
    
    file_name <- paste0("/mnt/NEOGENE4/projects/melih_2020/mus_musculus/mean_tajimasD_scores_per_transcript_GSE99791_",tissue,"_",p,".csv")
    write.csv(mean_tajimasD_subset, file = file_name)
    
    
    
  }

}
