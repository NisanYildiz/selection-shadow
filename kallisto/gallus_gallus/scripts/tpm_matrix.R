library(tools)

count_matrix <- data.frame(
  transcript_id = read.delim("/mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus/quant/SRR7127899/abundance.tsv")$target_id)

for(path in list.files(path = "/mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus/quant")){
  
  tpm <- read.delim(paste0("/mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus/quant/",path, "/abundance.tsv"))$tpm
  srr_name <- path
  count_matrix[[srr_name]] <- tpm
  
}

count_matrix$transcript_id <- gsub("\\..*","",count_matrix$transcript_id)

write.csv(count_matrix, file = "/mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus/GSE114129_tpm.csv")
