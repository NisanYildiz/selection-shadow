library(tools)

count_matrix <- data.frame(
  transcript_id = read.delim("/mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus/quant/TSRR5657253.fastq.gz/abundance.tsv")$target_id)

for(path in list.files(path = "/mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus/quant")){
  
  tpm <- read.delim(paste0("/mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus/quant/",path, "/abundance.tsv"))$tpm
  
  srr_name <- sub('.', '', path)
  srr_name <- sub(".fastq.gz", "", srr_name)
  
  count_matrix[[srr_name]] <- tpm
  
}

count_matrix$transcript_id <- gsub("\\..*","",count_matrix$transcript_id)

write.csv(count_matrix, file = "/mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus/GSE99791_tpm.csv")


