library(Biostrings)
library(magrittr)
library(optparse)


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input fasta file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output fasta file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Input fasta must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$out)){
  opt$out <- gsub(".fa", "_loi.fa", opt$file)
}


getLongestIsoform <- function(proteome_file = opt$file, output_file = opt$out){
  
  proteome <- readAAStringSet(proteome_file)
  headers <- names(proteome)
  
  gene_names <- sapply(headers, function(header){
    split_header <- strsplit(header, " ") %>% unlist()
    split_header[grepl("gene:", split_header)] %>% gsub("gene:", "", .)
  })
  
  transcript_names <- sapply(headers, function(header){
    split_header <- strsplit(header, " ") %>% unlist()
    split_header[grepl("transcript:", split_header)] %>% gsub("transcript:", "", .)
  })
  
  isoform_lengths <- proteome@ranges@width
  
  protein_df <- data.frame(gene_names, transcript_names, isoform_lengths)
  
  longestIsoformNames <- sapply(unique(gene_names), function(gene_name){
    
    subdf <- protein_df[protein_df$gene_names == gene_name,]
    
    #retriving the longest isoform name
    longestIsoform <- as.character(subdf$transcript_names[which.max(subdf$isoform_lengths)])
    
  })
  
  proteome_loI <- proteome[transcript_names %in% longestIsoformNames]
  names(proteome_loI) <- names(longestIsoformNames)
  
  writeXStringSet(proteome_loI, output_file)
}

getLongestIsoform()
