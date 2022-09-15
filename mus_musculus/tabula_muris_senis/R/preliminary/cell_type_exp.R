myPaths <- .libPaths()

myPaths <- c(myPaths, "/home/myildiz/RLibs/")

.libPaths(myPaths)  # add new path

library(Seurat, lib.loc = "/home/myildiz/RLibs/")


#library(tidyverse)
# data downloaded from: https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728
# metadata downloaded from: 
# https://s3.console.aws.amazon.com/s3/buckets/czb-tabula-muris-senis/Metadata/?region=us-west-2&tab=overview
lung = ReadH5AD("mus_musculus/tabula_muris_senis/R/data_files/tabula-muris-senis-facs-processed-official-annotations-Lung.h5ad")
liver = ReadH5AD("mus_musculus/tabula_muris_senis/R/data_files/tabula-muris-senis-facs-processed-official-annotations-Liver.h5ad")
muscle = 
  ReadH5AD("mus_musculus/tabula_muris_senis/R/data_files/tabula-muris-senis-facs-processed-official-annotations-Limb_Muscle.h5ad")
brain = 
  ReadH5AD("mus_musculus/tabula_muris_senis/R/data_files/tabula-muris-senis-facs-processed-official-annotations-Brain_Non-Myeloid.h5ad")
skin = ReadH5AD("mus_musculus/tabula_muris_senis/R/data_files/tabula-muris-senis-facs-processed-official-annotations-Skin.h5ad")
kidney = ReadH5AD("mus_musculus/tabula_muris_senis/R/data_files/tabula-muris-senis-facs-processed-official-annotations-Kidney.h5ad")


tis.exp = list(lung = Seurat::GetAssayData(lung, slot="data"),
               liver = Seurat::GetAssayData(liver, slot="data"),
               muscle = Seurat::GetAssayData(muscle, slot="data"),
               brain = Seurat::GetAssayData(brain, slot="data"),
               skin = Seurat::GetAssayData(skin, slot="data"),
               kidney = Seurat::GetAssayData(kidney, slot="data"))

metadata = read.csv("mus_musculus/tabula_muris_senis/R/data_files/metadata/tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids.csv")
metadata = list(lung = metadata[match(rownames(lung@meta.data),metadata$obs_names),],
                liver = metadata[match(rownames(liver@meta.data),metadata$obs_names),],
                muscle = metadata[match(rownames(muscle@meta.data),metadata$obs_names),],
                brain = metadata[match(rownames(brain@meta.data),metadata$obs_names),],
                skin = metadata[match(rownames(skin@meta.data),metadata$obs_names),],
                kidney = metadata[match(rownames(kidney@meta.data),metadata$obs_names),])

rm(lung,liver,muscle,brain, skin, kidney)
#########

for(x in 1:6){
  metadata[[x]]$mouse.id = as.character(metadata[[x]]$mouse.id)
}
for(x in 1:6){
  metadata[[x]]$obs_names = as.character(metadata[[x]]$obs_names)
}
for(x in 1:6){
  metadata[[x]]$cell_ontology_class = as.character(metadata[[x]]$cell_ontology_class)
}

# remove cell types that are identified in less than 15 cells in all time points. 
# apply only for lung and skin, other organs have high number of cells for each cell type:
rowSums(table(metadata$brain$cell_ontology_class,metadata$brain$mouse.id)) # min 17 cells
rowSums(table(metadata$muscle$cell_ontology_class,metadata$muscle$mouse.id)) # min 136 cells
rowSums(table(metadata$liver$cell_ontology_class,metadata$liver$mouse.id)) # min 30 cells
rowSums(table(metadata$lung$cell_ontology_class,metadata$lung$mouse.id))# min 4
rowSums(table(metadata$kidney$cell_ontology_class,metadata$kidney$mouse.id)) # min 16 cells
rowSums(table(metadata$skin$cell_ontology_class,metadata$skin$mouse.id)) # min 13 cells

rm.cc = rowSums(table(metadata$lung$cell_ontology_class,metadata$lung$mouse.id))
rm.cc = names(rm.cc[rm.cc<15]) # remove cell types with < 15 cells
metadata$lung = metadata$lung[!metadata$lung$cell_ontology_class%in%rm.cc,]
tis.exp$lung = tis.exp$lung[,colnames(tis.exp$lung)%in%metadata$lung$obs_names]

rm.cc = rowSums(table(metadata$skin$cell_ontology_class,metadata$skin$mouse.id))
rm.cc = names(rm.cc[rm.cc<15]) # remove cell types with < 15 cells
metadata$skin = metadata$skin[!metadata$skin$cell_ontology_class%in%rm.cc,]
tis.exp$skin = tis.exp$skin[,colnames(tis.exp$skin)%in%metadata$skin$obs_names]



# expression matrix columns and metadata cell rows are in same order:
for(i in 1:6) print(identical(metadata[[i]]$obs_names, colnames(tis.exp[[i]])))

# separate ages into different lists:
tis.exp.3m = sapply(names(tis.exp),function(x){
  tis.exp[[x]][, metadata[[x]]$age == "3m"]
})
tis.exp.18m = sapply(names(tis.exp),function(x){
  tis.exp[[x]][, metadata[[x]]$age == "18m"]
})
tis.exp.24m = sapply(names(tis.exp),function(x){
  tis.exp[[x]][, metadata[[x]]$age == "24m"]
})

# remove genes not expressed in any of the cells:
tis.exp.3m = sapply(names(tis.exp.3m), function(x){
  tis.exp.3m[[x]][!apply(tis.exp.3m[[x]]==0, 1, function(y) sum(y) == ncol(tis.exp.3m[[x]])), ]
})
tis.exp.18m = sapply(names(tis.exp.18m), function(x){
  tis.exp.18m[[x]][!apply(tis.exp.18m[[x]]==0, 1, function(y) sum(y) == ncol(tis.exp.18m[[x]])), ]
})
tis.exp.24m = sapply(names(tis.exp.24m), function(x){
  tis.exp.24m[[x]][!apply(tis.exp.24m[[x]]==0, 1, function(y) sum(y) == ncol(tis.exp.24m[[x]])), ]
})

# separate metadata into different lists:
metadata.3m = sapply(names(metadata), function(x){
  metadata[[x]][metadata[[x]]$age =="3m",]
},simplify = F)
metadata.18m = sapply(names(metadata), function(x){
  metadata[[x]][metadata[[x]]$age =="18m",]
},simplify = F)
metadata.24m = sapply(names(metadata), function(x){
  metadata[[x]][metadata[[x]]$age =="24m",]
},simplify = F)

#rm(metadata,tis.exp)

# expression matrix columns and metadata cell rows are in same order:
for(i in 1:6) print(identical(metadata.3m[[i]]$obs_names, colnames(tis.exp.3m[[i]])))
for(i in 1:6) print(identical(metadata.18m[[i]]$obs_names, colnames(tis.exp.18m[[i]])))
for(i in 1:6) print(identical(metadata.24m[[i]]$obs_names, colnames(tis.exp.24m[[i]])))


cell.type.exp.3m = sapply(names(tis.exp.3m),function(a){
  cell.class = unique(metadata.3m[[a]]$cell_ontology_class)#unique cell types in tissue "a"
  mouse.id = unique(metadata.3m[[a]]$mouse.id)#unique mouse id's for tissue "a" 
  cell.class.exp = sapply(1:length(cell.class), function(x){
    f = sapply(1:length(mouse.id), function(y){
      cell.names = metadata.3m[[a]][metadata.3m[[a]]$cell_ontology_class == cell.class[x] &
                                      metadata.3m[[a]]$mouse.id == mouse.id[y], 1]
      if(length(cell.names)!=0){
        rowMeans(as.matrix(tis.exp.3m[[a]][,cell.names]))
      }
    },simplify = F)
    f = f[!sapply(f,is.null)]
    rowMeans(do.call(cbind,f))
  })
  colnames(cell.class.exp) = cell.class
  return(cell.class.exp)
})

cell.type.exp.18m = sapply(names(tis.exp.18m),function(a){
  cell.class = unique(metadata.18m[[a]]$cell_ontology_class)
  mouse.id = unique(metadata.18m[[a]]$mouse.id)
  cell.class.exp = sapply(1:length(cell.class), function(x){
    f = sapply(1:length(mouse.id), function(y){
      cell.names = metadata.18m[[a]][metadata.18m[[a]]$cell_ontology_class == cell.class[x] &
                                       metadata.18m[[a]]$mouse.id == mouse.id[y], 1]
      if(length(cell.names)!=0){
        rowMeans(as.matrix(tis.exp.18m[[a]][,cell.names]))
      }
    },simplify = F)
    f = f[!sapply(f,is.null)]
    rowMeans(do.call(cbind,f))
  })
  colnames(cell.class.exp) = cell.class
  return(cell.class.exp)
})

cell.type.exp.24m = sapply(names(tis.exp.24m),function(a){
  cell.class = unique(metadata.24m[[a]]$cell_ontology_class)
  mouse.id = unique(metadata.24m[[a]]$mouse.id)
  cell.class.exp = sapply(1:length(cell.class), function(x){
    f = sapply(1:length(mouse.id), function(y){
      cell.names = metadata.24m[[a]][metadata.24m[[a]]$cell_ontology_class == cell.class[x] &
                                       metadata.24m[[a]]$mouse.id == mouse.id[y], 1]
      if(length(cell.names)!=0){
        rowMeans(as.matrix(tis.exp.24m[[a]][,cell.names]))
      }
    },simplify = F)
    f = f[!sapply(f,is.null)]
    rowMeans(do.call(cbind,f))
  })
  colnames(cell.class.exp) = cell.class
  return(cell.class.exp)
})

saveRDS(cell.type.exp.3m,file="mus_musculus/tabula_muris_senis/R/results/3m.celltype.rds")
saveRDS(cell.type.exp.18m,file="mus_musculus/tabula_muris_senis/R/results/18m.celltype.rds")
saveRDS(cell.type.exp.24m,file="mus_musculus/tabula_muris_senis/R/results/24m.celltype.rds")

#### calculate cell type expression per individual:

cell.type.ind.exp.3m = sapply(names(tis.exp.3m),function(a){
  cell.class = unique(metadata.3m[[a]]$cell_ontology_class)#unique cell types within tissue "a"
  mouse.id = unique(metadata.3m[[a]]$mouse.id)#unique mouse id's within tissue "a"
  cell.class.exp = sapply(1:length(cell.class), function(x){
    f = sapply(1:length(mouse.id), function(y){
      cell.names = metadata.3m[[a]][metadata.3m[[a]]$cell_ontology_class == cell.class[x] &
                                      metadata.3m[[a]]$mouse.id == mouse.id[y], 1]
      if(length(cell.names)!=0){
        rowMeans(as.matrix(tis.exp.3m[[a]][,cell.names]))
      }
    },simplify = F)
    rm_null_ind = sapply(f,is.null)
    f = f[!rm_null_ind]
    indX.exp = do.call(cbind,f)
    colnames(indX.exp) = mouse.id[!rm_null_ind]
    return(indX.exp)
  },simplify = F)
  names(cell.class.exp) = cell.class
  return(cell.class.exp)
})

cell.type.ind.exp.18m = sapply(names(tis.exp.18m),function(a){
  cell.class = unique(metadata.18m[[a]]$cell_ontology_class)
  mouse.id = unique(metadata.18m[[a]]$mouse.id)
  cell.class.exp = sapply(1:length(cell.class), function(x){
    f = sapply(1:length(mouse.id), function(y){
      cell.names = metadata.18m[[a]][metadata.18m[[a]]$cell_ontology_class == cell.class[x] &
                                       metadata.18m[[a]]$mouse.id == mouse.id[y], 1]
      if(length(cell.names)!=0){
        rowMeans(as.matrix(tis.exp.18m[[a]][,cell.names]))
      }
    },simplify = F)
    rm_null_ind = sapply(f,is.null)
    f = f[!rm_null_ind]
    indX.exp = do.call(cbind,f)
    colnames(indX.exp) = mouse.id[!rm_null_ind]
    return(indX.exp)
  },simplify = F)
  names(cell.class.exp) = cell.class
  return(cell.class.exp)
})

cell.type.ind.exp.24m = sapply(names(tis.exp.24m),function(a){
  cell.class = unique(metadata.24m[[a]]$cell_ontology_class)
  mouse.id = unique(metadata.24m[[a]]$mouse.id)
  cell.class.exp = sapply(1:length(cell.class), function(x){
    f = sapply(1:length(mouse.id), function(y){
      cell.names = metadata.24m[[a]][metadata.24m[[a]]$cell_ontology_class == cell.class[x] &
                                       metadata.24m[[a]]$mouse.id == mouse.id[y], 1]
      if(length(cell.names)!=0){
        rowMeans(as.matrix(tis.exp.24m[[a]][,cell.names]))
      }
    },simplify = F)
    rm_null_ind = sapply(f,is.null)
    f = f[!rm_null_ind]
    indX.exp = do.call(cbind,f)
    colnames(indX.exp) = mouse.id[!rm_null_ind]
    return(indX.exp)
  },simplify = F)
  names(cell.class.exp) = cell.class
  return(cell.class.exp)
})

saveRDS(cell.type.ind.exp.3m, file="mus_musculus/tabula_muris_senis/R/results/3m.celltype.per.ind.rds")
saveRDS(cell.type.ind.exp.18m, file="mus_musculus/tabula_muris_senis/R/results/18m.celltype.per.ind.rds")
saveRDS(cell.type.ind.exp.24m, file="mus_musculus/tabula_muris_senis/R/results/24m.celltype.per.ind.rds")
saveRDS(metadata, file="mus_musculus/tabula_muris_senis/R/results/metadata.rds")
