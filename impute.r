library(Seurat)

raw_tbl <- read.table(file = 'GrubmanData/scRNA_rawCounts.tsv', sep = '\t', header = TRUE)
rownames(raw_tbl) = raw_tbl[,1]
raw_tbl[,1] <- NULL

meta_tbl <- read.table(file = 'GrubmanData/scRNA_metadata.csv', sep = ',', header = TRUE)
rownames(meta_tbl) <- meta_tbl[,1]
meta_tbl[,1] <- NULL

tbl <- read.table(file = 'GrubmanData/impute.csv', sep = ',', header = TRUE)
rownames(tbl) = tbl[,1]
tbl[,1] <- NULL

mydata <- CreateSeuratObject(tbl)
mydata_norm <- NormalizeData(mydata)

mydata_norm <- ScaleData(mydata_norm)
mydata_norm <- FindVariableFeatures(mydata_norm)

mydata_norm <- RunPCA(mydata_norm, features=VariableFeatures(mydata_norm))
mydata_norm <- FindNeighbors(mydata_norm, dims = 1:10)
mydata_norm <- FindClusters(mydata_norm)

mydata_norm <- RunUMAP(mydata_norm, dims=1:10)
DimPlot(mydata_norm, reduction = "umap")

#data_to_write_out <- as.data.frame(as.matrix(mydata_norm@assays$RNA@data))
#fwrite(x = data_to_write_out, file = "outfile.csv")

#res <- VIPER(tbl[,1:500], num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
#             report = FALSE, outdir = NULL, prefix = NULL)
# gene.variances <- apply(res$imputed,1,var)

gene.variances <- apply(as.matrix(mydata_norm@assays$RNA@counts),1,var)

d <- density(mydata_norm@assays$RNA@data["LRP1B",])
plot(d)
d1 = density(mydata_norm@assays$RNA@data["LRP1B",][Idents(mydata_norm)[names(mydata_norm@assays$RNA@data["LRP1B",])] == 1])


s1 = Idents(mydata_norm)[names(mydata_norm@assays$RNA@data["LRP1B",])] == 2
s2 = as.vector(raw_tbl["LRP1B",][names(mydata_norm@assays$RNA@data["LRP1B",])] != 0)
s3 = meta_tbl[names(mydata_norm@assays$RNA@data["LRP1B",]),2] == "AD5"
d2 = density(mydata_norm@assays$RNA@data["LRP1B",][s1 & s2 & s3])
plot(d2)


get.variance <- function(gene, patient, cell_type) {
  s1 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),6] == cell_type
  #s2 <- as.vector(raw_tbl[gene,][names(mydata_norm@assays$RNA@data[gene,])] != 0)
  s3 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),2] == patient
  count.list <- mydata_norm@assays$RNA@data[gene,][s1 & s3]
  print(var(count.list))
  plot(density(count.list))
  return(count.list)
}

cell_pop_vector = c()
for (gene in rownames(tbl)[1:30])
{
  for (patient in c("AD1", "AD2", "AD3", "AD4", "AD5", "AD6", "AD-un", "Ct1", "Ct2", "Ct3", "Ct4", "Ct5", "Ct6", "Ct-un"))
  {
    for (cell_type in c("astro", "doublet", "endo", "mg", "neuron", "oligo", "OPC", "unID"))
    {
      s1 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),6] == cell_type
      s2 <- as.vector(raw_tbl[gene,][names(mydata_norm@assays$RNA@data[gene,])] != 0)
      s3 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),2] == patient
      cell_pop_vector <- c(cell_pop_vector, as.vector(mydata_norm@assays$RNA@data[gene,][s1 & s2 & s3]))
      print(names(cell_pop_vector))
      print()
      names(cell_pop_vector) <- c(names(cell_pop_vector), paste(gene,patient,cell_type,sep="_"))
    }
  }
}
