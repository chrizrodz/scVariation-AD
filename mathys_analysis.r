library(Seurat)
library(Matrix)
library(enrichR)
library(ggplot2)
library(reshape2)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015" , "ChEA_2016" ,"KEGG_2016")
patient_num = 1
tbl <- read.csv(paste("MathysData/patient",patient_num,".csv",sep=""), row.names = 1)
gene.names <- rownames(tbl)

# # LOAD GRUBMAN DATA
# raw_tbl <- read.table(file = 'GrubmanData/scRNA_rawCounts.tsv', sep = '\t', header = TRUE)
# rownames(raw_tbl) = raw_tbl[,1]
# raw_tbl[,1] <- NULL
# 
# meta_tbl <- read.table(file = 'GrubmanData/scRNA_metadata.csv', sep = ',', header = TRUE)
# rownames(meta_tbl) <- meta_tbl[,1]
# meta_tbl[,1] <- NULL
# 
# tbl <- read.table(file = 'GrubmanData/DeepImpute.csv', sep = ',', header = TRUE)
# rownames(tbl) = tbl[,1]
# tbl[,1] <- NULL

# LOAD MATHYS DATA
master_tbl <- read.csv("MathysData/patient10_impute.csv", row.names = 1)
master_sparse <- Matrix(as.matrix(master_tbl), sparse = TRUE)


master_tbl_not_impute <- read.csv("MathysData/patient10.csv", row.names = 1)
master_sparse_not_impute <- Matrix(as.matrix(master_tbl_not_impute, sparse = TRUE))


file_list <- list.files("MathysData/")
file_list <- file_list[grepl("patient",file_list) & grepl("impute", file_list)]

for (file_name in file_list[2:length(file_list)])
{
  full_path <- paste("MathysData/",file_name,sep="")
  tbl <- read.csv(full_path, row.names = 1)
  tbl_sparse <- Matrix(as.matrix(tbl), sparse = TRUE)
  #writeMM(tbl_sparse, paste0("MathysData/",file_name,".mtx"))
  master_sparse <- merge(master_sparse, tbl_sparse, by=0, all=TRUE)
}


cell_annos <- read.table("MathysData/filtered_column_metadata.txt", header = TRUE, row.names = 1, sep="\t")
gene.names = c()
for (patient_num in c(1:48))
{
  tbl <- read.csv(paste("MathysData/patient",patient_num,".csv",sep=""), row.names = 1)
  
  gene.names <- rownames(tbl)
  print(gene.names)
  for (cell_type in c("Oli","Per","Ex","In","Ast","End","Opc", "Mic"))
  {
    cell_names <- getType(cell_type, patient_num)
    tbl_sparse <- tbl_sparse <- Matrix(as.matrix(tbl[,cell_names]), sparse = TRUE)
    writeMM(tbl_sparse, paste0("MathysData/patient",patient_num,"_",cell_type,".mtx"))
  }
}


# WRITE PATIENT CELL =STATS
patient_stats = data.frame(patient_num=NA, cell_count=NA, Oli=NA, Per=NA, Ex=NA, In=NA, Ast=NA, End=NA, Opc=NA, Mic=NA)
for (patient_num in c(1:48))
{
  print(patient_num)
  counts = c(patient_num, length(which(sapply(strsplit(rownames(cell_annos),"\\."),`[[`, 2)==patient_num)))
  for (cell_type in c("Oli","Per","Ex","In","Ast","End","Opc", "Mic"))
  {
    counts = c(counts, length(getType(cell_type, patient_num)))
  }
  patient_stats = rbind(patient_stats, counts)
}
patient_stats = patient_stats[2:nrow(patient_stats),]
write.csv(patient_stats, "PatientCellStats.csv", quote=FALSE, row.names = FALSE)

getType <- function(type, patient)
{
  return(rownames(cell_annos)[grepl(type,cell_annos$broad.cell.type) & sapply(strsplit(rownames(cell_annos),"\\."),`[[`, 2)==patient])
}

getExpressionVals <- function(type, patient, gene)
{
  cell_names <- getType(type, patient)
  return(as.numeric(master_tbl[gene,cell_names]))
}


# VARIANCE ANALYSIS
gene.order.list = list()
cell_type="Oli"
for (cell_type in c("Ast", "Mic", "Ex","In","Oli",))
{
  df.list <- list()
  print("INITIALIZE")
  for (patient_num in c(1:48))
  {
    print(patient_num)
    #master_tbl <- read.csv(paste("MathysData/patient",patient_num,"_impute.csv",sep=""), row.names = 1)
    master_sparse <- readMM(paste("MathysData/patient",patient_num,"_",cell_type,"_impute.mtx",sep=""))
    master_tbl <- data.frame(master_sparse)
    rownames(master_tbl) <- gene.names
    scatter.df = data.frame(gene=character(), avg_exp=numeric(), var_exp=numeric(), coef_var=numeric(), stringsAsFactors = FALSE)
    for (gene in rownames(master_tbl))
    {
      express_vals = as.numeric(master_tbl[gene,])
      if (length(express_vals) == 0)
      {
        next()
      }
      avg_exp = mean(express_vals)
      var_exp = var(express_vals)
      coef_var = var_exp/avg_exp
      if (avg_exp > 1)
      {
       scatter.df[nrow(scatter.df) + 1,] = list(gene, avg_exp, var_exp, coef_var)       
      }
    }
    coef_order = rank(scatter.df$coef_var)
    scatter.df["coef_order"] = coef_order
    df.list = append(df.list, list(scatter.df))
    #scat_plot <- ggplot(scatter.df, aes(x=as.numeric(var_exp), y=as.numeric(avg_exp))) + geom_point() +labs(title="Patient 10 Excitatory Neurons", x ="Expression Variance", y = "Expression Mean")
    #scat_plot <- scat_plot +theme(axis.title=element_text(size=14,face="bold"), plot.title = element_text(size=16,face="bold", hjust=0.5))
  }
  
  patient.info <- read.csv("MathysData/PatientInfo.csv")
  id_maps <- read.csv("MathysData/id_mapping.csv")
  unique_ids <- unique(id_maps$Subject)
  patient.info <- patient.info.orig[,]
  
  count <- 1
  for (unique_val in unique(id_maps$Subject))
  {
    patient.info[patient.info$Subject==unique_val,"patient.num"] <- count
    count = count + 1
  }
  patient.info <- patient.info[order(patient.info$patient.num),]
  ctrl_pat_num_list <- patient.info[patient.info$pathologic.diagnosis.of.AD=="NO","patient.num"]
  ad_pat_num_list <- patient.info[patient.info$pathologic.diagnosis.of.AD=="YES","patient.num"]
    
  ## CONTROL PATIENTS
  print("CONTROL")
  gene.base = list()
  for (patient_num in ctrl_pat_num_list)
  {
    print(patient_num)
    stat.df <- df.list[[patient_num]]
    for (gene in stat.df$gene)
    {
      if (!gene %in% names(gene.base))
      {
        gene.base[[gene]] = c(stat.df[stat.df$gene == gene,"coef_var"])
      }
      else
      {
        gene.base[[gene]] = c(gene.base[[gene]], stat.df[stat.df$gene == gene,"coef_var"])
      }
    }
  }
  for (entry in names(gene.base))
  {
    gene.base[[entry]] = mean(gene.base[[entry]])
  }
    
  ## AD PATIENTS
  print("AD")
  gene.order <- list()
  for (patient_num in ad_pat_num_list)
  {
    print(patient_num)
    stat.df <- df.list[[patient_num]]
    if (nrow(stat.df) == 0)
    {
      next()
    }
    stat.df$diff_coef = NA
    for (gene in stat.df$gene)
    {
      if (!is.null(gene.base[[gene]]))
      {
        stat.df[stat.df$gene == gene, "diff_coef"] = stat.df[stat.df$gene == gene,"coef_var"] / gene.base[[gene]]
      }
      else
      {
        stat.df[stat.df$gene == gene, "diff_coef"] = 1
      }
    }
    diff_coef_order = rank(stat.df$diff_coef)
    stat.df["diff_coef_order"] = diff_coef_order
    
    for (gene in stat.df$gene)
    {
      if (!gene %in% names(gene.order))
      {
        gene.order[[gene]] = stat.df[stat.df$gene == gene,"diff_coef_order"]
      }
      else
      {
        gene.order[[gene]] = gene.order[[gene]] + stat.df[stat.df$gene == gene,"diff_coef_order"]
      }
    }
  }
  gene.order <- gene.order[order(unlist(gene.order), decreasing = TRUE)]
  gene.order.list[[cell_type]] = gene.order
  gene.order <- gene.order.list[[cell_type]]
  enriched <- enrichr(names(gene.order)[1:100], dbs)
  enrich.df <- enriched$GO_Biological_Process_2015
  #enrich.df <- enrich.df[order(enrich.df$Combined.Score, decreasing = TRUE),]
  enrich.df$Term <- sapply(strsplit(as.character(enrich.df$Term), "\\(GO"), `[`, 1)
  enrich.df$LogP <- -log(enrich.df$Adjusted.P.value)
  ggplot(data=enrich.df[1:5,], aes(x=reorder(Term,-LogP), y=LogP)) + geom_bar(stat="identity", fill="steelblue") + theme_minimal() + 
    theme(axis.text.x = element_text(angle = 20, hjust = 1), axis.text = element_text(size=12)) + labs(title=paste0(cell_type," Pathway Enrichment"), x="GO Term", y="-Log P Value") + 
    theme(axis.title=element_text(size=14,face="bold"), plot.margin = margin(1, 1, 1, 6, "cm"), plot.title = element_text(size=16,face="bold", hjust=0.5), axis.text.x = element_text(size=10))
  
}


#CROSS-CORRELATION ANALYSIS
for (cell_type in c("Ast", "Mic", "Ex","In","Oli",))
{
  patient_num <- 28
  gene.list <- names(gene.order.list[[cell_type]][1:10])
  gene.df <- as.data.frame(Matrix(,length(gene.list),length(gene.list)))
  rownames(gene.df) = gene.list
  colnames(gene.df) = gene.list
  master_sparse <- readMM(paste("MathysData/patient",patient_num,"_",cell_type,"_impute.mtx",sep=""))
  master_tbl <- data.frame(master_sparse)
  rownames(master_tbl) <- gene.names
  
  for (gene1 in gene.list)
  {
    for (gene2 in gene.list)
    {
      gene1.exp <- rank(as.numeric(master_tbl[gene1,]))
      gene2.exp <- rank(as.numeric(master_tbl[gene2,]))
      print(gene1.exp)
      corr <- cor.test(x=gene1.exp, y=gene2.exp, method = 'spearman')
      gene.df[gene1,gene2] = corr$estimate
    }
  }
  gene.mtx <- as.matrix(gene.df)
  gene.mtx[lower.tri(gene.mtx)]<-NA
  melted_cormat <- melt(gene.mtx, na.rm = TRUE)
  
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman's\nRank Correlation") +
    theme_minimal()+ # minimal theme
    labs(title=paste0(cell_type, " Rank Correlation"))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    theme(axis.text.y = element_text(vjust = 1, 
                                     size = 12, hjust = 1))+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size=16,face="bold", hjust=0.5))+
    coord_fixed()
  print(ggheatmap)
}








# Transcriptional Noise Analysis
expression.vector <- c(4,0,0,1,5,2)
names(expression.vector) <- c("A","B","C","D","E","F")
new.vector <- sample(c(1:length(expression.vector)), 4, replace=TRUE, prob=expression.vector)
for (i in 1:length(expression.vector)){
  expression.vector[[i]] <- sum(new.vector==i)
}

# Unlist the data to get an array that represents the bootstrap row.
site_sample <- unlist(site_sample)

new.expression.vector <- rep(0, length(expression.vector))
names(new.expression.vector) <- names(expression.vector)
for(gene.name in new.vector)
{
  new.expression.vector[[gene.name]] <- new.expression.vector[[gene.name]] + 1
}


#Start
patient_splits_impute <- list()
for (cell_type in c("Ast", "Mic","Ex","In"))
{
  patient_num <- 1
  master_sparse <- as.data.frame(readMM(paste("MathysData/patient", patient_num, "_",cell_type,"_impute.mtx",sep="")))
  rownames(master_sparse) <- gene.names
  colnames(master_sparse) <- paste0(rep(paste0("patient", patient_num, "_"),ncol(master_sparse)), 1:ncol(master_sparse))
  master_sparse <- master_sparse[,colSums(master_sparse)>500]
  cell_counts = c(ncol(master_sparse))
  for (patient_num in c(2:48))
  {
    print(patient_num)
    partial_sparse <- as.data.frame(readMM(paste("MathysData/patient", patient_num, "_",cell_type,"_impute.mtx",sep="")))
    partial_sparse <- partial_sparse[,colSums(partial_sparse)>500]
    if(is.null(ncol(partial_sparse)) || ncol(partial_sparse) == 0)
    {
      cell_counts <- c(cell_counts, 0)
      next
    }
    rownames(partial_sparse) <- gene.names
    colnames(partial_sparse) <- paste0(rep(paste0("patient", patient_num, "_"),ncol(partial_sparse)), 1:ncol(partial_sparse))
    cell_counts <- c(cell_counts, ncol(partial_sparse))
    print(ncol(partial_sparse))
    master_sparse <- merge(master_sparse, partial_sparse, by=0, all=TRUE)
    rownames(master_sparse) <- master_sparse$Row.names
    master_sparse <- master_sparse[,colnames(master_sparse) != "Row.names"]
  }
  
  
  min_cell_count <- min(cell_counts)
  min_expresssion_count <- min(colSums(master_sparse))
  
  for (patient_num in c(1:48))
  {
    master_tbl <- master_sparse[, startsWith(colnames(master_sparse),paste0("patient",patient_num))]
    for (cell.name in colnames(master_tbl))
    {
      expression.vector <- master_tbl[, cell.name]
      names(expression.vector) <- rownames(master_tbl)
      
      new.vector <- sample(c(1:length(expression.vector)), min_expresssion_count, replace=TRUE, prob=expression.vector)
      for (i in 1:length(expression.vector)){
        expression.vector[[i]] <- sum(new.vector==i)
      }
      
      master_sparse[, cell.name] <- expression.vector
    }
  }
  
  average.expression.vals <- rowMeans(master_sparse)
  average.expression.vals <- average.expression.vals[order(average.expression.vals, decreasing = TRUE)]
  middle_genes <- average.expression.vals[round(length(average.expression.vals)/10):round(length(average.expression.vals)-length(average.expression.vals)/10)]
  gene_split <- split(c(1:length(average.expression.vals)), sort(c(1:length(average.expression.vals))%%10))
  
  all.lowest.vars <- c()
  for (chunk in gene_split[2:length(gene_split)-1])
  {
    chunk_names <- names(average.expression.vals[chunk])
    coeff_var <- apply(master_sparse[chunk_names,], 1, sd)/rowMeans(master_sparse[chunk_names,])
    coeff_var <- coeff_var[order(coeff_var)]
    lowest_vars <- names(coeff_var)[1:round(length(coeff_var)/10)]
    all.lowest.vars <- c(all.lowest.vars, lowest_vars)
  }
  
  distance.list <- c()
  for (patient_num in c(1:48))
  {
    vector.list = list()
    average.vector = NA
    master_tbl <- sqrt(master_sparse[all.lowest.vars, startsWith(colnames(master_sparse),paste0("patient",patient_num))])
    if(is.null(ncol(master_tbl)) || ncol(master_tbl) == 0)
    {
      next
    }
    
    for(cell.num in 1:ncol(master_tbl))
    {
      vector.list[[cell.num]] <- as.vector(master_tbl[,cell.num])
      if (is.na(average.vector))
      {
        average.vector <- vector.list[[cell.num]]
      }
      else
      {
        average.vector <- average.vector + vector.list[[cell.num]]
      }
    }
    average.vector <- average.vector / ncol(master_tbl)
    print(average.vector)
    
    distances <- c()
    for (expression.vector in vector.list)
    {
      distances <- c(distances, dist(rbind(expression.vector, average.vector)))
    }
    distance.list[[paste0(patient_num,"-",cell_type)]] <- distances
  }
  
  
  patient.info.orig <- read.csv("MathysData/PatientInfo.csv")
  mean.list <- list()
  for (patient_num in c(1:48))
  {
    if(paste0(patient_num,"-",cell_type) %in% names(distance.list))
    {
      mean.list[[patient_num]] <- mean(distance.list[[paste0(patient_num,"-",cell_type)]])
    }
    else
    {
      mean.list[[patient_num]] <- NA
    }
  }
  print(cell_type)
  patient.info <- data.frame(patient.info.orig)
  patient.info$mean.var <- unlist(mean.list)
  patient_splits_impute[[cell_type]] <- patient.info[,]
  patient.info[cell_counts < 100, "mean.var"] = NA
  mean.list <- mean.list[!is.na(unlist(mean.list))]
  control_pop <- unlist(mean.list[1:24])
  ad_pop <- unlist(mean.list[25:48])
  if(length(control_pop) > 10 && length(ad_pop) > 10)
  {
    print(patient.info[!is.na(patient.info$mean.var),])
    p <- ggplot(patient.info[!is.na(patient.info$mean.var),], aes(x=pathologic.diagnosis.of.AD, y=mean.var, fill=pathologic.diagnosis.of.AD)) + geom_boxplot()
    print(p + scale_fill_manual(values=c("green","red"), name ="AD Diagnosis") + ggtitle(paste0(cell_type, " Transcriptional Noise")) + theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + ylab("Transcriptional Noise"))
    print(wilcox.test(x=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="NO","mean.var"]),
                      y=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="YES","mean.var"])))
  }
  save(patient_splits_impute, file="PatientNoiseData.Robj")
}



id_maps <- read.csv("MathysData/id_mapping.csv")
count <- 1
for (unique_val in unique(id_maps$Subject))
{
  patient.info[patient.info$Subject==unique_val,"patient.num"] <- count
  count = count + 1
}
patient.info <- patient.info[order(patient.info$patient.num),]

#LOAD DATA
load("PatientNoiseData.Robj")
for (cell_type in c("Oli","Ex","In","Ast","Opc", "Mic"))
{ 
  print(cell_type)
  patient.info <- patient_splits_impute[[cell_type]]
  print(patient.info[!is.na(patient.info$mean.var),])
  p <- ggplot(patient.info[!is.na(patient.info$mean.var),], aes(x=pathologic.diagnosis.of.AD, y=mean.var, fill=pathologic.diagnosis.of.AD)) + geom_boxplot()
  print(p + scale_fill_manual(values=c("green","red"), name ="AD Diagnosis") + ggtitle(paste0(cell_type, " Transcriptional Noise")) + theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + ylab("Transcriptional Noise"))
  print(wilcox.test(x=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="NO","mean.var"]),
                    y=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="YES","mean.var"])))
  
}



cell_type<-"Oli"
patient.info <- patient.info.orig[,]

count <- 1
for (unique_val in unique(id_maps$Subject))
{
  patient.info[patient.info$Subject==unique_val,"patient.num"] <- count
  count = count + 1
}
patient.info <- patient.info[order(patient.info$patient.num),]


#Oli
patient.info$mean.var <- c(10.414369,10.37273,10.362002,10.353396,9.734913,9.65723,9.61368,9.794884,9.768886,10.333224,10.523849,10.265127,10.337802,10.319725,10.313016,10.479568,10.189642,9.97422,10.222888,10.186376,10.370086,10.284193,10.392342,9.800375,10.078981,10.309585,10.369473,10.429319,10.478098,10.301264,10.462311,10.429068,10.411581,10.317995,10.101163,10.352719,10.195298,10.27636,10.256137,10.480797,10.136771,10.269222,10.429921,10.34631,10.3734,10.148304,10.315197,10.026855)
#In
patient.info$mean.var <- c(9.887978,9.910291,9.831943,9.97853,9.409386,9.286926,9.430934,9.59957,9.276963,9.927943,9.789092,9.64435,10.048589,9.865778,9.744501,9.930171,9.721217,10.026842,9.715261,9.833851,9.733851,10.1066,9.807982,9.563519,9.966459,9.823121,9.972719,9.728389,9.971399,9.46272,9.713612,9.919344,9.951561,9.924049,9.65428,9.859406,9.991742,10.126491,9.924184,10.098696,10.124375,10.061179,10.108612,10.103695,10.000568,10.040891,9.623821,9.584909)
#Ex
patient.info$mean.var <- c(9.734753,9.790219,9.67923,9.851276,9.099306,8.983766,8.990822,9.53634,9.086163,9.983817,9.590989,10.030465,9.734722,9.630982,9.458109,9.619216,9.386366,9.796523,9.571528,9.558022,9.64234,9.850922,9.382599,9.954678,9.904992,9.485747,10.025012,9.818956,9.758328,9.379309,9.425271,9.814263,9.736471,9.570549,9.918582,9.519908,9.73802,9.995929,9.654494,9.872005,10.030268,9.865651,9.927187,9.966441,9.755466,9.927242,9.699504,9.274082)
#Ast
patient.info$mean.var <- c(10.063862,10.147207,10.044828,9.986668,9.333961,9.385011,9.381227,9.578641,9.432678,9.744638,9.997193,10.026893,9.931712,9.750002,9.992337,9.765264,9.825713,9.933952,9.946177,10.031397,9.971886,10.218745,9.973554,9.6312,9.91783,9.943193,9.91047,10.027891,9.95496,9.569062,9.818381,10.056383,9.806774,9.876686,9.885494,9.889306,9.648344,10.099108,9.919887,10.226693,10.012761,9.921898,9.754396,9.91178,9.633009,9.908445,9.730456,9.462766)
#Mic
patient.info$mean.var <- c(9.754228,9.688024,9.829346,9.850825,9.51675,9.538507,9.546673,8.365324,9.194671,9.791699,8.172603,9.033636,9.32932,9.809502,9.833026,9.385552,9.52165,9.558805,9.381362,9.633723,9.570744,9.101186,9.134323,9.257766,9.643646,5.445457,9.857531,9.144464,9.192892,9.389986,9.583078,9.5729,9.702523,7.856503,9.270712,8.420753,9.625732,9.513657,9.68098,9.711305,9.613784,9.836485,9.412941,9.39439,9.442274,8.355736,9.730456,9.462766)



patient.info[patient.info$msex==0,"msex"] <- "FEMALE"
patient.info[patient.info$msex==1,"msex"] <- "MALE"
patient.info$pathologic.diagnosis.of.AD <- as.character(patient.info$pathologic.diagnosis.of.AD)
patient.info[patient.info$pathologic.diagnosis.of.AD=="YES","pathologic.diagnosis.of.AD"] <- "AD"
patient.info[patient.info$pathologic.diagnosis.of.AD=="NO","pathologic.diagnosis.of.AD"] <- "CTRL"
patient.info$pathologic.diagnosis.of.AD <- factor(patient.info$pathologic.diagnosis.of.AD, levels = c("CTRL","AD"),ordered=TRUE)
patient.info$combined <- paste(patient.info$msex, patient.info$pathologic.diagnosis.of.AD,sep = " ")
patient.info$combined <- factor(patient.info$combined, levels = c("FEMALE CTRL", "FEMALE AD", "MALE CTRL", "MALE AD"), ordered = TRUE)
p <- ggplot(patient.info[!is.na(patient.info$mean.var),], aes(x=combined, y=mean.var, fill=msex)) + geom_boxplot()
print(p + scale_fill_manual(values=c("pink","steelblue1"), name ="Sex") + ggtitle(paste0(cell_type, " Transcriptional Noise")) + theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + ylab("Transcriptional Noise"))
print(wilcox.test(x=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="YES" & patient.info$msex=="MALE","mean.var"]),
                  y=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="YES" & patient.info$msex=="FEMALE","mean.var"])))


p <- ggplot(patient.info[!is.na(patient.info$mean.var),], aes(x=pathologic.diagnosis.of.AD, y=mean.var, fill=pathologic.diagnosis.of.AD)) + geom_boxplot()
print(p + scale_fill_manual(values=c("green","red"), name ="AD Diagnosis") + ggtitle(paste0(cell_type, " Transcriptional Noise")) + theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + ylab("Transcriptional Noise"))
print(wilcox.test(x=as.numeric(patient.info[!is.na(patient.info$mean.var) & as.character(patient.info$pathologic.diagnosis.of.AD)=="AD","mean.var"]),
                  y=as.numeric(patient.info[!is.na(patient.info$mean.var) & as.character(patient.info$pathologic.diagnosis.of.AD)=="CTRL","mean.var"])))




patient.info$cogdx <- as.character(patient.info$cogdx)
p <- ggplot(patient.info[!is.na(patient.info$mean.var),], aes(x=cogdx, y=mean.var, fill=cogdx)) + geom_boxplot()
print(p + scale_fill_manual(values=c("pink","steelblue1"), name ="Sex") + ggtitle(paste0(cell_type, " Transcriptional Noise")) + theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + ylab("Transcriptional Noise"))
print(wilcox.test(x=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="YES" & patient.info$msex=="MALE","mean.var"]),
                  y=as.numeric(patient.info[!is.na(patient.info$mean.var) & patient.info$pathologic.diagnosis.of.AD=="YES" & patient.info$msex=="FEMALE","mean.var"])))














mydata <- CreateSeuratObject(master_tbl)
mydata_norm <- NormalizeData(mydata)

mydata_norm <- ScaleData(mydata_norm)
mydata_norm <- FindVariableFeatures(mydata_norm)

mydata_norm <- RunPCA(mydata_norm, features=VariableFeatures(mydata_norm))
mydata_norm <- FindNeighbors(mydata_norm, dims = 1:10)
mydata_norm <- FindClusters(mydata_norm)

mydata_norm2 <- RunUMAP(mydata_norm, dims=1:10)
DimPlot(mydata_norm2, reduction = "umap")
FeaturePlot(mydata_norm2, features=c("NRGN", "GAD1","VCAN","FLT1", "CD74", "AQP4"))



master_tbl_not_impute <- read.csv("MathysData/patient10.csv", row.names = 1)
mydata2 <- CreateSeuratObject(master_tbl_not_impute)
mydata_norm2 <- NormalizeData(mydata2)

mydata_norm2 <- ScaleData(mydata_norm2)
mydata_norm2 <- FindVariableFeatures(mydata_norm2)

mydata_norm2 <- RunPCA(mydata_norm2, features=VariableFeatures(mydata_norm2))
mydata_norm2 <- FindNeighbors(mydata_norm2, dims = 1:10)
mydata_norm2 <- FindClusters(mydata_norm2)

mydata_norm2 <- RunUMAP(mydata_norm2, dims=1:10)
DimPlot(mydata_norm2, reduction = "umap")
FeaturePlot(mydata_norm2, features=c("NRGN", "GAD1","VCAN","FLT1", "CD74", "AQP4"))


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


get.counts <- function(gene, patient, cell_type) {
  s1 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),6] == cell_type
  #s2 <- as.vector(raw_tbl[gene,][names(mydata_norm@assays$RNA@data[gene,])] != 0)
  s3 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),2] == patient
  count.list <- mydata_norm@assays$RNA@data[gene,][s1 & s3]
  print(var(count.list))
  plot(density(count.list))
  return(count.list)
}

compare.variance <- function(gene, cell_type) {
  AD.var.scores <- c()
  for (AD_patient in c("AD1", "AD2", "AD3", "AD4", "AD5", "AD6", "AD-un"))
  {
    print(AD_patient)
    AD.var.scores <- c(AD.var.scores, var(get.counts(gene, AD_patient, cell_type))) 
  }
  Ct.var.scores <- c()
  for (Ct_patient in c("Ct1", "Ct2", "Ct3", "Ct4", "Ct5", "Ct6", "Ct-un"))
  {
    print(Ct_patient)
    Ct.var.scores <- c(Ct.var.scores, var(get.counts(gene, Ct_patient, cell_type))) 
  }
  return(list(Ct.var.scores, AD.var.scores))
}








# cell_pop_vector = c()
# for (gene in rownames(tbl)[1:30])   
# {
#   for (patient in c("AD1", "AD2", "AD3", "AD4", "AD5", "AD6", "AD-un", "Ct1", "Ct2", "Ct3", "Ct4", "Ct5", "Ct6", "Ct-un"))
#   {
#     for (cell_type in c("astro", "doublet", "endo", "mg", "neuron", "oligo", "OPC", "unID"))
#     {
#       s1 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),6] == cell_type
#       s2 <- as.vector(raw_tbl[gene,][names(mydata_norm@assays$RNA@data[gene,])] != 0)
#       s3 <- meta_tbl[names(mydata_norm@assays$RNA@data[gene,]),2] == patient
#       cell_pop_vector <- c(cell_pop_vector, as.vector(mydata_norm@assays$RNA@data[gene,][s1 & s2 & s3]))
#       print(names(cell_pop_vector))
#       print()
#       names(cell_pop_vector) <- c(names(cell_pop_vector), paste(gene,patient,cell_type,sep="_"))
#     }
#   }
# }
# mtx=readMM("MathysData/filtered_count_matrix.mtx")
#> writeMM(mtx, "MathysData/count_matrix.csv")
#