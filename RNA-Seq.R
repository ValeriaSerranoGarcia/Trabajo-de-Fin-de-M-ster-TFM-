################################# RNA - Seq ####################################

######################### Introducción de libreria #############################
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(gplots)
library(DESeq2)
library(RColorBrewer)
library(corrplot)
library(ggfortify)
library(tidyverse)
library(pheatmap)
library(DEGreport)
library(ggrepel)
library(tibble)
library(openxlsx)
library(readxl)

################################################################################
######################### OBTENCIÓN DE TABLA DATOS #############################
################################################################################

# 1. Carga de directorio y creacion de variables.
directorio <- "C:/Users/valer/Documents/GSE123509_RAW"
archivos <- list.files(directorio, full.names = TRUE)
lista_datos <- list() #Lista para almacenar los datos. 


# 2. Cargar los datos del conteo.

# Leer y transformar cada archivo
for (archivo in archivos) {
  # Leer los datos del archivo
  datos <- read.table(archivo, header = FALSE)  # Suponiendo que los datos no tienen encabezado
  
  # Obtener el nombre del archivo sin extension para usar como identificador de muestra
  nombre_archivo <- gsub(".txt", "", basename(archivo))
  
  # Agregar los datos transformados a la lista
  lista_datos[[nombre_archivo]] <- datos
}

# 3. Crear un dataframe combinando los datos de cada archivo

genes <- unique(unlist(lapply(lista_datos, function(x) x[, 1]))) # Extraer los nombres de genes unicos de todos los archivos

# Crear un dataframe vacio con columnas para los genes y muestras
df_final <- data.frame(Gene = genes, stringsAsFactors = FALSE)

# Combinar los datos de cada archivo en el dataframe final
for (nombre_archivo in names(lista_datos)) {
  datos <- lista_datos[[nombre_archivo]]
  
  # Asignar los valores de expresiÓn a las muestras correspondientes en df_final
  df_final[nombre_archivo] <- datos$V2[match(df_final$Gene, datos$V1)]
}

# Guardamos el archivo en la carpeta con los datos. 
ruta_archivo <- file.path(directorio, "datos_procesados.xlsx")

# Guardar el dataframe en un archivo Excel
write.xlsx(df_final, file = ruta_archivo, rowNames = FALSE)

################################################################################
################################ PCA GLOBAL ####################################
################################################################################
# 1. Introducimos las variables a analizar.
setwd("C:/Users/valer/Documents/GSE123509_RAW")
# Variables
comparison <- "N_vs_T" #Comparamos las meustras mutantes con las muestras de ratones silvestres.
condition1 <- "F.M.N" #Muestras de mesencéfalo / prosencéfalo de ratones silvestres. 
condition2 <- "F.M.T" #Muestras de mesencéfalo / prosencéfalo de ratones mutantes.
condition3 <- "Cb.N"  #Muestras del cerebelo de ratones silvestres.
condition4 <- "Cb.T" #Muestras del cerebelo de ratones mutantes.
rep1 <- 3 # Introducimos tres réplicas por cada condición.
rep2 <- 3
rep3 <- 3
rep4 <- 3
pvt <- 0.05 #Establecemos el p-valor

# Condiciones
group <- factor(c(rep(condition1, rep1), rep(condition2, rep2), rep(condition3, rep3), rep(condition4, rep4)))
df_counts <-read.table("data.tsv", sep = "\t", header = TRUE, row.names = 1)

# Seleccion de las muestras que queremos analizar
keep <- grep(condition1,colnames(df_counts))
keep <- append(keep, grep(condition2,colnames(df_counts)))
keep <- append(keep, grep(condition3,colnames(df_counts)))
keep <- append(keep, grep(condition4,colnames(df_counts)))
df_counts <- df_counts[ ,keep]

# Directorio para los resultados
dir.create(comparison, showWarnings = FALSE)
dir.create(paste0(comparison, "/graphs"), showWarnings = FALSE)

#Limpieza de filas con NA
df_counts_clean <- df_counts[complete.cases(df_counts), ]

# Funcion
sampleCondition <- group
sampleTable <- data.frame(row.names = colnames(df_counts_clean), condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromMatrix(countData = df_counts_clean, colData = sampleTable, design = ~condition)
ddsHTSeq <- DESeq(ddsHTSeq)


# Transformacion rlog
ddsHTSeq.rld <- rlogTransformation(ddsHTSeq, blind = TRUE)

# PCA
png(file = paste0(comparison, "/graphs/PCA.png"), width = 7*300, height = 5*300, res = 300)
print(plotPCA(ddsHTSeq.rld, intgroup = "condition"))
dev.off()

################# Expresión diferencial genes del cerebelo (Cb) #################
# Variables
comparison <- "NCb_vs_TCb" #Comparamos muestras de cerebelo
condition1 <- "Cb.N" #Comparamos muestras del cerebelo de ratones silvestres (sin mutación).
condition2 <- "Cb.T" #Comparamos muestras del cerebelo de ratones con la mutación.

rep1 <- 3 #Realizamos el análisis teniendo en cuenta las tres réplicas por muestra. 
rep2 <- 3
pvt <- 0.05

# Condiciones
group <- factor(c(rep(condition1, rep1), rep(condition2, rep2))) #Realizamos la comapración de mutantes con silvestres. 
df_counts <- read.table("datos_procesados2.tsv", sep = "\t", header = TRUE, row.names = 1)

# Seleccion de las muestras que queremos analizar.
keep <- grep(condition1,colnames(df_counts))
keep <- append(keep, grep(condition2,colnames(df_counts)))
df_counts <- df_counts[ ,keep]

# Directorio para los resultados
dir.create(comparison, showWarnings = FALSE)
dir.create(paste0(comparison, "/graphs"), showWarnings = FALSE)

# Funcion
sampleCondition <- group
sampleTable <- data.frame(row.names = colnames(df_counts), condition = sampleCondition) #Tabla con información de las muestras.
ddsHTSeq <- DESeqDataSetFromMatrix(countData = df_counts, colData = sampleTable, design = ~condition) #Creamos las variables necesarias con los datos para realizar mi expresión diferencial.
ddsHTSeq <- DESeq(ddsHTSeq) #Realizamos la expresión diferencial. 

# Generamos los resultados del análisis de expresión diferencial para todas las comparaciones de genes entre las condiciones.
results_df <- results(ddsHTSeq, contrast=c("condition", condition2, condition1))
write.table(as.matrix(results_df), file = paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE, na = "")

# Normalizamos los conteos teniendo en cuenta cada muestra y cada gen.
normalized_counts <- counts(ddsHTSeq, normalized = TRUE)
write.table(as.matrix(normalized_counts), file = paste0(comparison, "/DESeq2_counts.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# Filtramos los genes que se expresan diferencialmente y se guardan.
ddsHTSeq.res <- results(ddsHTSeq, contrast=c("condition", condition2, condition1))
ddsHTSeq.res.fdr <- ddsHTSeq.res[!is.na(ddsHTSeq.res$padj), ]
ddsHTSeq.res.fdr <- ddsHTSeq.res.fdr[ddsHTSeq.res.fdr$padj < pvt, ]
write.table(as.matrix(ddsHTSeq.res.fdr), file = paste0(comparison, "/DESeq2_FDR.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# Transformacion rlog
ddsHTSeq.rld <- rlogTransformation(ddsHTSeq, blind = TRUE)

# Heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distRL <- dist(t(assay(ddsHTSeq.rld)))
mat <- as.matrix(distRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsHTSeq), paste0(condition))
png(file = paste0(comparison, "/graphs/heatmap.png"), width = 6*300, height = 6*300, res = 300)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13, 13))
dev.off()

# PCA
png(file = paste0(comparison, "/graphs/PCA.png"), width = 7*300, height = 5*300, res = 300)
print(plotPCA(ddsHTSeq.rld, intgroup = "condition"))
dev.off()

# Correlograma
png(file = paste0(comparison, "/graphs/corr_plot.png"), width = 7*300, height = 5*300, res = 300)
corrplot(cor(normalized_counts), method = "square", addCoef.col = "white")
dev.off()

# Volcano plot
data <- read.table(paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", header = TRUE, row.names = 1)

data <- data %>%
  mutate(
    Expression = case_when(
      log2FoldChange > 1.5 & pvalue < pvt ~ "Up",
      log2FoldChange < -1.5 & pvalue < pvt ~ "Down",
      TRUE ~ "No significance"),
    Significance = -log10(pvalue))
umbral_y <- -log10(pvt)

p <- ggplot(data, aes(x = log2FoldChange, y = Significance)) +
  geom_point(aes(color = Expression), alpha = 0.4, size = 1.6) +
  scale_color_manual(values = c("red3", "gray60", "green3")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  geom_vline(xintercept = -1, linetype = "dotdash", color = "gray25") +
  geom_vline(xintercept = 1, linetype = "dotdash", color = "gray25") +
  geom_hline(yintercept = umbral_y, linetype = "dotdash", color = "gray25") +
  coord_cartesian(xlim = c(-15, 15)) + 
  theme(panel.background = element_rect(fill = "gray96"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  xlab(expression(log[2] * "FC")) +
  ylab(expression("-log"[10] * "PValue"))

suppressWarnings(ggsave(paste0(comparison, "/graphs/volcano_plot.png"), plot = p, width = 7, height = 5))


####### Expresión diferencial genes del mesencéfalo / prosencéfalo #########
# Variables
comparison <- "NMF_vs_TMF"
condition1 <- "F.M.N"
condition2 <- "F.M.T"

rep1 <- 3
rep2 <- 3
pvt <- 0.05

# Condiciones. Muestras con 3 réplicas + cargar archivos. 
group <- factor(c(rep(condition1, rep1), rep(condition2, rep2)))
df_counts <- read.table("datos_procesados.tsv", sep = "\t", header = TRUE, row.names = 1)

# Selección de las muestras que queremos analizar
keep <- grep(condition1,colnames(df_counts))
keep <- append(keep, grep(condition2,colnames(df_counts)))
df_counts <- df_counts[ ,keep] #Seleccionamos las columnas que corresponden a las condiciones específicas y las añadimos junto con nuestro archivo de las muestras.

##############
### DESEQ2 ###
##############

# Directorio para los resultados
dir.create(comparison, showWarnings = FALSE)
dir.create(paste0(comparison, "/graphs"), showWarnings = FALSE)

# Función. Creamos las condiciones y creamos DESeq para el análisis. 
sampleCondition <- group
sampleTable <- data.frame(row.names = colnames(df_counts), condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromMatrix(countData = df_counts, colData = sampleTable, design = ~condition)
ddsHTSeq <- DESeq(ddsHTSeq)

# Se exponen los resultados del análisis de expresión diferencial. 
results_df <- results(ddsHTSeq, contrast=c("condition", condition2, condition1))
write.table(as.matrix(results_df), file = paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE, na = "")

# Normalizamos los conteos.
normalized_counts <- counts(ddsHTSeq, normalized = TRUE)
write.table(as.matrix(normalized_counts), file = paste0(comparison, "/DESeq2_counts.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# Obtenemos la expresión diferencial dada las condiciones estudiadas.
ddsHTSeq.res <- results(ddsHTSeq, contrast=c("condition", condition2, condition1))
ddsHTSeq.res.fdr <- ddsHTSeq.res[!is.na(ddsHTSeq.res$padj), ]
ddsHTSeq.res.fdr <- ddsHTSeq.res.fdr[ddsHTSeq.res.fdr$padj < pvt, ]
write.table(as.matrix(ddsHTSeq.res.fdr), file = paste0(comparison, "/DESeq2_FDR.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# Transformación rlog
ddsHTSeq.rld <- rlogTransformation(ddsHTSeq, blind = TRUE)

# Heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distRL <- dist(t(assay(ddsHTSeq.rld)))
mat <- as.matrix(distRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsHTSeq), paste0(condition))
png(file = paste0(comparison, "/graphs/heatmap.png"), width = 6*300, height = 6*300, res = 300)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13, 13))
dev.off()

# PCA
png(file = paste0(comparison, "/graphs/PCA.png"), width = 7*300, height = 5*300, res = 300)
print(plotPCA(ddsHTSeq.rld, intgroup = "condition"))
dev.off()

# Correlograma
png(file = paste0(comparison, "/graphs/corr_plot.png"), width = 7*300, height = 5*300, res = 300)
corrplot(cor(normalized_counts), method = "square", addCoef.col = "white")
dev.off()

# Volcano plot
data <- read.table(paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", header = TRUE, row.names = 1)

data <- data %>%
  mutate(
    Expression = case_when(
      log2FoldChange > 1 & pvalue < pvt ~ "Up",
      log2FoldChange < -1 & pvalue < pvt ~ "Down",
      TRUE ~ "No significance"),
    Significance = -log10(pvalue))
umbral_y <- -log10(pvt)

p <- ggplot(data, aes(x = log2FoldChange, y = Significance)) +
  geom_point(aes(color = Expression), alpha = 0.4, size = 1.6) +
  scale_color_manual(values = c("red3", "gray60", "green3")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  geom_vline(xintercept = -1, linetype = "dotdash", color = "gray25") +
  geom_vline(xintercept = 1, linetype = "dotdash", color = "gray25") +
  geom_hline(yintercept = umbral_y, linetype = "dotdash", color = "gray25") +
  coord_cartesian(xlim = c(-15, 15), ylim = c(0, max(data$Significance) + 10)) + 
  theme(panel.background = element_rect(fill = "gray96"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  xlab(expression(log[2] * "FC")) +
  ylab(expression("-log"[10] * "PValue"))

suppressWarnings(ggsave(paste0(comparison, "/graphs/volcano_plot.png"), plot = p, width = 7, height = 5))

################################################################################
############################# DIAGRAMA DE VENN #################################
################################################################################
setwd("D:/TFM results definitivo/RNA-Seq")
N_T = read.table("filtered_genesNT.tsv",header = TRUE, sep = "\t")
NCb_TCb = read.table("filtered_genesCb.tsv",header = TRUE, sep = "\t")

#Plot
venn.plot <- venn.diagram(
  x = list(C = N_T$X, D = NCb_TCb$X),
  category.names = c("M/F", "Cb"),
  fill = c("blue", "purple"),  # Especificar colores
  alpha = 0.5,  # Transparencia de los colores
  filename = NULL,
  output = TRUE
)

# Visualizar el diagrama
grid.draw(venn.plot)

# Guardar el diagrama de Venn en un archivo PNG
png("diagrama_venn.png", width = 2000, height = 2000, res = 300)
grid.draw(venn.plot)
dev.off()

#Genes comunes
genes_comunes <- intersect(N_T$X, NCb_TCb$X)
genes_comunes


################################################################################
########################### ENRIQUECIMEINTO FUNCIONAL ##########################
################################################################################

#Fijamos el organismo con el que vamos a trabajar
organism = "org.Mm.eg.db" #Mus musculus.
library(organism, character.only = TRUE)

#Instalamos directorio
setwd("C:/Users/valer/Documents/GSE123509_RAW/NMF_vs_TMF") #Si queremos que el análisis sea del cerebelo debemos de introducir "C:/../Documents/GSE123509_RAW/NCb_vs_TCb".

#Realizamos la operaci?n de lecturas de valores.
df = read.table("filtered_genes.tsv", header=TRUE)  

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#Formamos nuestra table
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ALIAS", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

# Calcular la matriz de similitud de términos
gse <- pairwise_termsim(gse)

#Dotplot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

#Mapa de enriquecimiento. 
emapplot(gse, showCategory = 10)

#Red categorica.
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

#Ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

#GSA Plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

write.csv2(as.data.frame(gse), file="MF") #Creamos el archivo que nos relacione las funciones enriquecidas con los genes. 

