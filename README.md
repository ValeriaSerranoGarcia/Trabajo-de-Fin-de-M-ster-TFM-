# Trabajo de Fin de Master (TFM).

A lo largo de este proyecto, se introducen los comandos originados para el análisis de las muestras GSE7768 (NCL1) y GSE123509 (NCL2). Se emplean diferentes procedimientos dependiendo de si la técnica empleada para la obtención de los datos ha sido Microarray o RNA-Seq. 

## Índice

1. [Análisis de NCL1, muestras obtenidas de Microarray](#Análisis-NCL1)
   
   1.1 [Normalización de los resultados](#Normalización-de-los-resultados)
   
   1.2 [Análisis de la expresión diferencial](#Análisis-de-la-expresión-diferencial)

  1.3 [Diagrama de Venn](#Diagrama-de-Venn)

  1.4 [Enriquecimiento funcional](#Enriquecimiento-funcional)

  
2. [Análisis de NCL2, muestras obtenidas de RNA-Seq](#Análisis-NCL2)
   
  2.1 [Obtención y ordenamiento de los datos](#Obtención-y-ordenamiento-de-los-datos)
   
   2.2 [Análisis de la expresión diferencial y filtrado de genes](#Análisis-de-la-expresión-diferencial-y-filtrado-de-genes)

  2.3 [Diagrama de Venn](#Diagrama-de-Venn)

  2.4 [Enriquecimiento funcional](#Enriquecimiento-funcional)

## 1. Análisis NCL1

### 1.1. Normalización de los resultados
Antes de la normalización de los resultados es necesario que carguemos las librerías que vamos a emplear a lo largo de este análisis. 
```r
library(oligo)
library(Biobase)
library(affy)
library(splitstackshape)
library(GEOquery)
library(limma)
library(affycoretools) 
library(mouse4302.db)
library(dplyr)
library(reshape2)
library(readxl)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyr)
library(VennDiagram)
library(readr)
library(dplyr)
```

Una vez introducido las librerías necesarias, a continuación vamos a cargar los resultados y a normalizarlos. Con ello conseguiremos modelar los archivos que contienen los datos de estudio para poder realizar el análisis de la expresión diferencial. 

```r
# 1. Introducimos directorios y cargamos los archivos. 
setwd("~/Documents")
celFiles<-list.celfiles() #Cargamos los archivos .cel que hemos descargado de nuestro directorio.
affyRaw<-read.celfiles(celFiles) #Leemos los archivos. 

# 2. Normalozación de los resultados. 
eset<-oligo::rma(affyRaw) #Normalizamos los archivos y calculamos la expresión para posteriormente guardarlo en la variable eset.
pData(eset) #Comprobamos los eset que se han cargado, esto nos debe proporcionar una información de exactamente 18 archivos que representan condiciones de mutantes y tiempos diferentes. 
geneSummaries <- annotateEset(eset,mouse4302.db)

# 3. Preparación de los resultados para el análisis de expresión diferencial.
target_table<-read.csv("target_table.csv",sep=";") #Introducimos nuestra matriz target_table que hemos creado recopilando los datos del estudio. 
groups<-target_table[,3] #Introducimos la información de mutante y KO de la target_table en mi variable groups. 
design<-model.matrix(~factor(groups)) #Introducidmos los grupos en dicha matriz.
colnames(design)<-c("WT","WTvsKO") #Cambiamos el título de la matriz que representa las condiciones que se van a comparar. 
design #Observamos nuestra matriz.
````

### 1.2. Análisis de la expresión diferencial
Para la realización del análisis de la expresión diferencial, debemos de establecer previamente las comparaciones que se quieren realizar. En este caso, al ser un análisis de expresión diferencial a lo largo del tiempo, debemos introducir las comparaciones a los tres, cinco y ocho meses. 

```r
#Introducimos las muestras a estudiar y las comparaciones a realizar.
lev<-c("KO.3mo","KO.3mo","KO.3mo","KO.5mo","KO.5mo","KO.5mo","KO.8mo","KO.8mo","KO.8mo","WT.3mo","WT.3mo","WT.3mo","WT.5mo","WT.5mo","WT.5mo","WT.8mo","WT.8mo","WT.8mo")
design <- model.matrix(~0+lev)
fit_t<-lmFit(geneSummaries,design) 
cont.dif <- makeContrasts(
  Dif3mo= (levKO.3mo)-(levWT.3mo),
  Dif5mo =(levKO.5mo)-(levWT.5mo),
  Dif8mo=(levKO.8mo)-(levWT.8mo),
  levels=design)


# Analisis de la expresion diferencial
fitime <- contrasts.fit(fit_t, cont.dif)
fittime<- eBayes(fitime)
log2_threshold <- log2(1.5) #Parámetro para realizar la expresión diferencial.
fit_T<-topTable(fittime, adjust="BH", number=Inf) #Expresion diferencial de los genes en función del tiempo.
res_time_def<-subset(fit_T,adj.P.Val<0.05 & (abs (Dif3mo)>log2_threshold | abs (Dif5mo)>log2_threshold | abs (Dif8mo)>log2_threshold))#Seleccionamos los genes que presenten un valor de p-value<0.05 y además estamos seleccionando aquellos cuyo valor absoluto sea mayor de 1. 
write.table(res_time_def,"differential_expressionTIME.txt",sep="\t") #Almacenamos datos expresión diferencial.

```

### 1.3. Diagrama de Venn
La realización de un Digrama de Venn nos permite poder identificar genes comunes a los tres tiempos estudiados. Además, nos permite poder observar que existen genes exclusivos de cada uno de los tiempos.

```r
############################ EXTRACCION DATOS ##############################
# 1. Cargamos los datos de expresión diferencial e información de la muestra.
exp_table<-read.table("differential_expressionTIME.tsv") #Datos de expresion diferencial.
target_table<-read.csv("target_table.csv",sep=";") #Datos con información de las muestras.

# 2. Modificacion de las tablas en formatos iguales y manejables.
# Poner la exp_table para tener las condiciones de tiempo en una sola columna
exp_table_melted <- melt(exp_table, id.vars = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME","P.Value"),
                         variable.name = "TIME", value.name = "Expression")

# Convertir la columna TIME a un formato que coincida con target_table
exp_table_melted$TIME <- gsub("Dif", "", exp_table_melted$TIME)

# Unir las tablas
result_table <- merge(target_table, exp_table_melted, by = "TIME")

# Seleccionar columnas deseadas
result_table <- result_table %>% select(GEO_ID, CONDITION, TIME, PROBEID, SYMBOL, Expression, P.Value)
print(result_table)

log2_threshold <- log2(1.5)
dif3_months= result_table %>% filter(TIME == "3mo" & Expression > abs(log2_threshold) & P.Value < 0.05)
#Vamos a descargar los genes de 5 meses que se expresan diferencialmente.
dif5_months= result_table %>% filter(TIME == "5mo" & Expression > abs(log2_threshold) & P.Value < 0.05)
#Vamos a descargar los genes de 8 meses que se expresan diferencialmente.
dif8_months= result_table %>% filter(TIME == "8mo" & Expression > abs(log2_threshold) & P.Value < 0.05)

write_tsv(dif3_months, "dif3mo.tsv")
write_tsv(dif5_months, "dif5mo.tsv")
write_tsv(dif8_months, "dif8mo.tsv")

################### Representación del Diagrama de Venn ####################
# 1. Eliminacion de valores NA para evitar errores en el resultado.
dif3_months = na.omit(dif3_months)
dif5_months = na.omit(dif5_months)
dif8_months = na.omit(dif8_months)

# 2. Creacion del Diagrama de Venn.
venn.plot <- venn.diagram(
  x = list(A = dif3_months$SYMBOL, B = dif5_months$SYMBOL, C = dif8_months$SYMBOL),
  category.names = c("3 meses", "5 meses", "8 meses"),
  fill = c("yellow", "orange", "green"),  # Especificar colores
  alpha = 0.5,  # Transparencia de los colores
  filename = NULL,
  output = TRUE
)

# 3. Visualizacion  del diagrama y guardado del mismo.
grid.draw(venn.plot)
png("diagrama_venn.png", width = 2000, height = 2000, res = 300)
grid.draw(venn.plot)
dev.off()

# 4. Extraccion de genes comunes y guardado. 
genes_comunes_35 <- intersect(dif3_months$SYMBOL, dif5_months$SYMBOL)
genes_comunes_38 <- intersect(dif3_months$SYMBOL, dif8_months$SYMBOL)
genes_comunes_58 <- intersect(dif5_months$SYMBOL, dif8_months$SYMBOL)
genes_comunes_358 <- Reduce(intersect, list(dif3_months$SYMBOL, dif5_months$SYMBOL, dif8_months$SYMBOL))

write.table(genes_comunes_35, file = "genes_comunes_35.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_comunes_38, file = "genes_comunes_38.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_comunes_58, file = "genes_comunes_58.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_comunes_358, file = "genes_comunes_358.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# 5. Busqueda de los genes exclusivos expresados en cada seccion temporal.
mo3 <- setdiff(dif3_months$SYMBOL, c(genes_comunes_358, genes_comunes_35, genes_comunes_38))
write.table(mo3, file = "mo3.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

mo5 <- setdiff(dif5_months$SYMBOL, c(genes_comunes_358, genes_comunes_35, genes_comunes_58))
write.table(mo5, file = "mo5.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

mo8 <- setdiff(dif8_months$SYMBOL, c(genes_comunes_358, genes_comunes_58, genes_comunes_38))
write.table(mo8, file = "mo8.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

### 1.4 Enriquecimiento funcional.
El análisis de enriquecimiento funcional permite establecer procesos biológicos que se encuentran enriquecidos en el conjunto de genes seleccionados. 

```r
# 1. Preparación de variables y datos que se necesitan. 

#Preparación y filtrado de los datos.
data = read.table("differential_expressionTIME.tsv")
mon3_dif <- data[, c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME", "Dif3mo", "AveExpr", "F", "P.Value", "adj.P.Val")]
write.table(mon3_dif,"diferencial_3.txt")

mon5_dif <- data[, c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME", "Dif5mo", "AveExpr", "F", "P.Value", "adj.P.Val")]
write.table(mon5_dif,"diferencial_5.txt")

mon8_dif <- data[, c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME", "Dif8mo", "AveExpr", "F", "P.Value", "adj.P.Val")]
write.table(mon8_dif,"diferencial_8.txt")

#Lectura de los datos y filtrado
df = read.table("diferencial_5.txt", header=TRUE) #filtered_genes.tsv <-- En RNA-Seq. || dif3mo.tsv --> En Microarray. 
log2_threshold <- log2(1.5)
df = df %>% filter(abs(Dif5mo) > abs(log2_threshold) & P.Value < 0.05)


#Fijamos el organismo con el que vamos a trabajar
organism = "org.Mm.eg.db" #Mus musculus.
library(organism, character.only = TRUE)

# Obtenemos el valor de la expresión para el analisis. 
original_gene_list <- df$Dif5mo
names(original_gene_list) <- df$SYMBOL

# Omitimos los valores de NA. 
gene_list<-na.omit(original_gene_list)

# Ordenamos la lista en orden decreciente (necesario para clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# 2. Análisis del enriquecimiento funcional. 
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

# 3. Graficos.

# Calcular la matriz de similitud de términos
gse <- pairwise_termsim(gse)

#Dotplot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

#Mapa de enriquecimiento. 
emapplot(gse, showCategory = 10)

#Red categorica.
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

# 4. Almacenamiento del enriquecimiento funcional
write.csv2(as.data.frame(gse), file="test5mo.csv")
```

## 2. Análisis de NCL2

### 2.1 Obtención y ordenamiento de los datos.
Antes de la realización de los diversos procesos que nos permitan obtener genes expresados diferencialmente y rutas enriquecidas, debemos de cargar las librerías que vamos a emplear. 

```r
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
```
Posteriormente, vamos a cargar los datos de RNA-Seq y a ordenarlos para hacer posible el análisis de los mismos.
```r
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
```

Realizamos un análisis de PCA global, para observar la distribución de los datos tanto del cerebelo como del mesencéfalo / prosencéfalo.
```r
setwd("C:/Users/valer/Documents/GSE123509_RAW")
# Variables
comparison <- "N_vs_T"
condition1 <- "F.M.N"
condition2 <- "F.M.T"
condition3 <- "Cb.N"
condition4 <- "Cb.T"
rep1 <- 3
rep2 <- 3
rep3 <- 3
rep4 <- 3
pvt <- 0.05

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
```

### 2.2 Análisis de la expresión diferencial y filtrado de genes.
En este caso, el análisis de la expresión diferencial se emplea mediante el paquete DESeq2. Hemos realizado por separado el análisis de la expresión diferencial del cerebelo (Cb) y del mesencéfalo / prosencéfalo (M/F)

```r
################# Expresión diferencial genes del cerebelo #################
# Variables
comparison <- "NCb_vs_TCb"
condition1 <- "Cb.N"
condition2 <- "Cb.T"

rep1 <- 3
rep2 <- 3
pvt <- 0.05

# Condiciones
group <- factor(c(rep(condition1, rep1), rep(condition2, rep2)))
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
sampleTable <- data.frame(row.names = colnames(df_counts), condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromMatrix(countData = df_counts, colData = sampleTable, design = ~condition)
ddsHTSeq <- DESeq(ddsHTSeq)

# Resultados
results_df <- results(ddsHTSeq, contrast=c("condition", condition2, condition1)) #Primero condición mutante condition2
write.table(as.matrix(results_df), file = paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE, na = "")

# Counts
normalized_counts <- counts(ddsHTSeq, normalized = TRUE)
write.table(as.matrix(normalized_counts), file = paste0(comparison, "/DESeq2_counts.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# FDR
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

# Resultados
results_df <- results(ddsHTSeq, contrast=c("condition", condition2, condition1))
write.table(as.matrix(results_df), file = paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE, na = "")

# Counts
normalized_counts <- counts(ddsHTSeq, normalized = TRUE)
write.table(as.matrix(normalized_counts), file = paste0(comparison, "/DESeq2_counts.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# FDR
ddsHTSeq.res <- results(ddsHTSeq, contrast=c("condition", condition2, condition1))
ddsHTSeq.res.fdr <- ddsHTSeq.res[!is.na(ddsHTSeq.res$padj), ]
ddsHTSeq.res.fdr <- ddsHTSeq.res.fdr[ddsHTSeq.res.fdr$padj < pvt, ]
write.table(as.matrix(ddsHTSeq.res.fdr), file = paste0(comparison, "/DESeq2_FDR.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

################
### GRÁFICOS ###
################

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

```

### 2.3 Diagrama de Venn.

```r
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
```

## 2.4 Enriquecimiento funcional.

```r

#Fijamos el organismo con el que vamos a trabajar
organism = "org.Mm.eg.db" #Mus musculus.
library(organism, character.only = TRUE)

#Instalamos directorio
setwd("C:/Users/valer/Documents/GSE123509_RAW/NMF_vs_TMF")

#Realizamos la operaci?n de lecturas de valores.
df = read.table("filtered_genes.tsv", header=TRUE) #filtered_genes.tsv <-- En RNA-Seq. || dif3mo.tsv --> En Microarray. 

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

#keytypes(org.Mm.eg.db)
#[1] "ACCNUM"       "ALIAS"        "ENSEMBL"     
#[4] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
#[7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
#[10] "GENENAME"     "GENETYPE"     "GO"          
#[13] "GOALL"        "IPI"          "MGI"         
#[16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"        
#[19] "PFAM"         "PMID"         "PROSITE"     
#[22] "REFSEQ"       "SYMBOL"       "UNIPROT"

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

write.csv2(as.data.frame(gse), file="MF")

```


