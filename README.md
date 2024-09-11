# Trabajo de Fin de Master (TFM)

A lo largo de este proyecto, se introducen los comandos originados para el análisis de las muestras GSE7768 (NCL1) y GSE123509 (NCL2). Se emplean diferentes procedimientos dependiendo de si la técnica empleada para la obtención de los datos ha sido Microarray o RNA-Seq. 

## Índice

1. Análisis de NCL1, muestras obtenidas de Microarray
   
   1.1 Normalización de los resultados
   
   1.2 Análisis de la expresión diferencial

   1.3 Diagrama de Venn

   1.4 Enriquecimiento funcional

  
2. Análisis de NCL2, muestras obtenidas de RNA-Seq
   
   2.1 Obtención y ordenamiento de los datos

   2.2 Análisis de la expresión diferencial y filtrado de genes

   2.3 Diagrama de Venn

   2.4 Enriquecimiento funcional


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
# 1. Introducción de directorios y cargar archivos. 
setwd("~/Documents") #Los archivos deben estar introducidos en el directorio de Documents. Si están en otra ubicación, se deberá cambiar este comando a: setwd("C:/ruta_absoluta_del_directorio")
celFiles<-list.celfiles() #Cargamos los archivos .cel que hemos descargado.
affyRaw<-read.celfiles(celFiles) #Leemos los archivos. 

# 2. Normalozación de los resultados. 
eset<-oligo::rma(affyRaw) #Normalizamos los archivos y calculamos la expresión para posteriormente guardarlo en la variable eset.
pData(eset) #Comprobamos los eset que se han cargado, esto nos debe proporcionar una información de exactamente 18 archivos que representan condiciones de mutantes y tiempos diferentes. 
geneSummaries <- annotateEset(eset,mouse4302.db)

# 3. Preparación de los resultados para el análisis de expresión diferencial.
target_table<-read.csv("target_table.csv",sep=";") #Introducimos nuestra matriz target_table que hemos creado recopilando los datos del estudio. Esta tabla se encuentra disponible en los archivos del Github.
groups<-target_table[,3] #Introducimos la información de mutante y KO de la target_table en la variable groups. 
design<-model.matrix(~factor(groups)) #Introducimos los grupos en dicha matriz.
colnames(design)<-c("WT","WTvsKO") #Cambiamos el título de la matriz que representa las condiciones que se van a comparar. 
design #Observamos nuestra matriz.
````

### 1.2. Análisis de la expresión diferencial
Para la realización del análisis de la expresión diferencial, debemos de establecer previamente las comparaciones que se quieren realizar. En este caso, al ser un análisis de expresión diferencial a lo largo del tiempo, debemos introducir las comparaciones a los tres, cinco y ocho meses. 

```r
#Introducción de muestras a comparar y definición de la comparación.
lev<-c("KO.3mo","KO.3mo","KO.3mo","KO.5mo","KO.5mo","KO.5mo","KO.8mo","KO.8mo","KO.8mo","WT.3mo","WT.3mo","WT.3mo","WT.5mo","WT.5mo","WT.5mo","WT.8mo","WT.8mo","WT.8mo")
design <- model.matrix(~0+lev)
fit_t<-lmFit(geneSummaries,design) # LmFit nos permite observar si existen diferencias en la expresión entre las condiciones estudiadas (WT vs KO) y los tiempos expuestos. 
cont.dif <- makeContrasts(
  Dif3mo= (levKO.3mo)-(levWT.3mo),
  Dif5mo =(levKO.5mo)-(levWT.5mo),
  Dif8mo=(levKO.8mo)-(levWT.8mo),
  levels=design)


# Analisis de la expresion diferencial
fitime <- contrasts.fit(fit_t, cont.dif)
fittime<- eBayes(fitime) #Identifica genes con expresión diferencial.
log2_threshold <- log2(1.5) #Parámetro para realizar la expresión diferencial.
fit_T<-topTable(fittime, adjust="BH", number=Inf) #Expresion diferencial de los genes en función del tiempo.
res_time_def<-subset(fit_T,adj.P.Val<0.05 & (abs (Dif3mo)>log2_threshold | abs (Dif5mo)>log2_threshold | abs (Dif8mo)>log2_threshold))#Seleccionamos los genes que presenten un valor de p-value<0.05 y además estamos seleccionando aquellos cuyo valor absoluto sea mayor de 1. 
write.table(res_time_def,"differential_expressionTIME.txt",sep="\t") #Almacenamos datos expresión diferencial.

```

### 1.3. Diagrama de Venn
La realización de un Digrama de Venn nos permite poder identificar genes comunes a los tres tiempos estudiados. Además, nos permite poder observar que existen genes exclusivos de cada uno de los tiempos.

```r

#1. Preparación y filtrado de los datos. Obtenemos tres archivos para los tres tiempos diferentes con los valores de expresión. 
data = read.table("differential_expressionTIME.tsv")
mon3_dif <- data[, c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME", "Dif3mo", "AveExpr", "F", "P.Value", "adj.P.Val")]
write.table(mon3_dif,"diferencial_3.txt")

mon5_dif <- data[, c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME", "Dif5mo", "AveExpr", "F", "P.Value", "adj.P.Val")]
write.table(mon5_dif,"diferencial_5.txt")

mon8_dif <- data[, c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME", "Dif8mo", "AveExpr", "F", "P.Value", "adj.P.Val")]
write.table(mon8_dif,"diferencial_8.txt")

# 2. Segundo filtrado para corroborar que la selección de genes diferenciales se ha hecho correctamente. 
log2_threshold <- log2(1.5)
dif3_months= mon3_dif %>% filter(abs(Dif3mo) > log2_threshold & P.Value < 0.05)
#Vamos a descargar los genes de 5 meses que se expresan diferencialmente.
dif5_months= mon5_dif %>% filter(abs(Dif5mo) > log2_threshold & P.Value < 0.05)
#Vamos a descargar los genes de 8 meses que se expresan diferencialmente.
dif8_months= mon8_dif %>% filter(abs(Dif8mo) > log2_threshold & P.Value < 0.05)

#Descargamos los archivos.
write_tsv(dif3_months, "dif3mo.tsv")
write_tsv(dif5_months, "dif5mo.tsv")
write_tsv(dif8_months, "dif8mo.tsv")


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
# 1. Lectura de los datos y filtrado
df = read.table("diferencial_5.txt", header=TRUE) #Sustituir por diferencial_3.txt o diferencial_8.txt en caso de que se desee realizar el análisis a los tres u ocho meses respectivamente.  
log2_threshold <- log2(1.5)
df = df %>% filter(abs(Dif5mo) > abs(log2_threshold) & P.Value < 0.05) #Realizamos un nuevo filtro para asegurarnos de que obtenemos los parámetros deseados. 

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

# 4. Almacenamiento de las funciones enriquecidas junto con los genes asociados. 
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
directorio <- "C:/Users/valer/Documents/GSE123509_RAW" #En este caso, el directorio se debe sustituir por el directorio donde se hayan almacenado los datos descargados. 
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
sampleCondition <- group #Asignamos las condiciones experimentales de las muestras a comparar.
sampleTable <- data.frame(row.names = colnames(df_counts_clean), condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromMatrix(countData = df_counts_clean, colData = sampleTable, design = ~condition)#Creamos las variables necesarias con los datos para realizar mi expresión diferencial.
ddsHTSeq <- DESeq(ddsHTSeq) #Cálculo de los genes expresados diferencialmente.


# Transformacion rlog de los datos.
ddsHTSeq.rld <- rlogTransformation(ddsHTSeq, blind = TRUE)

# Representación de PCA para las muestras
png(file = paste0(comparison, "/graphs/PCA.png"), width = 7*300, height = 5*300, res = 300)
print(plotPCA(ddsHTSeq.rld, intgroup = "condition"))
dev.off()
```

### 2.2 Análisis de la expresión diferencial y filtrado de genes.
En este caso, el análisis de la expresión diferencial se emplea mediante el paquete DESeq2. Hemos realizado por separado el análisis de la expresión diferencial del cerebelo (Cb) y del mesencéfalo / prosencéfalo (M/F)

```r
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
sampleCondition <- group #Asignamos las condiciones experimentales.
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

```

### 2.3 Diagrama de Venn.

```r
#Definición de las variables, en este caso definimos las variables asociadas a los genes que se expresan diferencialmente en el cerebelo y en el mesencéfalo / prosencéfalo. Se recomienda en este caso crear un directorio y poner los dos archivos filtered_genes en este. En este caso, hemos definido los archivos de forma diferente para identificar cual pertenece al cerebelo y cual pertenece al mesencéfalo / prosencéfalo.
N_T = read.table("filtered_genesNT.tsv",header = TRUE, sep = "\t")
NCb_TCb = read.table("filtered_genesCb.tsv",header = TRUE, sep = "\t")

#Representación del diagrama.
venn.plot <- venn.diagram(
  x = list(C = N_T$X, D = NCb_TCb$X),
  category.names = c("M/F", "Cb"),
  fill = c("blue", "purple"),  # Especificar colores
  alpha = 0.5,  # Transparencia de los colores
  filename = NULL,
  output = TRUE
)

# Visualización del diagrama
grid.draw(venn.plot)

# Guardar el diagrama de Venn en un archivo PNG
png("diagrama_venn.png", width = 2000, height = 2000, res = 300)
grid.draw(venn.plot)
dev.off()

#Obtención de los genes comunes resultantes de la intersección de la gráfica.
genes_comunes <- intersect(N_T$X, NCb_TCb$X)
genes_comunes
```

### 2.4 Enriquecimiento funcional.

```r
#Se fija el organismo con el que se va a trabajar, en el caso de este estudio, Mus Musculus,
organism = "org.Mm.eg.db" #Mus musculus.
library(organism, character.only = TRUE)

#Se introduce el directorio de trabajo.
setwd("C:/Users/valer/Documents/GSE123509_RAW/NMF_vs_TMF") #Si queremos que el análisis sea del cerebelo debemos de introducir "C:/../Documents/GSE123509_RAW/NCb_vs_TCb".

#Lectura de los genes expresados diferencialmente. 
df = read.table("filtered_genes.tsv", header=TRUE)  

# Se selecciona el valor de expresión de la lista de genes. 
original_gene_list <- df$log2FoldChange

# Se seleccionan los nombres de los genes asociados a nuestro archivo de genes expresados diferencialmente. 
names(original_gene_list) <- df$X

# Se eliminan los valores NA, para evitar problemas en el análisis. 
gene_list<-na.omit(original_gene_list)

# Se ordena la lista de genes, necesaria para realizar el análisis.
gene_list = sort(gene_list, decreasing = TRUE)

#Creación y ejecución del análisis de expresión diferencial. 
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

```


