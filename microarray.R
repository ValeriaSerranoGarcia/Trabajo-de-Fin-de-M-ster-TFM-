##################### Análisis de Microarray ###################################

######################### Introducción de libreria #############################
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

####################### Análisis de resultados Microarray ######################

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

#Introducimos las muestras a estudiar y las comparaciones a realizar.
lev<-c("KO.3mo","KO.3mo","KO.3mo","KO.5mo","KO.5mo","KO.5mo","KO.8mo","KO.8mo","KO.8mo","WT.3mo","WT.3mo","WT.3mo","WT.5mo","WT.5mo","WT.5mo","WT.8mo","WT.8mo","WT.8mo")
design <- model.matrix(~0+lev)
fit_t<-lmFit(geneSummaries,design) 
cont.dif <- makeContrasts(
  Dif3mo= (levKO.3mo)-(levWT.3mo),
  Dif5mo =(levKO.5mo)-(levWT.5mo),
  Dif8mo=(levKO.8mo)-(levWT.8mo),
  levels=design)


# 4. Analisis de la expresion diferencial
fitime <- contrasts.fit(fit_t, cont.dif)
fittime<- eBayes(fitime)
log2_threshold <- log2(1.5) #Parámetro para realizar la expresión diferencial.
fit_T<-topTable(fittime, adjust="BH", number=Inf) #Expresion diferencial de los genes en función del tiempo.
res_time_def<-subset(fit_T,adj.P.Val<0.05 & (abs (Dif3mo)>log2_threshold | abs (Dif5mo)>log2_threshold | abs (Dif8mo)>log2_threshold))#Seleccionamos los genes que presenten un valor de p-value<0.05 y además estamos seleccionando aquellos cuyo valor absoluto sea mayor de 1. 
write.table(res_time_def,"differential_expressionTIME.txt",sep="\t") #Almacenamos datos expresión diferencial.


################################################################################
################################## DIAGRAMA DE VENN ############################
###################################EXTRACCION DATOS ############################
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

################################################################################
################################## DIAGRAMA DE VENN ############################
################################################################################

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


################################################################################
########################### ENRIQUECIMIENTO FUNCIONAL ##########################
################################################################################

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