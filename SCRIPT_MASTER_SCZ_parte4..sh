# SCRIPT - TFM - MASTER BIOINFO - 

# ESTUDIO DE SELECCION NATURAL RECIENTE EN ESQUIZOFRENIA - PARTE 4

# Script creado para el análisis de la funcionalidad de los genes que contienen variantes asociadas a SCZ y sometidas a selección positiva repecto a las que no lo están.
# Diciembre 2019 - Javier González-Peñas.
# ==============================================
# 
# Descripción -
# 
# Mediante este script se hacen las siguientes procesos:
# -  Se preparan los archivos de variantes para seleccionar aquellas asociadas a esquizofrenia y sometidas a selección, y se mapean lo genes implicados
# -  Se analiza en FUMA, STRING para ver sus propiedades funcionales (plataformas online)
# -  Además, se analizan los patrones de expresión en células específicas a partir de datos de RNAseq single cell en cerebro.




# ==============================================
# 1- PREPARACIÓN DE ARCHIVOS PARA ANALISIS DE GENE SET
# ==============================================

SCZ_filtered_ancest = read.table("SCZ_filtered_ancest.txt", header = T, row.names = NULL)

# Cojo todas las variantes y ordeno por p de PCAdapt
qsort <- SCZ_filtered_ancest[order(SCZ_filtered_ancest$pvalues),]

# AÑADIMOS INFORMACIÓN DE PERTENENCIA AL TOP5% DE SELECCIÓN O AL RESTO
qsort5=head(qsort,11760)
qsort95=qsort[11761:235910,]
qsort5$SEL = "sel5"
qsort95$SEL = "sel95"
temp=rbind(qsort5,qsort95)
SCZ_filtered_ancest <- temp

# Filtramos top5% y p_SCZ < 0,05
temp = subset(SCZ_filtered_ancest, p < 0.05 | SEL == "sel5")

# Dividimos en riesgo y protección
prot = subset(temp, OR_NEW < 1) # 1344 SNPs
risk = subset(temp, OR_NEW > 1) # 1292 SNPs

# Seleccionamos las variables necesarias para IMPUT en FUMA
myvars <- c("SNP.x", "p")
temp_prot <- prot[myvars]
temp_risk <- risk[myvars]
colnames(temp_prot) <- c("rsID", "Pval")
colnames(temp_risk) <- c("rsID", "Pval")

SCZ_filtered_ancest_SEL5_prot <- temp_prot
SCZ_filtered_ancest_SEL5_risk <- temp_risk

rm(temp, temp_prot, temp_risk, prot, risk)

# Exportamos la tabla con información de pertenencia al top 5% de selección
write.table(SCZ_filtered_ancest_SEL5_prot, 'SCZ_filtered_ancest_SEL5_prot.txt', row.names = FALSE, quote = FALSE)
write.table(SCZ_filtered_ancest_SEL5_risk, 'SCZ_filtered_ancest_SEL5_risk.txt', row.names = FALSE, quote = FALSE)




# ==============================================
# 2 - ANÁLISIS EN FUMA
# ==============================================

# las variantes de cada set se analizan en fuma para
# Mapeo posicional, interación de cromatina (HiC) o eQTL
# Análisis de expresión en cerebro y en resto de tejidos
# Análisis de funciones biológicas (GO)
# (Plataforma online)


# ==============================================
# 3 - ANÁLISIS EN STRING
# ==============================================

# Subimos las listas de genes de protección y de riesgo de los mapeados y analizamos interacción PPI (Plataforma online)




# ==============================================
# 4 - ANÁLISIS DE EXPRESIÓN EN CÉLULAS ESPECÍFICAS
# ==============================================

### instalamos el paquete pSI y leemos datos de expresión RNAseq
install.packages("pSI")
library(pSI)
help(pSI)
data2<-read.table("FPKM_ZHANG_NEW3",header=T,row.names=1, sep="\t")
head(data2)

# convertimos matrix con logaritmo2 y normalizamos
data2_log=log2(data2)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(data2_log)
data2_log_norm=range01(data2_log)
write.table(data2_log_norm, file = "data3_log_norm", sep = "\t")

# Seleccionamos la matriz ya preparada y corremos el pSI
data_norm<-read.table("data3_log_norm",header=T,row.names=1, sep="\t")
head(data_norm)
pSI=specificity.index(data_norm, bts = 50, p_max = 1, e_min = 0.3, hist = FALSE, SI = FALSE)

# Cargamos las listas de genes de protección o riesgo de las listas de selección
RISK = read.table("GENES_SEL5_RISK_list")
PROT = read.table("GENES_SEL5_PROT_list")
RISK2=t(RISK)
PROT2=t(PROT)

# Corremos el anñálisis de enriquecimiento, a distintos umbrales de especificidad
Out_RISK = fisher.iteration(pSI, RISK2, p.adjust = TRUE)
Out_PROT = fisher.iteration(pSI, PROT2, p.adjust = TRUE)

# EXPORTAMOS resultados
write.table(Out_RISK, file = "Out_RISK", sep = "\t")
write.table(Out_PROT, file = "Out_PROT", sep = "\t")


## FIN PARTE 4






