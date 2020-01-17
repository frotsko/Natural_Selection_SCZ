# SCRIPT - TFM - MASTER BIOINFO - 

# ESTUDIO DE SELECCION NATURAL RECIENTE EN ESQUIZOFRENIA - PARTE 2

# Script creado para el cálculo de la correlación entre señales de selección natural reciente y de asociación con SCZ
# Diciembre 2019 - Javier González-Peñas.
# ==============================================
# 
# Descripción -
# 
# Mediante este script se hacen las siguientes procesos:
# -  Se analiza la correlación entre la magnitud de selección natural reciente y la de asociación a SCZ a lo largo de las 
# variantes independientes obtenidas en el apartado anterior
# -  Se comprueba que esa correlación obtenida no sea artefactual mediante análisis de aquellas variantes más asociadas con 
# SCZ (p < 0.01), mediante un análisis de compraración con la deistribución de correlaciones tras 10000 permutaciones y 
# mediante enriquecimiento del top 5% de variantes sometidas a selección a lo largo de distintos umbrales de asociación con SCZ




# ==============================================
# 1-  PREPARADO DE ARCHIVOS Y CÁLCULO DE CORRELACIÓN  SCZ - SELECCIÓN
# ==============================================

## A) EXTRAEMOS DE LOS 1KG LA INFORMACIÓN ACERCA DE LA FRECUENCIA ALÉLICA DE LAS VARIANTES EN POBLACIÓN EUROPEA (LA UTILIZAREMOS EN LOS ANÁLISIS DE CORRELACIÓN)
# Extraemos del archivo completo de los 1000G la frecuencia en Europeos, porque la necesitaremos para corregir los análisis de correlación
# SacAMOS las variantes del VCF
zcat ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | sed -e '1,252d' > 1000_Variants

# Extraigo la información y separo campos de frecuencia
awk 'NR==FNR {FILE1[$1]=$0; next} ($3 in FILE1) {print $0}' SNP_SCZ_rsID_pcadapt 1000_Variants > 1KG_SCZ
cat 1KG_SCZ | awk -v OFS="\t" '{print $3, $8}' | awk -F";" '$1=$1' OFS="\t" | awk -v OFS="\t" '{print $1, $10}' | awk -F"=" '$1=$1' OFS="\t" | awk -v OFS="\t" '{print $1, $3}' > 1KG_SCZ_EUR_MAF 

Rstudio
library(ppcor)
# Cargo archivo de MAF
MAF_SCZ = read.table("1KG_SCZ_EUR_MAF", header = F, row.names = NULL)
colnames(MAF_SCZ) <- c("SNP.x", "FREQ")
MAF_SCZ$MAF = ifelse(MAF_SCZ$FREQ < 0.5, MAF_SCZ$FREQ, 1-MAF_SCZ$FREQ)


## B) CARGAMOS LA INFORMACIÓN DE SELCCIÓN POR VARIANTE CALCULADA EN EL PASO ANTERIOR
# Cargo archivo de variantes y le añado el P valor de selección
rsID_SCZ = read.table("VCF_pcadapt_SCZ.bim", header = F, row.names = NULL)
pvalues = read.table("pvalues.txt", header = F, row.names = NULL)

head(rsID_SCZ)
temp = cbind(rsID_SCZ,pvalues)
temp$POS = paste(temp$V1, temp$V4, sep=":")
myvars <- c("POS", "V2", "pvalues")
temp2 <- temp[myvars]
colnames(temp2) <- c("POS", "SNP", "pvalues")


## C) CARGAMOS LA INFORMACIÓN DE ASOCIACIÓN A SCZ
# Cargo archivo de asociacion GWAS
SCZ = read.table("overlap_final_SCZ", header = T, row.names = NULL)
SCZ$POS = paste(SCZ$chr, SCZ$bp, sep=":")
SCZ$POS<-sub("chr", "", SCZ$POS)


## D) UNIMOS LAS TRES FUENTES DE INFORMACIÓN Y CALCULAMOS LA CORRELACIÓN (dependiente de frecuencia)
# HAgo unión de variante, selección, asociacion y frecuencia
SCZ_filtered1 = merge(SCZ,temp2, by = "POS")
SCZ_filtered = merge(SCZ_filtered1,MAF_SCZ, by = "SNP.x")
rm(SCZ_filtered1)

# calculo correlacion (Debemos coregir por frecuencia - análisis de correlación condicionada a frecuencia)
library(ppcor)
SCZ_sel = cor.test( ~ SCZ_filtered$p + SCZ_filtered$pvalues, data=SCZ_filtered, method = "spearman", continuity = FALSE, conf.level = 0.95)
SCZ_sel_corr = pcor.test(SCZ_filtered$p, SCZ_filtered$pvalues, SCZ_filtered$MAF, method = "spearman")

write.table(SCZ_sel, 'SCZ_sel.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(SCZ_sel_corr, 'SCZ_sel_corr.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(SCZ_filtered, 'SCZ_filtered.txt', row.names = FALSE, quote = FALSE)




# ==============================================
# 2 -  Análisis de los valores de correlación en el subset de variantes más asociadas (P_SCZ < 0.01; N_SNPs = 17121)
# ==============================================

dataPC_001=SCZ_filtered[SCZ_filtered$p < 0.01, ]
SCZ_MAF_PSCZ_01=pcor.test(dataPC_001$p, dataPC_001$pvalues, dataPC_001$MAF, method = "spearman")
write.table(SCZ_MAF_PSCZ_01, 'SCZ_sel_corr_001.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)



# ==============================================
# 3 -  Análisis por permutaciones. Comprobar si los valores de correlación en ambos casos (TOTAL de variantes y P_SCZ < 0.01) son superiores a la distribucion de 10.000 permutaciones al azar
# ==============================================

## Definimos función
cor_fun <- function(x, y)
	{
 	 cor <- cor.test(x, y, method = "spearman", continuity = FALSE, conf.level = 0.95)
 	 o=cbind(r=cor$estimate)
 	 return(o)
	}


## A) PERMUTACIONES CON TODAS LAS VARIANTES
## Definimos parámetros para la permutación
n <- length(SCZ_filtered$p)                 # tamaño de la muestra
n.simul <- 10000                            # numero de simulaciones
# Inicializo matrices que guardan por columna el resultado de cada remuestreado
b1.simul <- matrix(0, ncol=n.simul, nrow=n)
# Inicializo vector que guarda  cada experimento
pend <- numeric(n.simul)

set.seed(200785)                            # fijo semilla

## Inicio permutación simple (sin reemplazamiento)
for(i in 1:n.simul) 
	{
  	b1.simul[,i] <- sample(SCZ_filtered$p,n,replace=FALSE) 
	}

for(i in 1:n.simul) 
	{
  	pend[i] <- cor_fun (b1.simul[,i],SCZ_filtered$pvalues)
	}

## Exportamos distribución de correlaciones y calculamos P valor de permutación
write.table(pend, "Permutation_ALL")
hist(pend, border = "black", col = "grey", main = "Correlation distribution", xlab = "P_SCZ - P_PCAdapt Spearman correlation", xlim=c(-0.06,0.06), las = 3, breaks = 20)
yve = 0.049
abline (v = yve, col = 'red') 
P <- sum (pend >= yve)/(n.simul + 1) # calculate p value from permutaions
P


## B) PERMUTACIONES CON LAS VARIANTES P_SCZ < 0.01
## Definimos parámetros para la permutación
n <- length(dataPC_001$p)                   # tamaño de la muestra
n.simul <- 10000                            # numero de simulaciones
# Inicializo matrices que guardan por columna el resultado de cada remuestreado
b1.simul <- matrix(0, ncol=n.simul, nrow=n)
# Inicializo vector que guarda  cada experimento
pend <- numeric(n.simul)

set.seed(200785)                            # fijo semilla

## Inicio permutación simple (sin reemplazamiento)
for(i in 1:n.simul) 
	{
  	b1.simul[,i] <- sample(dataPC_001$p,n,replace=FALSE) 
	}

for(i in 1:n.simul) 
	{
  	pend[i] <- cor_fun (b1.simul[,i],dataPC_001$pvalues)
	}

## Exportamos distribución de correlaciones y calculamos P valor de permutación
write.table(pend, "Permutation_P_SCZ_001")
hist(pend, border = "black", col = "grey", main = "Correlation distribution", xlab = "P_SCZ - P_PCAdapt Spearman correlation (P_SCZ < 0.01)", xlim=c(-0.06,0.06), las = 3, breaks = 20)
yve = 0.024
abline (v = yve, col = 'red') 
P <- sum (pend >= yve)/(n.simul + 1) # calculate p value from permutaions
P



# ==============================================
# 4 -  Análisis de enriquecimiento del top 5% de variantes sometidas a selección sobre las variantes asociadas con SCZ
# ==============================================

# Cojo todas las variantes y ordeno por p de PCAdapt
qsort <- SCZ_filtered[order(SCZ_filtered$pvalues),]

# Cojo solo los que tienen p de PCA perteneciente al 5% top
qsort5=head(qsort,11760)
qsort95=qsort[11761:235910,]


thres=c(0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001)

dataRanges_asoc <- c() 
for (i in thres) 
	{
  	SCZ_P <- subset(qsort, p < i)   

 	 SCZ_P_seleccion_5 <- subset(qsort5, p < i)
  	No_SCZ_P_seleccion_5 <- subset(qsort5, p > i)
  
  	SCZ_P_No_seleccion_95 <- subset(qsort95, p < i)
  	No_SCZ_P_No_seleccion_95 <- subset(qsort95, p > i)

  	prop_SCZ_sel=dim(SCZ_P_seleccion_5)[1]/(dim(SCZ_P_seleccion_5)[1] + dim(SCZ_P_No_seleccion_95)[1])
  	prop_No_SCZ_sel=dim(No_SCZ_P_seleccion_5)[1]/(dim(No_SCZ_P_seleccion_5)[1] + dim(No_SCZ_P_No_seleccion_95)[1])
  	
  	x <- matrix(c(dim(SCZ_P_seleccion_5)[1], dim(No_SCZ_P_seleccion_5)[1], dim(SCZ_P_No_seleccion_95)[1], dim(No_SCZ_P_No_seleccion_95)[1]), byrow = TRUE, 2, 2)
  
  	df <- data.frame(Range = paste0( "P_SCZ_Thres_",i), prop_SCZ_sel = prop_SCZ_sel, N_SCZ_sel = dim(SCZ_P_seleccion_5)[1], N_SCZ_sel = dim(SCZ_P_seleccion_5)[1], prop_No_SCZ_sel = prop_No_SCZ_sel, N_No_SCZ_sel = dim(No_SCZ_P_seleccion_5)[1], Chisq_pval = chisq.test(x)$p.value)
  	dataRanges_asoc <- rbind(dataRanges_asoc, df)
	}

dataRanges_asoc
write.table(dataRanges_asoc, "PROP_seleccion_SCZ.txt", sep = "\t", dec = ",")


## CONSTRUíMOS TABLA CON RESULTADOS PREVIOS Y HACEMOS GRÁFICA
library(ggplot2)
medias <- read.table("PROP_seleccion_SCZ.txt",header=T, sep="\t", dec = ",")

level_order <- factor(medias$P_treshold, level = c('P < 0.5', 'P < 0.2', 'P < 0.1','P < 0.05', 'P < 0.01', 'P < 0.001','P < 1e-4', 'P < 1e-5', 'P < 1e-6','P < 1e-7', 'P < 5e-8'))
g <- ggplot(data=medias,aes(x=level_order,y=OR,label = label)) + geom_errorbar(aes(ymin=Cinf,ymax=Csup),width=0.1, size=1, stat = "identity") + labs(x="SCZ significance threshold", y = "OR") + theme_classic() + theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=12), strip.text.x = element_text(size = 12))+ geom_hline(aes(yintercept = 1), linetype="dashed") + ylim(0, 3.3)
+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(), panel.spacing = unit(0.8, "lines")) + geom_errorbar(aes(ymin=Cinf,ymax=Csup),width=0.2, size=1) + theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=12), strip.text.x = element_text(size = 12)) + geom_point(aes(shape=SAMPLE), size= 2.5) +  guides(shape=FALSE) + geom_text(data=medias, aes(x=SAMPLE, y=Csup, label=label), col='black', size=2.7, hjust=0.45 ,vjust=-0.7) + scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) + theme(panel.background = element_blank(),axis.line = element_line(colour = "black"), axis.ticks.y = element_line(colour = "black")) + geom_hline(aes(yintercept = 0), linetype="dashed")
g




### FINAL parte 2




