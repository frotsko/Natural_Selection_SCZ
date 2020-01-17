# SCRIPT - TFM - MASTER BIOINFO - 

# ESTUDIO DE SELECCION NATURAL RECIENTE EN ESQUIZOFRENIA - PARTE 5

# Script creado para el análisis del riesgo poligénico de las variantes sometidas a selección natural reciente respecto a las que no lo están sobre una muestra independiente
# Diciembre 2019 - Javier González-Peñas.
# ==============================================
# 
# Descripción -
# 
# Mediante este script se hacen las siguientes procesos:
# - Se preparan los datos de genotipado de la muestra descubrimiento (GWAS SCZ 2014) y de la muestra diana (cohorte independiente CIBERSAM)
# - Dividimos los datos de asociación en distintos umbrales de selección natural reciente (20 intervalos del 5% de variantes
# - Se compara el valor de la varanaza explicada por el PRS de SCZ en el umbral de mayor selección natural reciente frente al resto
# - Se compara este valor de varianza frente a 1000 permutaciones del mismo número de variantes, pero fuera de este umbral de selección.





# ==============================================
# 1 - PREPARADO DE MUESTRA DIANA
# ==============================================

## Eliminamos variantes en complejo mayor de histocompatibilidad (MHC)
rstudio
system2("./plink", args=c("--bfile PGC_imput_QC_filtered", "--exclude range MHC_filter.txt", "--make-bed", "--out tempMHC"))
dataTarget <- read.table('tempMHC.bim', header = FALSE, stringsAsFactors = FALSE)





# ==============================================
# 2 - LECTURA DE DATOS DE GWAS DE SCZ, DIVISIÓN ENTRE UMBRALES DE SELECCIÓN Y CÁLCULO DE PRS
# ==============================================

dataDiscovery_filtered <- read.table('SCZ_filtered.txt', header = TRUE, stringsAsFactors = FALSE)

# Ordenamos por p valor de selecciónm natural
df <- data.frame(dataDiscovery_filtered[order(dataDiscovery_filtered$pvalues),], stringsAsFactors = FALSE)
row.names(df) <- c()
N = 20 # definimos número de grupos.
varsperGroup <- floor(nrow(df)/N)

## calculamos el componente poligénico (PRS) sobre una muestra independiente (PGC_new_filtered) para 
for (i in 1:N)
	{
  	pos_i <- 1 + varsperGroup*(i-1)
  	pos_f <- varsperGroup*i
  	if (i == N)
  		{
    	pos_f <- nrow(df)
  		}
  	myData <- df[pos_i:pos_f,]
  
  	# SALVAMOS OR Y P VALOR SCZ
  	SCZ_SCORE <- myData[,c("POS", "a1", "or")]
  	SCZ_P <- myData[,c("POS", "p")]
  
  	write.table(SCZ_SCORE, paste0("SCZ_SCORE_",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  	write.table(SCZ_P, paste0("SCZ_P_", i, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  	# CALCULAMOS PRS
  	system(paste0("./plink --bfile tempMHC --score SCZ_SCORE_", i, ".txt --q-score-file SCZ_P_", i, ".txt --q-score-range Q_RANGES.txt --out SCZ_PRS_", i, "_SCORE"))
	}




# ==============================================
# 3 - GUARDAMOS LOS ARCHIVOS DE RIESGO POLIGÉNICO PARA CADA UMBNRAL EN UNA CARPETA DIFERENTE Y EXPORTAMOS A EXCEL
# ==============================================

# CREAMOS CARPETA PARA CADA 5% DE VARAINTES DE SELECCIÓN EN ORDEN DECRECIENTE
library(data.table)
library(xlsx)
folderNames <- c("0-5", "5-10", "10-15", "15-20","20-25", "25-30","30-35", "35-40","40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-85", "85-90","90-95", "95-100")
f_all <- 'SCZ_results_PRS'

# HACEMOS UNA LISTA CON LOS ARCHIVOS DE RIESGO CALCULADO EN CADA CASO (.PROFILE)
dataRange <- list()
for (i in 1:length(folderNames))
	{
	filesRange <- list.files(path = folderNames[i], pattern = ".profile", full.names = TRUE)
	dataGroup <- list()
  
	for (j in 1:length(filesRange))
		{
    	fin <- read.table(filesRange[j], header = TRUE)
    	scoreData <- data.table(fin$FID, fin$IID, fin$SCORE)
    	colnames(scoreData) <- c("FID", "IID", paste0("SCORE_", j-1))
    	setkey(scoreData, FID, IID)
    	dataGroup[[j]] <- scoreData
  		}
	dataRange[[i]] <- Reduce(function(x, y) merge(x, y, all=TRUE), list(dataGroup[[1]], dataGroup[[2]], dataGroup[[3]], dataGroup[[4]], dataGroup[[5]], dataGroup[[6]], dataGroup[[7]],dataGroup[[8]], dataGroup[[9]], dataGroup[[10]], dataGroup[[11]], dataGroup[[12]], dataGroup[[13]], dataGroup[[14]],dataGroup[[15]], dataGroup[[16]], dataGroup[[17]], dataGroup[[18]], dataGroup[[19]], dataGroup[[20]]))
	names(dataRange)[[i]]  <- paste0("range", i)
	}

for (i in 1:length(dataRange))
	{
  	# write.xlsx(dataRange[[i]], f_all, sheetName = names(dataRange)[[i]], col.names=TRUE, row.names = FALSE, append=TRUE)
  
  	write.xlsx(dataRange[[i]], paste0(f_all, i, '.xlsx'), sheetName = names(dataRange)[[i]], col.names=TRUE, row.names = FALSE, append=FALSE)
	}





# ==============================================
# 4 - CÁLCULO DE LA VARIANZA EXPLICADA DEL PRS DE SCZ ENTRE CASOS Y CONTROLES,  MEDIANTE MODELO DE REGRESIÓN LOGÍSTICA
# ==============================================

library(rms)

# CARGAMOS LA INFO DE LOS SUJETOIS DE LA MUESTRA (COVARIABLES DEL ANÁLISIS)
familyInfo <- read.table("PGC_new_filtered.fam", header=FALSE)
familyInfo <- data.frame (FID= familyInfo$V1, IID = familyInfo$V2, scode= familyInfo$V5, dcode = familyInfo$V6)
# scode = 1 (Males), scode = 2 (Females)
# dcode = 1 (Case),  dcode = 2 (Control), dcode = -9 (No Diagnosis Info)

familyInfo <- familyInfo[familyInfo$dcode %in% c(1, 2),] # Los pheno -9 están quitados. 

MDS       <- read.table("MDS.mds", header=T, stringsAsFactors = FALSE)
covars    <- merge(MDS, familyInfo, by = c("IID", "FID"))

folderNames <- c("0-5", "5-10", "10-15", "15-20","20-25", "25-30","30-35", "35-40","40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-85", "85-90","90-95", "95-100")

dataRange <- list()
i <- j <- 1
for (i in 1:length(folderNames)) 
	{
  	filesRange <- list.files(path = folderNames[i], pattern = ".profile", full.names = TRUE)
  	dataGroup <- c()
  
  	for (j in 1:length(filesRange)) 
  		{
    	fin <- read.table(filesRange[j], header = TRUE, stringsAsFactors =  FALSE)
    	variables_SCZ <- merge(fin, covars, by = c("IID", "FID"))
    
    	#Porcentaje de missing genotypes
    	N.snps <- summary(variables_SCZ$CNT)
    	variables_SCZ$percent_missing <- (N.snps[6] - variables_SCZ$CNT2) / N.snps[6] # N.snps[6] is Max number
    
    	#Variables_SCZ para regresión logística (check names!)
    	pheno <- as.numeric(variables_SCZ$PHENO)
    	score <- variables_SCZ$SCORE
   		miss  <- variables_SCZ$percent_missing
    
    	MDS1  <- variables_SCZ$C1
    	MDS2  <- variables_SCZ$C2
    	MDS3  <- variables_SCZ$C3
    	MDS4  <- variables_SCZ$C4
    	MDS5  <- variables_SCZ$C5
    	MDS6  <- variables_SCZ$C6
    	MDS7  <- variables_SCZ$C7
    	MDS8  <- variables_SCZ$C8
    	MDS9  <- variables_SCZ$C9
    	MDS10  <- variables_SCZ$C10
    	SEX   <- variables_SCZ$scode
    
    	# se calcula la varianza mediante la estimación de la R2 de Naggelkerke (diferencia de explicación entre un modelo solo con covariables y uno completo)
    	H0 <- lrm(pheno ~ miss + MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX)
   		print(H0)
   		H1 <- lrm(pheno ~ miss + MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX + score, scale = TRUE)
   		print(H1)
   		DELTA_R2 <- H1$stats[10] - H0$stats[10]
    	R2 <- DELTA_R2 * 100
    	print(R2)
    	DELTA_AUC <- H1$stats[6] - H0$stats[6]
    	AUC <- DELTA_AUC * 100
    	print(AUC)
    	
    	#Calculo del P valor del modelo
    	phenoN <- pheno - 1 #Coded for glm
    	model <- glm(phenoN ~ miss + MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX + score, family = binomial())
    	P <- summary(model)
    	print(P)
    
    	# Salvamos las varaibles de la regresión
    	df <- data.frame(R2 = R2, AUC = AUC, P = P$coefficients[14,4], Range = paste0( "Q_RANGE_S", j - 1, "_INTERVAL_", folderNames[i]))
    
    	dataGroup <- rbind(dataGroup, df)
  		}
  	dataRange[[i]] <- dataGroup
  	names(dataRange)[[i]]  <- paste0("Range_", folderNames[i])
	}

final <- rbind(dataRange[[1]], dataRange[[2]], dataRange[[3]], dataRange[[4]], dataRange[[5]], dataRange[[6]], dataRange[[7]], dataRange[[8]], dataRange[[9]], dataRange[[10]],dataRange[[11]], dataRange[[12]], dataRange[[13]], dataRange[[14]], dataRange[[15]], dataRange[[16]], dataRange[[17]], dataRange[[18]], dataRange[[19]], dataRange[[20]])

rownames(final) <- c() 
write.table(final, "output_PRS_selection_SCZ.txt", dec = ".", row.names = FALSE )





# ==============================================
# 5 - SELECCIONAMOS EL UMBRAL DE MAYOR EXPLICACIÓN EN SCZ (P < 0.05) Y REPRESENTAMOS UN HISTOGRAMA DE LAS R2 PARA CADA INTERVALO DE SELECCIÓN NATURAL RECIENTE
# ==============================================

table=read.table("output_PRS_selection_SCZ.txt", header = T, sep = "\t", dec = ",")
table3 <- table[grep("Q_RANGE_S4", table$Range), ]

head(table3)
library(ggplot2)
library(grid)
table3$Range <- as.character(table3$Range)
table3$Range <- factor(table3$Range, levels=unique(table3$Range))

library(RColorBrewer)

nb.cols <- 50
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(nb.cols)

ggplot(table3, aes(x=Range, y=R2, fill = Range)) + 
  geom_bar(stat="identity") +theme_minimal() + scale_fill_manual(values = mycolors)+
  theme_classic() +   labs(title="SCZ PGS in decreasing selection intervals ",  x="5% PCAdapt intervals", y = "pseudo-"~R^2) + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + theme(legend.position = "none")


## Analizamos la tendencia de explicación en los distintos umbrales de selecció natural

# Hacemos la regresión entre los R2 de los PGS y el ordinal de los veintiles
temp <- read.table("lrm_SCZ_20.txt", header = TRUE, stringsAsFactors = FALSE, dec=".")
model = lm(interval ~ R2, data = temp)
model$coefficients[2]

# Representamos recta
library(ggplot2)
ggplot(temp, aes(x=interval, y=R2)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm)


## Definimos parámetros para la permutación
n <- length(temp$R2)                 # tamaño de la muestra
n.simul <- 10000                            # numero de simulaciones
# Inicializo matrices que guardan por columna el resultado de cada remuestreado
b1.simul <- matrix(0, ncol=n.simul, nrow=n)
# Inicializo vector que guarda  cada experimento
pend <- numeric(n.simul)

set.seed(200785)                            # fijo semilla

## Inicio permutación simple (sin reemplazamiento)
for(i in 1:n.simul) 
  {
    b1.simul[,i] <- sample(temp$R2,n,replace=FALSE) 
  }

for(i in 1:n.simul) 
  {
    model = lm(temp$interval ~ b1.simul[,i])

    pend[i] = model$coefficients[2]
  }

## Exportamos distribución de correlaciones y calculamos P valor de permutación
write.table(pend, "Permutation_R2")
hist(pend, border = "black", col = "grey", main = "r2 real Vs 1K permutaciones", xlab = "r2 modelo regresión", xlim=c(-8,8), las = 3, breaks = 20)
yve = -6.57
abline (v = yve, col = 'red') 
P <- sum (pend <= yve)/(n.simul + 1) # calculate p value from permutaions
P



# ==============================================
# 6 - ANÁLISIS DE PERMUTACIONES PARA COMPARAR LA EXPLICACIÓN DE LAS VARIANTES TOP 5% DE SELECCIÓN CON EL MISMO NÚMERO DE VARIANTES DEL RESTO TOTAL.
# ==============================================

## Definimos número de permutaciones y set de 5% de selección y 95% restante
kNumsPerm = 1000 # número de permutaciones
dataDiscovery_filtered <- read.table("SCZ_filtered.txt", header = TRUE, stringsAsFactors = FALSE, dec=".")
df <- data.frame(dataDiscovery_filtered[order(dataDiscovery_filtered$pvalues),], stringsAsFactors = FALSE)
row.names(df) <- c()

vars5 <- floor(nrow(df)*0.05)
data_5 <- df[1:vars5,]
vars95 <- vars5+1
data_95 <- df[vars95:nrow(df),]
write.table(data_5, "SCZ_filtered_5_percent.txt", row.names = FALSE, col.names = TRUE)
write.table(data_95, "SCZ_filtered_95_percent.txt", row.names = FALSE, col.names = TRUE)

## Cálculo de PRS de mismo número de variantes de SEL 5% pero del grupo del 95% restante. calculamos knumsPerms PRS
data_95 <- read.table("SCZ_filtered_95_percent.txt", header = TRUE, stringsAsFactors = FALSE, dec=".")

for(k in 1:kNumsPerm)
  {
    sel <- sample(1:nrow(data_95), vars5)
    sel <- sel[sort(sel)]
    myData <- data_95[sel,]
  
    SCZ_SCORE <- myData[,c("POS", "a1", "or")]
    SCZ_SCORE$OR_EUR <- log(myData$OR_EUR)
    SCZ_P <- myData[,c("POS", "p")]
  
    write.table(SCZ_SCORE, "SCZ_SCORE_95.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(SCZ_P, "SCZ_P_95.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
    system("./plink --bfile PGC_new_filtered --score SCZ_SCORE_95.txt --q-score-file SCZ_P_95.txt --q-score-range Q_RANGES_S4.txt --out SCZ_PRS_SCORE_95")
  
    fin <- read.table("SCZ_PRS_SCORE_95.S4.profile", header = TRUE, stringsAsFactors =  FALSE)
    score <- data.frame(FID= fin$FID, IID= fin$IID, fin$SCORE)
    colnames(score)[3] <- paste0("Score95_", k)
    variables_SCZ <- merge(variables_SCZ, score, by = c("IID", "FID"))
    write.table(variables_SCZ, "data_SCZ_permutations.txt", dec= ".", row.names = FALSE)
  }

variables_SCZ <- read.table("data_SCZ_permutations.txt", header=TRUE, stringsAsFactors = FALSE)

## Variables_SCZ para regresión logística (check names!)
pheno <- as.numeric(variables_SCZ$PHENO)
MDS1  <- variables_SCZ$C1
MDS2  <- variables_SCZ$C2
MDS3  <- variables_SCZ$C3
MDS4  <- variables_SCZ$C4
MDS5  <- variables_SCZ$C5
MDS6  <- variables_SCZ$C6
MDS7  <- variables_SCZ$C7
MDS8  <- variables_SCZ$C8
MDS9  <- variables_SCZ$C9
MDS10 <- variables_SCZ$C10
SEX   <- variables_SCZ$scode
AGE   <- variables_SCZ$AGE

H0 <- lrm(PHENO ~ MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX + AGE, scale = TRUE)
# H0 <- lrm(pheno ~ MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX , scale = TRUE)
print(H0)

dataRanges <- c() 
for (j in 1:kNumsPerm+1) 
  {
    score <- variables_SCZ[,16+j]
    H1 <- lrm(PHENO ~ MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX + AGE + score, scale = TRUE)
    # H1 <- lrm(pheno ~ MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX  + score, scale = TRUE)
    # print(H1)
    DELTA_R2 <- H1$stats[10] - H0$stats[10]
    R2 <- DELTA_R2 * 100
    # print(R2)
    DELTA_AUC <- H1$stats[6] - H0$stats[6]
    AUC <- DELTA_AUC * 100
    # print(AUC)
  
    phenoN <- PHENO - 1 #Coded for glm
    model <- glm(phenoN ~ MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX + AGE + score, family = binomial())
    # model <- glm(phenoN ~ MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8 + MDS9 + MDS10 + SEX + score, family = binomial())
  
    P <- summary(model)
    # print(P)
    # Save the variables of the logistic regression
    df <- data.frame(R2 = R2, AUC = AUC, P = P$coefficients[13,4], Range = colnames(variables_SCZ)[16+j])
    dataRanges <- rbind(dataRanges, df)
  }

write.table(dataRanges, "data_SCZ_lrm_permutations.txt", dec= ".", row.names = FALSE)
rownames(dataRanges) <- c()


## Printamos el gráfico de permutaciones y calculamos el P valor de la comparación entre el valor de explicación delñ 5% frente a la distribución.

perm <- read.table("data_SCZ_lrm_permutations.txt", header = TRUE)
hist(perm$R2, border = "black", col = "grey", main = "SCZ_PRS (real Vs 1K permutations)", xlab = "SCZ_PRS", xlim=c(0,5), las = 3, breaks = 20)
yve = 3.69
abline (v = yve, col = 'red') 
P <- sum (perm$R2 >= yve)/(1000 + 1) 
P


## FIN PARTE 5





