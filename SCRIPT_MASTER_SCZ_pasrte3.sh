# SCRIPT - TFM - MASTER BIOINFO - 

# ESTUDIO DE SELECCION NATURAL RECIENTE EN ESQUIZOFRENIA - PARTE 3

# Script creado para el análisis de la direccionalidad de la selección natural reciente asociada con SCZ (variantes de riesgo o protección)
# Diciembre 2019 - Javier González-Peñas.
# ==============================================
# 
# Descripción -
# 
# Mediante este script se hacen las siguientes procesos:
# -  Se identifica el alelo ancestral y el derivado de todas las variantes analizadas a partir de la información de 1KG
# -  Se analiza la tendencia de que el alelo que surge en la población (derivado) y que está sometido a selección natural reciente
# tenga tendencia a ser de riesgo (selección positiva de variantes de SCZ) o de protección (selección positiva de variantes de protección)



# ==============================================
# 1- DESCRIPCIÓN DE ALELO ANCESTRAL Y ALELO DERIVADO, Y CÁLCULO DE OR DE SCZ PARA EL ALELO DERIVADO
# ==============================================

## A partir de la info de 1KG obtenemos el alelo ancedtral de cada variante
cat 1KG_SCZ | awk -v OFS="\t" '{print $3, $8}' | awk -F";" '$1=$1' OFS="\t" | awk -v OFS="\t" '{print $1, $12}' | sed 's/|||//g' |  awk -F"=" '$1=$1' OFS="\t" | awk -v OFS="\t" '{print $1, $3}' > 1KG_SCZ_ancestral 

# En R, solapamos la oinformacion con la de selección y asociación (SCZ_filtered)
rstudio

SCZ_ancest = read.table("1KG_SCZ_ancestral", header = F, row.names = NULL)
colnames(SCZ_ancest) <- c("SNP.x", "Al_ancest")

SCZ_filtered = read.table("SCZ_filtered.txt", header = T, row.names = NULL)
SCZ_filtered_ancest = merge(SCZ_filtered,SCZ_ancest, by = "SNP.x")
temp<-subset(SCZ_filtered_ancest, (Al_ancest %in% c("A", "C", "G", "T", "a", "c", "g", "t")))
head(temp)
SCZ_filtered_ancest <- temp
rm(temp)
## De las 235911 hay 231244 variantes con información de ancestralidad

# Como hay alelos en minuscula, los convierto a mayúscula
require(tidyverse)
temp = SCZ_filtered_ancest %>% mutate(Al_ancest = toupper(Al_ancest))
head(temp)
unique(temp$Al_ancest)
SCZ_filtered_ancest <- temp

## Adapto el valor de OR para el alelo derivado (el que aparece en la población y no es ancestral)
SCZ_filtered_ancest$OR_NEW <- ifelse(SCZ_filtered_ancest$a2 == SCZ_filtered_ancest$Al_ancest, SCZ_filtered_ancest$or, 1/(SCZ_filtered_ancest$or))
write.table(SCZ_filtered_ancest, 'SCZ_filtered_ancest.txt', row.names = FALSE, quote = FALSE)




# ==============================================
# 2 - ANÁLISIS DE LA TENDENCIA DE RIESGO PROTECCIÓN DE LOS ALELOS DERIVADOS SOMETIDOS A SELECCIÓN NATURAL RECIENTE
# ==============================================

SCZ_filtered_ancest = read.table("SCZ_filtered_ancest.txt", header = F, row.names = NULL)

# Cojo todas las variantes y ordeno por p de PCAdapt
qsort <- SCZ_filtered_ancest[order(SCZ_filtered_ancest$pvalues),]

# AÑADO INFORMACIÓN DE PERTENENCIA AL TOP5% DE SELECCIÓN O AL RESTO
qsort5=head(qsort,11760)
qsort95=qsort[11761:235910,]
qsort5$SEL = "sel5"
qsort95$SEL = "sel95"
temp=rbind(qsort5,qsort95)
SCZ_filtered_ancest <- temp


# HAGO UN LOOP PARA ANALIZAR ENRIQUECIMIENTO POR UMBRAL DE ASOCIACIÓN CON SCZ (test t)
thres=c(0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000005)

dataRanges_PCAdapt <- c() 
for (i in thres) {
  table_sel5 = subset(SCZ_filtered_ancest, SEL == "sel5")
  table_sel95 = subset(SCZ_filtered_ancest, SEL == "sel95")
  table_thres = subset(SCZ_filtered_ancest, p < i)
  
  table_SCZ_5 = subset(table_sel5, p < i)
  t_5 = t.test(table_SCZ_5$OR_NEW, conf.level = 0.95, mu = 1, paired = FALSE)
  table_SCZ_95 = subset(table_sel95, p < i)
  t_95 = t.test(table_SCZ_95$OR_NEW, conf.level = 0.95, mu = 1, paired = FALSE)
  t_twosamples = t.test(table_thres$OR_NEW ~ table_thres$SEL, conf.level = 0.95)
  df <- data.frame(P_sel5 = t_5$p.value, est_5 = t_5$estimate, CI_5_inf = t_5$conf.int[1], CI_5_sup = t_5$conf.int[2], P_sel95 = t_95$p.value, est_95 = t_95$estimate, CI_95_inf = t_95$conf.int[1], CI_95_sup = t_95$conf.int[2], P_compar = t_twosamples$p.value, Range = paste0( "", i))
  dataRanges_PCAdapt <- rbind(dataRanges_PCAdapt, df)
}
dataRanges_PCAdapt
write.table(dataRanges_PCAdapt, "dataRanges_PCA_105_095", sep = "\t", dec = ",")


# HAGO UN LOOP PARA ANALIZAR ENRIQUECIMIENTO POR UMBRAL DE ASOCIACIÓN CON SCZ (test de proporciones)
thres=c(0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000005)

dataRanges_PCA <- c() 
for (i in thres) {
  table_sel5 = subset(SCZ_filtered_ancest, SEL == "sel5")
  table_sel95 = subset(SCZ_filtered_ancest, SEL == "sel95")
  table_thres = subset(SCZ_filtered_ancest, p < i)
  table_NOthres = subset(SCZ_filtered_ancest, p > i)
  
  table_5_prot = subset(table_sel5, p < i & OR_NEW < 1)
  table_5_risk = subset(table_sel5, p < i & OR_NEW > 1)
  prop_5=dim(table_5_risk)[1]/(dim(table_5_prot)[1] + dim(table_5_risk)[1])
  
  table_95_prot = subset(table_sel95, p < i & OR_NEW < 1)
  table_95_risk = subset(table_sel95, p < i & OR_NEW > 1)
  prop_95=dim(table_95_risk)[1]/(dim(table_95_prot)[1] + dim(table_95_risk)[1])
  
  x <- matrix(c(dim(table_5_risk)[1], dim(table_95_risk)[1], dim(table_5_prot)[1], dim(table_95_prot)[1]), byrow = TRUE, 2, 2)
  
  df <- data.frame(Range = paste0( "", i), prop_5_risk = prop_5, prop_95_risk = prop_95, Chisq_pval = chisq.test(x)$p.value)
  dataRanges_PCA <- rbind(dataRanges_PCA, df)
}
dataRanges_PCA
write.table(dataRanges_PCA, "dataRanges_PROP_PCA_105_095", sep = "\t", dec = ",")



# ==============================================
# 3 - CON LA INFORMACIÓN DE PROPORCIÓ DE VARIANTES DE RIESGO, CONTRUYO GRÁFICOS DE LOS TEST EN CADA UMBRAL de SCZ
# ==============================================

table2=read.table("table_ttest", header = T, sep = "\t", dec = ",")
table3=read.table("table_prop", header = T, sep = "\t", dec = ",")
head(table2)
head(table3)
library(ggplot2)
library(grid)

table2$Range <- as.character(table2$Range)
table2$Range <- factor(table2$Range, levels=unique(table2$Range))
table3$Range <- as.character(table3$Range)
table3$Range <- factor(table3$Range, levels=unique(table3$Range))

p1<- ggplot(table2, aes(x=Range, y=est, group=VAR, color=VAR)) + geom_point() +
  geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2, position=position_dodge(0.05)) + 
   theme_classic() + scale_color_manual(values=c('#D55E00','#0072B2')) +  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.8) +
  labs(title="",  x="P_SCZ threshold", y = "OR_SCZ") + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + 
  theme(axis.text.x = element_text(size=8.5, angle=45, hjust = 1), axis.text.y = element_text(color="black", size=8, angle=0)) + 
  theme(legend.title = element_blank())

table3$prop_PROT = 1-table3$prop_risk
head(table3)

p2<- ggplot(data=table3, aes(x=Range, y=prop_risk, fill=VAR)) +
  geom_bar(width=0.6,stat="identity", position=position_dodge(width=0.6), color="black") + theme_minimal() + 
  scale_fill_manual(values=c('#D55E00','#0072B2' )) + geom_hline(yintercept=0.5, linetype="dashed", color = "black", size=0.8) +
  labs(title="",  x="P_SCZ threshold", y = "Proportion of risk alleles (OR > 1)") + theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + 
  theme(axis.text.x = element_text(size=9, angle=45, hjust = 0.8), axis.text.y = element_text(color="black", size=9, angle=0)) + 
  theme(legend.title = element_blank())

p1
p2


### FINAL parte 3

