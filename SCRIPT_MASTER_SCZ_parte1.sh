
# SCRIPT - TFM - MASTER BIOINFO - 

# ESTUDIO DE SELECCION NATURAL RECIENTE EN ESQUIZOFRENIQ - PARTE 1

# Script creado para el cálculo de seleccion natural reciente en poblaciones humanas sobre las variantes asociadas con esquizofrenia
# Diciembre 2019 - Javier González-Peñas.
# ==============================================
# 
# Descripción -
# 
# Mediante este script se hacen las siguientes procesos:
# - Se descargan los datos de genoma completo de la base de datos de los 1000 genomas (1KG) y se juntan en una base común de genotipos, 
# se modifican para poder trabajar con ellos y se excluyen a los individuos americanos por tener demasiada mezcla poblacional
# - Se utilizan los datos de GWAS de SCZ (Ripke et al., 2014) para filtrar los datos genotípicos anteriores en base a aquellas variantes 
# genotipadas en el GWAS, se realiza un proceso de podado (clumping) para quedarse con variantes independientes y se exporta a un formato adecuado (PLINK)
# - Se analiza la información genotípica para extraer valores de probabilidad de estar sometidos a selección natural reciente por variante.
# 
# ==============================================


# ==============================================
# 1-  PROCESADO DE DATOS DE 1KG 
# ==============================================

# Los archivos de secuenciación de genoma completo se descargan para trabajar con ellos.
for i in {1..22}
do
   wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/LL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done

# Los archivos de secuenciación por cromosoma se juntan en un único archivo final.
bcftools concat ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Oz -o VCF.recode_ALL.vcf.gz

# Se descargan la lista de individuos de distintas poblaciones (disponible en excel en https://www.internationalgenome.org/1000-genomes-browsers): europeas (EUR), Asiáticas-este (EAS), Asiçaticas sur (SAS), Americanas (AMR) y Africanas (AFR)

# Seleccionamos solo los datos de la muestra europea para hacer después después LD-clumping. Eliminamos duplicados, y ordenamos. Exportamos a formato PLINK
bcftools view -S EUR -Ov  VCF.recode_ALL.vcf.gz > VCF.recode_EUR.vcf
grep -v '^#' VCF.recode_EUR.vcf | cut -f 3 | sort | uniq -d > reference.dups
./plink --vcf VCF.recode_EUR.vcf --exclude reference.dups --make-bed --out VCF.recode_EUR.dedup

# Preparamos listas de individuos de klas poblaciones, y creamos un VCF con todos los individuos menos los americanos (los eliminamos por demasiada mezcla poblacional)
cat EUR AFR SAS EAS > All_noAMR
vcftools --gzvcf VCFtodos_ALL.recode.vcf.gz --keep All_noAMR --recode --recode-INFO-all --out VCF.recode_ALL.vcf
gzip VCF.recode_ALL.vcf.recode.vcf



# ==============================================
# 2 -  PROCESADO DE DATOS DE GWAS DE SCZ
# ==============================================

## descargamos y descomprimimos el summary statistics del GWAS más ereciente de SCZ (Ripke et al., 2019)
gunzip rall.gz > rall.txt

# Seleccionamos variantes de de 1KG (paso anterior) que estén en el archivo de summary stats de SCZ, y filtramos por calidad de imputación > 0.9
awk 'NR==FNR{a[$3];next} FNR==1 || ($2 in a)' VCF.recode_ALL.rsID rall.txt |  awk -v OFS="\t" '$6 > 0.9  {print ;}' > overlap_SCZ
cat overlap_Qual | wc -l  #total de 5220878 variantes

# LD-Clumping para quedarnos con variantes SNPs independientes (tag-SNPs), priorizando por P valores de SCZ y sobre estructura haplotipica de población europea (Igual que el GWAS)
./plink --bfile VCF.recode_EUR.dedup --clump overlap_SCZ --clump-field p --clump-p1 1 --clump-p2 1 --clump-r2 0.25 --clump-kb 250 --out overlap_SCZ_clump

# Filtramos para seleccionar solo variantes independientes del paso anterior
awk 'NR==FNR {FILE1[$3]=$0; next} ($2 in FILE1) {print $0}' overlap_SCZ_clump.clumped overlap_SCZ > overlap_final_SCZ
cat overlap_final_SCZ | wc -l  #total de 235911 variantes independientes

#Extraemos la fila de rsID de SNPs
cat overlap_final_SCZ | awk -v OFS="\t" '{print $2}' | sed '1d' > SNP_SCZ_rsID_pcadapt

#En vcftools, hacemos subset del archivo VCF (datos genotípicos individuales) completo de 1KG para aquellas variantes seleccionadas. Nos quedamos solo con SNPs bialélicos, y MAF > 5%. Convertimos a formato PLINK
vcftools --gzvcf VCF.recode_ALL.vcf.recode.vcf.gz --remove-indels --maf 0.05 --min-alleles 2 --max-alleles 2 --snps SNP_SCZ_rsID_pcadapt --recode --recode-INFO-all --out VCF_pcadapt_SCZ
./plink --vcf VCF_pcadapt_SCZ.recode.vcf --make-bed --out VCF_pcadapt_SCZ




# ==============================================
# 3 -  ANÁLISIS DE SELECCIÓN NATURAL RECIENTE EN VARIANTES INDEPENDIENTES
# ==============================================

## ANÁLISIS DE SELECCIÓN NATURAL RECIENTE

## Cargamos datos en PCAdapt 
rstudio
library(pcadapt)
path_to_file <- "~/Desktop/TRABAJO/JAVIP/RESULTADOS/RESULTADOS FINALES/trabajo solo con PCAdapt/ASD/VCF_pcadapt_SCZ.bed"
FILE1 <- read.pcadapt(path_to_file,type="bed")
x <- pcadapt(FILE1,K=10,min.maf = 0.05)
plot(x,option="screeplot") #Vemos que 3 componentes son las que explican variabilidad (FIGURA_1)
dev.copy(pdf,'CP.pdf')
dev.off()

	# Rehacemos el análisis para 3 componentes, que son las que más explican 
	x <- pcadapt(FILE1,K=3,min.maf = 0.05)

## Añadimos información de las poblaciones y subpoblaciones de 1KG y exportamos los gráficos de componentes principales
POP<-read.table("POP.txt",header = TRUE)
POP2<-t(POP)
poplist.names <- POP2
plot(x,option="scores",pop=poplist.names)
dev.copy(pdf,'POP_PC1-2.pdf')
dev.off()
plot(x,option="scores",i=1,j=3,pop=poplist.names)
dev.copy(pdf,'POP_PC1-3.pdf')
dev.off()

SUBPOP<-read.table("SUBPOP.txt",header = TRUE)
SUBPOP2<-t(SUBPOP)
poplist.names <- SUBPOP2
plot(x,option="scores",pop=poplist.names)
dev.copy(pdf,'SUBPOP_PC1-2.pdf')
dev.off()
plot(x,option="scores",i=1,j=3,pop=poplist.names)
dev.copy(pdf,'SUBPOP_PC1-3.pdf')
dev.off()

# Exportamos otros gráficos (manhattan plot, QQplot y grafico de histograma de frecuencia)
plot(x,option="manhattan")
dev.copy(pdf,'manhattan.pdf')
dev.off()
plot(x,option="qqplot",threshold=0.1)
dev.copy(pdf,'qqplot.pdf')
dev.off()
hist(x$pvalues,xlab="p-values",main=NULL,breaks=1000)
plot(x,option="p-values")
dev.copy(pdf,'HIST_1000_pvalues.pdf')
dev.off()

# Creamos un vector de valores P de selección de las variantes analizadas y exportamos 
pvalues<-x$pvalues
write.table(pvalues, 'pvalues.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)



### FINAL parte 1




