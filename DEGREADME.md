# Gene-Differential-Analysis Overview:
1)Reading in table of counts
2)Adding annotation
3)Filtering lowly expressed genes
4)Quality control
5)Normalisation for composition bias
6)Differential expression analysis
7)Testing relative to a threshold
8)Visualisation
9)Gene set testing

# Install required packages in R
#install.packages("xlsx")
library(xlsx)
#install.packages("readxl")
library(readxl)

# Read ovarian raw count data
data=read.xlsx2("C:/Users/shristi/Downloads/cancer_data/GSE212991_Raw_counts.xlsx", row.names=1,sheetIndex = 1  )

mat = data

# Normalization of raw count data as count matrix

for(i in 1:ncol(data)){
  mat[,i] = as.numeric(data[,i])
}

library(edgeR)
library(limma)
log_data = cpm(mat, log = TRUE)

zs = log_data

# Z score calculation

for (i in 1:nrow(log_data)){
  
  vec = as.numeric(log_data[i,])
  zs[i,1:ncol(zs)] = (vec-mean(vec))/sd(vec)
  
}
head(zs)
mat1 = log_data

# Heatmap
library(ComplexHeatmap)
library(circlize)
Heatmap(zs[1:30,],col=colorRamp2(c(-2,0,2),c("purple","white","light pink")))

# Variance is calculated by applying variance formula on log data

var1 = apply(mat1, 1, var)

varGenes = sort(var1, decreasing = TRUE)

# Top 50 genes extracted from sorted matrix

top_gene=varGenes[1:50]

top_gene

mat_data=mat1[names(top_gene),]

mat_data

names(top_gene)

# Differential gene expression
# t_test

saveRDS(log_data,file = 'log_ttest.RDS')

data1=readRDS('log_ttest.RDS')

#Creating matrix with 4 columns
matX= matrix(NA,ncol = 4,nrow = nrow(data1))
rownames(matX)=rownames(data1)
colnames(matX)=c('test','control','pval','log2FC')

for(i in 1:nrow(data1)){
vec1= as.numeric(data1[i,1:2])
vec2=as.numeric(data1[i,3:6])
res=t.test(vec1,vec2,paired=F,alternative='two.sided')
matX[i,1]=res$estimate[[1]]
matX[i,2]=res$estimate[[2]]
matX[i,3]=res$p.value
matX[i,4]=matX[i,1]-matX[i,2]
}
matX=as.data.frame(matX)
num=which(is.nan(matX$pval1))
matX[num,'pval']=1

# Volcano Plot
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(matX,lab = rownames(matX),x='log2FC',y='pval')
#if (!require("BiocManager", quietly = TRUE))
 #install.packages("BiocManager")

library(matrixStats)
