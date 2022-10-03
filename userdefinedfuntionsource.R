#install.packages("xlsx")
library(xlsx)
#install.packages("readxl")
library(readxl)

data=read.xlsx2("C:/Users/shristi/Downloads/cancer_data/GSE212991_Raw_counts.xlsx", row.names=1,sheetIndex = 1  )
data
mat = data

for(i in 1:ncol(data)){
  mat[,i] = as.numeric(data[,i])
}
logdata1=source('log2udf.r')
logcpm_data(mat)
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
heatmpdata=source('heatmapudf.r')
vec1= c("A1BG","A1BG-AS1","A1CF","A2ML1")
heatmp_fn(data1,vec1)

