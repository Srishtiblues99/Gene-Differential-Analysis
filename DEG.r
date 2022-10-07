#to identify which genes are differential in test vs control
data1=readRDS("log_ttest.rds")
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

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(matX,lab = rownames(matX),x='log2FC',y='pval')
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")