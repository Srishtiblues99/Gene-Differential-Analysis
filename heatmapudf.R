library(ComplexHeatmap)
library(circlize)
data1=readRDS('log_ttest.RDS')
vec1= c("A1BG","A1BG-AS1","A1CF","A2ML1")
heatmp_fn=function(data1,vec1)
  {
  
  data=data1[vec1,]
  matx= as.matrix(y)
  for(i in 1:nrow(matx)){
    vec1= as.numeric(matx[i,])
    
    matx[i, 1:ncol(matx)] = (vec1-mean(vec1))/sd(vec1)
 
  }
  head(matx)
  mat1 = data1
  Heatmap(matx,col=colorRamp2(c(-2,0,2),c("purple","white","light pink")))
}
dataheatmp=heatmp_fn(data1,vec1)

pdf('data1.pdf',width = 10,height = 10)

plot(1:5,pch=20)

dev.off()


