logcpm_data=function(y){
  mat_cpm=y
  
  for(i in 1:ncol(y)){
    mat_cpm[,i]=(y[,i]/sum(y[,i]))*1000000
  }
  logcpm= log2(mat_cpm+1)
  return(logcpm)
}
  datacpm=logcpm_data (y)
