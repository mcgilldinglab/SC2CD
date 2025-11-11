m_CARD<-function(clu,cell_dec){
  ms<-data.frame(cluster=c(unique(clu)),Var=NA,Sum_Var=NA,Spots_num=NA)
  for (ci in ms$cluster) {
    cell<-scale(cell_dec[clu==ci,],scale=F)
    variance<-apply(cell,2,function(x)sum(x^2))
    ms[ms$cluster==ci,]$Var<-sum(variance)/dim(cell)[1]
    ms[ms$cluster==ci,]$Sum_Var<-sum(variance)
    ms[ms$cluster==ci,]$Spots_num<-dim(cell)[1]
  }
  var_all<-sum(ms$Sum_Var)/length(clu)
  return(list(Var_info=ms,G=var_all))
}

Eva<-function(dec,W_spa,M_adj){
  W_spa<-as.matrix(W_spa)
  N<-dim(dec)[1]
  S0<-sum(W_spa)
  S1<-sum(M_adj)
  mean_v<-colMeans(dec)
  M_mat<-as.matrix(t(dec)-mean_v)
  denom<-sum(M_mat^2)
  M_I<-N/S0*sum(W_spa*(t(M_mat)%*%M_mat))/denom
  G_C<-(N-1)/(2*S1)*sum(M_adj*rdist(dec))/denom
  return(list(M_I=M_I,G_C=G_C))
}

