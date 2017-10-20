discretise.data<-function(Li,Zi,deltai,xi,M)
{
  #Li:truncation levels: Li<Ti, otherwise the datum is not registered
  #when there is no truncation Li=0
  #Zi: observation time; Zi=min(Ti,Ci)
  #deltai: censoring indicator
  #xi: grid of times. If xi=0, the grid is calculated automatically
  #M: length of the grid
  
  Zi<-sort(Zi);deltai<-deltai[order(Zi)];Li<-Li[order(Zi)]
  risk<-function(x,Li,Zi){sapply(1:length(x),function(i) length(which(Li<x[i] & x[i]<=Zi)))}
  
  n<-length(Zi)
  
  if(missing(xi)){
    if (missing(M)) M<-round(n/2)
    rango<-max(Zi)-min(Zi)
    d.grid<-rango/(M-1)
    xi<-seq(min(Zi),max(Zi),by=d.grid)
  }
  M<-length(xi)
  
  d.grid<-diff(xi)[1]
  
  Zi.times<-Zi[deltai==1];
  Oi<-hist(Zi.times,breaks=c(xi,xi[M]+d.grid),right=F,plot=F)$counts ##occurrences
  Yi<-risk(xi,Li,Zi) #subjets at risk
  Ei<-Yi*d.grid  #exposure
  return(data.frame(xi,Oi,Ei))
}
#################################################################################################
## by luz
# Zi<-c(3,6,7,7,8,10,11,11,11,12,13,13,14,16,20,20,22,32,34,36)
# Li<-rep(0,20)
# deltai<-double(20); ii<-c(1,3,4,8,9,13,14,15,16)
# deltai[ii]<-1; deltai<-1-deltai
# ti<-0;M<-6
# discretase.data(Li,Zi,deltai,ti,M)
#################################################################################################


