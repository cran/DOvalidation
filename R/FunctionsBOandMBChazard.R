# in this version we use functions from the first version of the package (no modification)
## LL hazard:
## hazard.LL(xi,Oi,Ei,x,b,K="epa",Ktype="symmetric",CI=FALSE) 
# for the rest we have modifcations since we allow different weighting (optional argument) 


## The multiplicative bias corrected estimator 

hazard.MBC<-function(xi,Oi,Ei,x,b,K="sextic",Ktype="symmetric") 
{
  #we have discretized data:
  # xi are the grid points (mid points of the considered intervals)
  # Oi is the vector with the number of occurrences into the interval
  # Ei is the vector with the exposure into the interval per unit of time,
  #     this is Yi*(x_i-x_{i-1})) with Yi being the total number of individual at risk
  #     at the beginning of the interval
  # x is the estimation point (or vector)
  # b is the bandwidth parameter (scalar)
  # K is the kernel function
  # Ktype is one of "symmetric", "left" or "right" so Ku is considered as
  #   the standard symmetric kernel or the one-sided versions
  
  #alpha.LL.xi is the local linear estimator evaluated at xi
  alpha.LL.xi<-as.vector(hazard.LL(xi,Oi,Ei,x=xi,b,K,Ktype,CI=FALSE)$hLL)
  #alpha.LL.x is the local linear estimator evaluated at x (it can be a vector)
  alpha.LL.x<-as.vector(hazard.LL(xi,Oi,Ei,x=x,b,K,Ktype,CI=FALSE)$hLL)
 
  M<-length(xi)
  ngrid<-length(x)
  if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
    u<-(x-xi)/b
  }
  # 'u' in general is a matrix with M rows and ngrid columns
  # if the evaluation is in only one single x then we have only a column vector
  if (K=="epa"){
    Ku<-function(u) K.epa(u)
    if (Ktype=="left") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u<0)  )
    if (Ktype=="right") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u>0)  )
  }
  if (K=="sextic"){
    Ku<-function(u) K.sextic(u)
    if (Ktype=="left") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u<0)  )
    if (Ktype=="right") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u>0)  )
  }
  
  # coefficients in the formula of K.bar
  # for the MBC estimator are calculated with weighting funcion W.tilde(s)=(alphaLL(t))^2
  Ku.u<-Ku(u)
  
  #coefficients in the formula of K.bar with the MBC weighting W.tilde(s)=(alphaLL(t))^2
  if (ngrid>1){
    aM.0<-b^(-1)*colSums(Ku.u*(alpha.LL.xi)^2*Ei,na.rm=T) 
    aM.1<-b^(-1)* colSums(Ku.u*(u*b)*(alpha.LL.xi)^2*Ei,na.rm=T)
    aM.2<-b^(-1)*colSums(Ku.u*((u*b)^2)*(alpha.LL.xi)^2*Ei,na.rm=T)
  }  else {
    aM.0<-b^(-1)*sum(Ku.u*(alpha.LL.xi)^2*Ei,na.rm=T) 
    aM.1<-b^(-1)*sum(Ku.u*(u*b)*(alpha.LL.xi)^2*Ei,na.rm=T)
    aM.2<-b^(-1)*sum(Ku.u*((u*b)^2)*(alpha.LL.xi)^2*Ei,na.rm=T)
  }
  
  Kbar.g<-function(u,Ku.u) 
  { 
    den<-(aM.0*aM.2-aM.1^2)
    correct<-(aM.2-aM.1*t((u*b)))/den  
    den.0<-which(abs(den)==0)
    correct[den.0,]<-NA
    Kbar<-t(correct)*Ku.u
    return(Kbar)
  }  
  
  if (ngrid>1) gM.x<-b^(-1)*colSums(Kbar.g(u,Ku.u)*(alpha.LL.xi)*Oi,na.rm=T) else {
    gM.x<-b^(-1)*sum(Kbar.g(u,Ku.u)*(alpha.LL.xi)*Oi,na.rm=T)} 
  gM.x[is.nan(gM.x)]<-1  ## if the MBC correction cannot be computed (infinity) then MBC=LL
  
  hMBC<- alpha.LL.x* gM.x
  
  hazard.res<-list(x=x , hMBC = hMBC)
  return(hazard.res)
}  



#########################################################
# BO-validation  
# if we pass the grid of bandwidths, it should be divided by C as we did for OSCV in
# the applications for the JRSSB
b.BO<-function(grid.b,nb,K="sextic",type.bo='Oi',xi,Oi,Ei,wei="same") 
{
  delta.M<-xi[2]-xi[1]
  M<-length(xi)
  
  # The rescaling constant for do-validation
  if (K=="epa") {Cval<-0.5371 }
  if (K=="sextic") {Cval<-0.5874 }
  M<-length(xi)
  if (missing(grid.b)){
    amp<-xi[M]-xi[1]
    b.min<-amp/(M+1)
    b.max<-amp/2
    if (missing(nb)) nb<-50
    grid.b<-seq(b.min,b.max,length=nb)
    grid.b<-grid.b/Cval
  }
  nb<-length(grid.b)
  # The hazard estimators with the best direction
  
  alphaLL.bo<-function(xi,Oi,Ei,x,b,K,type.bo='Oi')
  {
    M<-length(xi)
    
    ngrid<-length(x)
    if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
      u<-(x-xi)/b
    }
    uu<-u*b #is a matrix (ngrid x KK) same as u without scaling with b
    ### We calculate a vector with the same length of x indicating the type of kernel used
    #   to estimate at each component of x  
    #   Since the kernel is evaluated at u (with rows indicating the component of x):
    get.direct.i<-function(Oi,u,i)  
    {
      u<-matrix(u,M,ngrid)
      Oi.left<-sum(Oi[((abs(u[,i])<1) & (u[,i]<0))])
      Oi.right<-sum(Oi[((abs(u[,i])<1) & (u[,i]>0))])
      if (Oi.left<Oi.right) direct<-2 else direct<-1
      return(direct)
    }
    if (type.bo=='Ei') direct<-sapply(1:ngrid, get.direct.i,Oi=Ei,u=u) else direct<-sapply(1:ngrid, get.direct.i,Oi=Oi,u=u)
    ## direct is a vector with dimension ngrid
    
    Ku<-function(u,direct)
    {
      u<-matrix(u,M,ngrid)
      if (K=="epa"){  
        Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.epa(u[,i])*(u[,i]<0) else res<- 2*K.epa(u[,i])*(u[,i]>0)})
      } else {
        Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.sextic(u[,i])*(u[,i]<0) else res<- 2*K.sextic(u[,i])*(u[,i]>0)})
      }
      return(Ku.u)
    }    
    
    Ku.u<-Ku(u,direct)
    #coefficients in the formula of K.bar
    if (ngrid>1)
    {#a0.d<-b^(-1)*colSums(Ku.u*Ei,na.rm=T) #we don't need this
      a1.d<-b^(-1)* colSums(Ku.u*uu*Ei,na.rm=T)
      a2.d<-b^(-1)*colSums(Ku.u*(uu^2)*Ei,na.rm=T)
    }
    else 
    {#a0.d<-b^(-1)*sum(Ku.u*Ei,na.rm=T) #we don't need this
      a1.d<-b^(-1)*sum(Ku.u*uu*Ei,na.rm=T)
      a2.d<-b^(-1)*sum(Ku.u*(uu^2)*Ei,na.rm=T)
    }
    
    if (ngrid>1) {correct<-t(apply(uu,1,function(y){return(a2.d-y*a1.d)}))}
    else {correct<-(a2.d-a1.d*uu)}
    
    Kequiv.u<-correct*Ku.u
    
    if (ngrid>1)  OLL.d<-b^(-1)*colSums(Kequiv.u*Oi,na.rm=T) else {
      OLL.d<-b^(-1)*sum(Kequiv.u*Oi,na.rm=T)} 
    
    if (ngrid>1)  ELL.d<-b^(-1)*colSums(Kequiv.u*Ei,na.rm=T) else {
      ELL.d<-b^(-1)*sum(Kequiv.u*Ei,na.rm=T)} 
    
    hLL<-OLL.d/ELL.d
    ii<-which(is.nan(hLL))
    hLL[ii]<-0
    
    return(hLL)
  }
  ###
  
  cv.score<-function(b.ind)
  {
  # print(paste('b=',b.ind))
  cv1.d<-NA
  alpha.i<-as.vector(alphaLL.bo(xi=xi,Oi=Oi,Ei=Ei,x=xi,b=grid.b[b.ind],K=K,type.bo))
  
  if (wei=='exposure') cv1.d<-sum((alpha.i^2)*Ei,na.rm=TRUE)
  if (wei=='same') cv1.d<-sum((alpha.i^2)*delta.M,na.rm=TRUE)

  if (is.na(cv1.d)) cv.d<-NA else {
    # hi.xi.d is the estimator obtained after doing Oi=0
    O.i<-matrix(Oi,M,M);
    i0<-which(Oi==0)
    Oi0<-Oi-1; Oi0[i0]<-0
    diag(O.i)<-Oi0
    hi.xi.d<-sapply(1:M, function(i) {alphaLL.bo(xi=xi,Oi=O.i[,i],Ei=Ei,x=xi[i],
                                                 b=grid.b[b.ind],K,type.bo)})
    if (wei=='exposure') cv2.d<-sum(hi.xi.d*Oi,na.rm=TRUE)
    if (wei=='same') cv2.d<-sum(hi.xi.d*Oi*delta.M/Ei,na.rm=TRUE)
    cv.d<-(cv1.d-2*cv2.d)
  }
  cv.d[cv.d==0]<-NA
  
  return(cv.d)  # the CV-score (hatQ_0(b))
}

  cvbo.values<-as.vector(sapply(1:nb,cv.score)) 

  ind.bo<-which.min(cvbo.values)
  bbo<-grid.b[ind.bo]*Cval ## Rescale to the original grid of bandwidths
  if ((ind.bo==1)|(ind.bo==nb)) warning("The BO-validation score doesn't have a minumum in the grid of bandwidths")
  bo.res<-list(bbo=bbo,ind.bo=ind.bo,cvbo.values=cvbo.values,Cval=Cval,
                 grid.b=grid.b)
  return(bo.res)

}


b.BO.MBC<-function(grid.b,nb,K="sextic",type.bo='Oi',xi,Oi,Ei,wei="same")
{
  delta.M<-xi[2]-xi[1]
  M<-length(xi)
  
  # The rescaling constant for do-validation
  if (K=="epa") {Cval<-0.5947941}
  if (K=="sextic") {Cval<-0.6501} 
  
  if (missing(grid.b)){
    amp<-xi[M]-xi[1]
    b.min<-amp/(M+1)
    b.max<-amp/2
    if (missing(nb)) nb<-50
    grid.b<-seq(b.min,b.max,length=nb)
    grid.b<-grid.b/Cval
  }
  nb<-length(grid.b)
  
  ## a function to calculate the MBC hazard estimator with BO kernel
  alphaMBC.bo<-function(xi,Oi,Ei,x,b,K,type.bo='Oi')
  {
    M<-length(xi)
    
    ngrid<-length(x)
    if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
      u<-(x-xi)/b
    }
    
    ### We calculate a vector with the same length of x indicating the type of kernel used
    #   to estimate at each component of x  
    #   Since the kernel is evaluated at u (with rows indicating the component of x):
    get.direct.i<-function(Oi,u,i)  
    {
      u<-matrix(u,M,ngrid)
      Oi.left<-sum(Oi[((abs(u[,i])<1) & (u[,i]<0))])
      Oi.right<-sum(Oi[((abs(u[,i])<1) & (u[,i]>0))])
      if (Oi.left<Oi.right) direct<-2 else direct<-1
      return(direct)
    }
    if (type.bo=='Ei') direct<-sapply(1:ngrid, get.direct.i,Oi=Ei,u=u) else direct<-sapply(1:ngrid, get.direct.i,Oi=Oi,u=u)
    ## direct is a vector with dimension ngrid
    
    
    Ku<-function(u,direct)
    {
      u<-matrix(u,M,ngrid)
      if (K=="epa"){  
        Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.epa(u[,i])*(u[,i]<0) else res<- 2*K.epa(u[,i])*(u[,i]>0)})
      } else {
        Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.sextic(u[,i])*(u[,i]<0) else res<- 2*K.sextic(u[,i])*(u[,i]>0)})
      }
      return(Ku.u)
    }    
    
    Ku.u<-Ku(u,direct)
    
    ## first compute the LL for xi and then for x
    # The hazard estimators with the best direction
    
    alphaLL.bo<-function(xi,Oi,Ei,x,b,K,type.bo='Oi')
    {
      M<-length(xi)
      
      ngrid<-length(x)
      if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
        u<-(x-xi)/b
      }
      uu<-u*b #is a matrix (ngrid x KK) same as u without scaling with b
      ### We calculate a vector with the same length of x indicating the type of kernel used
      #   to estimate at each component of x  
      #   Since the kernel is evaluated at u (with rows indicating the component of x):
      get.direct.i<-function(Oi,u,i)  
      {
        u<-matrix(u,M,ngrid)
        Oi.left<-sum(Oi[((abs(u[,i])<1) & (u[,i]<0))])
        Oi.right<-sum(Oi[((abs(u[,i])<1) & (u[,i]>0))])
        if (Oi.left<Oi.right) direct<-2 else direct<-1
        return(direct)
      }
      if (type.bo=='Ei') direct<-sapply(1:ngrid, get.direct.i,Oi=Ei,u=u) else direct<-sapply(1:ngrid, get.direct.i,Oi=Oi,u=u)
      ## direct is a vector with dimension ngrid
      
      Ku<-function(u,direct)
      {
        u<-matrix(u,M,ngrid)
        if (K=="epa"){  
          Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.epa(u[,i])*(u[,i]<0) else res<- 2*K.epa(u[,i])*(u[,i]>0)})
        } else {
          Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.sextic(u[,i])*(u[,i]<0) else res<- 2*K.sextic(u[,i])*(u[,i]>0)})
        }
        return(Ku.u)
      }    
      
      Ku.u<-Ku(u,direct)
      #coefficients in the formula of K.bar
      if (ngrid>1)
      {#a0.d<-b^(-1)*colSums(Ku.u*Ei,na.rm=T) #we don't need this
        a1.d<-b^(-1)* colSums(Ku.u*uu*Ei,na.rm=T)
        a2.d<-b^(-1)*colSums(Ku.u*(uu^2)*Ei,na.rm=T)
      }
      else 
      {#a0.d<-b^(-1)*sum(Ku.u*Ei,na.rm=T) #we don't need this
        a1.d<-b^(-1)*sum(Ku.u*uu*Ei,na.rm=T)
        a2.d<-b^(-1)*sum(Ku.u*(uu^2)*Ei,na.rm=T)
      }
      
      if (ngrid>1) {correct<-t(apply(uu,1,function(y){return(a2.d-y*a1.d)}))}
      else {correct<-(a2.d-a1.d*uu)}
      
      Kequiv.u<-correct*Ku.u
      
      if (ngrid>1)  OLL.d<-b^(-1)*colSums(Kequiv.u*Oi,na.rm=T) else {
        OLL.d<-b^(-1)*sum(Kequiv.u*Oi,na.rm=T)} 
      
      if (ngrid>1)  ELL.d<-b^(-1)*colSums(Kequiv.u*Ei,na.rm=T) else {
        ELL.d<-b^(-1)*sum(Kequiv.u*Ei,na.rm=T)} 
      
      hLL<-OLL.d/ELL.d
      ii<-which(is.nan(hLL))
      hLL[ii]<-0
      
      return(hLL)
    }
    ###
    hLL.i<-alphaLL.bo(xi,Oi,Ei,x=xi,b,K,type.bo)
    
    if (length(x)==length(xi)){ 
      ind<-is.element(xi,x)
      hLL.x<-hLL.i[ind]
    } else hLL.x<-alphaLL.bo(xi,Oi,Ei,x=x,b,K,type.bo)
    ###
    
    #coefficients in the formula of K.bar with the MBC weighting W.tilde(s)=(alphaLL(t))^2
    if (ngrid>1){
      aM.0<-b^(-1)*colSums(Ku.u*(hLL.i)^2*Ei,na.rm=T) #we don't need this
      aM.1<-b^(-1)* colSums(Ku.u*(u*b)*(hLL.i)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*colSums(Ku.u*((u*b)^2)*(hLL.i)^2*Ei,na.rm=T)
    }  else {
      aM.0<-b^(-1)*sum(Ku.u*(hLL.i)^2*Ei,na.rm=T) #we don't need this
      aM.1<-b^(-1)*sum(Ku.u*(u*b)*(hLL.i)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*sum(Ku.u*((u*b)^2)*(hLL.i)^2*Ei,na.rm=T)
    }
    
    Kbar.g<-function(u,Ku.u) 
    { 
      den<-(aM.0*aM.2-aM.1^2)
      correct<-(aM.2-aM.1*t((u*b)))/den  
      den.0<-which(abs(den)==0)
      correct[den.0,]<-NA
      Kbar<-t(correct)*Ku.u
      return(Kbar)
    }  
    
    if (ngrid>1) gM.x<-b^(-1)*colSums(Kbar.g(u,Ku.u)*(hLL.i)*Oi,na.rm=T) else {
      gM.x<-b^(-1)*sum(Kbar.g(u,Ku.u)*(hLL.i)*Oi,na.rm=T)} 
    gM.x[is.nan(gM.x)]<-1  ## if the MBC correction cannot be computed (infinity) then MBC=LL
    
    hMBC<- hLL.x* gM.x
    #return(as.vector(hMBC.d))
    res<-list(hLL=hLL.i,hMBC=hMBC)
    return(res)
    #
  }
  ## loo version 
  alphaMBC.bo.loo<-function(xi,Oi,Ei,x,b,K,type.bo='Oi',hLL.i,hLL.x)
  {
    M<-length(xi)
    
    ngrid<-length(x)
    if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
      u<-(x-xi)/b
    }
    
    get.direct.i<-function(Oi,u,i)  
    {
      u<-matrix(u,M,ngrid)
      Oi.left<-sum(Oi[((abs(u[,i])<1) & (u[,i]<0))])
      Oi.right<-sum(Oi[((abs(u[,i])<1) & (u[,i]>0))])
      if (Oi.left<Oi.right) direct<-2 else direct<-1
      return(direct)
    }
    if (type.bo=='Ei') direct<-sapply(1:ngrid, get.direct.i,Oi=Ei,u=u) else direct<-sapply(1:ngrid, get.direct.i,Oi=Oi,u=u)
    ## direct is a vector with dimension ngrid
    
    
    Ku<-function(u,direct)
    {
      u<-matrix(u,M,ngrid)
      if (K=="epa"){  
        Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.epa(u[,i])*(u[,i]<0) else res<- 2*K.epa(u[,i])*(u[,i]>0)})
      } else {
        Ku.u<-sapply(1:ngrid, function(i) {if (direct[i]==1) res<- 2*K.sextic(u[,i])*(u[,i]<0) else res<- 2*K.sextic(u[,i])*(u[,i]>0)})
      }
      return(Ku.u)
    }    
    
    Ku.u<-Ku(u,direct)
    
    ## we do not do loo in the pilot LL 
    
    ## first compute the LL for xi and then for x
    
    #  hLL.i<-alphaLL.bo(xi,Oi,Ei,x=xi,b,type.bo)$hLL.d
    
    #  if (length(x)==length(xi)){ 
    #  ind<-is.element(xi,x)
    #  hLL.x<-hLL.i[ind]
    #} else hLL.x<-alphaLL.bo(xi,Oi,Ei,x=x,b,type.bo)$hLL.d
    ###
    ###
    
    #coefficients in the formula of K.bar with the MBC weighting W.tilde(s)=(alphaLL(t))^2
    if (ngrid>1){
      aM.0<-b^(-1)*colSums(Ku.u*(hLL.i)^2*Ei,na.rm=T) #we don't need this
      aM.1<-b^(-1)* colSums(Ku.u*(u*b)*(hLL.i)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*colSums(Ku.u*((u*b)^2)*(hLL.i)^2*Ei,na.rm=T)
    }  else {
      aM.0<-b^(-1)*sum(Ku.u*(hLL.i)^2*Ei,na.rm=T) #we don't need this
      aM.1<-b^(-1)*sum(Ku.u*(u*b)*(hLL.i)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*sum(Ku.u*((u*b)^2)*(hLL.i)^2*Ei,na.rm=T)
    }
    
    Kbar.g<-function(u,Ku.u) 
    { 
      den<-(aM.0*aM.2-aM.1^2)
      correct<-(aM.2-aM.1*t((u*b)))/den  
      den.0<-which(abs(den)==0)
      correct[den.0,]<-NA
      Kbar<-t(correct)*Ku.u
      return(Kbar)
    }  
    
    if (ngrid>1) gM.x<-b^(-1)*colSums(Kbar.g(u,Ku.u)*(hLL.i)*Oi,na.rm=T) else {
      gM.x<-b^(-1)*sum(Kbar.g(u,Ku.u)*(hLL.i)*Oi,na.rm=T)} 
    gM.x[is.nan(gM.x)]<-1  ## if the MBC correction cannot be computed (infinity) then MBC=LL
    
    hMBC<- hLL.x* gM.x
    return(as.vector(hMBC))
    
  }
  
  #  b.grid.do<-b.grid/Cval
  
  
  cv.score<-function(b.ind)
  {
    cv1.d<-NA
    res.alpha.i<-alphaMBC.bo(xi=xi,Oi=Oi,Ei=Ei,xi,b=grid.b[b.ind],K,type.bo)
    alpha.i<-res.alpha.i$hMBC
    hLL.i<-res.alpha.i$hLL
    
    if (wei=='exposure') cv1.d<-sum((alpha.i^2)*Ei,na.rm=TRUE)
    if (wei=='same') cv1.d<-sum((alpha.i^2)*delta.M,na.rm=TRUE)
    
    if (is.na(cv1.d)) cv.d<-NA else {
      # hi.xi.d is the estimator obtained after doing Oi=0
      O.i<-matrix(Oi,M,M);
      i0<-which(Oi==0)
      Oi0<-Oi-1; Oi0[i0]<-0
      diag(O.i)<-Oi0
      #hi.xi.d<-sapply(1:M, function(i) {alpha.bo(xi=xi,Oi=O.i[,i],Ei=Ei,x=xi[i],
      #        b=b.grid.do[b.ind],type.bo)})
     
      #using  alphaMBC.bo.loo<-function(xi,Oi,Ei,x,b,type.bo='Oi',hLL.i,hLL.x)
      #hLL.i<- as.vector(alphaLL.bo(xi=xi,Oi=Oi,Ei=Ei,xi,
      #           b=b.grid.do[b.ind],type.bo)$hLL.d)
      hi.xi.d<-sapply(1:M, function(i) {alphaMBC.bo.loo(xi=xi,Oi=O.i[,i],Ei=Ei,x=xi[i],
             b=grid.b[b.ind],K,type.bo,hLL.i=hLL.i,hLL.x=hLL.i[i])})
      
      if (wei=='exposure') cv2.d<-sum(hi.xi.d*Oi,na.rm=TRUE)
      if (wei=='same') cv2.d<-sum(hi.xi.d*Oi*delta.M/Ei,na.rm=TRUE)
      cv.d<-(cv1.d-2*cv2.d)
      cv.d[cv.d==0]<-NA
    }
    return(cv.d)  
    
  }
  cvbo.values<-as.vector(sapply(1:nb,cv.score)) 
  
  ind.bo<-which.min(cvbo.values)
  bbo<-grid.b[ind.bo]*Cval ## Rescale to the original grid of bandwidths
  if ((ind.bo==1)|(ind.bo==nb)) warning("The BO-validation score doesn't have a minumum in the grid of bandwidths")
  bo.res<-list(bbo=bbo,ind.bo=ind.bo,cvbo.values=cvbo.values,Cval=Cval,
               grid.b=grid.b)
  return(bo.res)# the CV-score (hatQ_0(b))
}

########################################################################
## The cross-validation method for MBC

b.CV.MBC<-function(grid.b,nb,K="sextic",xi,Oi,Ei,wei="same") 
{
  M<-length(xi)
  delta.M<-xi[2]-xi[1]
  if (missing(grid.b)){
    amp<-xi[M]-xi[1]
    b.min<-amp/(M+1)
    b.max<-amp/2
    if (missing(nb)) nb<-50
    grid.b<-seq(b.min,b.max,length=nb)
  }
  nb<-length(grid.b)

  hazard.MBC.info<-function(xi,Oi,Ei,x,b,K="sextic",Ktype="symmetric",alpha.LL.xi,alpha.LL.x) 
  {
    #alpha.LL.xi is the local linear estimator evaluated at xi
    #alpha.LL.x is the local linear estimator evaluated at x (it can be a vector)

    M<-length(xi)
    ngrid<-length(x)
    if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
      u<-(x-xi)/b
    }
    # 'u' in general is a matrix with M rows and ngrid columns
    # if the evaluation is in only one single x then we have only a column vector
    if (K=="epa"){
      Ku<-function(u) K.epa(u)
      if (Ktype=="left") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u<0)  )
      if (Ktype=="right") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u>0)  )
    }
    if (K=="sextic"){
      Ku<-function(u) K.sextic(u)
      if (Ktype=="left") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u<0)  )
      if (Ktype=="right") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u>0)  )
    }
    
    # coefficients in the formula of K.bar
    # for the MBC estimator are calculated with weighting funcion W.tilde(s)=(alphaLL(t))^2
    Ku.u<-Ku(u)
    
    #coefficients in the formula of K.bar with the MBC weighting W.tilde(s)=(alphaLL(t))^2
    if (ngrid>1){
      aM.0<-b^(-1)*colSums(Ku.u*(alpha.LL.xi)^2*Ei,na.rm=T) 
      aM.1<-b^(-1)* colSums(Ku.u*(u*b)*(alpha.LL.xi)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*colSums(Ku.u*((u*b)^2)*(alpha.LL.xi)^2*Ei,na.rm=T)
    }  else {
      aM.0<-b^(-1)*sum(Ku.u*(alpha.LL.xi)^2*Ei,na.rm=T) 
      aM.1<-b^(-1)*sum(Ku.u*(u*b)*(alpha.LL.xi)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*sum(Ku.u*((u*b)^2)*(alpha.LL.xi)^2*Ei,na.rm=T)
    }
    
    Kbar.g<-function(u,Ku.u) 
    { 
      den<-(aM.0*aM.2-aM.1^2)
      correct<-(aM.2-aM.1*t((u*b)))/den  
      den.0<-which(abs(den)==0)
      correct[den.0,]<-NA
      Kbar<-t(correct)*Ku.u
      return(Kbar)
    }  
    
    if (ngrid>1) gM.x<-b^(-1)*colSums(Kbar.g(u,Ku.u)*(alpha.LL.xi)*Oi,na.rm=T) else {
      gM.x<-b^(-1)*sum(Kbar.g(u,Ku.u)*(alpha.LL.xi)*Oi,na.rm=T)} 
    gM.x[is.nan(gM.x)]<-1  ## if the MBC correction cannot be computed (infinity) then MBC=LL
    
    hMBC<- alpha.LL.x* gM.x
    return(as.vector(hMBC))
  }  
  
  cv.score<-function(b.ind)
   {
     cv1.d<-NA
     b<-grid.b[b.ind]
     alpha.LL.i<-as.vector(hazard.LL(xi=xi,Oi=Oi,Ei=Ei,x=xi,b=b,K=K,Ktype="symmetric",CI=FALSE)$hLL)
     alpha.i<-as.vector(hazard.MBC.info(xi=xi,Oi=Oi,Ei=Ei,x=xi,b=b,K=K,Ktype="symmetric",
                                        alpha.LL.i,alpha.LL.i)) 
     
     
     if (wei=='exposure') cv1.d<-sum((alpha.i^2)*Ei,na.rm=TRUE)
     if (wei=='same') cv1.d<-sum((alpha.i^2)*delta.M,na.rm=TRUE)
         
     if (is.na(cv1.d)) cv.d<-NA else {
        # hi.xi.d is the estimator obtained after doing Oi=0
        O.i<-matrix(Oi,M,M); 
        i0<-which(Oi==0)
  	    Oi0<-Oi-1; Oi0[i0]<-0
        diag(O.i)<-Oi0 ; 
          
        hi.xi.d<-sapply(1:M, function(i) {hazard.MBC.info(xi=xi,Oi=O.i[,i],
                     Ei=Ei,x=xi[i], b=b,K=K,Ktype="symmetric",
                     alpha.LL.xi=alpha.LL.i, alpha.LL.x=alpha.LL.i[i])})
        if (wei=='exposure') cv2.d<-sum(hi.xi.d*Oi,na.rm=TRUE)
        if (wei=='same') cv2.d<-sum(hi.xi.d*Oi*delta.M/Ei,na.rm=TRUE)
        cv.d<-(cv1.d-2*cv2.d)
        cv.d[cv.d==0]<-NA
        
      }
     return(cv.d)  # the CV-score (hatQ_0(b))
  }
   cv.values<-sapply(1:nb,cv.score)
   ind.cv<-which.min(cv.values)
   if ((ind.cv==1)|(ind.cv==nb)) warning("The CV score doesn't have a minumum in the grid of bandwidths")
   bcv<-grid.b[ind.cv]
   cv.res<-list(bcv=bcv,ind.cv=ind.cv,cv.values=cv.values,grid.b=grid.b)
   return(cv.res)
}


b.OSCV.MBC<-function(grid.b,nb,K="sextic",Ktype="left",xi,Oi,Ei,wei="same") 
{
  # The rescaling constant for do-validation
  if (K=="epa") {Cval<-0.5947941}
  if (K=="sextic") {Cval<-0.6501}  
  
  M<-length(xi)
  delta.M<-xi[2]-xi[1]
  
  if (missing(grid.b)){
    amp<-xi[M]-xi[1]
    b.min<-amp/(M+1)
    b.max<-amp/2
    if (missing(nb)) nb<-50
    grid.b<-seq(b.min,b.max,length=nb)
    grid.b<-grid.b/Cval
  }
  nb<-length(grid.b)
  
  
  
  hazard.MBC.info<-function(xi,Oi,Ei,x,b,K="sextic",Ktype="symmetric",alpha.LL.xi,alpha.LL.x) 
  {
    #alpha.LL.xi is the local linear estimator evaluated at xi
    #alpha.LL.x is the local linear estimator evaluated at x (it can be a vector)
    
    M<-length(xi)
    ngrid<-length(x)
    if (ngrid>1) {u<-sapply(1:M,function(i) (x-xi[i])/b); u<-t(u)} else {
      u<-(x-xi)/b
    }
    # 'u' in general is a matrix with M rows and ngrid columns
    # if the evaluation is in only one single x then we have only a column vector
    if (K=="epa"){
      Ku<-function(u) K.epa(u)
      if (Ktype=="left") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u<0)  )
      if (Ktype=="right") Ku<-function(u) K.epa(u)*2*( (abs(u)<1) & (u>0)  )
    }
    if (K=="sextic"){
      Ku<-function(u) K.sextic(u)
      if (Ktype=="left") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u<0)  )
      if (Ktype=="right") Ku<-function(u) K.sextic(u)*2*( (abs(u)<1) & (u>0)  )
    }
    
    # coefficients in the formula of K.bar
    # for the MBC estimator are calculated with weighting funcion W.tilde(s)=(alphaLL(t))^2
    Ku.u<-Ku(u)
    
    #coefficients in the formula of K.bar with the MBC weighting W.tilde(s)=(alphaLL(t))^2
    if (ngrid>1){
      aM.0<-b^(-1)*colSums(Ku.u*(alpha.LL.xi)^2*Ei,na.rm=T) 
      aM.1<-b^(-1)* colSums(Ku.u*(u*b)*(alpha.LL.xi)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*colSums(Ku.u*((u*b)^2)*(alpha.LL.xi)^2*Ei,na.rm=T)
    }  else {
      aM.0<-b^(-1)*sum(Ku.u*(alpha.LL.xi)^2*Ei,na.rm=T) 
      aM.1<-b^(-1)*sum(Ku.u*(u*b)*(alpha.LL.xi)^2*Ei,na.rm=T)
      aM.2<-b^(-1)*sum(Ku.u*((u*b)^2)*(alpha.LL.xi)^2*Ei,na.rm=T)
    }
    
    Kbar.g<-function(u,Ku.u) 
    { 
      den<-(aM.0*aM.2-aM.1^2)
      correct<-(aM.2-aM.1*t((u*b)))/den  
      den.0<-which(abs(den)==0)
      correct[den.0,]<-NA
      Kbar<-t(correct)*Ku.u
      return(Kbar)
    }  
    
    if (ngrid>1) gM.x<-b^(-1)*colSums(Kbar.g(u,Ku.u)*(alpha.LL.xi)*Oi,na.rm=T) else {
      gM.x<-b^(-1)*sum(Kbar.g(u,Ku.u)*(alpha.LL.xi)*Oi,na.rm=T)} 
    gM.x[is.nan(gM.x)]<-1  ## if the MBC correction cannot be computed (infinity) then MBC=LL
    
    hMBC<- alpha.LL.x* gM.x
    return(as.vector(hMBC))
  }  
  
  cv.score<-function(b.ind)
  {
    cv1.d<-NA
    alphaLL.i<-hazard.LL(xi=xi,Oi=Oi,Ei=Ei,xi,b=grid.b[b.ind],K,Ktype=Ktype)$hLL
    alpha.i<-hazard.MBC.info(xi,Oi,Ei,x=xi,
                          b=grid.b[b.ind],K=K,Ktype=Ktype,
                          alpha.LL.xi=alphaLL.i,alpha.LL.x=alphaLL.i) 
      
    if (wei=='exposure') cv1.d<-sum((alpha.i^2)*Ei,na.rm=TRUE)
    if (wei=='same') cv1.d<-sum((alpha.i^2)*delta.M,na.rm=TRUE)
  
    if (is.na(cv1.d)) cv.d<-NA else {
      # hi.xi.d is the estimator obtained after doing Oi=0
      O.i<-matrix(Oi,M,M);
      i0<-which(Oi==0)
	    Oi0<-Oi-1; Oi0[i0]<-0
      diag(O.i)<-Oi0
      
      hi.xi.d<-sapply(1:M, function(i) {
          hazard.MBC.info(xi=xi,Oi=O.i[,i],Ei=Ei,x=xi[i],
           b=grid.b[b.ind],K=K,Ktype=Ktype,
           alpha.LL.xi=alphaLL.i,alpha.LL.x=alphaLL.i[i])})
      } 
      
      if (wei=='exposure') cv2.d<-sum(hi.xi.d*Oi,na.rm=TRUE)
      if (wei=='same') cv2.d<-sum(hi.xi.d*Oi*delta.M/Ei,na.rm=TRUE)
  
      cv.d<-(cv1.d-2*cv2.d)
      cv.d[cv.d==0]<-NA
      return(cv.d)  # the CV-score (hatQ_0(b))
  }
  
  oscv.values<-sapply(1:nb,cv.score) 
  
  ind.oscv<-which.min(oscv.values)
  boscv<-grid.b[ind.oscv]*Cval ## Rescale to the original grid of bandwidths
  if ((ind.oscv==1)|(ind.oscv==nb)) warning("The OSCV score doesn't have a minumum in the grid of bandwidths")
  oscv.res<-list(boscv=boscv,ind.oscv=ind.oscv,oscv.values=oscv.values,Cval=Cval,
                 grid.b=grid.b)
  return(oscv.res)
}




