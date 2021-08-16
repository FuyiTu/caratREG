#Parametric Estimators
samplevar<-function(A,S,R,strt,strt_num){
  r0 = rep(0,strt_num)
  r1 = rep(0,strt_num)
  nk0 = rep(0,strt_num)
  nk1 = rep(0,strt_num)
  sr0 = rep(0,strt_num)
  sr1 = rep(0,strt_num)
  for(i in 1:strt_num){
    inds = which(S == strt[i])
    nk1[i] = sum(A[inds])
    nk0[i] = sum(1 - A[inds])
    if(nk0[i] == 0){
        r0[i] = 0
        sr0[i] = 0
        r1[i] = sum(R[inds]*A[inds]/nk1[i])
        sr1[i] = sum((R[inds]-r1[i])^2*A[inds])/nk1[i]
    }
    else if(nk1[i] == 0){
        r1[i] = 0
        sr1[i] = 0
        r0[i] = sum(R[inds]*(1-A[inds])/nk0[i])
        sr0[i] = sum((R[inds]-r0[i])^2*(1-A[inds]))/nk0[i]
    }
    else{
      r0[i] = sum(R[inds]*(1-A[inds])/nk0[i])
      r1[i] = sum(R[inds]*A[inds]/nk1[i])
      sr0[i] = sum((R[inds]-r0[i])^2*(1-A[inds]))/nk0[i]
      sr1[i] = sum((R[inds]-r1[i])^2*A[inds])/nk1[i]
    }
  }
  pn = (nk0+nk1)/sum(nk0+nk1)
  tau = r1 - r0
  r1 = r1 - sum(A*R)/sum(nk1)
  r0 = r0 - sum((1-A)*R)/sum(nk0)
  return(cbind(pn,sr0,sr1,r0,r1,nk1/(nk0+nk1),tau))
}

samplevarr<-function(pn,sr0,sr1,pi){
  return(sum(pn*sr1)/pi+sum(pn*sr0)/(1-pi))
}

samplevarHr<-function(pn,r0,r1){
  return(sum(pn*(r1-r0)^2))
}

samplevarAr<-function(pn,sq,r0,r1,pi){
  return(sum(pn*sq*(r1/pi+r0/(1-pi))^2))
}

samplevarpr<-function(pn,sq,r0,r1,pi){
  return(((1-2*pi)/(pi*(1-pi)))^2*sum(pn*sq*(r1-r0)^2))
}

DME_var<-function(A,S,Y,sq,strt,strt_num,pi){
  n = length(A)
  sv = samplevar(A,S,Y,strt,strt_num)
  s = samplevarr(sv[,1],sv[,2],sv[,3],pi)+
    samplevarHr(sv[,1],sv[,4],sv[,5])+
    samplevarAr(sv[,1],sq,sv[,4],sv[,5],pi)
  return(sqrt(s/n))
}

FEE_var<-function(A,S,Y,sq,strt,strt_num,pi){
  n = length(A)
  sv = samplevar(A,S,Y,strt,strt_num)
  s = samplevarr(sv[,1],sv[,2],sv[,3],pi)+
    samplevarHr(sv[,1],sv[,4],sv[,5])+
    samplevarpr(sv[,1],sq,sv[,4],sv[,5],pi)
  return(sqrt(s/n))
}

SE_var<-function(A,S,Y,strt,strt_num,pi){
  n = length(A)
  sv = samplevar(A,S,Y,strt,strt_num)
  s = samplevarr(sv[,1],sv[,2],sv[,3],pi)+
    samplevarHr(sv[,1],sv[,4],sv[,5])
  return(sqrt(s/n))
}

regc<-function(A,Y,X){
  ind1 = which(A == 1)
  ind0 = which(A == 0)
  X1 = X[ind1,]
  X0 = X[ind0,]
  Y1 = Y[ind1]
  Y0 = Y[ind0]
  Y1 = Y1 - mean(Y1)
  Y0 = Y0 - mean(Y0)
  X1 = X1 - rep(colMeans(X1),rep.int(nrow(X1),ncol(X1)))
  X0 = X0 - rep(colMeans(X0),rep.int(nrow(X0),ncol(X0)))
  return(list(beta0 = MASS::ginv(t(X0)%*%X0)%*%t(X0)%*%Y0,
              beta1 = MASS::ginv(t(X1)%*%X1)%*%t(X1)%*%Y1))
}


regct<-function(A,S,Y,X,strt,strt_num){
  ind1 = which(A == 1)
  ind0 = which(A == 0)
  X1 = X[ind1,]
  X0 = X[ind0,]
  Y1 = Y[ind1]
  Y0 = Y[ind0]
  S1 = S[ind1]
  S0 = S[ind0]
  p = ncol(X)
  for(i in 1:strt_num){
    inds0 = which(S0 == strt[i])
    inds1 = which(S1 == strt[i])
    if(is.null(inds0)){
      if(length(inds1) == 1){
        X1[inds1,] = rep(0,ncol(X1))
        Y1[inds1] = 0
      }
      else{
        X1[inds1,] = X1[inds1,] - rep(colMeans(X1[inds1,]),rep.int(length(inds1),ncol(X1)))
        Y1[inds1] = Y1[inds1] - mean(Y1[inds1])
      }
    }
    else if(is.null(inds1)){
      if(length(inds0) == 1){
        X0[inds0,] = rep(0,ncol(X0))
        Y0[inds0] = 0
      }
      else{
        X0[inds0,] = X0[inds0,] - rep(colMeans(X0[inds0,]),rep.int(length(inds0),ncol(X0)))
        Y0[inds0] = Y0[inds0] - mean(Y0[inds0])
      }
    }
    else if(length(inds0) == 1 & length(inds1) == 1){
      X0[inds0,] = rep(0,ncol(X0))
      Y0[inds0] = 0
      X1[inds1,] = rep(0,ncol(X1))
      Y1[inds1] = 0
    }
    else if(length(inds0) == 1){
      X0[inds0,] = rep(0,ncol(X0))
      Y0[inds0] = 0
      X1[inds1,] = X1[inds1,] - rep(colMeans(X1[inds1,]),rep.int(length(inds1),ncol(X1)))
      Y1[inds1] = Y1[inds1] - mean(Y1[inds1])
    }
    else if(length(inds1) == 1){
      X0[inds0,] = X0[inds0,] - rep(colMeans(X0[inds0,]),rep.int(length(inds0),ncol(X0)))
      Y0[inds0] = Y0[inds0] - mean(Y0[inds0])
      X1[inds1,] = rep(0,ncol(X1))
      Y1[inds1] = 0
    }
    else{
      X0[inds0,] = X0[inds0,] - rep(colMeans(X0[inds0,]),rep.int(length(inds0),ncol(X0)))
      Y0[inds0] = Y0[inds0] - mean(Y0[inds0])
      X1[inds1,] = X1[inds1,] - rep(colMeans(X1[inds1,]),rep.int(length(inds1),ncol(X1)))
      Y1[inds1] = Y1[inds1] - mean(Y1[inds1])
    }
  }
  return(list(beta0 = MASS::ginv(t(X0)%*%X0)%*%t(X0)%*%Y0,
              beta1 = MASS::ginv(t(X1)%*%X1)%*%t(X1)%*%Y1))
}

SR_var<-function(A,S,Y,X,sq,strt,strt_num,pi){
  n = nrow(X)
  reg = regc(A,Y,X)
  Y = Y - X%*%(pi*reg$beta1 + (1-pi)*reg$beta0)
  sv = samplevar(A,S,Y,strt,strt_num)
  s = samplevarr(sv[,1],sv[,2],sv[,3],pi)+
    samplevarHr(sv[,1],sv[,4],sv[,5])+
    samplevarAr(sv[,1],sq,sv[,4],sv[,5],pi)
  return(sqrt(s/n))
}

AN_var<-function(A,S,Y,X,sq,strt,strt_num,pi){
  n = nrow(X)
  reg = regct(A,S,Y,X,strt,strt_num)
  Y = Y - X%*%(pi*reg$beta1 + (1-pi)*reg$beta0)
  sv = samplevar(A,S,Y,strt,strt_num)
  s = samplevarr(sv[,1],sv[,2],sv[,3],pi)+
    samplevarHr(sv[,1],sv[,4],sv[,5])+
    samplevarpr(sv[,1],sq,sv[,4],sv[,5],pi)
  return(sqrt(s/n))
}

IR_var<-function(A,S,Y,X,strt,strt_num,pi){
  n = nrow(X)
  reg = regct(A,S,Y,X,strt,strt_num)
  Y = Y - X%*%((1-pi)*reg$beta1 + pi*reg$beta0)
  sv = samplevar(A,S,Y,strt,strt_num)
  s = samplevarr(sv[,1],sv[,2],sv[,3],pi)+
    samplevarHr(sv[,1],sv[,4],sv[,5])
  return(sqrt(s/n))
}