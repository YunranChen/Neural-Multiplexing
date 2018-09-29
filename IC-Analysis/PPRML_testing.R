library("statmod")
library(purrr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(rbenchmark)
library(pROC)
library(stringr)
#library(truncnorm)
#rm(list=ls())


## logsum and logmean
logsum=function(x){
  #to calculate log(sum(exp(x)))
  res=max(x)+log(sum(exp(x-max(x))))
  return(res)
}
logmean=function(x){
  #to calculate log(sum(exp(x)))
  res=max(x)+log(sum(exp(x-max(x))))-log(length(x))
  return(res)
}

## decide the range for x according to rule-of-thumb


x_range=function(xs_vec,n_mu,l,u){
  n=length(xs_vec)
  bw_cut=sd(xs_vec)*(4/3/n)^(1/5)*3
  ub=max(xs_vec)+bw_cut
  lb=min(xs_vec)-bw_cut
  xxs=seq(max(lb,l),min(ub,u),length.out = n_mu)
  return(xxs)
}



rtgamma=function(size,shape_,rate_,a,b){
  u=runif(n = size)
  c_inv=pgamma(q = b,shape = shape_,rate = rate_)-pgamma(q = a,shape = shape_,rate = rate_)
  x=qgamma(p = u*c_inv+pgamma(q=a,shape=shape_,rate=rate_),shape=shape_,rate = rate_)
  return(x)
}
dtgamma=function(x_,shape_,rate_,a,b){
  c_inv=pgamma(q = b,shape = shape_,rate = rate_)-pgamma(q = a,shape = shape_,rate = rate_)
  x=dgamma(x = x_,shape=shape_,rate = rate_)/c_inv
  return(x)
}

#PPRML for single trial
PPRML_int=function(xs_bn,nGQ=20,n_mu=100,nP=100,alpha=0.5){
  #return log(m_x)
  #permutation
  n=length(xs_bn)
  ind=sapply(1:nP,function(x){sample.int(n = n,size = n,replace = FALSE)})
  #guass.quad
  out=gauss.quad(nGQ)
  ## weight seq
  w=1/(1:n+1)
  ##Recursion
    Xs_bn_p=sapply(1:nP,function(i){xs_bn[ind[,i]]}) #each col is a permutation
    lq=quantile(xs_bn,0.25)
    uq=quantile(xs_bn,0.75)
    a=lq-alpha*(uq-lq)
    b=uq+alpha*(uq-lq)
    names(a)=names(b)=NULL
    xs_gq=out$nodes*(b-a)/2+(a+b)/2 #change interval
    theta=seq(a,b,length.out = n_mu)
    ## Initial Guess: Uniform on [0,10]
    f=rep(1/(b-a),n_mu)
    f[1]=0
    f[n_mu]=0
    f_gq=rep(1/(b-a),nGQ)
    res_t=map(1:nP,function(ix){
      m_x=0
      for (j in 1:n){
        poi_x=map_dbl(theta,~dpois(x = Xs_bn_p[j,ix],lambda = .x))
        px_grid=map_dbl(xs_gq,~dpois(x = Xs_bn_p[j,ix],lambda = .x))
        mi_x=((b-a)/2*sum(out$weights*px_grid*f_gq))
        f_gq=(1-w[j])*f_gq+w[j]*px_grid*f_gq/mi_x
        f=(1-w[j])*f+w[j]*poi_x*f/mi_x
        m_x=log(mi_x)+m_x
      }
      return(list(f=f,m_x=m_x))})
    pf=sapply(res_t,function(res_t){res_t$f})%>%apply(., 1,mean)
    res_py=logmean(map_dbl(res_t,~.x$m_x))
  #return(list(res_f=pf,res_m=res_py))
  return(res_m=res_py)
}
CML=function(xs_bn,alpha=0.5){
    yi=xs_bn
    alpha=0.5
    lq=quantile(yi,0.25)
    uq=quantile(yi,0.75)
    a=lq-alpha*(uq-lq)
    b=uq+alpha*(uq-lq)
    names(a)=names(b)=NULL
    res=-log(b-a)-sum(lgamma(yi+1))+lgamma(sum(yi)+1)-(sum(yi)+1)*log(length(yi))+log(pgamma(q = b,shape = sum(yi)+1,rate = length(yi))-pgamma(q = a,shape = sum(yi)+1,rate = length(yi)))
  return(res)
}
PPRML.poi=function(xs_bn,nGQ=20,n_mu=100,nP=100,alpha=0.5){
  res=exp(CML(xs_bn,alpha)-PPRML_int(xs_bn,nGQ,n_mu,nP,alpha))
  return(res)
}

#PPRML for multiplexing
PPRML_outA_lp=function(xs_bn,xs_a,r_a=0.5,s_a=2e-10,mu_l=0,nGQ=20,n_mu=100,nP=100){
  
  r_A=sum(xs_a)+r_a
  s_A=length(xs_a)+s_a
  imu_a=mean(xs_a)
  
  n=length(xs_bn)

  indx=sapply(1:nP,function(x){sample.int(n = n,size = n,replace = FALSE)})
  
  #guass.quad
  out=gauss.quad(nGQ)
  
    allf=function(phi){
      a=0
      b=1
      xs_gq=out$nodes*(b-a)/2+(a+b)/2 #change interval
      #ind=ceiling((xs_gq-a)/((b-a)/n_mu))
      z=seq(a,b,length.out = n_mu)
      
      ## Initial Guess:
      
      f=rep(1,n_mu) #uniform
      f[n_mu]=0
      d_f=d_fgq=0
      f_gq=rep(1,nGQ)
      ## weight seq
      w=1/(1:n+1)
      
      ##Recursion
      Xs_bn_p=sapply(1:nP,function(i){xs_bn[indx[,i]]}) #each col is a permutation
      res_t=map(1:nP,function(ix){
        m_x=0
        grad=0
        for (j in 1:n){
          #step1:prml
          px_grid=map_dbl(xs_gq,~dpois(x = Xs_bn_p[j,ix],lambda = .x*(phi-mu_l)+mu_l))
          mi_x=((b-a)/2*sum(out$weights*px_grid*f_gq))
          poi_x=map_dbl(z,~dpois(x = Xs_bn_p[j,ix],lambda = .x*(phi-mu_l)+mu_l))
          #step2:gradient
          d_poi_x=poi_x*(Xs_bn_p[j,ix]*z/(mu_l+z*(phi-mu_l))-z)
          d_px_grid=px_grid*(Xs_bn_p[j,ix]*xs_gq/(mu_l+xs_gq*(phi-mu_l))-xs_gq)
          G1=d_poi_x*f+poi_x*d_f
          G1_gq=d_px_grid*f_gq+px_grid*d_fgq
          G2=poi_x*f
          G2_gq=px_grid*f_gq
          d_logmi_x=(b-a)/2*sum(out$weights*G1_gq)/mi_x
          #step3:update f,gr_f
          f=(1-w[j])*f+w[j]*poi_x*f/mi_x
          f_gq=(1-w[j])*f_gq+w[j]*px_grid*f_gq/mi_x
          d_f=(1-w[j])*d_f+w[j]*(G1-G2*d_logmi_x)/mi_x
          d_fgq=(1-w[j])*d_fgq+w[j]*(G1_gq-G2_gq*d_logmi_x)/mi_x
          #store the result
          m_x=log(mi_x)+m_x
          grad=grad+d_logmi_x
        }
        return(list(f=f,m_x=m_x,grad=grad))
      })
      pf=sapply(res_t,function(res_t){res_t$f})%>%apply(., 1,mean)
      pm_x=map_dbl(res_t,~.x$m_x)%>%logmean(.)
      pgrad=map_dbl(res_t,function(res_t){res_t$grad})%>%logmean(.)
      return(list(f=pf,m_x=pm_x,grad=pgrad))
    }
    obj=function(theta){
      phi=exp(theta)+mu_l
      p_mua=dgamma(x = phi,shape = r_A,rate = s_A,log = TRUE)
      return(-(allf(phi)$m_x+p_mua))
    }
    grad=function(theta){
      phi=exp(theta)+mu_l
      g_mua=((r_A-1)/phi-s_A)
      prml=allf(phi)$grad
      res=-(prml+g_mua)*exp(theta)
      return(res)
    }
    res=optim(par = log(imu_a-mu_l),fn = obj,gr = grad,method = "BFGS",hessian = TRUE)
    #optim is for minimization
    H=(-res$hessian-grad(res$par))*exp(-2*res$par)
    logmx=0.5*log(2*pi)+0.5*log((-H)^(-1))-res$value
  res_py=logmx
  return(res_py)
}
PPRML_outB_lp=function(xs_bn,xs_b,r_b=0.5,s_b=2e-10,mu_u=180,nGQ=20,n_mu=100,nP=100){
  
  r_B=sum(xs_b)+r_b
  s_B=length(xs_b)+s_b
  imu_b=mean(xs_b)
  
  n=length(xs_bn)
  indx=sapply(1:nP,function(x){sample.int(n = n,size = n,replace = FALSE)})
  
  #guass.quad
  out=gauss.quad(nGQ)
  
    allf=function(phi){
      a=0
      b=1
      xs_gq=out$nodes*(b-a)/2+(a+b)/2 #change interval
      #ind=ceiling((xs_gq-a)/((b-a)/n_mu))
      z=seq(a,b,length.out = n_mu)
      
      ## Initial Guess:
      
      f=rep(1,n_mu) #uniform
      f[1]=0
      d_f=d_fgq=0
      f_gq=rep(1,nGQ)
      
      ## weight seq
      w=1/(1:n+1)
      
      ##Recursion
      Xs_bn_p=sapply(1:nP,function(i){xs_bn[indx[,i]]}) #each col is a permutation
      res_t=map(1:nP,function(ix){
        
        m_x=0
        grad=0
        for (j in 1:n){
          #step1:prml
          px_grid=map_dbl(xs_gq,~dpois(x = Xs_bn_p[j,ix],lambda = .x*(mu_u-phi)+phi))
          mi_x=((b-a)/2*sum(out$weights*px_grid*f_gq))
          poi_x=map_dbl(z,~dpois(x = Xs_bn_p[j,ix],lambda = .x*(mu_u-phi)+phi))
          #step2:gradient
          d_poi_x=poi_x*(Xs_bn_p[j,ix]*(1-z)/(phi+z*(mu_u-phi))-(1-z))
          d_px_grid=px_grid*(Xs_bn_p[j,ix]*(1-xs_gq)/(phi+xs_gq*(mu_u-phi))-(1-xs_gq))
          G1=d_poi_x*f+poi_x*d_f
          G1_gq=d_px_grid*f_gq+px_grid*d_fgq
          G2=poi_x*f
          G2_gq=px_grid*f_gq
          d_logmi_x=(b-a)/2*sum(out$weights*G1_gq)/mi_x
          #step3:update f,gr_f
          f=(1-w[j])*f+w[j]*poi_x*f/mi_x
          f_gq=(1-w[j])*f_gq+w[j]*px_grid*f_gq/mi_x
          d_f=(1-w[j])*d_f+w[j]*(G1-G2*d_logmi_x)/mi_x
          d_fgq=(1-w[j])*d_fgq+w[j]*(G1_gq-G2_gq*d_logmi_x)/mi_x
          #store the result
          m_x=log(mi_x)+m_x
          grad=grad+d_logmi_x
        }
        return(list(f=f,m_x=m_x,grad=grad))
      })
      pf=sapply(res_t,function(res_t){res_t$f})%>%apply(., 1,mean)
      pm_x=map_dbl(res_t,~.x$m_x)%>%logmean(.)
      pgrad=map_dbl(res_t,function(res_t){res_t$grad})%>%logmean(.)
      return(list(f=pf,m_x=pm_x,grad=pgrad))
    }
    
    obj=function(theta){
      phi=exp(theta)*mu_u/(1+exp(theta))
      p_mua=dgamma(x = phi,shape = r_B,rate = s_B,log = TRUE)
      return(-(allf(phi)$m_x+p_mua))
    }
    grad=function(theta){
      phi=exp(theta)*mu_u/(1+exp(theta))
      g_mua=((r_B-1)/phi-s_B)
      prml=allf(phi)$grad
      res=-(prml+g_mua)*phi*(1-phi)/mu_u
      return(res)
    }
    res=optim(par = log(imu_b/(mu_u-imu_b)),fn = obj,gr = grad,method = "BFGS",hessian = TRUE)
    #optim is for minimization
    phi_r=exp(res$par)*mu_u/(1+exp(res$par))
    H=(-res$hessian+grad(res$par)*(-1/phi_r+1/(1-phi_r)))*mu_u/(phi_r-phi_r^2)
    logmx=0.5*log(2*pi)+0.5*log((-H)^(-1))-res$value

  res_py=logmx
  return(res_py)
}
PPRML_sin_lp=function(xs_bn,xs_a,r_a=0.5,s_a=2e-10,nGQ=20){
  
  r_A=sum(xs_a)+r_a
  s_A=length(xs_a)+s_a
  imu_a=mean(xs_a)
  
  n=length(xs_bn)
  

    obj=function(theta){
      phi=exp(theta)
      m_x=dpois(x = xs_bn,lambda = phi,log = TRUE)%>%sum(.)
      p_mua=dgamma(x = phi,shape = r_A,rate = s_A,log = TRUE)
      return(-(m_x+p_mua))
    }
    grad=function(theta){
      phi=exp(theta)
      grad_x=sum(xs_bn/phi-1)
      g_mua=(r_A-1)/phi-s_A
      res=-(grad_x+g_mua)*phi
      return(res)
    }
    res=optim(par = log(imu_a),fn = obj,gr = grad,method = "BFGS",hessian = TRUE)
    H=(-res$hessian-grad(res$par))*exp(-2*res$par)
    logmx=0.5*log(2*pi)+0.5*log((-H)^(-1))-res$value

  res_py=logmx
  return(res_py)
}
PPRML_int_lp=function(xs_bn,xs_a,xs_b,e=0,r_a=0.5,s_a=2e-10,r_b=0.5,s_b=2e-10,nGQ=20,n_mu=100,nP=100){
  
  r_A=sum(xs_a)+r_a
  s_A=length(xs_a)+s_a
  r_B=sum(xs_b)+r_b
  s_B=length(xs_b)+s_b
  imu_a=mean(xs_a)
  imu_b=mean(xs_b)
  
  n=length(xs_bn)

  indx=sapply(1:nP,function(x){sample.int(n = n,size = n,replace = FALSE)})
  
  #guass.quad
  out=gauss.quad(nGQ)
  
    allf=function(phi){
      a=0
      b=1
      xs_gq=out$nodes*(b-a)/2+(a+b)/2 #change interval
      ind=ceiling((xs_gq-a)/((b-a)/n_mu))
      z=seq(a,b,length.out = n_mu)
      eind=ceiling((e-a)/((b-a)/n_mu))
      
      ## Initial Guess:
      if (e==0){
        f=rep(1,n_mu) #uniform
        f[n_mu]=0
        f[1]=0
        f_gq=rep(1,nGQ)
      }else{
        f=rep(1/(1-2*e),n_mu)
        f[1:eind]=0
        f[(n_mu-eind+1):n_mu]=0
        f_gq=rep(1/(1-2*e),n_mu)
        f_gq[xs_gq<=e]=0
        f_gq[xs_gq>(1-e)]=0
      }
      
      d_f=matrix(0,nrow=2,ncol=n_mu)
      d_fgq=matrix(0,nrow=2,ncol=n_mu)
      ## weight seq
      w=1/(1:n+1)
      
      ##Recursion
      Xs_bn_p=sapply(1:nP,function(i){xs_bn[indx[,i]]}) #each col is a permutation
      res_t=map(1:nP,function(ix){
        m_x=0
        grad=c(0,0)
        for (j in 1:n){
          #step1:prml
          px_grid=map_dbl(xs_gq,~dpois(x = Xs_bn_p[j,ix],lambda = phi[1]+.x*(phi[2]-phi[1])))
          mi_x=((b-a)/2*sum(out$weights*px_grid*f_gq))
          poi_x=map_dbl(z,~dpois(x = Xs_bn_p[j,ix],lambda = phi[1]+.x*(phi[2]-phi[1])))
          #step2:gradient
          d_poi_x_a=poi_x*(Xs_bn_p[j,ix]*(1-z)/(phi[1]+z*(phi[2]-phi[1]))-(1-z))
          d_poi_x_b=poi_x*(Xs_bn_p[j,ix]*z/(phi[1]+z*(phi[2]-phi[1]))-z)
          d_px_grid_a=px_grid*(Xs_bn_p[j,ix]*(1-xs_gq)/(phi[1]+xs_gq*(phi[2]-phi[1]))-(1-xs_gq))
          d_px_grid_b=px_grid*(Xs_bn_p[j,ix]*xs_gq/(phi[1]+xs_gq*(phi[2]-phi[1]))-xs_gq)
          G_a=d_poi_x_a*f+poi_x*d_f[1,]
          G_b=d_poi_x_b*f+poi_x*d_f[2,]
          G_a_gq=d_px_grid_a*f_gq+px_grid*d_fgq[1,]
          G_b_gq=d_px_grid_b*f_gq+px_grid*d_fgq[2,]
          G2=poi_x*f
          G2_gq=px_grid*f_gq
          d_logmi_x_a=(b-a)/2*sum(out$weights*G_a_gq)/mi_x
          d_logmi_x_b=(b-a)/2*sum(out$weights*G_b_gq)/mi_x
          d_logmi_x=c(d_logmi_x_a,d_logmi_x_b)
          #step3:update f,gr_f
          d_f[1,]=(1-w[j])*d_f[1,]+w[j]*(G_a-G2*d_logmi_x_a)/mi_x
          d_f[2,]=(1-w[j])*d_f[2,]+w[j]*(G_b-G2*d_logmi_x_b)/mi_x
          f=(1-w[j])*f+w[j]*poi_x*f/mi_x
          d_fgq[1,]=(1-w[j])*d_fgq[1,]+w[j]*(G_a_gq-G2_gq*d_logmi_x_a)/mi_x
          d_fgq[2,]=(1-w[j])*d_fgq[2,]+w[j]*(G_b_gq-G2_gq*d_logmi_x_b)/mi_x
          f_gq=(1-w[j])*f_gq+w[j]*px_grid*f_gq/mi_x
          #store the result
          m_x=log(mi_x)+m_x
          grad=grad+d_logmi_x
        }
        return(list(f=f,m_x=m_x,grad=grad))
      })
      pf=sapply(res_t,function(res_t){res_t$f})%>%apply(., 1,mean)
      pm_x=map_dbl(res_t,~.x$m_x)%>%logmean(.)
      pgrad=sapply(res_t,function(res_t){res_t$grad})%>%apply(., 1,logmean)
      return(list(f=pf,m_x=pm_x,grad=pgrad))}
    
    obj=function(theta){
      phi=exp(sort(theta))
      p_mua=dgamma(x = phi[1],shape = r_A,rate = s_A,log = TRUE)
      p_mub=dgamma(x = phi[2],shape = r_B,rate = s_B,log = TRUE)
      return(-(allf(phi)$m_x+p_mua+p_mub))
    }
    grad=function(theta){
      phi=exp(sort(theta))
      g_mua=(r_A-1)/phi[1]-s_A
      g_mub=(r_B-1)/phi[2]-s_B
      prml=allf(phi)$grad
      res=-(prml+c(g_mua,g_mub))*phi
      return(res[order(theta)])
    }
    res=optim(par = c(log(imu_a),log(imu_b)),fn = obj,gr = grad,method = "BFGS",hessian = TRUE)
    #optim is for minimization
    H=matrix(0,nrow=2,ncol=2)
    res_par=res$par%>%sort(.)
    H[1,1]=(-res$hessian[1,1]-grad(res_par)[1])*exp(-2*res_par[1])
    H[2,2]=(-res$hessian[2,2]-grad(res_par)[2])*exp(-2*res_par[2])
    H[1,2]=-res$hessian[1,2]*exp(-sum(res_par))
    H[2,1]=-res$hessian[2,1]*exp(-sum(res_par))
    
    logmx=log(2*pi)-0.5*log(det(-H))-res$value
  res_py=logmx
  return(res_py)
}
PPRML_mix_lp=function(xs_bn,xs_a,xs_b,e=0,r_a=0.5,s_a=2e-10,r_b=0.5,s_b=2e-10,nGQ=20,n_mu=100,nP=100){
  
  r_A=sum(xs_a)+r_a
  s_A=length(xs_a)+s_a
  r_B=sum(xs_b)+r_b
  s_B=length(xs_b)+s_b
  imu_a=mean(xs_a)
  imu_b=mean(xs_b)
  
  n=length(xs_bn)

  indx=sapply(1:nP,function(x){sample.int(n = n,size = n,replace = FALSE)})
  
  
    allf=function(phi){
      a=0
      b=1
      
      z=c(a+e,b-e)
      
      ## Initial Guess:
      f=c(0.5,0.5)
      d_f=matrix(0,nrow=2,ncol=2)
      
      ## weight seq
      w=1/(1:n+1)
      
      ##Recursion
      Xs_bn_p=sapply(1:nP,function(i){xs_bn[indx[,i]]}) #each col is a permutation
      res_t=map(1:nP,function(ix){
        m_x=0
        grad=c(0,0)
        for (j in 1:n){
          #step1:prml
          poi_x=map_dbl(z,~dpois(x = Xs_bn_p[j,ix],lambda = phi[1]+.x*(phi[2]-phi[1])))
          mi_x=sum(poi_x*f)
          
          #step2:gradient
          d_poi_x_a=poi_x*(Xs_bn_p[j,ix]*(1-z)/(phi[1]+z*(phi[2]-phi[1]))-(1-z))
          d_poi_x_b=poi_x*(Xs_bn_p[j,ix]*z/(phi[1]+z*(phi[2]-phi[1]))-z)
          G_a=d_poi_x_a*f+poi_x*d_f[1,]
          G_b=d_poi_x_b*f+poi_x*d_f[2,]
          G2=poi_x*f
          d_logmi_x_a=sum(G_a)/mi_x
          d_logmi_x_b=sum(G_b)/mi_x
          d_logmi_x=c(d_logmi_x_a,d_logmi_x_b)
          #step3:update f,gr_f
          d_f[1,]=(1-w[j])*d_f[1,]+w[j]*(G_a-G2*d_logmi_x_a)/mi_x
          d_f[2,]=(1-w[j])*d_f[2,]+w[j]*(G_b-G2*d_logmi_x_b)/mi_x
          f=(1-w[j])*f+w[j]*poi_x*f/mi_x
          
          #store the result
          m_x=log(mi_x)+m_x
          grad=grad+d_logmi_x
        }
        return(list(f=f,m_x=m_x,grad=grad))
      })
      pf=sapply(res_t,function(res_t){res_t$f})%>%apply(., 1,mean)
      pm_x=map_dbl(res_t,~.x$m_x)%>%logmean(.)
      pgrad=sapply(res_t,function(res_t){res_t$grad})%>%apply(., 1,logmean)
      return(list(f=pf,m_x=pm_x,grad=pgrad))}
    
    obj=function(theta){
      phi=exp(sort(theta))
      p_mua=dgamma(x = phi[1],shape = r_A,rate = s_A,log = TRUE)
      p_mub=dgamma(x = phi[2],shape = r_B,rate = s_B,log = TRUE)
      return(-(allf(phi)$m_x+p_mua+p_mub))
    }
    grad=function(theta){
      phi=exp(sort(theta))
      g_mua=(r_A-1)/phi[1]-s_A
      g_mub=(r_B-1)/phi[2]-s_B
      prml=allf(phi)$grad
      res=-(prml+c(g_mua,g_mub))*phi
      return(res[order(theta)])
    }
    res=optim(par = c(log(imu_a),log(imu_b)),fn = obj,gr = grad,method = "BFGS",hessian = TRUE)
    #optim is for minimization
    H=matrix(0,nrow=2,ncol=2)
    res_par=res$par%>%sort(.)
    H[1,1]=(-res$hessian[1,1]-grad(res_par)[1])*exp(-2*res_par[1])
    H[2,2]=(-res$hessian[2,2]-grad(res_par)[2])*exp(-2*res_par[2])
    H[1,2]=-res$hessian[1,2]*exp(-sum(res_par))
    H[2,1]=-res$hessian[2,1]*exp(-sum(res_par))
    
    logmx=log(2*pi)-0.5*log(det(-H))-res$value
  res_py=logmx
  return(res_py)
}

PPRML.poi.mix=function(xs_bn,xs_a,xs_b,mu_l=0,mu_u=180,e=0,
                       r_a=0.5,s_a=2e-10,r_b=0.5,s_b=2e-10,nGQ=20,n_mu=100,nP=100){
  outA=PPRML_outA_lp(xs_bn,xs_a,r_a,s_a,mu_l,nGQ,n_mu,nP)
  outB=PPRML_outB_lp(xs_bn,xs_b,r_b,s_b,mu_u,nGQ,n_mu,nP)
  sin=  PPRML_sin_lp(xs_bn,xs_a,r_a,s_a,nGQ)
  int=  PPRML_int_lp(xs_bn,xs_a,xs_b,e,r_a,s_a,r_b,s_b,nGQ,n_mu,nP)
  mix=  PPRML_mix_lp(xs_bn,xs_a,xs_b,e,r_a,s_a,r_b,s_b,nGQ,n_mu,nP)
  mxs=c(outA,outB,sin,int,mix)
  logsum_mxs=logsum(mxs)
  post.prob=exp(mxs-logsum_mxs)
  #models = c("OutA","OutB","Sing", "Int","Mix")
  win.model = which.max(post.prob)
  out <- list(post.prob=post.prob,
              win.model=win.model)
  return(out)
}
#PPRML.poi.mix(xs_bn,xs_a,xs_b)


log.pm <- function(x, a, b){
  a.x <- a + sum(x)
  b.x <- b + length(x)
  mu.map <- a.x / b.x
  return(sum(dpois(x, mu.map, log = TRUE)) - diff(dgamma(mu.map, c(a, a.x), c(b, b.x), log = TRUE)))
}

## PPRML.test for all
PPRML.tests <- function(xA, xB, xAB, labels = c("A", "B", "AB"), remove.zeros = FALSE,mu_l=0,mu_u=180,e=0,gamma.pars = c(0.5, 2e-10),nGQ=20,n_mu=100,nP=100,alpha=0.5){
  
  a=r_a=r_b= gamma.pars[1]
  b=s_a=s_b= gamma.pars[2]
  
  if(remove.zeros){
    xA <- xA[xA != 0]
    xB <- xB[xB != 0]
    xAB <- xAB[xAB != 0]
  }
  
  nA <- length(xA)
  nB <- length(xB)
  nAB <- length(xAB)
  if(nA == 0 | nB == 0) stop("not enough data in single sound")
  
  ## how different are the two pure trials?
  
  two.poi.ibf <- Vectorize(function(i, j) return(log.pm(xA[-i],xA[i]+a,1+b) + log.pm(xB[-j],xB[j]+a,1+b) - log.pm(c(xA[-i],xB[-j]),xA[i]+xB[j]+a,2+b)))
  lbf.pure <- mean(c(outer(1:length(xA), 1:length(xB), two.poi.ibf)))
  
  pvls=c(PPRML.poi(xA,nGQ,n_mu,nP,alpha),PPRML.poi(xB,nGQ,n_mu,nP,alpha))
  poimix.res=PPRML.poi.mix(xAB,xA,xB,mu_l,mu_u,e,r_a,s_a,r_b,s_b,nGQ,n_mu,nP)
  out <- list(separation.logBF = lbf.pure,
              post.prob = poimix.res$post.prob,
              win.model = poimix.res$win.model,
              pois.pvalue = pvls,
              samp.sizes = c(nA, nB, nAB))
    return(out)
}
#PPRML.tests(Acounts,Bcounts,ABcounts)

## Poisson analysis customized for spike train triplets

esti.PPRML <- function(trials, spiketimes, frq = c(1100, 742),
                       pos = c(24, -6), on.reward = TRUE, start.time = 0, end.time = 600,
                       match.level = FALSE, AB.eqlevel = FALSE, go.by.soff = FALSE, ...){
  
  require(neuromplex)
  
  attach(trials)
  attach(spiketimes)
  
  timestamps <- split(TIMES, TRIAL2)
  ntrials <- length(timestamps)
  trial.id <- as.numeric(names(timestamps)) ## same as unique(TRIAL2)
  
  ix1 <- TASKID == 8 & A_FREQ == frq[1] & XA == pos[1]
  ix2 <- TASKID == 8 & A_FREQ == frq[2] & XA == pos[2]
  ix3 <- TASKID == 12 & (A_FREQ == frq[1] & B_FREQ == frq[2] & XA == pos[1] & XB == pos[2]) | (A_FREQ == frq[2] & B_FREQ == frq[1] & XA == pos[2] & XB == pos[1])
  
  if(on.reward){
    ix1 <- ix1 & REWARD == 1
    ix2 <- ix2 & REWARD == 1
    ix3 <- ix3 & REWARD == 1
  } 
  
  blev <- sort(unique(B_LEVEL[ix3]))
  targ.lev <- blev[blev > 0]
  lev <- "*"
  if(match.level) {
    if(length(targ.lev) > 1) {
      targ.lev <- max(targ.lev)
      warning("Multiple single sound levels, choosing largest one")
    }
    ix1 <- ix1 & A_LEVEL == targ.lev
    ix2 <- ix2 & A_LEVEL == targ.lev
    lev <- as.character(targ.lev)
  }
  
  if(AB.eqlevel) ix3 <- ix3 & (A_LEVEL == B_LEVEL)
  
  sing1 <- trials[ix1, 1]
  sing2 <- trials[ix2, 1]
  duplx <- trials[ix3, 1]    
  success <- REWARD[ix3]
  
  if(go.by.soff) end.time <- min(SOFF[ix1 | ix2 | ix3])
  
  spike.counter <- function(jj){
    jj1 <- match(jj, trial.id)
    spks <- timestamps[[jj1]]
    return(sum(spks > start.time & spks < end.time))
  }
  
  Acounts <- sapply(sing1, spike.counter)
  Bcounts <- sapply(sing2, spike.counter)
  ABcounts <- sapply(duplx, spike.counter)
  
  detach(trials)
  detach(spiketimes)
  
  
  s1 <- paste(frq[1],"Hz ",pos[1],"deg ",lev,"Db ", sep = "")
  s2 <- paste(frq[2],"Hz ",pos[2],"deg ",lev,"Db ", sep = "")
  dp <- paste("Duplex: ", lev, "Db ", sep = "")
  res=PPRML.tests(Acounts, Bcounts, ABcounts, labels = c(s1, s2, dp), mu_l=0,mu_u=300*end.time/1000,e=0,...)
  return(res)
}
#esti.PPRML(trials,spiketimes)

poi.PRML.VC <- function(fname, cell, data.path = "Data", on.reward = TRUE,
                   match.level = FALSE, AB.eqlevel = FALSE, outfile = "",
                   start = 0, end = 600, go.by.soff = TRUE, remove.zeros = FALSE){
  
  infile1 <- paste(data.path, "/VC/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))
  
  infile2 <- paste(data.path, "/VC/", fname, "_cell", cell, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  FREQS <- unique(c(trials$A_FREQ, trials$B_FREQ))
  alt.freq <- sort(FREQS[FREQS > 0 & FREQS != 742])  
  alt.pos <- c(-24, -6, 6, 24)
  
  
  
  for(fr in alt.freq){
    for(po in alt.pos){
      try({lbf <-esti.PPRML(trials, spiketimes, c(fr, 742), c(po, -144/po), on.reward, start, end, match.level, AB.eqlevel, go.by.soff, remove.zeros)%>%unlist(.);
      cat(fname, cell, c(fr, po, lbf), "\n", file = outfile, append = TRUE)});
    }
  }
}


poi.PRML.JA <- function(fname, on.reward = TRUE, data.path = "Data",
                   match.level = FALSE, AB.eqlevel = FALSE, outfile = "",
                   start = 0, end = 600, remove.zeros = FALSE){
  
  infile1 <- paste(data.path, "/JA/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL"))
  
  infile2 <- paste(data.path, "/JA/", fname, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  FREQS <- unique(c(trials$A_FREQ, trials$B_FREQ))
  alt.freq <- sort(FREQS[FREQS > 0 & FREQS != 742])
  alt.pos <- c(-24, -6, 6, 24)
  
  par(mfcol = c(length(alt.freq),2), mar = c(2,3,3,0) + .1)
  
  for(fr in alt.freq){
    for(po in alt.pos){
      try({lbf <-esti.PPRML(trials, spiketimes, c(fr, 742), c(po, -144/po), on.reward, start, end, match.level, AB.eqlevel, go.by.soff=FALSE, remove.zeros)%>%unlist(.);
      cat(fname, c(fr, po, lbf), "\n", file = outfile, append = TRUE);
      })
    }  
  }
}
set.seed(123)
tic=proc.time()
vc.files.cell1.ls=list.files(path = "./Data/VC/", pattern = "*_cell1_*", all.files = FALSE,
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
vc.files.cell2.ls=list.files(path = "./Data/VC/", pattern = "*_cell2_*", all.files = FALSE,
                             full.names = FALSE, recursive = FALSE,
                             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#ja.files.ls=list.files(path = "./Data/JA/", pattern = "*d.txt", all.files = FALSE,
#                       full.names = FALSE, recursive = FALSE,
#                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
vc.fnames.cell1=str_split(vc.files.cell1.ls,"_cell1_",simplify = TRUE)[,1]
vc.fnames.cell2=str_split(vc.files.cell2.ls,"_cell2_",simplify = TRUE)[,1]
#ja.fnames=str_split(ja.files.ls,".txt",simplify = TRUE)[,1]
ja.fnames=read.table(file = "JA230_list_randomly_interleaved.txt",stringsAsFactors = FALSE)%>%.[,1]
for (vc.name in vc.fnames.cell1){
  poi.PRML.VC(vc.name, 1, data.path = "Data", on.reward = TRUE,
                          match.level = FALSE, AB.eqlevel = FALSE, outfile = "VC_cell1.txt",
                          start = 0, end = 600, go.by.soff = TRUE, remove.zeros = FALSE)
}
for (vc.name in vc.fnames.cell2){
  poi.PRML.VC(vc.name, 2, data.path = "Data", on.reward = TRUE,
              match.level = FALSE, AB.eqlevel = FALSE, outfile = "VC_cell2.txt",
              start = 0, end = 600, go.by.soff = TRUE, remove.zeros = FALSE)
}
for (ja.name in ja.fnames){
  poi.PRML.JA(vc.name, data.path = "Data", on.reward = TRUE,
              match.level = FALSE, AB.eqlevel = FALSE, outfile = "JA_.txt",
              start = 0, end = 600, remove.zeros = FALSE)
}

toc=proc.time()
tictoc=toc-tic