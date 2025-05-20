# Simulation Study for the paper

# Checking adequacy of variance function in nonparametric
# regression with unknown mean function


# Data Generating

start_time=Sys.time()

set.seed(9876)
lend=-1
rend=1
total=500

for(a in c(0.5,0.8,1,1.2,1.5))
 {
  power1=power2=matrix(0,nrow=6,ncol=5)
  mtht1=matrix(0,nrow=6,ncol=5)
  mtht2=matrix(0,nrow=6,ncol=5)
  mse=matrix(0,nrow=6,ncol=5)
  nr=1
  for(modl in c(0,1,2,3,4,5))
   {  
    nc=1
    for(n in c(100,300,500,800,1000))
     {
      tht1=rep(0,total)
      tht2=rep(0,total)
      TnCn1=TnCn2=rep(0,total)
      for(k in seq(total))
       {
        x1=runif(n,lend,rend)
        x2=runif(n,lend,rend)
        e=rnorm(n,0,1)
        vx0=0.5+0.5*x1^2+0.5*x2^2
        mx=1+2*x1+3*x2^2
      
        h=a*n^(-1/6)
        
        # Alternative: H1
        vx1=exp(x1+x2+1) 
      
        # Alternative: H2
        dex=(x1^2+x2^2+0.1)^(-1)
        vx2=dex 
      
        # Alternative: H3
        vx3=vx0+dex/sqrt(n)
      
        # Alternative: H4
        vx4=vx0+dex/sqrt(n*h)  
      
        # Alternative: H5
        vx5=vx0+dex/(log(n))  
      
        vx=vx0*(modl==0)+vx1*(modl==1)+vx2*(modl==2)+
           vx3*(modl==3)+vx4*(modl==4)+vx5*(modl==5)
        y=mx+sqrt(vx)*e
      
        w=(log(n)/n)^(1/6)
      
        Kh=function(u)
         {0.75*(1-(u/h)^2)*(abs(u/h)<=1)/h}
        Kw=function(u)
         {0.75*(1-(u/w)^2)*(abs(u/w)<=1)/w}
      
      
        # Estimation of mean function m(x)
      
        x1dif=kronecker(x1,x1,'-')
        x2dif=kronecker(x2,x2,'-')
        Kxx=matrix(Kh(x1dif)*Kh(x2dif),nrow=n)
        fx=rowMeans(Kxx)
        fxy=(Kxx%*%y)/n
        mhat=fxy/fx
      
        # Integration Limits
      
        intlim=function(u1,u2)
         {
          if(((u2-h)>(u1+h))|((u1-h)>(u2+h)))
           {
            return(c(999,999))
           }  else
           {
            intl=max(u1-h,u2-h)
            intr=min(u1+h,u2+h)
            return(c(intl,intr))
           }
         }
      
        n2=n*n
        A1=matrix(0,nrow=n2,ncol=2)
        A2=matrix(0,nrow=n2,ncol=2)
        
        kk=1
        for(ii in seq(n))
         {
          for(jj in seq(n))
           {
            A1[kk,1]=intlim(x1[ii],x1[jj])[1]
            A1[kk,2]=intlim(x1[ii],x1[jj])[2]
            A2[kk,1]=intlim(x2[ii],x2[jj])[1]
            A2[kk,2]=intlim(x2[ii],x2[jj])[2]
            kk=kk+1
           }
         }   
      
        A1[abs(A1)==999]=0 
        A2[abs(A2)==999]=0 
      
        cC1=A1[,1]
        dD1=A1[,2]
      
        cC2=A2[,1]
        dD2=A2[,2]
      
      
        # Estimate theta
      
        x1s=x1^2   
        xx=kronecker(x1,x1,"*")
        xxs=kronecker(x1s,x1s,"*")
        xax=kronecker(x1,x1,"+")
        xaxs=kronecker(x1s,x1s,"+")
      
        AA1=matrix((0.75^2/h^2)*(1-xaxs/h^2+xxs/h^4)*(dD1-cC1)+
                  (0.75^2/h^2)*(2*xax/h^2-2*xx*xax/h^4)*(dD1^2-cC1^2)/2-
                  (0.75^2/h^2)*(2/h^2-(xax^2+2*xx)/h^4)*(dD1^3-cC1^3)/3-
                  (0.75^2/h^2)*(2/h^4)*xax*(dD1^4-cC1^4)/4+
                  (0.75^2/h^6)*(dD1^5-cC1^5)/5,nrow=n,ncol=n)
      
        x2s=x2^2   
        xx=kronecker(x2,x2,"*")
        xxs=kronecker(x2s,x2s,"*")
        xax=kronecker(x2,x2,"+")
        xaxs=kronecker(x2s,x2s,"+")
      
        AA2=matrix((0.75^2/h^2)*(1-xaxs/h^2+xxs/h^4)*(dD2-cC2)+
                   (0.75^2/h^2)*(2*xax/h^2-2*xx*xax/h^4)*(dD2^2-cC2^2)/2-
                   (0.75^2/h^2)*(2/h^2-(xax^2+2*xx)/h^4)*(dD2^3-cC2^3)/3-
                   (0.75^2/h^2)*(2/h^4)*xax*(dD2^4-cC2^4)/4+
                   (0.75^2/h^6)*(dD2^5-cC2^5)/5,nrow=n,ncol=n)
      
        zz=(y-mhat)^2-0.5
      
        X=cbind(x1s*fx^2,x2s*fx^2)
        W=AA1*AA2
        Y=zz*fx^2
        temp=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y
      
        tht1[k]=temp[1]
        tht2[k]=temp[2]
        #cat(tht1[k],tht2[k],"\n")
      
        xi=(y-mhat)^2-0.5-tht1[k]*x1s-tht2[k]*x2s
        Tn=t(xi*fx^2)%*%W%*%(xi*fx^2)/n^2
        Cn=t((xi^2)*(fx^4))%*%diag(W)/n^2
        diag(W)=0
        Gn=(h^2/n^2)*t(xi^2*fx^4)%*%(W^2)%*%(xi^2*fx^4)
        TnCn1[k]<-n*h*(Tn-Cn)/sqrt(2*Gn)
        TnCn2[k]<-n*h*abs(Tn-Cn)/sqrt(2*Gn)
        #cat("TnCn1=",TnCn1[k],"\n\n")
      }  
    
      power1[nr,nc]=round(sum(TnCn1>=qnorm(0.95))/total,3)
      power2[nr,nc]=round(sum(TnCn2>=qnorm(0.975))/total,3)
      mtht1[nr,nc]=round(mean(tht1),3)
      mtht2[nr,nc]=round(mean(tht2),3)
      mse[nr,nc]=round(mean((tht1-0.5)^2+(tht2-0.5)^2),3)
      cat("nr=",nr,"nc=",nc,"\n")
      nc=nc+1
     }
     nr=nr+1
    }
  cat("a=",a,"\n")
  print(power1)
  print(power2)
  #mtht
  #mse
}

end_time=Sys.time()
end_time-start_time