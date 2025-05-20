# Simulation Study for the paper

# Checking adequacy of variance function in nonparametric
# regression with unknown mean function


# Data Gegenerating


  set.seed(9876)
  lend=0
  rend=1
  total=500
  power1=power2=matrix(0,nrow=3,ncol=3)
  mtht=matrix(0,nrow=3,ncol=3)
  mse=matrix(0,nrow=3,ncol=3)
  nr=1
  for(modl in c(0,1,2))
  {  
    nc=1
    for(n in c(50,100,200))
    {
      tht=rep(0,total)
      TnCn1=TnCn2=rep(0,total)
      for(k in seq(total))
      {
        x=runif(n,lend,rend)
        e=rnorm(n,0,1)
        
        mx=1+sin(x)
       
        cc=0*(modl==0)+0.5*(modl==1)+1*(modl==2)
        
        vx=0.5*exp(cc*x)
        
        y=mx+sqrt(vx)*e
        
        h=n^(-1/5)
        
        Kh=function(u)
         {
          0.75*(1-(u/h)^2)*(abs(u/h)<=1)/h
         }
       
        # Estimation of mean function m(x)
        
        xdif=kronecker(x,x,'-')
        Kxx=matrix(Kh(xdif),nrow=n)
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
        A=matrix(0,nrow=n2,ncol=2)
        kk=1
        for(ii in seq(n))
         {
          for(jj in seq(n))
           {
            A[kk,1]=intlim(x[ii],x[jj])[1]
            A[kk,2]=intlim(x[ii],x[jj])[2]
            kk=kk+1
           }
         }   
        
        A[abs(A)==999]=0 
        dD=A[,2]
        cC=A[,1]
        
        ind2=which(abs(cC-dD)>0)
        
        for(ii in ind2)
        {
          if(1<=cC[ii]) 
          {
            cC[ii]=0 
            dD[ii]=0
          } else
            if((0<cC[ii])&(cC[ii]<1)&(1<dD[ii]))
            {
              cC[ii]=cC[ii]
              dD[ii]=1
            } else
              if((0<cC[ii])&(dD[ii]<1))
              {
                cC[ii]=cC[ii]
                dD[ii]=dD[ii]
              } else
                if((cC[ii]<0)&(1<dD[ii]))  
                {
                  cC[ii]=0
                  dD[ii]=1
                }  else
                  if((cC[ii]<0)&(0<dD[ii])&(dD[ii]<1))
                  {
                    cC[ii]=0
                    dD[ii]=dD[ii]
                  }else
                  {
                    cC[ii]=0
                    dD[ii]=0
                  }   
         }
        
        
        # Estimate theta
        
        x2=x^2   
        xx=kronecker(x,x,"*")
        xx2=kronecker(x2,x2,"*")
        xax=kronecker(x,x,"+")
        xax2=kronecker(x2,x2,"+")
        
        AA=matrix((0.75^2/h^2)*(1-xax2/h^2+xx2/h^4)*(dD-cC)+
                    (0.75^2/h^2)*(2*xax/h^2-2*xx*xax/h^4)*(dD^2-cC^2)/2-
                    (0.75^2/h^2)*(2/h^2-(xax^2+2*xx)/h^4)*(dD^3-cC^3)/3-
                    (0.75^2/h^2)*(2/h^4)*xax*(dD^4-cC^4)/4+
                    (0.75^2/h^6)*(dD^5-cC^5)/5,nrow=n,ncol=n)
        
        zz=(y-mhat)^2
        tht[k]=t(zz*fx^2)%*%AA%*%(fx^2)/(t(fx^2)%*%AA%*%(fx^2))
        cat(tht[k],"\n")
        
        xi=(y-mhat)^2-tht[k]
        Tn=t(xi*fx^2)%*%AA%*%(xi*fx^2)/n^2
        Cn=t((xi^2)*(fx^4))%*%diag(AA)/n^2
        diag(AA)=0
        Gn=(h/n^2)*t(xi^2*fx^4)%*%(AA^2)%*%(xi^2*fx^4)
        TnCn1[k]<-n*h^(0.5)*(Tn-Cn)/sqrt(2*Gn)
        TnCn2[k]<-n*h^(0.5)*abs(Tn-Cn)/sqrt(2*Gn)
      }  
      
      power1[nr,nc]=round(sum(TnCn1>=qnorm(0.95))/total,3)
      power2[nr,nc]=round(sum(TnCn2>=qnorm(0.975))/total,3)
      mtht[nr,nc]=round(mean(tht),3)
      mse[nr,nc]=round(mean((tht-0.5)^2),3)
      #cat("nr=",nr,"nc=",nc,"\n")
      nc=nc+1
    }
    nr=nr+1
  }
  print(power1)
  print(power2)
  # mtht
  # mse

