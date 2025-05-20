# Simulation Study for the paper

# Checking adequacy of variance function in nonparametric
# regression with unknown mean function

# Simulation Studies for Table 2: Weight function
#     dG(x)=I[-2,2], dphi(x)=f_w^{-6}(x)I[-2,2]dx 


# Data Gegenerating

Thet=matrix(0,nrow=500,ncol=5)


for(a in c(0.5,0.8,1,1.2,1.5))
{
  set.seed(9876)
  lend=-1
  rend=1
  total=500
  power1=power2=matrix(0,nrow=6,ncol=5)
  mtht=matrix(0,nrow=6,ncol=5)
  mse=matrix(0,nrow=6,ncol=5)
  nr=1
  for(modl in c(0))#,1,2,3,4,5))
  {  
    nc=1
    for(n in c(100,200,300,400,500))
    {
      tht=rep(0,total)
      TnCn1=TnCn2=rep(0,total)
      for(k in seq(total))
      {
        x=runif(n,lend,rend)
        e=rnorm(n,0,1)
        vx0=1+0.5*x^2
        mx=1+2*x+3*x^2
        
        # Alternative: H1
        vx1=exp(x+1) 
        
        # Alternative: H2
        dex=(x^2+0.1)^(-1)
        vx2=dex 
        
        # Alternative: H3
        vx3=vx0+dex/sqrt(n)
        
        # Alternative: H4
        h=a*n^(-1/5)
        vx4=vx0+dex/sqrt(n*h^(0.5))  
        
        # Alternative: H5
        vx5=vx0+dex/(log(n))  
        
        vx=vx0*(modl==0)+vx1*(modl==1)+vx2*(modl==2)+
          vx3*(modl==3)+vx4*(modl==4)+vx5*(modl==5)
        y=mx+sqrt(vx)*e
        
        h=a*n^(-1/5)
        w=(log(n)/n)^(1/5)
        
        Kh=function(u)
        {
          0.75*(1-(u/h)^2)*(abs(u/h)<=1)/h
        }
        Kw=function(u)
        {
          0.75*(1-(u/w)^2)*(abs(u/w)<=1)/w
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
            intv = intlim(x[ii],x[jj])
            
            if(((-2)>intv[2])|(intv[1]>2))
            {
              intv[1]=999
              intv[2]=999
            }  else
            {
              intv[1]=max(-2,intv[1])
              intv[2]=min(2,intv[2])
            }
            
            A[kk,1]=intv[1]
            A[kk,2]=intv[2]
            kk=kk+1
          }
        }   
        
        cC=A[,1]
        dD=A[,2]
        
        
        
        # The computation of A for general weight  w(x)
       
        ind2=which((cC!=999))
        A0=A1=A2=A3=A4=rep(0,n2)
        for(jj in ind2)
        {
          byv=(dD[jj]-cC[jj])/50
          xseq=seq(cC[jj],dD[jj],by=byv)
          xX=kronecker(xseq,x,"-")
          KxX=matrix(Kw(xX),nrow=n)
          Wx=(colMeans(KxX))^(-6)
          A0[jj]=sum(Wx*byv)
          A1[jj]=sum(xseq*Wx)*byv
          A2[jj]=sum(xseq^2*Wx)*byv
          A3[jj]=sum(xseq^3*Wx)*byv
          A4[jj]=sum(xseq^4*Wx)*byv
        }
        
        # Estimate theta
        
        xx=kronecker(x,x,"*")
        x2=x^2   
        xx2=kronecker(x2,x2,"*")
        xax=kronecker(x,x,"+")
        xax2=kronecker(x2,x2,"+")
        
        AA=matrix((0.75^2/h^2)*(1-xax2/h^2+xx2/h^4)*A0+
                    (0.75^2/h^2)*(2*xax/h^2-2*xx*xax/h^4)*A1-
                    (0.75^2/h^2)*(2/h^2-(xax^2+2*xx)/h^4)*A2-
                    (0.75^2/h^2)*(2/h^4)*xax*A3+
                    (0.75^2/h^6)*A4,nrow=n,ncol=n)
        
        
        zz=(y-mhat)^2-1
        tht[k]=t(zz*fx^2)%*%AA%*%(x^2*fx^2)/(t(x^2*fx^2)%*%AA%*%(x^2*fx^2))
        cat(tht[k],"\n")
        
        xi=(y-mhat)^2-1-tht[k]*x2
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
      if(a==0.8)
      {
        Thet[,nc] = tht 
      }
      nc=nc+1
    }
    nr=nr+1
  }
  cat("a=",a,"\n")
  #print(power1)
  #print(power2)
  print(mtht[1,])
  print(mse[1,])
}












