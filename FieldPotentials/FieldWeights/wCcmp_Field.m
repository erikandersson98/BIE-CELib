function wCcmp = wCcmp_Field(a,b,c,ztg,zsc,wzpsc,ifleft)
% wCcmp = wCcmp_Field(a,b,c,ztg,zsc,wzpsc,ifleft)
% Returns: wCcmp - Compensation weight for Cauchy-singularity
  Ng=length(zsc);
  Nt=length(ztg);
  cc=(b-a)/2;
  ztgtr=(ztg-(b+a)/2)/cc;
  zsctr=(zsc-(b+a)/2)/cc;
  
  gam=-1i*ones(Nt,1);
  loggam=-0.5i*pi*ones(Nt,1); % log(gam)
  gam(ifleft)=1i;
  loggam(ifleft)=0.5i*pi;

  P=zeros(Nt,Ng);
  P(:,1)=loggam+log((1-ztgtr)./(gam.*(-1-ztgtr))); 
  for k=1:Ng-1
    P(:,k+1)=ztgtr.*P(:,k)+c(k);
  end
  V = zsctr.^(0:Ng-1);
  wCcmp=P/V-(wzpsc.'./(zsc.'-ztg));
end