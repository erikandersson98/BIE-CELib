function wLcorr = wLcorr_Field(a,b,c,ztg,zsc,awzpsc,ifleft,nz)
% wLcorr = wLcorr_Field(a,b,c,ztg,zsc,awzpsc,ifleft,nz)
% Returns: wLcorr - Correction weight for field evaluation, Logarithmic singularity
  Ng=length(zsc);
  Nt=length(ztg);
  cc=(b-a)/2;
  ztgtr=(ztg-(b+a)/2)/cc;
  zsctr=(zsc-(b+a)/2)/cc;

  gam=-1i*ones(Nt,1);
  gam(ifleft)=1i;
  loggam=-1i*pi/2*ones(Nt,1); % log(gam)
  loggam(ifleft)=1i*pi/2;

  P = zeros(Nt,Ng+1);
  Q = zeros(Nt,Ng);
  P(:,1)=loggam+log((1-ztgtr)./(gam.*(-1-ztgtr)));
  for k=1:Ng
    P(:,k+1)=ztgtr.*P(:,k)+c(k);
    if mod(k,2) == 0
        Q(:,k) = (-P(:,k+1) + P(:,1))/k;
    elseif mod(k,2) == 1
        Q(:,k) = (-P(:,k+1) + log((1-ztgtr).*(-1-ztgtr)))/k;
    end
  end
  V=zsctr.^(0:Ng-1); %Vandermonde matrix
  wLcorr=real(cc*(Q/V)./(1i*nz).'./awzpsc.')-log(abs(zsc.'-ztg)/abs(cc));
end