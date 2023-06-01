function [T,W]=GaussTW_gen(N)
% [T,W]=GaussTW_gen(N)
% Returns: T - N-point Gauss-Legendre quadrature nodes
%          W - N-point Gauss-Legendre quadrature weights 
% *** Written by Greg von Winckel - 02/25/2004
  N=N-1;
  N1=N+1; 
  N2=N+2;
  xu=linspace(-1,1,N1)';
  y=cos(pi/(2*N+2)*(1:2:2*N+1))'+0.27/N1*sin(pi*xu*N/N2);
  L=zeros(N1,N2);
  y0=2;
  while max(abs(y-y0))>eps
    L(:,1)=1;
    L(:,2)=y;
    for k=2:N1
      L(:,k+1)=((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
    end        
    Lp=N2*(L(:,N1)-y.*L(:,N2))./(1-y.^2);           
    y0=y;
    y=y0-L(:,N2)./Lp;        
  end
  T=flipud(y);      
  W=2*(N2/N1)^2./((1-y.^2).*Lp.^2);
end