function NP=NP_init(z,zp,zpp,w,wzp,N)
% NP = NP_init(z,zp,zpp,w,wzp,N)
% Returns: NP - Neumann-Poincare operator as NxN matrix
% Neumann--Poincaré operator = 2*double-layer operator  
  NP = imag(wzp.'./(z-z.'));
  NP(1:N+1:N^2) = -w.*imag(zpp./zp)/2;
  NP = NP/pi;
end