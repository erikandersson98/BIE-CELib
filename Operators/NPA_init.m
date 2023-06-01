function NPA=NPA_init(z,zp,zpp,w,awzp,N,nz)
% NPA = NPA_init(z,zp,zpp,w,awzp,N,nz)
% Returns: NPA - Adjoint NP operator as NxN matrix
  NPA = real(nz./(z.'-z)).*awzp';
  NPA(1:N+1:N^2) = -w.*imag(zpp./zp)/2;
  NPA = NPA/pi;
end