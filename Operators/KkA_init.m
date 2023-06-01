function KkA = KkA_init(z,zp,zpp,w,wzp,nz,N,LogC,kadz)
% KkA = KkA_init(z,zp,zpp,w,wzp,nz,N,LogC,kadz)
% Returns: KkA - Adjoint Kk operator on boundary as NxN matrix
    temp1 = kadz.*imag(-nz.*conj(nz).'.*wzp.'./(z-z.'));
    
    KkAz = 1i/2*besselh(1,kadz).*temp1;

    Kk0z = -w/(2*pi).*imag(zpp./zp);
    KkAz(1:N+1:N^2) = Kk0z;

    KkALz = -1/pi*besselj(1,kadz).*temp1;
    KkALz(1:N+1:N^2) = 0;

    KkA = KkAz+(KkALz.*LogC);
end