function Kk = Kk_init(z,zp,zpp,w,wzp,N,LogC,kadz)
% Kk = Kk_init(z,zp,zpp,w,wzp,N,LogC,kadz)
% Returns: Kk - Kk operator on boundary as NxN matrix
    temp1 = kadz.*imag(wzp.'./(z-z.'));
    
    Kkz = 1i/2*besselh(1,kadz).*temp1;

    Kk0z = -w/(2*pi).*imag(zpp./zp);
    Kkz(1:N+1:N^2) = Kk0z;

    KkLz = -1/pi*besselj(1,kadz).*temp1;
    KkLz(1:N+1:N^2) = 0;

    Kk = Kkz+(KkLz.*LogC);
end