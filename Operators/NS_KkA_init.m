function KkA = NS_KkA_init(z,wzp,nz,N,LogC,CauC,kadz)
% KkA = NS_KkA_init(z,wzp,nz,N,LogC,CauC,kadz)
% Returns: KkA - Adjoint Kk operator on boundary as NxN matrix
%---
% Variant of KkA_init which works on non-smooth boundaries
    temp1 = kadz.*imag(-nz.*conj(nz).'.*wzp.'./(z-z.'));
    
    KkAz = 1i/2*besselh(1,kadz).*temp1;

    Kk0z = 0;
    KkAz(1:N+1:N^2) = Kk0z;

    KkALz = -1/pi*besselj(1,kadz).*temp1;
    KkALz(1:N+1:N^2) = 0;

    KkACz = 1/pi*nz.*nz';
    
    KkA = KkAz+(KkALz.*LogC) + imag(KkACz.*CauC);
end