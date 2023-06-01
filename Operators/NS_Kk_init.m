function Kk = NS_Kk_init(z,wzp,N,LogC,CauC,kadz)
% Kk = NS_Kk_init(z,wzp,N,LogC,CauC,kadz)
% Returns: Kk - Kk operator on boundary as NxN matrix
% ---
% Variant of Kk_init which works on non-smooth boundaries
    temp1 = kadz.*imag(wzp.'./(z-z.'));
    
    Kkz = 1i/2*besselh(1,kadz).*temp1;
    Kkz(1:N+1:N^2) = 0;

    KkLz = -1/pi*besselj(1,kadz).*temp1;
    KkLz(1:N+1:N^2) = 0;

    KkCz = -1/pi;

    Kk = Kkz+(KkLz.*LogC)+imag(KkCz.*CauC);
end