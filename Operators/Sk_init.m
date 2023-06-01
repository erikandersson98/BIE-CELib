function Sk = Sk_init(awzp,N,LogC,k,kadz)
% Sk = Sk_init(awzp,N,LogC,k,kadz)
% Returns: Sk - Sk operator as NxN matrix
    Skz = 1i/2*besselh(0,kadz);

    Sk0 = 1i/2-1/pi*(log(k/2)-psi(1));
    Skz(1:N+1:N^2) = Sk0;

    SkLz = -1/pi*besselj(0,kadz);
    SkLz(1:N+1:N^2) = -1/pi;

    Sk = (Skz+SkLz.*LogC).*awzp.';
end