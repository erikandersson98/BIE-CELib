function Tk = Tk_init(z,wzp,awzp,nz,N,LogC,HypC,k,kadz)
% Tk = Tk_init(z,wzp,awzp,nz,N,LogC,HypC,k,kadz)
% Returns: Tk - Tk operator on boundary as NxN matrix
    temp1 = k*imag(conj(nz).*wzp.')./abs(z-z.');
    temp2 = kadz.^2.*real(-nz./(z-z.')).*imag(wzp.'./(z-z.'));
    Tkz = 1i/2*(besselh(1,kadz).*temp1 + besselh(2,kadz).*temp2);
    Tk0z = (1i/4*k^2-1/(4*pi)*k^2*(2*log(k/2)-2*psi(1)-1)).*awzp;
    Tkz(1:N+1:N^2) = Tk0z;

    TkLz = -1/pi*(besselj(1,kadz).*temp1 + besselj(2,kadz).*temp2);
    TkLz(1:N+1:N^2) = -k^2/(2*pi)*imag(conj(nz).*wzp);

    TkHz = -nz/pi;

    Tk = Tkz + TkLz.*LogC + imag(TkHz.*HypC);
end