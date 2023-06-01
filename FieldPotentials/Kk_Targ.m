function KkTg = Kk_Targ(ztg,zsc,wzp,kadz)
% KkTg = Kk_Targ(ztg,zsc,wzp,kadz)
% Returns: KkTg - Standard quadrature term for Kk operator, field eval
    KkTg = 1i/2*kadz.*besselh(1,kadz).*imag(wzp.'./(ztg-zsc.'));
end