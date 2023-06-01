function SkTg = Sk_Targ(awzp,kadz)
% SkTg = Sk_Targ(awzp,kadz)
% Returns: SkTg - Standard quadrature term for Sk operator, field eval
    Skz = 1i/2*besselh(0,kadz);
    SkTg = Skz.*awzp.';
end