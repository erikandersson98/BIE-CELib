function KkCmpTg = Kk_CmpTarg(ztg,zsc,wzpsc,kadz,wLcorr,wCcmp)
% KkCmpTg = Kk_CmpTarg(ztg,zsc,wzpsc,kadz,wLcorr,wCcmp)
% Returns: KkCmpTg - Compensation term for Kk operator, field eval.
    KkLz = -1/pi*kadz.*besselj(1,kadz).*imag(wzpsc.'./(ztg-zsc.'));
    KkCz = -1/pi;
    KkCmpTg = KkLz.*wLcorr + imag(KkCz.*wCcmp);
end