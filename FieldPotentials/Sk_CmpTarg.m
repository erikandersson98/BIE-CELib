function SkCmpTg = Sk_CmpTarg(awzpsc,kadz,wLcorr)
% SkCmpTg = Sk_CmpTarg(awzpsc,kadz,wLcorr)
% Returns: SkCmpTg - Compensation term for Sk operator, field eval.
    SkLz = -1/pi*besselj(0,kadz);
    SkCmpTg = SkLz.*wLcorr.*awzpsc.';
end