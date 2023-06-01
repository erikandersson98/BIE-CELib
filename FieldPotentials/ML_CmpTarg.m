function MLCmpTg = ML_CmpTarg(awzpsc,wLcorr)
% MLCmpTg=ML_CmpTarg(awzpsc,wLcorr)
% Returns: MLCmpTg - Compensation term for ML operator, field eval
    MLCmpTg = wLcorr.*awzpsc.';
end