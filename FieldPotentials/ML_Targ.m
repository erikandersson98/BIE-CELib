function MLTg = ML_Targ(ztg,zsc,awzp)
% MLTg = ML_Targ(ztg,zsc,awzp)
% Returns: MLTg - Standard quadrature term for ML operator, field eval
    MLTg= log(abs(zsc.'-ztg)).*awzp.'; 
end