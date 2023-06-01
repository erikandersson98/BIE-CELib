function MCCmpTg = MC_CmpTarg(wCcmp)
% MCCmpTg = MC_CmpTarg(wCcmp)
% Returns: MCCmpTg - Compensation term for MC operator, field eval.
    MCCmpTg = 1/(1i*pi)*wCcmp;
end