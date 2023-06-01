function MCTg = MC_Targ(ztg,zsc,wzpsc)
% MCTg = MC_Targ(ztg,zsc,wzpsc)
% Returns: MCTg - Standard quadrature term for MC operator, field eval
    MCTg = 1/(1i*pi)*wzpsc.'./(zsc.'-ztg);
end