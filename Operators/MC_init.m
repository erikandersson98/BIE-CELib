function MC = MC_init(z,wzp,N,CauC)
% MC = MC_init(z,wzp,N,CauC)
% Returns: MC - Cauchy operator as NxN matrix
    MC = wzp.'./(z.'-z);
    MC(1:N+1:N^2) = 0;
    MC = MC+CauC;
    MC = MC/(1i*pi);
end