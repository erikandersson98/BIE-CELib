function MH = MH_init(z,wzp,N,HypC)
% MH = MH_init(z,wzp,N,HypC)
% Returns: MH - Hypersingluar operator as NxN matrix
    MH = wzp.'./(z.'-z).^2;
    MH(1:N+1:N^2) = 0;
    MH = MH+HypC;
    MH = MH/(1i*pi);
end