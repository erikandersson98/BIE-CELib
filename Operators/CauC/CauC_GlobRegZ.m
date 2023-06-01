function CauC = CauC_GlobRegZ(z,zp,w,wzp,zPan,ngl,N,nPan)
% CauC = CauC_GlobRegZ(z,zp,w,wzp,zPan,ngl,N,nPan)
% Returns: CauC - Cauchy compensation term as NxN matrix
% ---
% Uses Global regularization with interpolation in complex variable

    CauC = zp.'./(z.'-z);
    CauC(1:N+1:N^2) = 0;
    CauC = -diag(CauC*w);
    D = diag(1:ngl-1,1); %Derivative matrix
    %Loop through target panels
    for i = 1:nPan
        a = zPan(i);
        b = zPan(i+1);
        cc=(b-a)/2;
        Tr = @(z) (z-(b+a)/2)/cc;
        tgInd = (ngl*(i-1)+1):ngl*i;
        ztg = z(tgInd);
        ztgtr = Tr(ztg);
        wzptg = wzp(tgInd);
        V=ztgtr.^(0:ngl-1);
        CauC(tgInd,tgInd) = CauC(tgInd,tgInd) + V*D/V.*wzptg/cc;
    end
    CauC = CauC + eye(N)*(1i*pi);
end