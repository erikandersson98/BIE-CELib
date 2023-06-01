function HypC = HypC_ProdIntZ(z,wzp,zPan,ngl,N,nPan)
% HypC = HypC_ProdIntZ(z,wzp,zPan,ngl,N,nPan)
% Returns: HypC - Hypersingular compensation term as NxN matrix
% ---
% Uses panelwise product integration in complex variable
    HypC = zeros(N);
    for i = 1:nPan
        tgInd = ngl*(i-1)+1:ngl*i;
        sbiPan = mod(i-2:i+1,nPan)+1;
        for j = 1:3
            pa = sbiPan(j);
            scInd = ((pa-1)*ngl+1):ngl*pa;
            b = zPan(sbiPan(j+1));
            a = zPan(pa);
            wHcorr = wHinitZ(z(tgInd),z(scInd),wzp(scInd),a,b,j,ngl);
            HypC(tgInd,scInd) = wHcorr;
        end
    end
end

function wHcmp = wHinitZ(ztg,zsc,wzpsc,a,b,j,ngl)
%Calculates local compensation term wHcmp
    c=(1-(-1).^(1:ngl))./(1:ngl);
    cc=(b-a)/2;
    Tr = @(z) (z-(b+a)/2)/cc;
    ztgtr = Tr(ztg);
    zsctr = Tr(zsc);
    if j == 2
        closetg = (1:ngl).';
    else
        closetg = find(abs(ztgtr) < 2);
    end
    Nc = length(closetg);
    ztgtrc = ztgtr(closetg);
    wHcmpTemp = wzpsc.'./(zsc.'-ztg(closetg)).^2;
    P = zeros(Nc,ngl+1);
    R = zeros(Nc,ngl);
    if j == 2
        wHcmpTemp(1:ngl+1:ngl^2) = 0;
        sgn = ones(Nc,1);
        sgn(imag(ztgtrc) < 0) = -1;
        argAdd = -sgn*pi*1i;
        P(:,1)=argAdd+log((1-ztgtrc)./(-1-ztgtrc));
    else
        P(:,1)=log((1-ztgtrc)./(-1-ztgtrc));
    end
    R(:,1) = -1./(1-ztgtrc)+1./(-1-ztgtrc);
    for k=1:ngl-1
        P(:,k+1)=ztgtrc.*P(:,k)+c(k);
        R(:,k+1) = -1./(1-ztgtrc)+(-1)^k./(-1-ztgtrc) + k*P(:,k);
    end
    V=zsctr.^(0:ngl-1);
    wHcmp = zeros(ngl);
    wHcmp(closetg,:) = (R/V)/cc-wHcmpTemp;
end