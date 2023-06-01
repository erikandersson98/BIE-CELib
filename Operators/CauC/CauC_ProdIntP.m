function CauC = CauC_ProdIntP(zp,zpp,GP,GW,w,dp,ngl,N,nPan)
% CauC = CauC_ProdIntP(zp,zpp,GP,GW,w,dp,ngl,N,nPan)
% Returns: CauC - Cauchy compensation term as NxN matrix
% ---
% Uses panelwise product integration in parameter variable

    CauC = zeros(N);
    for i = 1:nPan
        tgInd = ngl*(i-1)+1:ngl*i;
        sbiPan = mod(i-2:i+1,nPan)+1;
        for j = 1:3
            pa = sbiPan(j);
            scInd = ((pa-1)*ngl+1):ngl*pa;
            dpj = dp(pa);
            dp2 = dp(sbiPan(2));
            psctr = GP;
            ptgtr = (GP+(2-j))*dpj/dp2+2-j;
            wCcmp = wCinitP(ptgtr,psctr,GW,j,ngl);
            CauC(tgInd,scInd) = wCcmp;
        end
    end
    CauC(1:N+1:N^2) = CauC(1:N+1:N^2).' + (w.*zpp)./(2*zp);
end

function wCcmp = wCinitP(ptgtr,psctr,GW,j,ngl)
%Calculates local compensation term wCcmp
    c=(1-(-1).^(1:ngl))./(1:ngl);
    if j == 2
        closetg = (1:ngl).';
    else
        closetg = find(abs(ptgtr) < 2);
    end
    ptgtrc = ptgtr(closetg);
    wCcmpTemp = -GW.'./(psctr.'-ptgtrc);
    if j == 2
        wCcmpTemp(1:ngl+1:ngl^2) = 0;
    end
    P = zeros(length(closetg),ngl);
    P(:,1)=log(abs(1-ptgtrc)./abs(-1-ptgtrc));
    for k=1:ngl-1
        P(:,k+1)=ptgtrc.*P(:,k)+c(k);
    end
    V=psctr.^(0:ngl-1);
    wCcmp = zeros(ngl);
    wCcmp(closetg,:) = P/V+wCcmpTemp;
end