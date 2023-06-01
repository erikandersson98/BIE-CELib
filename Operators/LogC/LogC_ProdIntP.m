function LogC = LogC_ProdIntP(zp,awzp,GP,dp,ngl,N,nPan)
% LogC = LogC_ProdIntP(zp,awzp,GP,dp,ngl,N,nPan)
% Returns: LogC - Logarithmic correction term as NxN matrix
% ---
% Uses panelwise product integration in parameter variable
    LogC = zeros(N);
    for i = 1:nPan
        tgInd = ngl*(i-1)+1:ngl*i;
        sbiPan = mod(i-2:i+1,nPan)+1;
        for j = 1:3
            pa = sbiPan(j);
            scInd = ((pa-1)*ngl+1):ngl*pa;
            dpj = dp(pa);
            dp2 = dp(sbiPan(2));
            cc = dpj/2;
            psctr = GP;
            ptgtr = (GP+(2-j))*dpj/dp2+2-j;
            wLcorr = wLinitP(ptgtr,psctr,zp(scInd),zp(tgInd),awzp(scInd),j,ngl,cc);
            LogC(tgInd,scInd) = wLcorr;
        end
    end
end

function wLcorr = wLinitP(ptgtr,psctr,zpsc,zptg,awzpsc,j,ngl,cc)
%Calculates local correction term wLcorr
    c=(1-(-1).^(1:ngl))./(1:ngl);
    if j == 2
        closetg = (1:ngl).';
    else
        closetg = find(abs(ptgtr) < 2);
    end
    ptgtrc = ptgtr(closetg);
    wLcorrTemp = -log(abs(psctr.'-ptgtrc));
    if j == 2
        wLcorrTemp(1:ngl+1:ngl^2) = log(abs(cc*zptg));
    end
    P = zeros(length(closetg),ngl+1);
    Q = zeros(length(closetg),ngl);
    P(:,1)=log(abs(1-ptgtrc)./abs(-1-ptgtrc));
    for k=1:ngl
        P(:,k+1)=ptgtrc.*P(:,k)+c(k);
        if mod(k,2) == 0
            Q(:,k) = (-P(:,k+1) + P(:,1))/k;
        elseif mod(k,2) == 1
            Q(:,k) = (-P(:,k+1) + log(abs( (1-ptgtrc).*(-1-ptgtrc) )))/k;
        end
    end
    V=psctr.^(0:ngl-1);
    wLcorr = zeros(ngl);
    wLcorr(closetg,:) = cc*(Q/V).*abs(zpsc).'./awzpsc.'+wLcorrTemp;
end