function LogC = LogC_ProdIntZ(z,awzp,nz,zPan,ngl,N,nPan)
% LogC = LogC_ProdIntZ(z,awzp,nz,zPan,ngl,N,nPan)
% Returns: LogC - Logarithmic correction term as NxN matrix
% ---
% Uses panelwise product integration in complex variable
    LogC = zeros(N);
    for i = 1:nPan
        tgInd = ngl*(i-1)+1:ngl*i;
        sbiPan = mod(i-2:i+1,nPan)+1;
        for j = 1:3
            pa = sbiPan(j);
            scInd = ((pa-1)*ngl+1):ngl*pa;
            b = zPan(sbiPan(j+1));
            a = zPan(pa);
            ztg = z(tgInd);
            zsc = z(scInd);
            wLcorr = wLinitZ(ztg,zsc,awzp(scInd),nz(scInd),a,b,j,ngl);
            LogC(tgInd,scInd) = wLcorr;
        end
    end
end

function wLcorr = wLinitZ(ztg,zsc,awzpsc,nzsc,a,b,j,ngl)
%Calculates local correction term wLcorr
    c=(1-(-1).^(1:ngl))./(1:ngl);
    cc=(b-a)/2;
    Tr = @(z) (z-(b+a)/2)/cc;
    ztgtr = Tr(ztg);
    zsctr = Tr(zsc);
    if j == 2
        closetg = (1:ngl).';
    else
        closetg = find(abs(ztgtr)<2);
    end
    Nc = length(closetg);
    ztgtrc = ztgtr(closetg);
    wLcorrTemp = log(abs(zsc.'-ztg(closetg)));
    P = zeros(Nc,ngl+1);
    Q = zeros(Nc,ngl);
    if j == 2
        wLcorrTemp(1:ngl+1:ngl^2) = 0;
        sgn = ones(Nc,1);
        sgn(imag(ztgtrc) < 0) = -1;
        argAdd = -sgn*pi*1i;
        P(:,1)=argAdd+log((1-ztgtrc)./(-1-ztgtrc));
    else
        P(:,1)=log((1-ztgtrc)./(-1-ztgtrc));
    end
    for k=1:ngl
        P(:,k+1)=ztgtrc.*P(:,k)+c(k);
        if mod(k,2) == 0
            Q(:,k) = (-P(:,k+1) + P(:,1))/k;
        elseif mod(k,2) == 1
            Q(:,k) = (-P(:,k+1) + ...
                    log((1-ztgtrc).*(-1-ztgtrc)))/k;
        end
    end
    V=zsctr.^(0:ngl-1);
    wLcorr = zeros(ngl);
    wLcorr(closetg,:) = imag(cc*(Q/V).*conj(nzsc).')./awzpsc.'+log(abs(cc))-wLcorrTemp;
end