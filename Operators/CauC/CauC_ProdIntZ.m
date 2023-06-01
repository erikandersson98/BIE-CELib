function CauC = CauC_ProdIntZ(z,wzp,zPan,ngl,N,nPan)
% CauC = CauC_ProdIntZ(z,wzp,zPan,ngl,N,nPan)
% Returns: CauC - Cauchy compensation term as NxN matrix
% ---
% Uses panelwise product integration in complex variable

    CauC = zeros(N);
    %Loop through target panels
    for i = 1:nPan
        tgInd = ngl*(i-1)+1:ngl*i;
        sbiPan = mod(i-2:i+1,nPan)+1; %Index of all sb panels
        %Loop through source panels
        for j = 1:3
            pa = sbiPan(j);
            scInd = ((pa-1)*ngl+1):ngl*pa;
            b = zPan(sbiPan(j+1));
            a = zPan(pa);
            wCcorr = wCinitZ(z(tgInd),z(scInd),wzp(scInd),a,b,j,ngl);
            CauC(tgInd,scInd) = wCcorr;
        end
    end
end

function wCcmp = wCinitZ(ztg,zsc,wzpsc,a,b,j,ngl)
%Calculates local compensation term wCcmp
    c=(1-(-1).^(1:ngl))./(1:ngl);
    cc=(b-a)/2;
    Tr = @(z) (z-(b+a)/2)/cc; %Transformation function
    zsctr = Tr(zsc);
    ztgtr = Tr(ztg);
    if j == 2
        closetg = (1:ngl).';
    else
        closetg = find(abs(ztgtr)<2);
    end
    ztgtrc = ztgtr(closetg);
    wCcmpTemp = wzpsc.'./(zsc.'-ztg(closetg));
    Nc = length(closetg);
    P = zeros(Nc,ngl);
    if j == 2
        wCcmpTemp(1:ngl+1:ngl^2) = 0;
        sgn = ones(Nc,1);
        sgn(imag(ztgtrc) < 0) = -1;
        argAdd = -sgn*pi*1i;
        P(:,1)=argAdd+log((1-ztgtrc)./(-1-ztgtrc)); 
    else
        P(:,1)=log((1-ztgtrc)./(-1-ztgtrc));
    end
    for k=1:ngl-1
        P(:,k+1)=ztgtrc.*P(:,k)+c(k);
    end
    V=zsctr.^(0:ngl-1); %Vandermonde matrix
    wCcorrc = P/V;
    wCcmp = zeros(ngl);
    wCcmp(closetg,:) = wCcorrc-wCcmpTemp;
end