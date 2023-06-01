function HypC = HypC_LocRegP(z,zp,zpp,w,wzp,GP,dp,zPan,ngl,N,nPan)
% HypC = HypC_LocRegP(z,zp,zpp,w,wzp,GP,dp,zPan,ngl,N,nPan)
% Returns: HypC - Hypersingular compensation term as NxN matrix
% ---
% Uses local regularization with interpolation in parameter variable

    HypC = zeros(N);
    diagVect = zeros(N,1);
    D = diag(1:ngl-1,1); %Derivative matrix
    loopInd = [(N-ngl+1):(N),1:N, 1:ngl]';
    %Loop through target panels
    for i = 1:nPan
        temp = dp(i)/2;
        sba = zPan(mod(i-2,nPan)+1); %Starting point of sb in plane
        sbb = zPan(mod(i+1,nPan)+1); %End point of sb in plane
        sbcc = (sbb-sba)/2;
        Trsb = @(z) (z-(sbb+sba)/2)/sbcc;
        tgInd = (ngl*(i-1)+1):ngl*i; %Target points index
        sbInd = loopInd((ngl*(i-1)+1):ngl*(i+2)); %Short boundary index
        zsb = z(sbInd);
        ztg = z(tgInd);
        ztgtrsb = Trsb(ztg);
        zptg = zp(tgInd);
        zpptg = zpp(tgInd);
        
        sgn = ones(ngl,1);
        sgn(imag(ztgtrsb) < 0) = -1;
        gamma = sgn*1i;
        argAdd = -sgn*pi/2*1i;
        
        Hyptemp = zp(sbInd).'./(zsb.'-ztg).^2;
        derivTemp = zp(sbInd).'./(zsb.'-ztg);
        Hyptemp(ngl^2+1:ngl+1:2*ngl^2) = 0;
        derivTemp(ngl^2+1:ngl+1:2*ngl^2) = 0;
        derivVect = (argAdd + log((1-ztgtrsb)./((-1-ztgtrsb).*gamma))-derivTemp*w(sbInd));
        diagVect(tgInd) = diagVect(tgInd) - 1/sbcc*(1./(1-ztgtrsb)-1./(-1-ztgtrsb));
        diagVect(tgInd) = diagVect(tgInd) - Hyptemp*w(sbInd);
        
        
        V=GP.^(0:ngl-1); %Vandermonde matrix
        derivF = V*D/V/temp./zptg;
        dderivF = V*D^2/V/temp^2./zptg.^2-zpptg.*derivF./zptg.^2;
        HypC(tgInd,tgInd) = HypC(tgInd,tgInd) + derivF.*derivVect;
        HypC(tgInd,tgInd) = HypC(tgInd,tgInd) + dderivF.*wzp(tgInd)/2;
    end
    HypC(1:N+1:N^2) = HypC(1:N+1:N^2) + diagVect.';
end