function HypC = HypC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan)
% HypC = HypC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan)
% Returns: HypC - Hypersingular compensation term as NxN matrix
% ---
% Uses local regularization with interpolation in complex variable

    HypC = zeros(N);
    diagVect = zeros(N,1);
    D = diag(1:ngl-1,1); %Derivative matrix
    loopInd = [(N-ngl+1):(N),1:N, 1:ngl]';
    %Loop through target panels
    for i = 1:nPan
        sba = zPan(mod(i-2,nPan)+1); %Starting point of sb in plane
        sbb = zPan(mod(i+1,nPan)+1); %End point of sb in plane
        a = zPan(i);
        b = zPan(i+1); 
        cc=(b-a)/2;
        sbcc = (sbb-sba)/2;
        Tr = @(z) (z-(b+a)/2)/cc; %Transformation function
        Trsb = @(z) (z-(sbb+sba)/2)/sbcc;
        tgInd = (ngl*(i-1)+1):ngl*i; %Target points index
        sbInd = loopInd((ngl*(i-1)+1):ngl*(i+2)); %Short boundary index
        zsb = z(sbInd);
        ztg = z(tgInd);
        ztgtr = Tr(ztg);
        ztgtrsb = Trsb(ztg);

        sgn = ones(ngl,1);
        sgn(imag(ztgtrsb) < 0) = -1;
        gamma = sgn*1i;
        argAdd = -sgn*pi/2*1i;

        Hyptemp = zp(sbInd).'./(zsb.'-ztg).^2;
        derivTemp = zp(sbInd).'./(zsb.'-ztg);
        Hyptemp(ngl^2+1:ngl+1:2*ngl^2) = 0;
        derivTemp(ngl^2+1:ngl+1:2*ngl^2) = 0;
        derivVect = argAdd + log((1-ztgtrsb)./((-1-ztgtrsb).*gamma))-derivTemp*w(sbInd);
        diagVect(tgInd) = diagVect(tgInd) - 1/sbcc*(1./(1-ztgtrsb)-1./(-1-ztgtrsb));
        diagVect(tgInd) = diagVect(tgInd) - Hyptemp*w(sbInd);
        
        
        V=ztgtr.^(0:ngl-1); %Vandermonde matrix
        derivF = V*D/V/cc;
        dderivF = V*D^2/V/cc^2;
        HypC(tgInd,tgInd) = HypC(tgInd,tgInd) + derivF.*derivVect;
        HypC(tgInd,tgInd) = HypC(tgInd,tgInd) + dderivF.*wzp(tgInd)/2;
    end
    HypC(1:N+1:N^2) = HypC(1:N+1:N^2) + diagVect.';
end