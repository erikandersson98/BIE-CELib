function CauC = CauC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan)
% CauC = CauC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan)
% Returns: CauC - Cauchy compensation term as NxN matrix
% ---
% Uses Local regularization with interpolation in complex variable
    
    CauC = zeros(N);
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
        Trsb = @(z) (z-(sbb+sba)/2)/sbcc;
        Tr = @(z) (z-(b+a)/2)/cc; %Transformation function
        tgInd = (ngl*(i-1)+1):ngl*i; %Target points index
        sbInd = loopInd((ngl*(i-1)+1):ngl*(i+2)); %Short boundary index
        ztg = z(tgInd);
        ztgtr = Tr(ztg);
        ztgtrsb = Trsb(ztg);
        
        sgn = ones(ngl,1);
        sgn(imag(ztgtrsb) < 0) = -1;
        gamma = sgn*1i;
        argAdd = -sgn*pi/2*1i;
        
        CauCtemp = zp(sbInd).'./(z(sbInd).'-z(tgInd));
        CauCtemp(ngl^2+1:ngl+1:2*ngl^2) = 0;
        diagVect(tgInd) = diagVect(tgInd)- CauCtemp*w(sbInd);
        diagVect(tgInd) = diagVect(tgInd) + argAdd + ...
        log((1-ztgtrsb)./((-1-ztgtrsb).*gamma));
        V=ztgtr.^(0:ngl-1); %Vandermonde matrix
        CauC(tgInd,tgInd) = CauC(tgInd,tgInd) + V*D/V.*wzp(tgInd)/cc;
    end
    CauC(1:N+1:N^2) = CauC(1:N+1:N^2) + diagVect.';
end