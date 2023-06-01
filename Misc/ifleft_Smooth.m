function ifleft = ifleft_Smooth(zf,zpf,psc,ztg,maxIt,anglim)
% ifleft = ifleft_Smooth(zf,zpf,psc,ztg,maxIt,anglim)
% Inputs:
% - zf, zpf - curve parametrization and derivative
% - psc - paremeter points for first pass, expects zf(psc(1)) =
%       zf(psc(end))
% - maxIt - max iterations for refining
% - anglim - upper limit of angle between nz and zdiff for resolving point
%       Recommended anglim: pi/8
% ---
% ifleft for smooth curve by refining curve locally until normal vector
% of closest point points almost straight at target point

    nlim = 1e8; %So Matlab doesn't run out of memory if large n
    coslim = cos(anglim);
    if coslim < 0 || coslim >=1
        error('cos(anglim) in interval (0 1) expected')
    end

    Ntg = length(ztg);
    Nsc = length(psc)-1;
    zsc = zf(psc(1:end-1));
    zpsc = zpf(psc(1:end-1));
    nzsc =-1i*zpsc./abs(zpsc);
    ifleft = false(Ntg,1);
    resolved = false(Ntg,1);
    minIndztg = zeros(Ntg,1);
    cosztg = zeros(Ntg,1);
    
    isMin = false(Nsc,1);
    %First pass
    for i = 1:Ntg
        zdiff=ztg(i)-zsc;
        d2 = abs2(zdiff);
        [~,minInd]=min(d2);
        minIndztg(i) = minInd;
        isMin(minInd) = true;
        tempdotprod = real( (ztg(i)-zsc(minInd))*conj( nzsc(minInd)) );
        ifleft(i) = tempdotprod < 0;
        cosztg(i) = tempdotprod/abs(ztg(i)-zsc(minInd));
        resolved(i) = abs(cosztg(i))>coslim;
    end
    j = 0;
    len = Nsc;
    while ~all(resolved) && j < maxIt && 2*len < nlim
        len = 2*len;
        pscNew = zeros(len+1,1);
        pscNew(1:2:end) = psc;
        pscNew(2:2:end-1) = (psc(2:end)+psc(1:end-1))/2;
        psc = pscNew;
        minIndztg(~resolved) = 2*minIndztg(~resolved)-1;
        zsc = zf(psc(1:end-1));
        zpsc = zpf(psc(1:end-1));
        nzsc =-1i*zpsc./abs(zpsc);
        isMin = false(len,1);
        for i = 1:Ntg
            if ~resolved(i)
                if minIndztg(i) == 1
                    itemp = [len;1;2];
                else
                    itemp = (minIndztg(i)-1):(minIndztg(i)+1);
                end
                zsctemp = zsc(itemp);
                nzsctemp = nzsc(itemp);
                zdiff=ztg(i)-zsctemp;
                d2 = abs2(zdiff);
                [~,minInd]=min(d2);
                minIndztg(i) = itemp(minInd);
                isMin(itemp(minInd)) = true;
                tempdotprod = real( (ztg(i)-zsctemp(minInd))*conj(nzsctemp(minInd)) );
                ifleft(i) = tempdotprod < 0;
                cosztg(i) = tempdotprod/abs(ztg(i)-zsctemp(minInd));
                resolved(i) = abs(cosztg(i))>coslim;
            end
        end
        %Removing points which are no longer relevant, change array indx
        tempBool = isMin;
        tempBool = any([ tempBool([len,(1:len-1)]),tempBool,tempBool([2:len,1]) ],2);
        translateindx = cumsum(~tempBool);
        minIndztg(~resolved) = minIndztg(~resolved) - translateindx(minIndztg(~resolved));
        psc = psc([tempBool;true]);
        len = length(psc)-1;
        j = j+1;
    end
    if j == maxIt && ~all(resolved)
        disp('Max iteration reached, all points not resolved')
    elseif 2*len >= nlim && ~all(resolved)
        disp('Process aborted for memory concerns, all points not resolved')
    end
end

function d2 = abs2(z)
    xd=real(z);
    yd=imag(z);
    d2=xd.*xd+yd.*yd;
end