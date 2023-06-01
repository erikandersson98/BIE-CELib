function closeCell = Screening_dlim(z,awzp,u,nPan,ngl,dlim)
% closeCell = Screening_dlim(z,awzp,u,nPan,ngl,dlim)
% Returns: closeCell - Cell array containing the indices of close panels
    Nt = length(u);
    closeCell = cell(nPan,1);
    S = sum(reshape(awzp,[ngl,nPan])).';
    S2 = S.^2;
    dlim2 = dlim^2;
    for i = 1:Nt
        zdiff=u(i)-z;
        xd=real(zdiff);
        yd=imag(zdiff);
        d2=xd.*xd+yd.*yd;
        d2=reshape(d2,ngl,nPan);
        d=min(d2)'./S2;
        closePans=d<dlim2;
        if nnz(closePans)>0
            closeInd=find(closePans)';
            for i2=closeInd
                closeCell{i2}=[closeCell{i2} i];
            end
        end
    end

end