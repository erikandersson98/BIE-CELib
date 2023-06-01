function closeCell = Screening_Ellipse(awzp,u,zPan,nPan,ngl,C)
%closeCell = Screening_Ellipse(awzp,u,zPan,nPan,ngl,C)
%Returns: closeCell - Cell array containing the indices of close panels
    S = sum(reshape(awzp,[ngl,nPan])); %arc lengths
    closeCell = cell(nPan,1);
    for i = 1:nPan
        closePts = abs(u-zPan(i))+abs(u-zPan(i+1)) < C*S(i);
        closeCell{i} = find(closePts);
    end
end