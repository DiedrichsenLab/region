function region_plotcon(SPM,R,numc,varargin);
if (isempty(SPM))
    SPM=evalin('caller','SPM');
end;
if (isempty(R))
    hReg=evalin('caller','hReg');
    loc=spm_XYZreg('GetCoords',hReg)
    R={region('point',loc',sprintf('%2.0f %2.0f %2.0f',loc(1),loc(2),loc(3)))}; 
end;
if (isstruct(R))
    R={R}; 
end;
numreg=length(R); 
Y=region_getdata(vertcat(SPM.xCon(numc).Vcon),R); 
for i=1:numreg
    subplot(numreg,1,i); 
    barplot(numc',mean(Y{i},2)); 
    title(R{i}.name); 
end; 
