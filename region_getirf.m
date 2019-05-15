function [D,t] = region_getirf(SPM,regions,varargin);
% function [D,t] = region_getirf(SPM,regions,varargin);
% gets the evoked respones from all regions 
% for all onsets for all events 
% SPM: estimated SPM model
% Regions: Cell array of region structures 
% VARARGIN: 
%   pre: How many images Pre event 
%   Post: How many images post event 
%   fig: Print figure? 
% Out: D-structure of individual 
pre=3; 
post=8; 
method='adjusted_mean';
fig=0;

vararginoptions(varargin,{'pre','post','fig'}); 

if (ischar(SPM))
    load(SPM);
end; 

nscans = SPM.nscan;
t=[-pre:post]*SPM.xY.RT; 

% 1. Get ROI Data
[y_raw,y_adj,y_hat,y_res, B] = region_getts(SPM,regions,varargin);

% 2. Get onset information for all blocks
[D,session_onset] = spmj_get_ons_struct(SPM);

for i=1:length(D.block) 
    switch (method) 
        case 'adjusted_mean'
            D.y_adj(i,:)=cut(y_adj,pre,round(D.ons(i)),post,'padding','nan'); 
            D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan'); 
            D.y_res(i,:)=cut(y_res,pre,round(D.ons(i)),post,'padding','nan');             
        case 'evoked' 
    end; 
end; 

if (fig==1)
    traceplot(t,D.y_adj,'linestyle','-'); hold on;
    traceplot(t,D.y_hat,'linestyle',':');hold off;
end;