function [y_raw,y_adj,y_hat,y_res, B,y_filt] = region_getts(spm_file,regions,varargin)
% [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM_FILE,REGIONS,varargin);
% 
% Extracts raw time series for a number of regions
% Refits the model on the mean or first principle component of region
% OUTPUT: 
% y_raw:    Raw time series for each region, a N*num_regions matrix 
% y_adj:    The time series adjusted for nuiscance regressors 
% y_hat:    The prediced time series from regressors of interest
% y_res:    Residual time series 
% B:        The Beta weights for the time regions 
% INPUT: 
% spm_file: Model SPM.mat 
% regions:  Cell array for region strutures, or single region structure 
% VARARGIN: 
%   'stats':  Sufficient statistic used for the area 
%               mean:       Mean of the voxels 
%               whitemean:  Mean of the voxel, spatially weighted by noise
%               pca:        first principle component 
%               any other:  function name of function 
%   'reg_interest',ind  : Regressors that are considered to be regressors
%                           of interest, all other taken out 
%                       default: xX.iH, xX.iC (all indicators + covariates)  
% Written by Tim Verstynen & Joern Diedrichsen (2007)
% -----------------------------------------------------------

% Defaults 
stats='mean'; 
reg_interest=[]; 

vararginoptions(varargin,{'stats','reg_interest'}); 
if nargin < 2 | isempty(spm_file);
    spm_file = spm_get(1,'*.mat','Get SPM File');
end;

if nargin < 1 | isempty(regions)
    regions = spmj_getregions; 
end

if (~iscell(regions))
    dummy{1}=regions; 
    regions=dummy; 
end;
num_regions=length(regions);

% 2. Get the time series
if (ischar(spm_file))
    load(spm_file)
elseif (isstruct(spm_file))
    SPM=spm_file; 
end;

D = region_getdata(SPM.xY.VY,regions);
% R = region_getdata(SPM.VResMS,regions); 

% 3. Calculate sufficient stats on all regions 
for r=1:num_regions 
    switch(stats) 
        case 'mean' 
            y_raw(:,r)=mean(D{r},2); 
        case 'whitemean' 
            white_data=bsxfun(@rdivide,D{r},sqrt(R{r})); 
            y_raw(:,r)=mean(white_data,2); 
        case 'pca' 
            error('not implemented yet.');
            
        otherwise 
            try 
                y_raw(:,r)=eval([stats '(D{r}']); 
            catch 
                error(['Unknown stats function: ' stats]); 
            end; 
    end;
end;

% 3. Get the Design Matrix
if (isempty(reg_interest))
    reg_interest=[SPM.xX.iH SPM.xX.iC]; 
end;

% 5. Do the regression 
if (~isfield(SPM.xX,'pKX')) % SPM not run - 
    y_filt = spm_filter(SPM.xX.K,y_raw); 
    SPM.xX.xKXs.X = spm_filter(SPM.xX.K,SPM.xX.X);
    SPM.xX.pKX = pinv(SPM.xX.xKXs.X); 
    B = SPM.xX.pKX*y_filt;                              %-Parameter estimates
    y_res   = y_filt - SPM.xX.xKXs.X*B;             %-Residuals
    y_hat = SPM.xX.xKXs.X(:,reg_interest)*B(reg_interest,:); %- predicted values 
else
    y_filt = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
    B = SPM.xX.pKX*y_filt;                              %-Parameter estimates
    y_res   = spm_sp('r',SPM.xX.xKXs,y_filt);             %-Residuals
    y_hat = SPM.xX.xKXs.X(:,reg_interest)*B(reg_interest,:); %- predicted values 
end;
y_adj  = y_hat + y_res; 
    
