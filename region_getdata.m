function D=region_getdata(V,regions,varargin)
% function D=region_getdata(V,regions)
% Gets data from a number of Volumes 
% For a number of regions 
% INPUT: 
%   V: List of N memory-mapped volumes 
%   regions: (Cell Array) of P region structures, already calculated 
% OPTIONS: 
%   'interp',0-x: Interpolation: 0 nearest neighbor, 1: trilinear, etc 
% OUTPUT: 
%   D: (cell Array - P) of N*numvoxel data matrices 
% 
if (~iscell(regions));
    regions={regions}; 
    oneregion=1;
else
    oneregion=0;
end;
num_images=length(V); 
num_regions=length(regions); 
verbose =0; 
interp=0; %  Interpolation: 0 nearest neighbour, 1: trilinear interpol 

vararginoptions(varargin,{'interp','verbose'}); 

% See if regions are already calculated, if not, do so 
for r=1:num_regions
    if (~isfield(regions{r},'data'))
        regions{r}=region_calcregions(regions{r}); 
    end;
end;

% Loop over images and extract the data 
% Do this for all regions at the same time to increase speed 
mat=zeros(4,4); 
for i=1:num_images
    if (any(V(i).mat(:)~=mat(:)))   % aligment changed or end 
        X=[];Y=[];Z=[];
        for r=1:num_regions
            if (~isempty(regions{r}))
                [x,y,z]=spmj_affine_transform(regions{r}.data(:,1),regions{r}.data(:,2),regions{r}.data(:,3),inv(V(i).mat));
                from(r)=size(X,1)+1;
                X=[X;x];Y=[Y;y];Z=[Z;z];
                to(r)=size(X,1);
            else 
                from(r)=size(X,1)+1;
                to(r)=size(X,1);
            end; 
        end;
    end;
    A(:,i)=spm_sample_vol(V(i),X,Y,Z,interp); % Sampling of volume
    if(verbose & mod(i,10)==0)
        fprintf('.'); 
    end; 
    if(verbose & mod(i,100)==0)
        fprintf('\n'); 
    end; 
end; 

% Sort out the results according to the region 
for r=1:num_regions 
    D{r}=A(from(r):to(r),:)'; 
end;
if (oneregion)
    D=D{1}; 
end