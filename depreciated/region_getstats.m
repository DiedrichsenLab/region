function [T,refstruct]=region_getstats(varargin)
%  function T=spmj_regionstat('images',image_names,'regions',regions,...);
%       Extracts the values of the images over the regions
%       defined in the volumes regionname
%       Extracts the values of images/metric over the regions or points
% 'images': A cell array of cell arrays/char arrays 
%           of file names or Volume structures 
% 'deformations': Deformations from atlas space into original space. For
%           caret, this can be deformed_fiducial metric file, if not given,
%           no deformation is assumed
% 'regions': cell array of region-structures (see region_getregions)
% 'ignore_zeros', 1/0: Flag to ignore zeros within regions
% OUTPUT: 
%   A data structure T with fields: 
%       T.region     : number of region in structure    
%       T.regionname : region name (if present) 
%       T.image      : Image number 
%       T.image_name : File name of Image 
%       T.mean       : mean over the region 
T.region=[];
T.region_name={};
T.image=[];
T.image_name={};
T.reg_size=[];
T.mean=[];

images=[];
regions=[];
deformations=[];
mask=[];
ignore_zeros=0;

vararginoptions(varargin,{'images','regions','deformations','mask','ignore_zeros'});

if (ignore_zeros)
    T.non_zero=[];
end;

% -------------------------------------
% Prompt with GUI if arguments are empty
if (isempty(regions))
    regions=region_getregions();
end;
if (~iscell(regions))
    regions={regions};
end;

if (isempty(images))
    images=spm_get(inf,'.img',sprintf('Select images',s));
end;
if (iscell(images)) 
    images=char(images);
end;

if (isstruct(images))
    V=images;
    images=vertcat(V.fname);
else
    V=spm_vol(images);
end;


% -----------------------------------------
% APPLY ACTIVITY-MASK TO ROIS?
if (~isempty(mask))
    MAV=spm_vol(mask);
    MASK=[];
    T.in_mask=[];
    T.mask={};
    ma=0;num_masks=size(MAV,1);
end;

% -----------------------------------------
% [prepare regions],
num_images=size(images,1);
num_regions=length(regions);
for r=1:num_regions
    if (~isfield(regions{r},'data'))
        regions{r}=region_calcregions(regions{r});
    end;
end;

D=region_getdata(V,regions);
for i=1:num_images 
    for r=1:num_regions
        T.image(end+1,1)=i;
        T.image_name{end+1,1}=deblank(images(i,:));
        T.region(end+1,1)=r;
        if (isfield(regions{r},'name'));
            T.region_name{end+1,1}=regions{r}.name;
        end;
        T.reg_size(end+1,1)=size(regions{r}.data,1); 
        T.mean(end+1,1)=mean(D{r}(i,:),2);
    end;
end;
