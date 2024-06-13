function mycifti=region_deform_cifti(R,cifti,varargin)
% function new_cifti=region_deform_cifti(R,cifti)
% Deforms data stored in a cifti structure (returned by cifti_read)
% into the original group space 
% Comments: 
% This function relies on cifti-matlab 
% https://github.com/Washington-University/cifti-matlab
% Extension: should make this work with surface regions as well 
% INPUT: 
%       R: Regions 
%       cfiti: Cifti file information (returned by cift-read)
% OPTIONAL: 
%       vol: Volume file information 
% OUPUT: 
%       new_cifti: deformed cifti
% 2023 joern.diedrichsen@googlemail.com 

vol =[]; 
vararginoptions(varargin,{'vol'}); 

if ~iscell(R)
    R={R}; 
end
num_regions = length(R); 

if length(cifti.diminfo{1}.models) ~=num_regions
    error('number of regions must match number of brainmodels in cifti file'); 
end

mycifti=cifti; 
mycifti.metadata = cifti_metadata_set(mycifti.metadata, 'Provenance', 'region_deform_cifti');

if isempty(vol)
    spm_vol(R{1}.file); 
end

mycifti.cdata = []; 
% make the brainmodel axis
mycifti.diminfo{1}.vol.dims = vol.dim;  % Image dimension 
mycifti.diminfo{1}.vol.sform = vol.mat; % Image affine 

% Now loop over regions and make them into brain models
c_old = 1; % Count old 
c_new = 1; % New count 
for r=1:length(R)
    if size(R{r}.data,1) ~= cifti.diminfo{1}.models{r}.count
        error('Number voxels in brainmodel does not fit region');
    end
    bm.start=c_new; 
    bm.count = size(R{r}.original,1); 
    bm.struct = cifti.diminfo{1}.models{r}.struct;  % Unfortunately, we can't set the structure name to anything arbitary - need to 
    bm.type = cifti.diminfo{1}.models{r}.type;

    % Map data into new space 
    if bm.type=='vox'
        [i,j,k] = spmj_affine_transform(R{r}.original(:,1),R{r}.original(:,2),R{r}.original(:,3),inv(vol.mat));
        bm.voxlist = int16([i j k]'); 
        data = nan(bm.count,size(cifti.cdata,2));
        idx = R{r}.map>0; 
        data(idx,:) = cifti.cdata(R{r}.map(idx),:); % Get the orginal data
    else
        error('Surface ROI not implemented yet'); 
    end 
    c_new = c_new + bm.count; 
    c_old = c_old + cifti.diminfo{1}.models{r}.count; 
    mycifti.cdata = [mycifti.cdata;data];
    mycifti.diminfo{1}.models{r}=bm; 
end    
mycifti.diminfo{1}.length=size(mycifti.cdata,1);