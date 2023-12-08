function R=region_deformation(R,deffile,varargin)
% function R=spmj_deformation(R,deffile,varargin)
% Deforms regions of interest from atlas into individual space, retaining
% the mapping
% INPUT:
%      R: Single region or cell array of regions
%      deffile: Deformation info
%           For elastic (old): 'paramfile'
%           For dartel:{'u_a_<file>','Affine_file'}
%           in general: {Def,mat} (from spmdefs)
% VARARGIN:
%      'vol':  Volume or file with dim,mat to define voxelspace
%      'mask': Volume or file to define mask - if mask is given, it is also
%               used to define the voxel space 
% RETURN: 
%       region structure in individual space
%       
% joern.diedrichsen@googlemail.com
mask = [];
vol = [];
vararginoptions(varargin,{'mask','vol'});


if (~iscell(R))
    Rnew{1}=R;
    R=Rnew;
end;

% Get the deformation
if iscell(deffile)
    if ischar(deffile{1})
        [Def,mat]=spmdefs_get_dartel(deffile{1},deffile{2});
    elseif (isnumeric(deffile{1}))
        Def = deffile{1};
        mat = deffile{2};
    else
        error('deformation cell array needs to contain filenames or matrices');
    end
elseif ischar(deffile)
    [Def,mat]=spmdefs_get_sn2def(paramfile);
else
    error('deformation cell array needs to cell array or a single filename (for all normalization');
end

% Apply the deformation:
for r=1:length(R)
    R{r}.original = R{r}.data;
    [R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3)]=spmdefs_transform(Def,mat,R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3));
end

% Resample region into new space if vol/mask is given
if (~isempty(vol) && ~isempty(mask))
    error('Only provide vol or mask - not both');
elseif (~isempty(vol))
    mimg = vol;
elseif (~isempty(mask))
    mimg = mask;
else
    mimg = [];
end

if (~isempty(mimg))
    
    % Get the image
    if isstr(mimg)
        mimg = spm_vol(mimg);
    end
    
    % For each region
    for r=1:length(R)
        % Find the available voxels
        [i,j,k]=spmj_affine_transform(R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3),inv(mimg.mat));
        vox = [i j k];
        vox = round(vox);
        good = ~any(vox<[1 1 1] | (vox > mimg.dim),2);

        % If mask, only use voxels within mask
        if (~isempty(mask))
            good(good) = spm_sample_vol(mimg,vox(good,1),vox(good,2),vox(good,3),1)>0;
        end
        
        % Boil down to the unique voxel, retain mapping 
        R{r}.deform = R{r}.data;
        [R{r}.data,~,iC]=unique(vox(good,:),'rows');
        R{r}.map = zeros(size(good,1),1,'uint32');
        R{r}.map(good) = iC;
   
        % Return the coordinates to world coordinates
        [R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3)]=spmj_affine_transform(R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3),mimg.mat);
    
    end
end

% If single region - restore 
if (exist('Rnew','var'))
    R=R{1};
end
