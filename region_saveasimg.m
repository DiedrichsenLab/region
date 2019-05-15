function R=region_saveasimg(R,Vol,varargin)
% function R=region_saveasimg(R,Vol,varargin)
% Saves a region field as a *.nii file in the same
% resolution and orientation as Vol
% VARGINOPTION
% 'name',name     : Output name of the Volum
% 'z',z           : varaible that gives brightness to voxels
% Works only with SPM5
% j.diedrichsen@bangor.ac.uk
name=[];
z=[];
vararginoptions(varargin,{'name','z'});

% Get the sample image for orientation and resolution
if (nargin<2 || isempty(Vol))
    Vol=spm_select(1,'image','Select volume for determining resolution and orientation');
end;
if ischar(Vol)
    Vol=spm_vol(Vol);
end;

% Make the data matrix for the ROI.
% Right now I use only the nearest neighbor approach:
% Every voxel which is nearest neighbor to sampled locations is set to 1
X=zeros(Vol.dim);
numCoords=size(R.data,1);
Coords=[R.data';ones(1,numCoords)];
vox=round((inv(Vol.mat)*Coords)');     % Go into voxel space
j=find(vox(:,1)>0 & vox(:,1)<=Vol.dim(1) & ...
    vox(:,2)>0 & vox(:,2)<=Vol.dim(2) & ...
    vox(:,3)>0 & vox(:,3)<=Vol.dim(3));

if (~isempty(z))
    for i=1:length(j)
        X(round(vox(j(i),1)),round(vox(j(i),2)),round(vox(j(i),3)))=z(j(i));
    end;
else
    for i=1:length(j)
        X(round(vox(j(i),1)),round(vox(j(i),2)),round(vox(j(i),3)))=1;
    end;
end;

% Now Make the new file name and set data type and save
[a,b,ext,d]=spm_fileparts(Vol.fname);
if (isempty(name))
    Vol.fname=[R.name '.nii'];
else
    Vol.fname=name;
end;
Vol.dt(1)=2;                  % data type is uint8
% Vol.descip=R.name;
Vol.pinfo=[1 0 0]';
Vol=spm_create_vol(Vol);
for z=1:size(X,3)
    spm_write_plane(Vol,squeeze(X(:,:,z)),z);
end;

% recast the old region into the new voxel space
a=find(X>0);
[i,j,k]=ind2sub(size(X),a);
[x,y,z]=spmj_affine_transform(i,j,k,Vol.mat);
R.data=[x y z];