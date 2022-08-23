function mycifti=region_make_cifti(R,Vol,varargin)
% function R=region_make_cifti(R,Data,varargin)
% Creates a dense cifti file (with BrainModel axis) for the data for a set
% of regions. 
% By default all are made into volume-based parcels. 
% Surface functionality to be added 
% INPUT: 
%       R: Cell array of regions 
%       Vol: Volume information (spm_vol): all need to refer to the same volume space
% VARGINOPTION
%       'data', data: Cell array of data matrices (voxels x measures)
%                     If no data is given, it just makes the first
%                     brainmodel axis. 
%       'dtype', {'scalars','series','labels'}: Type of axis for the data
%       'dnames', cell: Cell array of names for dtype=scalar 

% August 23,22 joern.diedrichsen@googlemail.com 

data=[];
dtype='scalars'; 
dnames={}; 
struct={'CEREBELLUM','BRAIN_STEM','OTHER','THALAMUS_LEFT','THALAMUS_RIGHT'};
vararginoptions(varargin,{'data','dtype','dnames','struct'});

mycifti.metadata = []; 
mycifti.metadata = cifti_metadata_set(mycifti.metadata, 'Provenance', 'region_make_cifti');
mycifti.cdata = []; 
% make the brainmodel axis
bmaxis.type = 'dense'; % First 
bmaxis.vol.dims = Vol.dim;  % Image dimension 
bmaxis.vol.sform = Vol.mat; % Image affine 

% Now loop over regions and make them into brain models
start = 1; 
for r=1:length(R)
    bm.start=start; 
    bm.count = size(R{r}.data,1); 
    bm.struct = struct{r};  % Unfortunately, we can't set the structure name to anything arbitary.
    bm.type = 'vox'; 
    [i,j,k] = spmj_affine_transform(R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3),inv(Vol.mat)); 
    bm.voxlist = [i j k]'-1;  % voxlist is zero-based indices 
    start = start + bm.count; 
    bmaxis.models{r}=bm;
    if ~isempty(data)
        if size(data{r},1)~=bm.count
            error('Number of columns does not correspond to region %d',r); 
        end
        mycifti.cdata=[mycifti.cdata;data{r}];
    end
end
bmaxis.length=size(mycifti.cdata,1);

% Generate the correct data axis 
if ~isempty(data)
    daxis.type=dtype; 
    daxis.length = size(mycifti.cdata,2);
    if (isempty(dnames) && ~strcmp(dtype,'series')) 
        for i=1:daxis.length
            dnames{i}=sprintf('row %d',i); 
        end
    end
    switch(dtype)
        case 'series'
            daxis.seriesStart=0; 
            daxis.seriesStep =1.0; 
            daxis.seriesUNIT='SECOND';
        case 'scalars'
            for i=1:daxis.length
                daxis.maps(i).name = dnames{i}; 
                daxis.maps(i).metadata=[]; 
            end
        case 'labels'
            for i=1:daxis.length
                daxis.maps(i).name = dnames{i}; 
                daxis.maps(i).metadata=[]; 
                % for c=0:max(data)
                % daxis.maps(i).table(c).name=sprintf('label %c'); 
                % daxis.maps(i).table(c).key=c
                % daxis.maps(1).table(c).rgba = []
            end 
        otherwise 
            error('Unknown dtype: {scalars,series,labels} allowed'); 
    end
else 
    daxis =[]; 
end

% Assemble mycifti and return 
mycifti.diminfo={bmaxis,daxis}; 

    


