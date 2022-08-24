# Region: Lightweight and simple Region toolbox for Matlab/SPM
Contains functions for basic region of interest analysis 
regions can be defined based on simple geometrical forms or images 
and are kept in data structures, rather than as image
All references to locations are in (mm), such that they are independent 
from a certain voxel geometry. 

Authors:  Joern Diedrichsen (joern.diedrichsen@googlemail.com)
Timothy Vertstynen 

* region			  - create region structure 
* region_calcregions - Calculates the locations in a regions and stores them as R.data  
* region_getdata     - Gets the values from a series of Image files for a series of Regions 
* region_getts       - Gets Raw, predicted, adjusted and residual time  series from an SPM fr a series of regions 
* region_getirf      - Extracts the ts(using getts) and gets the evoked response for all events
* region_saveasimg   - Saves a certain region as an image
* region_deformation - Deforms regions into individual space over a non-linear transformation 


### Examples of usage
```
% generate the regions (or make structures by hand)  
R{1}=region('roi_image',file,value,name)
R{2}=region('surf_nodes',nodes,white,pial,linedef,vol,name)
....

% Calculate the region 
R=region_calcregions(R); 

% Deform the regions to the space of one individual 
R1=region_deformation(R,'mysubj_sn.mat'); 

% Save the region 1 as an image for checking or illustration purposes 
region_saveasimg(R{1},'myimage.img');   % Saves the region

% Extracts the time-series from an SPM and reestimates the betas 
[y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R1); 

% OR just get the raw betas 
[beta] = region_getdata(SPM.xVbeta,R1); 

% Save extracted data as CIFTI file:  
cii = region_make_cifti(R1,Vol,'data',beta,'dtype','scalars');


```

### Dependencies: 
* SPM
* https://github.com/jdiedrichsen/dataframe
* https://github.com/Washington-University/cifti-matlab
