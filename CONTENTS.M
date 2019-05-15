% Lightweight and simple Region toolbox for SPM5
% Contains functions for basic region of interest analysis 
% regions can be defined based on simple geometrical forms or images 
% and are kept in data structures, rather than as image
% All references to locations are in (mm), such that they are independent 
% from a certain voxel geometry. 
% Version 1.0 9/08/08
% Joern Diedrichsen (j.diedrichsen@bangor.ac.uk)
% Timothy Vertstynen 
% --------------------------------------------------------------
% 1. Data Structure functions 
%  region_getregions  - GUI for creating ROIs
%  region_calcregions - Calculates the locations in a regions and stores them as R.data  
%  region_getdata     - Gets the values from a series of Image files for a series of Regions 
%  region_getts       - Gets Raw, predicted, adjusted and residual time
%                       series from an SPM fr a series of regions 
%  region_getirf      - Extracts the ts(using getts) and gets the evoked
%                       response for all events
%  region_saveasimg   - Saves a certain region as an image
%  region_deformation - Deforms regions into individual space over a
%                       linear or non-linear transformation 
% 
% EXAMPLE OF USAGE 
% % get the regions through the GUI, or hardcode 
% R=region_getregions; 
% 
% % Calculate the region 
% R=region_calcregions(R); 
% 
% % Deform the regions to the space of one individual 
% R1=region_deformation(R,'mysubj_sn.mat'); 
% 
% % Extracts the time-series from an SPM and reestimates the betas 
% [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R1); 
% 
% % OR just get the raw betas 
% [beta] = region_getdata(SPM.xVbeta,R1); 
% 
% % Save the region 1 as an image for checking or illustration purposes 
% region_saveasimg(R{1},'myimage.img');   % Saves the region
