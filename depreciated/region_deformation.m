function R=region_deformation(R,paramfile,varargin)
% function R=spmj_deformation(R,paramfile,varargin)
% Deforms regions of interest from atlas into inidividual space 
% j.diedrichsen@bangor.ac.uk 
if (~iscell(R))
    Rnew{1}=R;
    R=Rnew;
end;

if ( nargin<2 || isempty(paramfile))
    paramfile=spm_select(1,'.mat','Select deformation map'); 
end; 

[Def,mat]=spmdefs_get_sn2def(paramfile); 

for r=1:length(R)
    [R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3)]=spmdefs_transform(Def,mat,R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3)); 
end;

if (exist('Rnew')) 
    R=R{1};
end;
