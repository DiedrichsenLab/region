function spmj_roi_brodman(name)
AAL=spm_vol('c:/projects/Atlas_templates/colin/aal.img');
SPM=spm_vol('c:/projects/Atlas_templates/spm/T1.img');
[Num,Name,Z]=textread('aal.txt','%d%s%d');
if (nargin<1 | isempty(strmatch(name,Name))
    fprintf('Pick region:\n');
    fprintf('%s\n',Name);
return;
end;
    region=strmatch(name,Name);
M=inv(AAL.mat)*SPM.mat;
for z=1:SPM.dim(3)
    Mz=M;
    Mz(3,4)=z*M(3,3)+M(3,4);
    AALv(:,:,z)=spm_slice_vol(AAL,Mz,[SPM.dim(1) SPM.dim(2)],0);
end;
ROI=SPM;ROI.fname=name;ROI.private.hdr.dime.pixdim(1)=2;
X=(AALv==region);
spm_write_vol(ROI,X);
