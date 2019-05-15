 function R=region_getregions(varargin) 
% function R=region_getregions(varargin) 
% GUI interface for creating regions 
% OUTPUT: a cell array of structures, one for each region 
%   Point: Just a single location 
%           R.type='point';
%           R.location(s)=[x,y,z];
%           R.name 
%           --------------------
%   Sphere: A sphere of a certain radius around a center location 
%           R.type='sphere'
%           R.location(s)=[x,y,z];
%           R.radius=x;
%           -------------------- 
%   Box: rectangular box around a certain center location 
%           R.type='box'
%           R.location(s)=[x,y,z];
%           R.size=[width,depth,height];
%           --------------------
%   Image: An Image file, with the voxels of interest set to > 0
%           R.type='image'
%           R.file(s)=filenames 
%           R.theshold: Threshold applied to image to get mask 
%           --------------------
%   ROIImage: An Image file, with voxel values indicating the number of the
%            region. 
%           R.type='image'
%           R.file=filenames 
%           R.value=value(s) in ROIImage.  
%           --------------------
%   Cluster: A group of connected super-threshold voxels 
%           R.type='cluster'
%           R.location(s)=Seed location 
%           R.file=Contrast file 
%           R.threshold=Threshold applied to contrast file  
%           --------------------
%   Paint: Paint file (caret), one region per entry is created 
%           R.type='paint'
%           R.location=[filename];
%           --------------------
%    Node: Node(s) on a caret Surface  
%           R.type='node'
%           R.location=nodenum
% 
%  When multiple locations, files, valus etc are given, multiple regions
%  are generated 
% 
% j.diedrichsen, 2007 

types={'done','point','sphere','box','image','roi_image','cluster','paint','vertex','vertexcircle'};
R={};p='x';r=1;
while (~strcmp(p,'done'))
    [p,YPos] = spm_input('Select region type',1,'m',types,types,1);
    type=p{1};
    switch (type)
        case 'point'
            location=spm_input('Location',2,'r','[0 0 0]',[inf 3]);
            for i=1:size(location,1) 
                R{r}.type=type; 
                R{r}.location=location(i,:);
                R{r}.name=spm_input(sprintf('Name %d',i),2+i,'s',sprintf('region_%d',r));
                r=r+1; 
            end;            
        case 'sphere'
            location=spm_input('Center(s)',2,'r','[0 0 0]',[inf 3]);           
            radius=spm_input('Radius',3,'r','3',[1]);           
            for i=1:size(location,1) 
                R{r}.type=type; 
                R{r}.location=location(i,:);
                R{r}.radius=radius; 
                R{r}.name=spm_input(sprintf('Name %d',i),2+i,'s',sprintf('region_%d',r));
                r=r+1; 
            end;
        case 'box'
            location=spm_input('Center(s)',2,'r','[0 0 0]',[inf 3]);           
            radius=spm_input('Size',3,'r','[5 5 5]',[inf 3]);           
            for i=1:size(location,1) 
                R{r}.type=type; 
                R{r}.location=location(i,:);
                R{r}.radius=radius; 
                R{r}.name=spm_input(sprintf('Name %d',i),2+i,'s',sprintf('region_%d',r));
                r=r+1; 
            end;
        case 'image'
            files=spm_select(inf,'image','Select ROI images');
            threshold=spm_input('Threshold',3,'r','0.1',[1]);
            for i=1:size(files,1) 
                R{r}.type=type; 
                R{r}.file=files(i,:);
                R{r}.threshold=threshold; 
                R{r}.name=spm_input(sprintf('Name %d',i),2+i,'s',R{r}.file);
                r=r+1; 
            end;
        case 'roi_image'
            file=spm_select(1,'image','Select ROI images');
            values=spm_input('values',2,'r','1',[inf]);
            for i=1:size(values,1) 
                R{r}.type=type; 
                R{r}.file=file;
                R{r}.value=values(i); 
                R{r}.name=spm_input(sprintf('Name %d',i),2+i,'s',sprintf('region_%d',values(i)));
                r=r+1; 
            end;
        case 'cluster'
            files=spm_select(1,'image','Select Contrast Image');
            location=spm_input('Center',2,'r','[0 0 0]',[inf 3]);
            threshold=spm_input('Threshold',3,'r','0',[1]);
            for i=1:size(location,1) 
                R{r}.type=type; 
                R{r}.file=files;
                R{r}.location=location(i,:);
                R{r}.threshold=threshold; 
                R{r}.name=spm_input(sprintf('Name %d',i),2+i,'s',sprintf('region_%d',r));
                r=r+1; 
            end;
        case 'paint'
            R{r}.type=type;
            R{r}.file=spm_get(1,'*.paint','Select ROI images');
            r=r+1;
        case 'vertex'
            R{r}.type=type;
            R{r}.location=spm_input('Node Number','+1','e');
            r=r+1;
    end;
end;
