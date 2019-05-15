function R=region(varargin)
% function R=region(varargin)
%   region constructor 
% Region constructor
% region('point',[x,y,z],'name');
% region('sphere',[x,y,z],4);
% region('surf_circle',nodenum,radius,white,pial,topo,varargin); 
%       if Radius is [mm P], then it uss maximally a radius of mm, but
%       until it has found P voxels 
% region('surf_nodes',vertices,white,pial,linedef,image); 
%       Selected nodes from a surface - voxels touching between white and
%       pial surface from image are included
% region('roi_image',file,value,name); 
switch(varargin{1})
    case 'point'
        R.type='point';
        R.location=varargin{2};
        if (nargin>2)
            R.name=varargin{3};
        else 
            R.name=sprintf('%2.0f %2.0f %2.0f',varargin{2}); 
        end
    case 'sphere'
        R.type='sphere';
        R.location=varargin{2};
        R.radius=varargin{3};
        if (nargin>3)
            R.name=varargin{4};
        end
    case 'surf_circle' 
        R.type='surf_circle' 
        R.location=varargin{2}; 
        R.radius=varargin{3}; 
        R.white=varargin{4}; 
        R.pial=varargin{5}; 
        R.topo=varargin{6}; 
    case 'surf_nodes' 
        R.type='surf_nodes' 
        R.location=varargin{2}; 
        R.white=varargin{3}; 
        R.pial=varargin{4}; 
        R.linedef = varargin{5}; 
        R.image = varargin{6}; 
    case 'roi_image'
        R.type='roi_image';
        R.file=varargin{2}; 
        R.value=varargin{3};         
        if (nargin>3)
            R.name=varargin{4};
        end; 
    case 'image'
        R.type='image';
        R.file=varargin{2}; 
        R.threshold=varargin{3};         
    otherwise
        error('unknown region type');
end;
