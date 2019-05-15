function R=region_calcregions(R,varargin)
% function R=spmj_calcregions(R,varargin)
% (Re-) calculates the locations for regions of interest
% adds .data field to reach region structure for locations
% VARGINOPTIONS
% 'gridwidth',2             : in mm the distance for different locations to
%                             sample(set to voxel size)
% 'voxelspace',file/volume  : Calculate the regions in the pixel space of a particular image
% 'exclude','all'        : Make all regions mutual exclusive
% 'exclude',[1 2;3 2]..  : List the numbers of the ROIs which should be
%                             made exclusive
% 'exlcude_thres',0.8
% If the region has a structure in field 'mask', the region will be masked
% by that image
% j.diedrichsen@bangor.ac.uk
% 31/1/2011: Fixed a number of issues with the masking of surface-based ROI
% 25/03/2017: Maedbh King - added check for .gii surface files in 'surf_nodes' and
% 'surf_circle'

gridwidth=2;        % 2mm grid width for sample_image
voxelspace=[];
mask=[];
exclude=[];       % Define the ROIs in a way that is mutual exclusive
exclude_thres=1;

vararginoptions(varargin,{'gridwidth','voxelspace','exclude','exclude_thres'});

if (~iscell(R))
    Rnew{1}=R;
    R=Rnew;
end;

if (~isempty(voxelspace))
    if (isstruct(voxelspace))
        Vpixel=voxelspace;
    elseif (ischar(voxelspace))
        Vpixel=spm_vol(voxelspace);
    else
        error('voxelspace must be either volume of file name');
    end;
end;

for c=1:length(R)
    
    %Printing progress
    if ~mod(c,10); fprintf('%d out of %d \n',c,length(R)); end
    
    if (~isempty(R{c}))
        switch  (R{c}.type)
            case 'point'
                for i=1:size(R{c}.location,1);
                    R{c}.data=R{c}.location(i,:);
                end;
            case 'sphere'
                % make sphere
                if (isempty(voxelspace))
                    left=-R{c}.radius.*[1 1 1];
                    right=+R{c}.radius.*[1 1 1];
                    s=R{c}.radius/(ceil(R{c}.radius/gridwidth));
                    [x,y,z]=meshgrid(left(1):s:right(1),left(2):s:right(2),left(3):s:right(3));
                    x=x(:);y=y(:);z=z(:);
                    rad=sqrt(x.^2+y.^2+z.^2);
                    indx=find(rad<=R{c}.radius);
                    x=x(indx);y=y(indx);z=z(indx);
                    % add for all locations
                    cent=R{c}.location;
                    R{c}.data=[x+cent(1) y+cent(2) z+cent(3)];
                else
                    [x,y,z]=meshgrid(1:Vpixel.dim(1),1:Vpixel.dim(2),1:Vpixel.dim(3));
                    [x,y,z]=spmj_affine_transform(x(:),y(:),z(:),Vpixel.mat);
                    cent=R{c}.location;
                    rad=sqrt((x-cent(1)).^2+(y-cent(2)).^2+(z-cent(3)).^2);
                    indx=find(rad<=R{c}.radius);
                    x=x(indx);y=y(indx);z=z(indx);
                    % add for all locations
                    R{c}.data=[x y z];
                    
                end;
            case 'sphere_voxel' % Tobias addition: Radius is defined in Voxels!
                
                %----get the center of the volume
                volume_center= round([Vpixel.dim(1), Vpixel.dim(2), Vpixel.dim(3)]/2);
                %----calculate the distance from the volume center to all the coordinates in the volume
                rad=sqrt((x(:)-volume_center(1)).^2+(y(:)-volume_center(2)).^2+(z(:)-volume_center(3)).^2);
                %----transform the volume coordinates in the image space
                [x,y,z]=spmj_affine_transform(x(:),y(:),z(:),Vpixel.mat);
                %----find the coords within the radius of the sphere
                indx=find(rad<=R{c}.radius);
                %----select the coordinates accordingly and move the sphere to the origin
                xyz=[x y z];
                xyz=[x(indx,:) y(indx,:) z(indx,:)]-repmat(xyz(find(rad==0),:), size(indx,1), 1);
                %----move the sphere to the expected location in the image
                xyz= xyz+repmat(R{c}.location, size(indx,1), 1);
                %----add location to the output sphere
                R{c}.data=[xyz(:,1) xyz(:,2) xyz(:,3)];
                
                %
            case 'box'
                % make Box
                left=-R{c}.size/2;
                right=+R{c}.size/2;
                [x,y,z]=meshgrid(left(1):gridwidth:right(1),left(2):gridwidth:right(2),left(3):gridwidth:right(3));
                x=x(:);y=y(:);z=z(:);
                cent=R{c}.location;
                R{c}.data=[x+cent(1) y+cent(2) z+cent(3)];
            case 'image'
                % Load images
                V=spm_vol(R{c}.file);
                X=spm_read_vols(V);
                [x,y,z]=ind2sub(size(X),find(X>R{c}.threshold));
                [x,y,z]=spmj_affine_transform(x,y,z,V.mat);
                R{c}.data=[x y z];
            case 'roi_image'
                % Load images
                V=spm_vol(R{c}.file);
                X=spm_read_vols(V);
                if (isfield(R{c},'image'))                                           % Tobias addition: Mask defined?
                    V.mask=spm_read_vols(spm_vol(R{c}.image));                       % Keep backwards compatible
                    [x,y,z]=ind2sub(size(X),find(round(X)==R{c}.value & V.mask~=0)); %
                else
                    [x,y,z]=ind2sub(size(X),find(round(X)==R{c}.value)); % Keep
                end;
                [x,y,z]=spmj_affine_transform(x,y,z,V.mat);
                R{c}.data=[x y z];
            case 'cluster'
                % Load image
                V=spm_vol(R{c}.file);
                X=spm_read_vols(V);
                Xindx=find(X>R{c}.threshold);
                [x,y,z]=ind2sub(size(X),Xindx);
                A=spm_clusters([x y z]');
                cntr=round(inv(V.mat)*[R{c}.location';1]);
                i=findrow([x y z],cntr(1:3)');
                if (isempty(i))
                    R{c}.data=[];
                    warning(sprintf('Region %s does not contain any data',R{c}.name));
                else
                    clindx=find(A==A(i));
                    [x,y,z]=spmj_affine_transform(x(clindx),y(clindx),z(clindx),V.mat);
                    R{c}.data=[x y z];
                end;
            case 'surf_nodes'  % Set of surface nodes
                % open with gifti or caret_load
                if strcmp(R{c}.white(end-3:end),'.gii'), % MK
                    c1=gifti(R{c}.white);
                    c1_data=c1.vertices;
                else
                    c1=caret_load(R{c}.white); 
                    c1_data=c1.data;
                end
                if strcmp(R{c}.pial(end-3:end),'.gii'),
                    c2=gifti(R{c}.pial);
                    c2_data=c2.vertices;
                else
                    c2=caret_load(R{c}.pial);
                    c2_data=c2.data;
                end
                
                V=spm_vol(R{c}.image);
                V.mask=spm_read_vols(V)~=0; % Define mask image: This is fixed 31/1
                
                % Find the coordinates between surfaces and all touching
                % voxels
                allcoords=surfing_nodeidxs2coords(c1_data(R{c}.location,:)',c2_data(R{c}.location,:)',[],R{c}.linedef);
                alllinvoxidxs=surfing_coords2linvoxelidxs(allcoords,V);
                indices=unique(alllinvoxidxs(~isnan(alllinvoxidxs)));
                N=length(indices);
                R{c}.linvoxidxs=indices;
                
                % calculate the Eucledian coordinates of voxel centers
                [i,j,k]=ind2sub(V.dim,double(indices));
                [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
                R{c}.data=[x,y,z];
                
                % Track the weight of each voxels: how many vertices are in
                % the voxel?
                weight=zeros(size(indices));
                for i=1:N
                    weight(i)=sum(sum(alllinvoxidxs==indices(i)));
                end;
                R{c}.weight=weight;
                
                % If flat surface given: Compute projection for each voxel onto
                % the surface and determine the cortical depth. Express the
                % projection both in barycentric coordinates and in average
                % flat coordinates
                if (isfield(R{c},'flat'))
                    if strcmp(R{c}.flat(end-3:end),'.gii'), % MK
                        c3=gifti(R{c}.flat);
                        c3_data=c3.vertices;
                    else
                        c3=caret_load(R{c}.flat);
                        c3_data=c3.data;
                    end
                    % Load topology for accurate mapping
                    if strcmp(R{c}.topo(end-3:end),'.gii'), % MK
                        topo=gifti(R{c}.topo);
                        topo_data=topo.vertices;
                    else
                        topo=caret_load(R{c}.topo);
                        topo_data=topo.data;
                    end
                    ca=(c1_data+c2_data)/2; % Average coordinateds
                    TR=triangulation(topo_data,ca(:,1),ca(:,2),ca(:,3));                % Create triangulation
                    
                    % Preallocate
                    R{c}.flatcoord=zeros(N,3)*NaN;              % Flat X,y coordinates
                    R{c}.tile=zeros(N,1)*NaN;                   % Triangle
                    R{c}.barycentric=zeros(N,3)*NaN;            % barycentric weights
                    R{c}.depth=zeros(N,1)*NaN;                  % Depth 0:pial 1: white
                    
                    numCols = size(alllinvoxidxs,2);
                    
                    % Loop over all the voxels that you found
                    for i=1:N
                        p0        = R{c}.data(i,:);                             % XYZ of voxel
                        indx      = find(alllinvoxidxs==indices(i));            % Find the places where the voxel was assigned
                        [row,col] = ind2sub(size(alllinvoxidxs),indx);          % Get the row (vertex) and column (depth)
                        nodes     = R{c}.location(unique(row));                 % Find the nodes that used this voxel
                        tri       = find(any(ismember(topo_data,nodes),2));     % Find the triangles that are involved
                        numTri    = length(tri);
                        baryc     = cartesianToBarycentric(TR,tri,repmat(R{c}.data(i,:),numTri,1));
                        baryc(baryc<0)=0;
                        baryc     = bsxfun(@rdivide,baryc,sum(baryc,2));        % Force the points to lie on the triangle
                        cart      = barycentricToCartesian(TR,tri,baryc);       % Find the closest point
                        sqdist    = sum(bsxfun(@minus,cart,p0).^2,2);           % Calcualte squared distance
                        [~,best]  = min(sqdist);
                        
                        % Now determine average coordinates
                        bestTri   = topo_data(tri(best),:);
                        p1        = baryc(best,:)*c1_data(bestTri',:);                    % XYZ coords on white
                        p2        = baryc(best,:)*c2_data(bestTri',:);                    % XYZ coords on pial
                        v2        = p0-p2;
                        v1        = p1-p2;
                        R{c}.depth(i,:) = (v2*v1')/(v1*v1');        % Take the mean of the depth: 0 at the pial (C2) and 1 at the white surface(C1)
                        R{c}.flatcoord(i,:) = baryc(best,:)*c3_data(bestTri',:);          % Coordinates on the flat surface
                    end;
                end;
            case 'surf_circle' % Does a surface-based circular ROI at a certain node
                % open with gifti or caret_load
                if strcmp(R{c}.white),
                    c1=gifti(R{c}.white);
                    c1_data=c1.vertices;
                else
                    c1=caret_load(R{c}.white);
                    c1_data=c1.data;
                end
                if strcmp(R{c}.pial),
                    c2=gifti(R{c}.pial);
                    c2_data=c2.vertices;
                else
                    c2=caret_load(R{c}.pial);
                    c2_data=c2.data;
                end
                if strcmp(R{c}.topo),
                    to=gifti(R{c}.topo);
                    to_data=to.vertices;
                else
                    to=caret_load(R{c}.topo);
                    to_data=to.data;
                end
                V=spm_vol(R{c}.image);
                V.mask=spm_read_vols(V)~=0;
                [indices,VMIN,VMAX,vORr]=surfing_voxelselection(c1_data',c2_data',to_data',R{c}.radius,V,R{c}.location);
                [i,j,k]=ind2sub(V.dim,double(indices{1}'));
                [x,y,z]=spmj_affine_transform(i,j,k,V.mat);
                R{c}.data=[x,y,z];
        end;
    end;
    % If necessary, apply mask.
    if (isfield(R{c},'mask'))
        A=R{c}.data;
        Vmask=R{c}.mask;
        [i,j,k]=spmj_affine_transform(A(:,1),A(:,2),A(:,3),inv(Vmask.mat));
        D=spm_sample_vol(Vmask,i,j,k,0);
        indx=find(D>0);
        R{c}.data=A(indx,:);
    end;
end;

% Deal with making ROIs mutual exclusive
if (~isempty(exclude))
    
    % Start with including all voxels
    for i=1:length(R)
        R{i}.excl=zeros(size(R{i}.data,1),1);
    end;
    
    % Deal with 'all" option
    if strcmp(exclude,'all')
        exclude=[];
        for i=1:length(R)
            for j=i:length(R)
                exclude=[exclude;i j];
            end;
        end;
    end;
    
    % Make voxels mutually exclusive
    for i=1:size(exclude,1)
        j=exclude(i,1);
        k=exclude(i,2);
        EQ=bsxfun(@eq,R{j}.linvoxidxs,R{k}.linvoxidxs');
        rows=find(sum(EQ,2)>0);
        cols=find(sum(EQ,1)>0)';
        for v=1:length(rows)
            wj=R{j}.weight(rows(v));
            wk=R{k}.weight(cols(v));
            if (wj/(wj+wk)>exclude_thres)       % Kill from region k
                R{k}.excl(cols(v))=1;
            elseif (wk/(wj+wk)>exclude_thres)   % Kill from region j
                R{j}.excl(rows(v))=1;
            else                                % Kill from both
                R{k}.excl(cols(v))=1;
                R{j}.excl(rows(v))=1;
            end;
        end;
    end;
    
    % Now exclude those voxels
    for i=1:length(R)
        R{i}.data=R{i}.data(~R{i}.excl,:);
        R{i}.numexcl=sum(R{i}.excl);
    end;
    
end;

if (exist('Rnew'))
    R=R{1};
end;
