function mesh = make_square_mesh(h,type)
% MAKE_SQUARE_MESH creates a square mesh
% INPUT
%   ne1d: number of element in each dimension
%   type (optional): options are 'structured' and 'unstructured'
% OUTPUT
%   mesh: mesh structure. The boundary edge groups are as follows:
%         1 = left; 2 = right; 3 = bottom; 4 = top
if (nargin < 1)
    h = 1;
end
if (nargin < 2)
    type = 'structured';
end

switch type
    case 'structured'
        ne1d = ceil(1.0/h);
        x1 = linspace(0,1,ne1d+1);
        [xx1,xx2] = ndgrid(x1,x1);
        coord = [xx1(:), xx2(:)];
        tri = delaunay(coord);
        
        %ind = reshape(1:(ne1d+1)^2,[ne1d+1,ne1d+1]);
        %quad = zeros(ne1d^2,4);
        %quad(:,1) = reshape(ind(1:end-1,1:end-1),[],1);
        %quad(:,2) = reshape(ind(2:end,1:end-1),[],1);
        %quad(:,3) = reshape(ind(1:end-1,2:end),[],1);
        %quad(:,4) = reshape(ind(2:end,2:end),[],1);
        %tri = [quad(:,[1,2,3]); quad(:,[4,3,2])];
    case 'unstructured'
        figh = figure;
        fd=@(p) drectangle(p,0,1,0,1);
        [coord,tri]=distmesh2d(fd,@huniform,h,[0,0;1,1],[0,0;0,1;1,0;1,1]);
        close(figh);
    otherwise
        error('unsupported type');
end

mesh.coord = coord;
mesh.tri = tri;

% create boundary edge groups
edge = [tri(:,2), tri(:,3)
        tri(:,3), tri(:,1)
        tri(:,1), tri(:,2)];
tol = 1e-6;

% mid edge coordinate
xe = reshape(mean(reshape(coord(edge(:),:),[size(edge),2]),2),size(edge));

% left
ii = abs(xe(:,1) - 0.0) < tol;
mesh.bgrp{1} = sortrows(edge(ii,:),1);

% right
ii = abs(xe(:,1) - 1.0) < tol;
mesh.bgrp{2} = sortrows(edge(ii,:),1);

% bottom 
ii = abs(xe(:,2) - 0.0) < tol;
mesh.bgrp{3} = sortrows(edge(ii,:),1);

% top
ii = abs(xe(:,2) - 1.0) < tol;
mesh.bgrp{4} = sortrows(edge(ii,:),1);

end